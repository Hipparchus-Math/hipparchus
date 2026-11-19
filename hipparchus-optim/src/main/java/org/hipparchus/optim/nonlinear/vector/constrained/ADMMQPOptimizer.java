/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.LocalizedOptimFormats;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Alternating Direction Method of Multipliers Quadratic Programming Optimizer.
 * \[
 *  min \frac{1}{2} X^T Q X + G X a\\
 *  A  X    = B_1\\
 *  B  X    \ge B_2\\
 *  l_b \le C X \le u_b
 * \]
 * Algorithm based on paper:"An Operator Splitting Solver for Quadratic Programs(Bartolomeo Stellato, Goran Banjac, Paul Goulart, Alberto Bemporad, Stephen Boyd,February 13 2020)"
 * @since 3.1
 */

public class ADMMQPOptimizer extends QPOptimizer {

    /** Algorithm settings. */
    private ADMMQPOption settings;

    /** Equality constraint (may be null). */
    private LinearEqualityConstraint eqConstraint;

    /** Inequality constraint (may be null). */
    private LinearInequalityConstraint iqConstraint;

    /** Boundary constraint (may be null). */
    private LinearBoundedConstraint bqConstraint;

    /** Objective function. */
    private QuadraticFunction function;

    /** Problem solver. */
    private final ADMMQPKKT solver;

    /** Problem convergence checker. */
    private ADMMQPConvergenceChecker checker;

    /** Convergence indicator. */
    private boolean converged;

    /** Current step size. */
    private double rho;

    /** Simple constructor.
     * <p>
     * This constructor sets all {@link ADMMQPOption options} to their default values
     * </p>
     */
    public ADMMQPOptimizer() {
        settings  = new ADMMQPOption();
        solver    = new ADMMQPKKT();
        converged = false;
        rho       = 0.1;
    }

    /** {@inheritDoc} */
    @Override
    public ConvergenceChecker<LagrangeSolution> getConvergenceChecker() {
        return checker;
    }

    /** {@inheritDoc} */
    @Override
    public LagrangeSolution optimize(OptimizationData... optData) {
        return super.optimize(optData);
    }

    /** {@inheritDoc} */
    @Override
    protected void parseOptimizationData(OptimizationData... optData) {
        super.parseOptimizationData(optData);
        for (OptimizationData data: optData) {

             if (data instanceof ObjectiveFunction) {
                function = (QuadraticFunction) ((ObjectiveFunction) data).getObjectiveFunction();
                continue;
            }

            if (data instanceof LinearEqualityConstraint) {
                eqConstraint = (LinearEqualityConstraint) data;
                continue;
            }
            if (data instanceof LinearInequalityConstraint) {
                iqConstraint = (LinearInequalityConstraint) data;
                continue;
            }

            if (data instanceof LinearBoundedConstraint) {
                bqConstraint = (LinearBoundedConstraint) data;
                continue;
            }

            if (data instanceof ADMMQPOption) {
                settings = (ADMMQPOption) data;
            }

            if (data instanceof MatrixDecompositionTolerance) {
                final MatrixDecompositionTolerance tolerance = (MatrixDecompositionTolerance) data;
                solver.updateDecompositionEpsilon(tolerance.getEpsMatrixDecomposition());
            }

        }
        // if we got here, convexObjective exists
        int n = function.dim();
        if (eqConstraint != null) {
            int nDual = eqConstraint.dimY();
            if (nDual >= n) {
                throw new MathIllegalArgumentException(LocalizedOptimFormats.CONSTRAINTS_RANK, nDual, n);
            }
            int nTest = eqConstraint.getA().getColumnDimension();
            if (nDual == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_NOT_ALLOWED);
            }
            MathUtils.checkDimension(nTest, n);
        }

    }

    /** {@inheritDoc} */
    @Override
    public LagrangeSolution doOptimize() {
        final int n = function.dim();
        int me = 0;
        int mi = 0;
        int mb = 0;
        int rhoUpdateCount = 0;

        //PHASE 1 First Solution


        //QUADRATIC TERM
        RealMatrix H = function.getP();
        //GRADIENT
        RealVector q = function.getQ();


        //EQUALITY CONSTRAINT
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
        }
        //INEQUALITY CONSTRAINT
        if (iqConstraint != null) {
            mi = iqConstraint.dimY();
        }
        //BOUNDED CONSTRAINT
        if (bqConstraint != null) {
            mb = bqConstraint.dimY();
        }

        RealVector lb = new ArrayRealVector(me + mi + mb);
        RealVector ub = new ArrayRealVector(me + mi + mb);

        //COMPOSE A MATRIX AND LOWER AND UPPER BOUND
        RealMatrix A = new Array2DRowRealMatrix(me + mi + mb, n);
        if (eqConstraint != null) {
            A.setSubMatrix(eqConstraint.jacobian(null).getData(), 0, 0);
            lb.setSubVector(0,eqConstraint.getLowerBound());
            ub.setSubVector(0,eqConstraint.getUpperBound());
        }
        if (iqConstraint != null) {
            A.setSubMatrix(iqConstraint.jacobian(null).getData(), me, 0);
            ub.setSubVector(me,iqConstraint.getUpperBound());
            lb.setSubVector(me,iqConstraint.getLowerBound());
        }

        if (mb > 0) {
            A.setSubMatrix(bqConstraint.jacobian(null).getData(), me + mi, 0);
            ub.setSubVector(me + mi,bqConstraint.getUpperBound());
            lb.setSubVector(me + mi,bqConstraint.getLowerBound());
        }

        checker = new ADMMQPConvergenceChecker(H, A, q, settings.getEps(), settings.getEps());

        //SETUP WORKING MATRIX
        RealMatrix Hw = H.copy();
        RealMatrix Aw = A.copy();
        RealVector qw = q.copy();
        RealVector ubw = ub.copy();
        RealVector lbw = lb.copy();
        RealVector x;
        if (getStartPoint() != null) {
            x = new ArrayRealVector(getStartPoint());
        } else {
            x = new ArrayRealVector(function.dim());
        }

        ADMMQPModifiedRuizEquilibrium dec = new ADMMQPModifiedRuizEquilibrium(H, A,q);

        if (settings.isScaling()) {
           //
            dec.normalize(settings.getEps(), settings.getScaleMaxIteration());
            Hw = dec.getScaledH();
            Aw = dec.getScaledA();
            qw = dec.getScaledQ();
            lbw = dec.getScaledLUb(lb);
            ubw = dec.getScaledLUb(ub);

            x = dec.scaleX(x.copy());

        }

        final ADMMQPConvergenceChecker checkerRho = new ADMMQPConvergenceChecker(Hw, Aw, qw, settings.getEps(), settings.getEps());
        //SETUP VECTOR SOLUTION

        RealVector z = Aw.operate(x);
        RealVector y = new ArrayRealVector(me + mi + mb);

        solver.initialize(Hw, Aw, qw, me, lbw, ubw, rho, settings.getSigma(), settings.getAlpha());
        RealVector xstar = null;
        RealVector ystar = null;
        RealVector zstar;

        while (iterations.getCount() <= iterations.getMaximalCount()) {
            ADMMQPSolution sol = solver.iterate(x, y, z);
            x = sol.getX();
            y = sol.getLambda();
            z = sol.getZ();
            //new ArrayRealVector(me + mi + mb);
            if (rhoUpdateCount < settings.getMaxRhoIteration()) {
                double rp       = checkerRho.residualPrime(x, z);
                double rd       = checkerRho.residualDual(x, y);
                double maxP     = checkerRho.maxPrimal(x, z);
                double maxD     = checkerRho.maxDual(x, y);
                boolean updated = manageRho(me, rp, rd, maxP, maxD);

                if (updated) {
                    ++rhoUpdateCount;
                }
            }


            if (settings.isScaling()) {

                xstar = dec.unscaleX(x);
                ystar = dec.unscaleY(y);
                zstar = dec.unscaleZ(z);

            } else {

                xstar = x.copy();
                ystar = y.copy();
                zstar = z.copy();

            }

            double rp        = checker.residualPrime(xstar, zstar);
            double rd        = checker.residualDual(xstar, ystar);
            double maxPrimal = checker.maxPrimal(xstar, zstar);
            double maxDual   = checker.maxDual(xstar, ystar);

            if (checker.converged(rp, rd, maxPrimal, maxDual)) {
                converged = true;
                break;
            }
            iterations.increment();

        }

        //SOLUTION POLISHING
        if (settings.isPolishing()) {
            ADMMQPSolution finalSol = polish(Hw, Aw, qw, lbw, ubw, x, y, z);
            if (settings.isScaling()) {
                xstar = dec.unscaleX(finalSol.getX());
                ystar = dec.unscaleY(finalSol.getLambda());
            } else {
                xstar = finalSol.getX();
                ystar = finalSol.getLambda();
            }
        }
        for (int i = 0; i < me + mi; i++) {
            ystar.setEntry(i, -ystar.getEntry(i));
        }

        return new LagrangeSolution(xstar, ystar, function.value(xstar));

    }

    /** Check if convergence has been reached.
     * @return true if convergence has been reached
     */
    public boolean isConverged() {
        return converged;
    }

    /** Polish solution.
     * @param H quadratic term matrix
     * @param A constraint coefficients matrix
     * @param q linear term matrix
     * @param lb lower bound
     * @param ub upper bound
     * @param x primal problem solution
     * @param y dual problem solution
     * @param z auxiliary variable
     * @return polished solution
     */
    private ADMMQPSolution polish(RealMatrix H, RealMatrix A, RealVector q, RealVector lb, RealVector ub,
                                  RealVector x, RealVector y, RealVector z) {

        List<double[]> Aentry    = new ArrayList<>();
        List<Double>  lubEntry   = new ArrayList<>();
        List<Double>  yEntry     = new ArrayList<>();

        // FIND ACTIVE ON LOWER BAND
        for (int j = 0; j < A.getRowDimension(); j++) {
            if (z.getEntry(j) - lb.getEntry(j) < -y.getEntry(j)) {  // lower-active

                Aentry.add(A.getRow(j));
                lubEntry.add(lb.getEntry(j));
                yEntry.add(y.getEntry(j));

            }
        }
        //FIND ACTIVE ON UPPER BAND
        for (int j = 0; j < A.getRowDimension(); j++) {
            if (-z.getEntry(j) + ub.getEntry(j) < y.getEntry(j)) { // lower-active

                Aentry.add(A.getRow(j));
                lubEntry.add(ub.getEntry(j));
                yEntry.add(y.getEntry(j));

            }

        }
        RealMatrix Aactive;
        RealVector lub;

        RealVector ystar;
        RealVector xstar = x.copy();
        //!Aentry.isEmpty()
        if (!Aentry.isEmpty()) {

            Aactive = new Array2DRowRealMatrix(Aentry.toArray(new double[0][]));
            lub = new ArrayRealVector(lubEntry.toArray(new Double[0]));
            ystar = new ArrayRealVector(yEntry.toArray(new Double[0]));
            solver.initialize(H, Aactive, q, 0, lub, lub,
                              settings.getSigma(), settings.getSigma(), settings.getAlpha());

            for (int i = 0; i < settings.getPolishIteration(); i++) {
                RealVector kttx = (H.operate(xstar)).add(Aactive.transpose().operate(ystar));
                RealVector ktty = Aactive.operate(xstar);
                RealVector b1 = q.mapMultiply(-1.0).subtract(kttx);
                RealVector b2 = lub.mapMultiply(1.0).subtract(ktty);
                ADMMQPSolution dxz = solver.solve(b1,b2);
                xstar = xstar.add(dxz.getX());
                ystar = ystar.add(dxz.getV());
            }

            return new ADMMQPSolution(xstar, null, y, A.operate(xstar));

        } else {
            return new ADMMQPSolution(x, null, y, z);
        }
    }

    /** Manage step size.
     * @param me number of equality constraints
     * @param rp primal residual
     * @param rd dual residual
     * @param maxPrimal primal vectors max
     * @param maxDual dual vectors max
     * @return true if rho has been updated
     */
    private boolean manageRho(int me, double rp, double rd, double maxPrimal, double maxDual) {
        boolean updated = false;
        if (settings.updateRho()) {

            // estimate new step size
            double rhonew = FastMath.min(FastMath.max(rho * FastMath.sqrt((rp * maxDual) / (rd * maxPrimal)),
                                                      settings.getRhoMin()),
                                         settings.getRhoMax());

            if ((rhonew > rho * 5.0) || (rhonew < rho / 5.0)) {

                rho = rhonew;
                updated = true;

                solver.updateSigmaRho(settings.getSigma(), me, rho);
            }
        }
        return updated;
    }

}
