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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Sequential Quadratic Programming Optimizer (extended version).
 *
 * <p>
 * Minimizes a nonlinear objective function subject to equality, inequality, and
 * box constraints using a Sequential Quadratic Programming method. This
 * implementation is inspired by the algorithm described in: "On the convergence
 * of a sequential quadratic programming method" by Klaus Schittkowski (1982).
 * </p>
 *
 * <p>
 * Supports: equality constraints, inequality constraints,
 * penalty function updates, line search strategies, BFGS Hessian update, and
 * augmented QP formulation using a relaxation variable.</p>
 *
 * @since 3.1
 */
public class SQPOptimizerS2 extends AbstractSQPOptimizer2 {

    /** Logger. */
    private final SQPLogger formatter = SQPLogger.defaultLogger();

    /**
     * Value of the equality constraints.
     */
    private RealVector eqEval;

    /**
     * Value of the inequality constraints.
     */
    private RealVector ineqEval;

    /**
     * Gradient of the objective function.
     */
    private RealVector J;

    /**
     * Hessian approximation.
     */
    private RealMatrix H;

    /**
     * Jacobian of the inequality constraints.
     */
    private RealMatrix JI;

    /**
     * Jacobian of the equality constraints.
     */
    private RealMatrix JE;

    /** Evaluation of the objective function. */
    private double functionEval;

    /** Current point. */
    private RealVector x;

    /**
     * {@inheritDoc}
     */
    @Override
    public LagrangeSolution doOptimize() {
        formatter.setEps(getSettings().getEps());
        formatter.logHeader();
        int me = 0;
        int mi = 0;

        //EQUALITY CONSTRAINT
        if (this.getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        //INEQUALITY CONSTRAINT
        if (this.getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
        }
        final int m = me + mi;

        double alpha;
        double rho = getSettings().getRhoCons(); //get initial value of slack variable from SQPOption

        if (this.getStartPoint() != null) {
            x = new ArrayRealVector(this.getStartPoint());
        } else {
            x = new ArrayRealVector(this.getObj().dim());
        }

        RealVector y = new ArrayRealVector(me + mi, 0.0);

        //all the function and constraint evaluation will be performed inside the penalty function
        MeritFunctionL2 penalty = new MeritFunctionL2(this.getObj(), this.getEqConstraint(), this.getIqConstraint(), x);

        LineSearch lineSearch = new LineSearch(getSettings().getEps(), 5, getSettings().getMu(), getSettings().getB(),
                                               getSettings().getMaxLineSearchIteration(), 2);
        H = MatrixUtils.createRealIdentityMatrix(x.getDimension());

        BFGSUpdater bfgs = new BFGSUpdater(H, 1.0e-11, getMatrixDecompositionTolerance().getEpsMatrixDecomposition());
        //INITIAL VALUES

        functionEval = penalty.getObjEval();
        if (this.getEqConstraint() != null) {
            eqEval = penalty.getEqEval();
        }
        if (this.getIqConstraint() != null) {
            ineqEval = penalty.getIqEval();
        }

        computeGradients();

        RealVector dx = new ArrayRealVector(x.getDimension());
        RealVector u = new ArrayRealVector(y.getDimension());
        penalty.update(J, JE, JI, x, y, dx, u);
        boolean augmented = true;


        boolean crit0 = false;
        boolean crit1 = false;
        boolean crit2 = false;
        boolean crit3 = false;
        for (int i = 0; i < this.getMaxIterations(); i++) {
            iterations.increment();

            LagrangeSolution sol1 = null;

            int qpLoop = 0;
            double sigma = getSettings().getSigmaMax() * 10.0;

            //log("SQP:hessian" + H);
            //LOOP TO FIND SOLUTION WITH SIGMA<SIGMA THRESHOLD

            while ((sigma > getSettings().getSigmaMax() || sigma < - Precision.EPSILON) &&
                   qpLoop < getSettings().getQpMaxLoop()) {
                sol1 = augmented ? solveAugmentedQP(y, rho) : solveQP();
                sigma = (sol1.getX().getDimension() == 0) ? getSettings().getSigmaMax() * 10.0 : sol1.getValue();

                if (sigma > getSettings().getSigmaMax() || sigma < - Precision.EPSILON) {
                    rho *= 10.0; //increment slack variable by factor 10.0
                    qpLoop += 1;
                    augmented = true;
                }

            }
            //IF SIGMA>SIGMA THRESHOLD ASSIGN DIRECTION FROM PENALTY GRADIENT

            if (qpLoop == getSettings().getQpMaxLoop()) {

                dx = (MatrixUtils.inverse(H).operate(penalty.gradX())).mapMultiply(-1.0);
                u = y.subtract(penalty.gradY());
                sigma = 0;
                augmented = true;

            } else {
                dx = sol1.getX();
                u = sol1.getLambda();
                sigma = sol1.getValue();
                penalty.updateRj(H, y, dx, u, sigma, iterations.getCount());
                //switch to normal QP if additional variable is small enough
                if (FastMath.abs(sigma) < getSettings().getEps()) {
                    augmented = false;
                }

            }

            penalty.update(J, JE, JI, x, y, dx, u);

            //if penalty gradient is >= 0 skip line search and try again with augmented QP
            if (penalty.getGradient() < 0) {

                rho = updateRho(dx, u, sigma);

                //LINE SEARCH
                alpha = lineSearch.search(penalty);

                //STORE OLD DERIVATIVE VALUE
                RealVector JOLD = J;
                RealMatrix JIOLD = JI;
                RealMatrix JEOLD = JE;

                //OLD LAGRANGIAN GRADIENT UPDATE WITH NEW MULTIPLIER
                if (m > 0) {
                    y = y.add(u.subtract(y).mapMultiply(alpha));
                }
                RealVector lagOld = lagrangianGradX(JOLD, JEOLD, JIOLD, x, u);
                //UPDATE ALL VARIABLE FOR THE NEXT STEP
                x = x.add(dx.mapMultiply(alpha));
                //PENALTY STORE ALL VARIABLE CALCULATED WITH THE LAST STEP
                //penalty function memorize last calculation done in the line search
                functionEval = penalty.getObjEval();
                eqEval = penalty.getEqEval();
                ineqEval = penalty.getIqEval();
                //calculate new Gradients
                computeGradients();
                //internalGradientCentered(x);
                //NEW LAGRANGIAN GRADIENT UPDATE WITH NEW MULTIPLIER AND NEW VALUES
                RealVector lagnew = lagrangianGradX(J, JE, JI, x, u);

                //CONVERGENCE CHECK
                //STATIONARITY
                double kkt = lagnew.getNorm();
                double step1 = dx.mapMultiply(alpha).dotProduct(H.operate(dx.mapMultiply(alpha)));
                double step2 = alpha * dx.getNorm();
                double violation = constraintViolation();
                crit0 = kkt <= FastMath.sqrt(getSettings().getEps());
                crit1 = step1 <= getSettings().getEps() * getSettings().getEps();
                crit2 = step2 <= getSettings().getEps() * (1.0 + x.getNorm());
                crit3 = violation <= FastMath.sqrt(getSettings().getEps());
                formatter.logRow(iterations.getCount(),
                                 alpha, lineSearch.getIteration(),
                                 step2, step1, kkt, violation, sigma,
                                 penalty.getPenaltyEval(), functionEval);

                if ((crit0 || (crit1 && crit2)) && crit3) {
                   break;
                }

                //HESSIAN UPDATE WITH THE LOGIC OF LINE SEARCH AND WITH THE INTERNAL LOGIC(DUMPING)
                if (lineSearch.isBadStepFailed()) {
                    //reset hessian and initialize for augmented QP solution with multiplier to zero
                    // bfgs.resetHessian(FastMath.min(gamma,10e6));
                    bfgs.resetHessian();

                    H         = bfgs.getHessian();
                    augmented = true;
                    rho = getSettings().getRhoCons(); //get initial value of slack variable
                    penalty.resetRj();
                    //y.set(0.0);
                    lineSearch.resetBadStepCount();

                } else if (lineSearch.isBadStepDetected()) {
                    //maintain the same Hessian

                    H = bfgs.getHessian();
                    augmented = true;// in case of bad step solve augmented QP

                } else {
                    // good step detected proceed with hessian update
                    bfgs.update(dx.mapMultiply(alpha), lagnew.subtract(lagOld));
                    H = bfgs.getHessian();
                }

            } else {
                augmented = true;
                rho *= 10.0; //increment slack variable by factor 10.0
            }
        }

        formatter.logRow(crit2, crit1, crit0, crit3);

        return new LagrangeSolution(x, y, functionEval);

    }

    /** Compute gradients.
     */
    private void computeGradients() {
        switch (getSettings().getGradientMode()) {
            case EXTERNAL :
                externalGradient();
                break;
            case FORWARD :
                forwardGradient();
                break;
            default :
                centralGradient();
                break;
        }
    }

    /** Compute constraints violations.
     * @return constraints violations
     */
    private double constraintViolation() {

        double crit = 0;
        if (this.getEqConstraint() != null) {
            crit = crit + eqEval.subtract(this.getEqConstraint().getLowerBound()).getL1Norm();
        }
        if (this.getIqConstraint() != null) {
            RealVector violated = ineqEval.subtract(this.getIqConstraint().getLowerBound());
            for (int k = 0; k < violated.getDimension(); k++) {
                violated.setEntry(k, FastMath.min(0.0, violated.getEntry(k)));
            }
            crit = crit + violated.getL1Norm();
        }
        return crit;
    }

    private double updateRho(final RealVector dx, final RealVector dy, final double additionalVariable) {
        int me = JE != null ? JE.getRowDimension() : 0;
        int mi = JI != null ? JI.getRowDimension() : 0;
        RealMatrix JAC;
        if (me + mi > 0) {
            JAC = new Array2DRowRealMatrix(me + mi, x.getDimension());
                if (JE != null) {
                    JAC.setSubMatrix(JE.getData(), 0, 0);
                }
                if (JI != null) {
                    JAC.setSubMatrix(JI.getData(), me, 0);
                }

            double num = 10.0 * FastMath.pow(dx.dotProduct(JAC.preMultiply(dy)), 2);
            double den = (1.0 - additionalVariable) * (1.0 - additionalVariable) * dx.dotProduct(H.operate(dx));
            //double den = (1.0 - additionalVariable) * dx.dotProduct(H.operate(dx));

            return FastMath.max(10.0, num / den);

        }
        return 0;
    }

    /** Solve augmented problem.
     * @param y Lagrange multipliers
     * @param rho rho
     * @return problem solution
     */
    private LagrangeSolution solveAugmentedQP(final RealVector y, final double rho) {

        RealVector g = J;

        int me = 0;
        int mi = 0;
        int add = 0;
        boolean violated = false;
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        if (getIqConstraint() != null) {

            mi = getIqConstraint().dimY();
            violated = ineqEval.subtract(getIqConstraint().getLowerBound()).getMinValue() <= getSettings().getEps() ||
                       y.getMaxValue() > 0;

        }
        // violated = true;
        if (me > 0 || violated) {
            add = 1;
        }

        RealMatrix H1 = new Array2DRowRealMatrix(H.getRowDimension() + add, H.getRowDimension() + add);
        H1.setSubMatrix(H.getData(), 0, 0);
        if (add == 1) {
            H1.setEntry(H.getRowDimension(), H.getRowDimension(), rho);
        }

        RealVector g1 = new ArrayRealVector(g.getDimension() + add);
        g1.setSubVector(0, g);

        LinearEqualityConstraint eqc = null;
        RealVector conditioneq;
        if (getEqConstraint() != null) {
            RealMatrix eqJacob = JE;
            RealMatrix Ae = new Array2DRowRealMatrix(me, x.getDimension() + add);
            RealVector be = new ArrayRealVector(me);
            Ae.setSubMatrix(eqJacob.getData(), 0, 0);
            conditioneq = this.eqEval.subtract(getEqConstraint().getLowerBound());
            Ae.setColumnVector(x.getDimension(), conditioneq.mapMultiply(-1.0));

            be.setSubVector(0, getEqConstraint().getLowerBound().subtract(this.eqEval));
            eqc = new LinearEqualityConstraint(Ae, be);

        }
        LinearInequalityConstraint iqc = null;

        if (getIqConstraint() != null) {

            RealMatrix iqJacob = JI;
            RealMatrix Ai = new Array2DRowRealMatrix(mi, x.getDimension() + add);
            RealVector bi = new ArrayRealVector(mi);
            Ai.setSubMatrix(iqJacob.getData(), 0, 0);

            RealVector conditioniq = this.ineqEval.subtract(getIqConstraint().getLowerBound());

            if (add == 1) {

                for (int i = 0; i < conditioniq.getDimension(); i++) {
                    if (!(conditioniq.getEntry(i) <= getSettings().getEps() || y.getEntry(me + i) > 0)) {
                        conditioniq.setEntry(i, 0);
                    }
                }

                Ai.setColumnVector(x.getDimension(), conditioniq.mapMultiply(-1.0));

            }
            bi.setSubVector(0, getIqConstraint().getLowerBound().subtract(this.ineqEval));
            iqc = new LinearInequalityConstraint(Ai, bi);

        }
        int box = 0;
        if (getBoxConstraint() != null) {
            box = getBoxConstraint().dimY();
        }
        //this.log("MI:" + box);
        LinearBoundedConstraint bc = null;

        if (add == 1) {

            RealMatrix sigmaA = new Array2DRowRealMatrix(1 + box, x.getDimension() + 1);
            sigmaA.setEntry(0, x.getDimension(), 1.0);

            ArrayRealVector lb = new ArrayRealVector(1 + box, 0.0);
            ArrayRealVector ub = new ArrayRealVector(1 + box, 1.0);
            bc = new LinearBoundedConstraint(sigmaA, lb, ub);

        }

        QuadraticFunction q = new QuadraticFunction(H1, g1, 0);

        LagrangeSolution sol = this.getQPSolver().optimize(new ObjectiveFunction(q), iqc, eqc, bc);

        // Solve the QP problem
        if (sol.getX().getDimension() == 0) {
            return sol;
        }
        double sigma;
        if (add == 1) {
            sigma = sol.getX().getEntry(x.getDimension());
        }
        else {
            sigma = 0;
        }

        return (me + mi == 0) ?
               new LagrangeSolution(sol.getX().getSubVector(0, x.getDimension()), null, sigma) :
               new LagrangeSolution(sol.getX().getSubVector(0, x.getDimension()), sol.getLambda().getSubVector(0, me + mi), sigma);

    }

    /**
     * Solves the Quadratic Programming (QP) subproblem in the current SQP iteration.
     * @return a {@link LagrangeSolution} representing the QP solution, or {@code null} if the QP failed
     */
    private LagrangeSolution solveQP() {

        final QuadraticFunction q  = new QuadraticFunction(this.H, this.J, 0);
        int                     n  = x.getDimension();
        int                     me = 0;
        int                     mi = 0;

        // Equality constraints
        LinearEqualityConstraint eqc = null;
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
            RealMatrix Ae = new Array2DRowRealMatrix(me, n);
            RealVector be = getEqConstraint().getLowerBound().subtract(eqEval);
            Ae.setSubMatrix(JE.getData(), 0, 0);
            eqc = new LinearEqualityConstraint(Ae, be);
        }

        // Inequality constraints
        LinearInequalityConstraint iqc = null;
        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
            RealMatrix Ai = new Array2DRowRealMatrix(mi, n);
            RealVector bi = getIqConstraint().getLowerBound().subtract(ineqEval);
            Ai.setSubMatrix(JI.getData(), 0, 0);
            iqc = new LinearInequalityConstraint(Ai, bi);
        }

        // Solve the QP problem
        LagrangeSolution sol = getQPSolver().optimize(new ObjectiveFunction(q), iqc, eqc);
        if (sol.getX().getDimension() == 0) {
            return sol;
        }

        // Extract primal and dual components
        RealVector solutionX      = sol.getX().getSubVector(0, n);
        RealVector solutionLambda = (me + mi > 0) ? sol.getLambda().getSubVector(0, me + mi) : new ArrayRealVector(0, 0);

        return new LagrangeSolution(solutionX, solutionLambda, 0.0);
    }

    /** Computes the gradient of the Lagrangian function with respect to the primal variable {@code x}.
     * <p>
     * The Lagrangian is defined as:
     * </p>
     * <pre>
     *     L(x, y) = f(x) - yₑᵗ·cₑ(x) - yᵢᵗ·cᵢ(x)
     * </pre>
     * <p>
     * where:
     * </p>
     * <ul>
     *   <li>{@code f(x)} is the objective function</li>
     *   <li>{@code cₑ(x)} are the equality constraints</li>
     *   <li>{@code cᵢ(x)} are the inequality constraints</li>
     *   <li>{@code y = [yₑ; yᵢ]} is the vector of Lagrange multipliers</li>
     * </ul>
     * <p>
     * The gradient with respect to {@code x} is given by:
     * </p>
     * <pre>
     *     ∇ₓ L(x, y) = ∇f(x) - JEᵗ·yₑ - JIᵗ·yᵢ
     * </pre>
     * @param otherJ  the gradient of the objective function {@code ∇f(x)}, length {@code n}
     * @param otherJE the Jacobian of the equality constraints, shape {@code [me x n]} (nullable)
     * @param otherJI the Jacobian of the inequality constraints, shape {@code [mi x n]} (nullable)
     * @param ignoredX  the current point in the primal space (not used directly, included for API symmetry)
     * @param y  the stacked Lagrange multipliers {@code [yₑ; yᵢ]}, length {@code me + mi}
     * @return the gradient of the Lagrangian with respect to {@code x}, length {@code n}
     */
    public RealVector lagrangianGradX(final RealVector otherJ, final RealMatrix otherJE, final RealMatrix otherJI,
                                      final RealVector ignoredX, final RealVector y) {

        RealVector gradL  = new ArrayRealVector(otherJ);
        int        offset = 0;

        // Subtract JEᵗ · yₑ if equality constraints exist
        if (otherJE != null) {
            int        me     = otherJE.getRowDimension();
            RealVector yEq    = y.getSubVector(0, me);
            RealVector termEq = otherJE.preMultiply(yEq);
            gradL = gradL.subtract(termEq);
            offset += me;
        }

        // Subtract JIᵗ · yᵢ if inequality constraints exist
        if (otherJI != null) {
            int        mi     = otherJI.getRowDimension();
            RealVector yIq    = y.getSubVector(offset, mi);
            RealVector termIq = otherJI.preMultiply(yIq);
            gradL = gradL.subtract(termIq);
        }

        return gradL;
    }

    /** Compute gradient directly.
     */
    private void externalGradient()
    {
        J = this.getObj().gradient(x);
        if (this.getEqConstraint() != null) {
            JE = this.getEqConstraint().jacobian(x);
        }
        if (this.getIqConstraint() != null) {
            JI = this.getIqConstraint().jacobian(x);
        }

    }

    /** Computes the gradient of the objective function and the Jacobians of the constraints
     * using forward finite differences (first-order accurate).
     * <p>
     * Each variable is perturbed independently by a small step size proportional to
     * the square root of machine precision, and partial derivatives are approximated
     * using forward differencing.
     * </p>
     */
    private void forwardGradient() {

        int    n       = x.getDimension();
        double sqrtEps = FastMath.sqrt(Precision.EPSILON);

        double     fRef  = this.functionEval;
        RealVector eqRef = this.eqEval;
        RealVector iqRef = this.ineqEval;

        RealVector gradF  = new ArrayRealVector(n);
        RealMatrix gradEq = (getEqConstraint() != null) ? new Array2DRowRealMatrix(eqRef.getDimension(), n) : null;
        RealMatrix gradIq = (getIqConstraint() != null) ? new Array2DRowRealMatrix(iqRef.getDimension(), n) : null;

        for (int i = 0; i < n; i++) {
            double xi = x.getEntry(i);
            double h  = sqrtEps * FastMath.max(1.0, FastMath.abs(xi));

            RealVector xPerturbed = new ArrayRealVector(x);
            xPerturbed.setEntry(i, xi + h);

            double fPerturbed = getObj().value(xPerturbed);
            gradF.setEntry(i, (fPerturbed - fRef) / h);

            if (gradEq != null) {
                RealVector eqPerturbed = getEqConstraint().value(xPerturbed);
                RealVector diffEq      = eqPerturbed.subtract(eqRef).mapMultiply(1.0 / h);
                gradEq.setColumnVector(i, diffEq);
            }

            if (gradIq != null) {
                RealVector iqPerturbed = getIqConstraint().value(xPerturbed);
                RealVector diffIq      = iqPerturbed.subtract(iqRef).mapMultiply(1.0 / h);
                gradIq.setColumnVector(i, diffIq);
            }
        }

        this.J  = gradF;
        this.JE = gradEq;
        this.JI = gradIq;
    }

    /**
     * Computes the gradient of the objective function and the Jacobians of the constraints
     * using centered finite differences (second-order accurate).
     */
    private void centralGradient() {

        int    n     = x.getDimension();
        double hBase = FastMath.cbrt(Precision.EPSILON);

        double fPlus;
        double fMinus;
        RealVector gradF  = new ArrayRealVector(n);
        RealMatrix gradEq = (getEqConstraint() != null) ? new Array2DRowRealMatrix(eqEval.getDimension(), n) : null;
        RealMatrix gradIq = (getIqConstraint() != null) ? new Array2DRowRealMatrix(ineqEval.getDimension(), n) : null;

        for (int i = 0; i < n; i++) {
            double xi = x.getEntry(i);
            double h  = hBase * FastMath.max(1.0, FastMath.abs(xi));

            RealVector xPlus  = new ArrayRealVector(x);
            RealVector xMinus = new ArrayRealVector(x);
            xPlus.addToEntry(i, h);
            xMinus.addToEntry(i, -h);

            fPlus  = getObj().value(xPlus);
            fMinus = getObj().value(xMinus);
            gradF.setEntry(i, (fPlus - fMinus) / (2.0 * h));

            if (gradEq != null) {
                RealVector eqPlus  = getEqConstraint().value(xPlus);
                RealVector eqMinus = getEqConstraint().value(xMinus);
                RealVector dEq     = eqPlus.subtract(eqMinus).mapDivide(2.0 * h);
                gradEq.setColumnVector(i, dEq);
            }

            if (gradIq != null) {
                RealVector iqPlus  = getIqConstraint().value(xPlus);
                RealVector iqMinus = getIqConstraint().value(xMinus);
                RealVector dIq     = iqPlus.subtract(iqMinus).mapDivide(2.0 * h);
                gradIq.setColumnVector(i, dIq);
            }
        }

        this.J  = gradF;
        this.JE = gradEq;
        this.JI = gradIq;
    }

    /** Set debug printer.
     * @param printer debug printer
     */
    public void setDebugPrinter(final DebugPrinter printer) {
        formatter.setDebugPrinter(printer);
    }

}
