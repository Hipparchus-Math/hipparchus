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

import java.util.Arrays;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
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
 * Supports: equality constraints, inequality constraints, penalty function
 * updates, line search strategies, BFGS Hessian update, and augmented QP
 * formulation using a relaxation variable.</p>
 *
 * @since 3.1
 */
public class SQPOptimizerS2 extends AbstractSQPOptimizer2 {

    private enum QPMode {
        /**
         * Solving the standard QP subproblem (no slack variable).
         */
        QP_STARDARD,
        /**
         * Solving the augmented QP subproblem (with slack variable).
         */
        QP_AUGMENTED
    }

    /**
     * Logger.
     */
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
     * Value of the bounds .
     */
    private RealVector bEval;

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

    /**
     * Jacobian of the bounds .
     */
    private RealMatrix JB;

    /**
     * Evaluation of the objective function.
     */
    private double functionEval;

    /**
     * old objective function.
     */
    private double functionEvalOld;

    /**
     * Current point.
     */
    private RealVector x;

    /**
     * Bounds in form o inequality constraints.
     */
    private LinearInequalityConstraint bounds;

    /**
     * Lower bound.
     */
    private ArrayRealVector LB;

    /**
     * Upper bound.
     */
    private ArrayRealVector UB;

    /**
     * {@inheritDoc}
     */
    @Override
    public LagrangeSolution doOptimize() {
        formatter.setEps(getSettings().getEps());
        formatter.logHeader();
        int me = 0;
        int mi = 0;
        int mb = 0;

        //EQUALITY CONSTRAINT
        if (this.getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        //INEQUALITY CONSTRAINT
        if (this.getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
        }

        mb = buildBoundsAsInequalities();

        final int m = me + mi + mb;

        double rho = getSettings().getRhoCons(); //get initial value of slack variable from SQPOption
        x = (this.getStartPoint() != null) ? new ArrayRealVector(this.getStartPoint()) : initGuess();

        RealVector y = new ArrayRealVector(me + mi + mb, 0.0);

        //all the function and constraint evaluation will be performed inside the penalty function
        MeritFunctionL2 penalty = new MeritFunctionL2(this.getObj(), this.getEqConstraint(), this.getIqConstraint(), this.bounds, x);

        LineSearch lineSearch = new LineSearch(getSettings().getEps(), 5, getSettings().getMu(), getSettings().getB(),
                getSettings().getMaxLineSearchIteration(), 2);

        //INITIAL VALUES
        functionEval = penalty.getObjEval();
        if (this.getEqConstraint() != null) {
            eqEval = penalty.getEqEval();
        }
        if (this.getIqConstraint() != null) {
            ineqEval = penalty.getIqEval();
        }
        if (this.bounds != null) {
            bEval = penalty.getBEval();
        }
        double EPS=this.getSettings().getEps();
        double EPS2=EPS*EPS;
        double sqrtEPS=FastMath.sqrt(EPS);
        computeGradients();
        double gamma = 1.0;
        H = MatrixUtils.createRealIdentityMatrix(x.getDimension()).scalarMultiply(FastMath.sqrt(gamma));

        BFGSUpdater bfgs = new BFGSUpdater(H, EPS, true, getMatrixDecompositionTolerance().getEpsMatrixDecomposition());

        RealVector dx = new ArrayRealVector(x.getDimension());
        RealVector u = new ArrayRealVector(y.getDimension());
        penalty.update(J, JE, JI, x, y, dx, u);
        QPMode QPMODE = QPMode.QP_AUGMENTED;

        boolean crit0 = false;
        boolean crit1 = false;
        boolean crit2 = false;
        boolean crit3 = false;
        boolean crit4 = false;
        double sigma = 0.0;
        double alpha = 0.0;
        RealVector lagOld = null;
        RealVector lagNew = null;

        LagrangeSolution qpSolution = null;

        boolean GRADFAIL = false;
        boolean FALLBACK = false;
        boolean RECOVERYMODE=false;
        double VIOLATION = constraintViolation();
        double KKT = Double.POSITIVE_INFINITY;
        double DHD = 0.0;
        double XNORM = 0.0;
        double FUNDIFF=0.0;
        int BFGSUPDATE=0;
        for (int i = 0; i < this.getMaxIterations(); i++) {

            sigma = getSettings().getSigmaMax() * 10.0;
            while ((sigma > getSettings().getSigmaMax() || sigma < 0.0) && rho < 1.0e9 && !FALLBACK) {

                qpSolution = (QPMODE == QPMode.QP_AUGMENTED) ? solveAugmentedQP(y, rho) : solveQP();
                sigma = (qpSolution == null || qpSolution.getX().getDimension() == 0) ? getSettings().getSigmaMax() * 10.0 : qpSolution.getValue();
                if ((sigma > getSettings().getSigmaMax() || sigma < 0.0)) {
                    if (QPMODE == QPMode.QP_AUGMENTED) rho = rhoUp(rho);
                    QPMODE = QPMode.QP_AUGMENTED;
                }

            }
            //IF SIGMA>SIGMA THRESHOLD AFTER SEVERAL ATTEMPT ASSIGN DIRECTION FROM PENALTY GRADIENT
            if (rho >= 1.0e9 || FALLBACK) {

//




                qpSolution = solveQPFallBack(penalty.gradX());

                if (qpSolution == null || qpSolution.getX().getDimension() == 0 )
                {
                    break;//infesible
                }

                dx = qpSolution.getX();
                projectDirectionInPlace(x, dx, true);
                //estimation of multiplier from penalty grad y
                if (m > 0) {
                    u = y.subtract(penalty.gradY());
                }
                //estimation of bounds multiplier from QP
                if (mb > 0) {
                    RealVector db = qpSolution.getLambda();
                    u.setSubVector(mi + me, db);
                }
                QPMODE = QPMode.QP_AUGMENTED;
                sigma = 0.0;
                rho = getSettings().getRhoCons();;

//

            } else {
                dx = qpSolution.getX();
                u = qpSolution.getLambda();

                sigma = qpSolution.getValue();
                projectDirectionInPlace(x, dx, true);
                penalty.updateRj(H, y, dx, u, sigma, iterations.getCount());
                //switch to normal QP if additional variable is small enough
                if (QPMODE == QPMode.QP_AUGMENTED && FastMath.abs(sigma) < getSettings().getEps() && !GRADFAIL) {

                    QPMODE = QPMode.QP_STARDARD;
                    rho = getSettings().getRhoCons();

                }

            }

            penalty.update(J, JE, JI, x, y, dx, u);
            double rmax = 0.0;
            while (penalty.getGradient() >= 0.0 && rmax < penalty.getRmax()) {
                penalty.rUp();
                rmax = penalty.getR().getMaxValue();
            }

             if(penalty.getGradient()>=0 )
            {
                if (QPMODE == QPMode.QP_AUGMENTED)  rho = rhoUp(rho);
                QPMODE = QPMode.QP_AUGMENTED;
                GRADFAIL = true;

            }
            //if penalty gradient is > 0 skip line search and try again with augmented QP
             else  {
                iterations.increment();
                GRADFAIL = false;
                FALLBACK = false;

                if(QPMODE == QPMode.QP_AUGMENTED) rho = updateRho(dx, u, sigma);
                lagOld = lagrangianGradX(J, JE, JI, x, u);

                //LINE SEARCH
                alpha = lineSearch.search(penalty);




                if (m > 0) {
                    y = y.add(u.subtract(y).mapMultiply(alpha));
                }

                x = x.add(dx.mapMultiply(alpha));


                //penalty function memorize last calculation done in the line search
                functionEvalOld = functionEval;
                functionEval = penalty.getObjEval();
                eqEval = penalty.getEqEval();
                ineqEval = penalty.getIqEval();
                bEval = penalty.getBEval();

                computeGradients();
                RealVector lagnew = lagrangianGradX(J, JE, JI, x, u);

//              CONVERGENCE CHECK
                KKT = lagnew.getNorm();
                VIOLATION = constraintViolation();
                //STEP LENGHT
                DHD = dx.mapMultiply(alpha).dotProduct(H.operate(dx.mapMultiply(alpha)));
                XNORM = alpha * dx.getNorm();
                FUNDIFF=FastMath.abs(functionEval - functionEvalOld);              
                crit0 = KKT <= sqrtEPS;
                crit1 = DHD <= EPS2;
                crit2 = XNORM <= EPS * (1.0 + x.getNorm());
                crit3 = VIOLATION <= sqrtEPS;
                crit4 = FUNDIFF < EPS* (1.0 + FastMath.abs(functionEval));
                
                

                if ((crit0 || (crit1 && crit2)) && crit3) {
                    break;
                }

                //HESSIAN UPDATE WITH THE LOGIC OF LINE SEARCH AND WITH THE INTERNAL LOGIC(DUMPING)
                if (lineSearch.isBadStepFailed()) {
                    FALLBACK = false;
                    //reset hessian and initialize for augmented QP
                    bfgs.resetHessian();
                    H = bfgs.getHessian();
                    QPMODE = QPMode.QP_AUGMENTED;
                    penalty.resetRj();
                    if (m > 0) {
                        y = u.copy();

                    }
                    penalty.update(J, JE, JI, x, y, dx, u);
                    lineSearch.resetBadStepCount();
                    BFGSUPDATE=-2;

                } else if (lineSearch.isBadStepDetected()) {
                    //maintain the same Hessian
                    BFGSUPDATE=-1;
                    H = bfgs.getHessian();
                } else {
                    BFGSUPDATE = bfgs.update(dx.mapMultiply(alpha), lagnew.subtract(lagOld));
                    //
                    if(BFGSUPDATE==1 && XNORM==0)
                    {
                        FALLBACK = true;
                        bfgs.resetHessian();
                        QPMODE = QPMode.QP_AUGMENTED;
                        penalty.resetRj();
                        if (m > 0) {
                            y = u.copy();

                            }
                        penalty.update(J, JE, JI, x, y, dx, u);
                        
                    }
                    H = bfgs.getHessian();
                    
                    

                }
                
                formatter.logRow(iterations.getCount(),
                        alpha, lineSearch.getIteration(),
                        XNORM, DHD, KKT, VIOLATION, sigma,
                        penalty.getPenaltyEval(), functionEval,FUNDIFF,BFGSUPDATE);

            }
        }
        formatter.logRow(iterations.getCount(),
                        alpha, lineSearch.getIteration(),
                        XNORM, DHD, KKT, VIOLATION, sigma,
                        penalty.getPenaltyEval(), functionEval,FUNDIFF,-99);
        formatter.logRow(crit2, crit1, crit0, crit3,crit4);

        return new LagrangeSolution(x, y, functionEval);

    }

    /**
     * Compute gradients.
     */
    private void computeGradients() {
        switch (getSettings().getGradientMode()) {
            case EXTERNAL:
                externalGradient();
                break;
            case FORWARD:
                forwardGradient();
                break;
            default:
                centralGradient();
                break;
        }
    }

    /**
     * Compute constraints violations.
     *
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

    private double rhoUp(double rho) {
        return rho * 10.0;
    }

    private double updateRho(final RealVector dx, final RealVector dy, final double additionalVariable) {
        int me = JE != null ? JE.getRowDimension() : 0;
        int mi = JI != null ? JI.getRowDimension() : 0;
        int mb = JB != null ? JB.getRowDimension() : 0;
        RealMatrix JAC;
        if (me + mi + mb > 0) {
            JAC = new Array2DRowRealMatrix(me + mi + mb, x.getDimension());
            if (JE != null) {
                JAC.setSubMatrix(JE.getData(), 0, 0);
            }
            if (JI != null) {
                JAC.setSubMatrix(JI.getData(), me, 0);
            }

            if (JB != null) {
                JAC.setSubMatrix(JB.getData(), me + mi, 0);
            }

            double num = 10.0 * FastMath.pow(dx.dotProduct(JAC.preMultiply(dy)), 2);
            double den = (1.0 - additionalVariable) * (1.0 - additionalVariable) * dx.dotProduct(H.operate(dx));
            if (den < Precision.SAFE_MIN) {
                den = Precision.SAFE_MIN;
            }

            return FastMath.max(10.0, num / den);

        }
        return 1.0;

    }

    /**
     * Solve augmented problem.
     *
     * @param y Lagrange multipliers
     * @param rho rho
     * @return problem solution
     */
    private LagrangeSolution solveAugmentedQP(final RealVector y, final double rho) {

        RealVector g = J;

        int me = 0;
        int mi = 0;
        int mb = 0;
        int mc = 0;
        int add = 0;
        boolean violated = false;
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        if (getIqConstraint() != null) {

            mi = getIqConstraint().dimY();
            violated = ineqEval.subtract(getIqConstraint().getLowerBound()).getMinValue() <= getSettings().getEps()
                    || y.getSubVector(me, mi).getMaxValue() > 0;

        }

        if (bounds != null) {
            mb = bounds.dimY();
        }
        mc = mi + mb;
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
        RealMatrix A = null;
        RealVector B = null;
        if (mc > 0) {
            A = new Array2DRowRealMatrix(mc, x.getDimension() + add);
            B = new ArrayRealVector(mc);
        }
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
            A.setSubMatrix(Ai.getData(), 0, 0);

            B.setSubVector(0, bi);

        }
        if (bounds != null) {
            A.setSubMatrix(JB.getData(), mi, 0);
            B.setSubVector(mi, bounds.getLowerBound().subtract(this.bEval));
        }
        if (mc > 0) {
            iqc = new LinearInequalityConstraint(A, B);
        }

        LinearBoundedConstraint bc = null;

        if (add == 1) {

            RealMatrix sigmaA = new Array2DRowRealMatrix(1, x.getDimension() + 1);
            sigmaA.setEntry(0, x.getDimension(), 1.0);

            ArrayRealVector lb = new ArrayRealVector(1, 0.0);
            ArrayRealVector ub = new ArrayRealVector(1, 1.0);
            bc = new LinearBoundedConstraint(sigmaA, lb, ub);

        }

        QuadraticFunction q = new QuadraticFunction(H1, g1, 0);

        LagrangeSolution sol = this.getQPSolver().optimize(new ObjectiveFunction(q), iqc, eqc, bc);

        // Solve the QP problem
        if (sol == null || sol.getX().getDimension() == 0) {
            return sol;
        }
        double sigma;
        if (add == 1) {
            sigma = sol.getX().getEntry(x.getDimension());
        } else {
            sigma = 0;
        }
       
        return (me + mi + mb == 0)
                ? new LagrangeSolution(sol.getX().getSubVector(0, x.getDimension()), null, sigma)
                : new LagrangeSolution(sol.getX().getSubVector(0, x.getDimension()), sol.getLambda().getSubVector(0, me + mi + mb), sigma);

    }

    /**
     * Solves the Quadratic Programming (QP) subproblem in the current SQP
     * iteration.
     *
     * @return a {@link LagrangeSolution} representing the QP solution, or
     * {@code null} if the QP failed
     */
    private LagrangeSolution solveQP() {

        final QuadraticFunction q = new QuadraticFunction(this.H, this.J, 0);
        int n = x.getDimension();
        int me = 0;
        int mi = 0;
        int mb = 0;

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
        RealMatrix A = null;
        RealVector B = null;
        int mc = 0;
        if (getIqConstraint() != null) {
            mc += getIqConstraint().dimY();
        }
        if (bounds != null) {
            mc += bounds.dimY();
        }
        if (mc > 0) {
            A = new Array2DRowRealMatrix(mc, n);
            B = new ArrayRealVector(mc);
        }
        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();

            RealVector bi = getIqConstraint().getLowerBound().subtract(ineqEval);
            A.setSubMatrix(JI.getData(), 0, 0);

            B.setSubVector(0, bi);
        }

        if (bounds != null) {
            mb = bounds.dimY();

            RealVector bb = bounds.getLowerBound().subtract(bEval);
            A.setSubMatrix(JB.getData(), mi, 0);
            B.setSubVector(mi, bb);
            // iqc = new LinearInequalityConstraint(Ab, bb);
        }
        if (mc > 0) {
            iqc = new LinearInequalityConstraint(A, B);
        }

        // Solve the QP problem
        LagrangeSolution sol = getQPSolver().optimize(new ObjectiveFunction(q), iqc, eqc);
        // Solve the QP problem
        if (sol == null || sol.getX().getDimension() == 0) {
            return sol;
        }

        // Extract primal and dual components
        RealVector solutionX = sol.getX();
        
        RealVector solutionLambda = (me + mi + mb > 0) ? sol.getLambda() : new ArrayRealVector(0, 0);

        return new LagrangeSolution(solutionX, solutionLambda, 0.0);
    }

    /**
     * Solves the Quadratic Programming (QP) subproblem in the current SQP
     * iteration if penalty gradient is >=0
     *
     * @return a {@link LagrangeSolution} representing the QP solution, or
     * {@code null} if the QP failed
     */
    private LagrangeSolution solveQPFallBack(RealVector penaltyGradx) {

        final QuadraticFunction q = new QuadraticFunction(this.H, penaltyGradx, 0);
        int n = x.getDimension();

        int mb = 0;

        // Inequality constraints
        LinearInequalityConstraint iqc = null;

        if (bounds != null) {
            mb = bounds.dimY();
            RealMatrix A = new Array2DRowRealMatrix(JB.getData());
            RealVector B = bounds.getLowerBound().subtract(bEval);
            iqc = new LinearInequalityConstraint(A, B);
        }

        // Solve the QP problem
        LagrangeSolution sol = getQPSolver().optimize(new ObjectiveFunction(q), iqc);
        // Solve the QP problem
        if (sol == null || sol.getX().getDimension() == 0) {
            return sol;
        }

        // Extract primal and dual components
        RealVector solutionX = sol.getX();
        RealVector solutionLambda = (mb > 0) ? sol.getLambda() : new ArrayRealVector(0, 0);

        return new LagrangeSolution(solutionX, solutionLambda, 0.0);
    }

    /**
     * Computes the gradient of the Lagrangian function with respect to the
     * primal variable {@code x}.
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
     * <li>{@code f(x)} is the objective function</li>
     * <li>{@code cₑ(x)} are the equality constraints</li>
     * <li>{@code cᵢ(x)} are the inequality constraints</li>
     * <li>{@code y = [yₑ; yᵢ]} is the vector of Lagrange multipliers</li>
     * </ul>
     * <p>
     * The gradient with respect to {@code x} is given by:
     * </p>
     * <pre>
     *     ∇ₓ L(x, y) = ∇f(x) - JEᵗ·yₑ - JIᵗ·yᵢ
     * </pre>
     *
     * @param otherJ the gradient of the objective function {@code ∇f(x)},
     * length {@code n}
     * @param otherJE the Jacobian of the equality constraints, shape
     * {@code [me x n]} (nullable)
     * @param otherJI the Jacobian of the inequality constraints, shape
     * {@code [mi x n]} (nullable)
     * @param ignoredX the current point in the primal space (not used directly,
     * included for API symmetry)
     * @param y the stacked Lagrange multipliers {@code [yₑ; yᵢ]}, length
     * {@code me + mi}
     * @return the gradient of the Lagrangian with respect to {@code x}, length
     * {@code n}
     */
    public RealVector lagrangianGradX(final RealVector otherJ, final RealMatrix otherJE, final RealMatrix otherJI,
            final RealVector ignoredX, final RealVector y) {

        RealVector gradL = new ArrayRealVector(otherJ);
        int offset = 0;

        // Subtract JEᵗ · yₑ if equality constraints exist
        if (otherJE != null) {
            int me = otherJE.getRowDimension();
            RealVector yEq = y.getSubVector(0, me);
            RealVector termEq = otherJE.preMultiply(yEq);
            gradL = gradL.subtract(termEq);
            offset += me;
        }

        // Subtract JIᵗ · yᵢ if inequality constraints exist
        if (otherJI != null) {
            int mi = otherJI.getRowDimension();
            RealVector yIq = y.getSubVector(offset, mi);
            RealVector termIq = otherJI.preMultiply(yIq);
            gradL = gradL.subtract(termIq);
            offset += mi;
        }

        // Subtract JBᵗ · yᵢ if bounds exist
        if (JB != null) {
            int mb = JB.getRowDimension();
            RealVector yIq = y.getSubVector(offset, mb);
            RealVector termIq = JB.preMultiply(yIq);
            gradL = gradL.subtract(termIq);
        }

        return gradL;
    }

    /**
     * Compute gradient directly.
     */
    private void externalGradient() {
        J = this.getObj().gradient(x);
        if (this.getEqConstraint() != null) {
            JE = this.getEqConstraint().jacobian(x);
        }
        if (this.getIqConstraint() != null) {
            JI = this.getIqConstraint().jacobian(x);
        }

    }

    /**
     * Computes the gradient of the objective function and the Jacobians of the
     * constraints using forward finite differences (first-order accurate).
     * <p>
     * Each variable is perturbed independently by a small step size
     * proportional to the square root of machine precision, and partial
     * derivatives are approximated using forward differencing.
     * </p>
     */
    private void forwardGradient() {

        int n = x.getDimension();
        double sqrtEps = FastMath.sqrt(Precision.EPSILON);

        double fRef = this.functionEval;
        RealVector eqRef = this.eqEval;
        RealVector iqRef = this.ineqEval;

        RealVector gradF = new ArrayRealVector(n);
        RealMatrix gradEq = (getEqConstraint() != null) ? new Array2DRowRealMatrix(eqRef.getDimension(), n) : null;
        RealMatrix gradIq = (getIqConstraint() != null) ? new Array2DRowRealMatrix(iqRef.getDimension(), n) : null;

        for (int i = 0; i < n; i++) {
            double xi = x.getEntry(i);
            double h = sqrtEps * FastMath.max(1.0, FastMath.abs(xi));

            RealVector xPerturbed = new ArrayRealVector(x);
            xPerturbed.setEntry(i, xi + h);

            double fPerturbed = getObj().value(xPerturbed);
            gradF.setEntry(i, (fPerturbed - fRef) / h);

            if (gradEq != null) {
                RealVector eqPerturbed = getEqConstraint().value(xPerturbed);
                RealVector diffEq = eqPerturbed.subtract(eqRef).mapMultiply(1.0 / h);
                gradEq.setColumnVector(i, diffEq);
            }

            if (gradIq != null) {
                RealVector iqPerturbed = getIqConstraint().value(xPerturbed);
                RealVector diffIq = iqPerturbed.subtract(iqRef).mapMultiply(1.0 / h);
                gradIq.setColumnVector(i, diffIq);
            }
        }

        this.J = gradF;
        this.JE = gradEq;
        this.JI = gradIq;
    }

    /**
     * Computes the gradient of the objective function and the Jacobians of the
     * constraints using centered finite differences (second-order accurate).
     */
    private void centralGradient() {

        int n = x.getDimension();
        double hBase = FastMath.cbrt(Precision.EPSILON);

        double fPlus;
        double fMinus;
        RealVector gradF = new ArrayRealVector(n);
        RealMatrix gradEq = (getEqConstraint() != null) ? new Array2DRowRealMatrix(eqEval.getDimension(), n) : null;
        RealMatrix gradIq = (getIqConstraint() != null) ? new Array2DRowRealMatrix(ineqEval.getDimension(), n) : null;

        for (int i = 0; i < n; i++) {
            double xi = x.getEntry(i);
            double h = hBase * FastMath.max(1.0, FastMath.abs(xi));

            RealVector xPlus = new ArrayRealVector(x);
            RealVector xMinus = new ArrayRealVector(x);
            xPlus.addToEntry(i, h);
            xMinus.addToEntry(i, -h);

            fPlus = getObj().value(xPlus);
            fMinus = getObj().value(xMinus);
            gradF.setEntry(i, (fPlus - fMinus) / (2.0 * h));

            if (gradEq != null) {
                RealVector eqPlus = getEqConstraint().value(xPlus);
                RealVector eqMinus = getEqConstraint().value(xMinus);
                RealVector dEq = eqPlus.subtract(eqMinus).mapDivide(2.0 * h);
                gradEq.setColumnVector(i, dEq);
            }

            if (gradIq != null) {
                RealVector iqPlus = getIqConstraint().value(xPlus);
                RealVector iqMinus = getIqConstraint().value(xMinus);
                RealVector dIq = iqPlus.subtract(iqMinus).mapDivide(2.0 * h);
                gradIq.setColumnVector(i, dIq);
            }
        }

        this.J = gradF;
        this.JE = gradEq;
        this.JI = gradIq;
    }

    /**
     * Set debug printer.
     *
     * @param printer debug printer
     */
    public void setDebugPrinter(final DebugPrinter printer) {
        formatter.setDebugPrinter(printer);
    }

    /**
     * Build Bound bounds constraint in the form of inequaliy: x_i ≥ lb_i e -x_i
     * ≥ -ub_i
     *
     * @return m, dimension of bound constraints
     */
    private int buildBoundsAsInequalities() {
        if (getSimpleBounds() == null) {
            int n = getObj().dim();
            this.bounds = null;
            LB = new ArrayRealVector(n, Double.NEGATIVE_INFINITY);
            UB = new ArrayRealVector(n, Double.POSITIVE_INFINITY);
            this.JB = null;
            return 0;
        }
        int n = getObj().dim();
        LB = (getSimpleBounds().getLower() != null) ? new ArrayRealVector(getSimpleBounds().getLower()) : new ArrayRealVector(n, Double.NEGATIVE_INFINITY);
        UB = (getSimpleBounds().getUpper() != null) ? new ArrayRealVector(getSimpleBounds().getUpper()) : new ArrayRealVector(n, Double.POSITIVE_INFINITY);

        // Raccogli righe e termini noti
        java.util.ArrayList<double[]> rows = new java.util.ArrayList<>();
        java.util.ArrayList<Double> rhs = new java.util.ArrayList<>();

        for (int i = 0; i < n; i++) {
            // lower:  e_i^T x ≥ lb_i
            if (Double.isFinite(LB.getEntry(i))) {
                double[] row = new double[n];
                row[i] = 1.0;
                rows.add(row);
                rhs.add(LB.getEntry(i));
            }
            // upper: -e_i^T x ≥ -ub_i  (equivalente a x_i ≤ ub_i)
            if (Double.isFinite(UB.getEntry(i))) {
                double[] row = new double[n];
                row[i] = -1.0;
                rows.add(row);
                rhs.add(-UB.getEntry(i));
            }
        }

        if (rows.isEmpty()) {
            this.bounds = null;
            this.JB = null;
            return 0;
        }

        // Costruisci A e b
        double[][] amat = rows.toArray(new double[rows.size()][]);
        double[] bvec = new double[rhs.size()];
        for (int k = 0; k < rhs.size(); k++) {
            bvec[k] = rhs.get(k);
        }

        RealMatrix A = new Array2DRowRealMatrix(amat, false); // no copy
        RealVector b = new ArrayRealVector(bvec, false);

        this.bounds = new LinearInequalityConstraint(A, b);
        this.JB = this.bounds.jacobian(null);
        return JB.getRowDimension();
    }

    /**
     * Initial guess build initial guess taking account the bounds
     *
     * @return initial guess vector
     */
    private RealVector initGuess() {
        int n = this.getObj().dim();       
        SimpleBounds sb = this.getSimpleBounds();
        double[] lb = (sb != null) ? sb.getLower() : null;
        double[] ub = (sb != null) ? sb.getUpper() : null;

        double[] xArr = hasStartPoint() ? Arrays.copyOf(getStartPoint(), n) : new double[n];

        for (int i = 0; i < n; i++) {
            double li = (lb != null && i < lb.length && Double.isFinite(lb[i])) ? lb[i] : Double.NEGATIVE_INFINITY;
            double ui = (ub != null && i < ub.length && Double.isFinite(ub[i])) ? ub[i] : Double.POSITIVE_INFINITY;

            // create a initial guess in the bounds
            if (Double.isFinite(li) && Double.isFinite(ui)) {
                xArr[i] = 0.5 * (li + ui);     // midpoint
            } else if (Double.isFinite(ui)) {
                xArr[i] = ui - 1.0;
            } else if (Double.isFinite(li)) {
                xArr[i] = li + 1.0;
            } else {
                xArr[i] = 0.0;
            }

        }
        return new ArrayRealVector(xArr, false);
    }

    /**
     * Project direction adjust QP solution direction to be sure respect the
     * bounds;
     */
    private void projectDirectionInPlace(RealVector x,
            RealVector d, boolean onDirection) {
        final int n = x.getDimension();
        if (onDirection) {
            for (int i = 0; i < n; i++) {
                double xi = x.getEntry(i);
                double minStep = LB.getEntry(i) - xi;  // può essere -∞
                double maxStep = UB.getEntry(i) - xi;  // può essere +∞
                double di = d.getEntry(i);
                if (di < minStep) {
                    di = minStep;
                }
                if (di > maxStep) {
                    di = maxStep;
                }

                d.setEntry(i, di);
            }
        } else {
            for (int i = 0; i < n; i++) {
                if (x.getEntry(i) < LB.getEntry(i)) {
                    x.setEntry(i, LB.getEntry(i));
                }
                if (x.getEntry(i) > UB.getEntry(i)) {
                    x.setEntry(i, UB.getEntry(i));
                }
            }
        }
    }

}
