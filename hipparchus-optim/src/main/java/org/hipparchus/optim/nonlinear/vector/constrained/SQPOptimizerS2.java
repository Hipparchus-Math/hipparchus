/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
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
import static org.hipparchus.optim.nonlinear.vector.constrained.GradientMode.EXTERNAL;
import static org.hipparchus.optim.nonlinear.vector.constrained.GradientMode.FORWARD;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Pair;
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
 * Supports: equality constraints, inequality constraints,bounds constraints, penalty function
 * updates, line search strategies, BFGS Hessian update, and augmented QP
 * formulation using a relaxation variable.
 * </p>
 *
 * @since 3.1
 */
public class SQPOptimizerS2 extends AbstractSQPOptimizer2 {

    private RealMatrix L;

    private enum QPMode {
        /**
         * Solving the standard QP subproblem (no slack variable).
         */
        QP_STARDARD,
        /**
         * Solving the augmented QP subproblem (with slack variable).
         */
        QP_AUGMENTED;
    }

    /**
     * Logger.
     */
    private final SQPLogger formatter = SQPLogger.defaultLogger();

    /**
     * Value of the equality constraints (residual v = f(x) - lower).
     */
    private RealVector eqEval;

    /**
     * Value of the inequality constraints (residual v = f(x) - lower).
     */
    private RealVector ineqEval;

    /**
     * Value of the bounds (residual v = Ax - lower).
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
     * Jacobian of the bounds.
     */
    private RealMatrix JB;

    /**
     * Evaluation of the objective function.
     */
    private double functionEval;

    /**
     * Old objective function evaluation.
     */
    private double functionEvalOld;

    /**
     * Current point.
     */
    private RealVector x;

    /**
     * Bounds in form of inequality constraints.
     */
    private LinearInequalityConstraint BIQ;

    /**
     * Lower bound.
     */
    private ArrayRealVector LB;

    /**
     * Upper bound.
     */
    private ArrayRealVector UB;

    /** Simple constructor.
     */
    public SQPOptimizerS2() {
        // nothing to do
    }

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

        // Equality constraints configuration
        if (this.getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        // Inequality constraints configuration
        if (this.getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
        }

        mb = buildBoundsAsInequalities();

        final int m = me + mi + mb;

        // Get initial value of slack variable from SQPOption
        double rho = getSettings().getRhoCons(); 
        x = (this.getStartPoint() != null) ? new ArrayRealVector(this.getStartPoint()) : initGuess();

        RealVector y = new ArrayRealVector(me + mi + mb, 0.0);

        // All the function and constraint evaluations will be performed inside the penalty function
        MeritFunctionL2 penalty = new MeritFunctionL2(this.getObj(), this.getEqConstraint(), this.getIqConstraint(), this.BIQ, x);

        LineSearch lineSearch = new LineSearch(getSettings().getAlphaMin(), getSettings().getHistory(), getSettings().getMu(), getSettings().getB(),
                getSettings().getMaxLineSearchIteration(), 0);

        // Initial values evaluation (Penalty function stores residuals v = eval - lower)
        functionEval = penalty.getObjEval();
        if (this.getEqConstraint() != null) {
            eqEval = penalty.getEqEval();
        }
        if (this.getIqConstraint() != null) {
            ineqEval = penalty.getIqEval();
        }
        if (this.BIQ != null) {
            bEval = penalty.getBEval();
        }
        
        double EPS = this.getSettings().getEps();
        double EPS2 = EPS * EPS;
        double sqrtEPS = FastMath.sqrt(EPS);
        
        computeGradients();
        
        double gamma = 1.0;
        H = MatrixUtils.createRealIdentityMatrix(x.getDimension()).scalarMultiply(FastMath.sqrt(gamma));

        BFGSUpdater bfgs = new BFGSUpdater(H, EPS, true, getMatrixDecompositionTolerance().getEpsMatrixDecomposition());
        
        L = bfgs.getL();
        RealVector dx = new ArrayRealVector(x.getDimension());
        RealVector u = new ArrayRealVector(y.getDimension());
        penalty.update(J, JE, JI, x, y, dx, u);
        
        QPMode QPMODE = QPMode.QP_STARDARD;

        boolean crit0 = false;
        boolean crit1 = false;
        boolean crit2 = false;
        boolean crit3 = false;
        boolean crit4 = false;
        boolean crit5 = false;
        double sigma = 0.0;
        double alpha = 0.0;
        RealVector lagOld = null;
        RealVector lagNew = null;

        LagrangeSolution qpSolution = null;
        
        boolean GRADFAIL = false;
        boolean FALLBACK = false;
        double VIOLATION = constraintViolation();
        double COMPLEMENTARY = Double.POSITIVE_INFINITY;
        double KKT = Double.POSITIVE_INFINITY;
        double DHD = 0.0;
        double XNORM = 0.0;
        double FUNDIFF = 0.0;
        int BFGSUPDATE = 0;
        int FUNEVALCOUNT=0;
        SQPCustomCriterion customCriteria = getSettings().getConvergenceFunction();  
        for (int i = 0; i < getSettings().getMaxIteration(); i++) {
            
            sigma = getSettings().getSigmaMax() * 10.0;
            while ((sigma > getSettings().getSigmaMax() || sigma < 0.0) && rho < 1.0e9 && !FALLBACK) {
                
                qpSolution = (QPMODE == QPMode.QP_AUGMENTED) ? solveAugmentedQP(y, rho) : solveQP();
                sigma = (qpSolution == null || qpSolution.getX().getDimension() == 0 || Double.isInfinite(qpSolution.getX().getMaxValue()) ||  Double.isInfinite(qpSolution.getX().getMinValue())) ? getSettings().getSigmaMax() * 10.0 : qpSolution.getValue();
                if ((sigma > getSettings().getSigmaMax() || sigma < 0.0)) {
                     if (QPMODE == QPMode.QP_AUGMENTED) rho = rhoUp(rho);
                     QPMODE = QPMode.QP_AUGMENTED;
                }
            }
            
            // If sigma > sigma threshold after several attempts, assign direction from penalty gradient
            if (rho >= 1.0e9 || FALLBACK) {
                boolean RESET = false;
                qpSolution = solveQPFallBack(penalty.gradX());

                    if (qpSolution == null || qpSolution.getX().getDimension() == 0 ) {
                    
                        bfgs.resetHessian();
                        H=bfgs.getHessian();
                        L=bfgs.getL();
                        qpSolution = solveQPFallBack(penalty.gradX());
                        if (qpSolution == null || qpSolution.getX().getDimension() == 0 ) {
                            penalty.resetRj();
                            if (m > 0) {
                                u.set(0.0);
                                y.set(0.0);
                            }
                            penalty.update(J, JE, JI, x, y, new ArrayRealVector(x.getDimension()), u);
                            qpSolution = solveQPFallBack(penalty.gradX());
                            RESET=true;

                            if (qpSolution == null || qpSolution.getX().getDimension() == 0 ) break;
//                        if (mb > 0) {
//                         RealVector db = qpSolution.getLambda();
//                         u.setSubVector(mi + me, db);
//                              }
                        
                        }  // Infeasible
                    // Estimation of multiplier from penalty grad y
                
                
                    }
                if (m > 0 ) {
                    u = y.subtract(penalty.gradY());
                }
                
                
                dx = qpSolution.getX();
                projectDirectionInPlace(x, dx, true);
                
                
                QPMODE = QPMode.QP_AUGMENTED;
                sigma = 0.0;
                rho = getSettings().getRhoCons();

            } else {
                dx = qpSolution.getX();
                u = qpSolution.getLambda();
                
                sigma = qpSolution.getValue();
                projectDirectionInPlace(x, dx, true);
                
                penalty.updateRj(H, y, dx, u, sigma, iterations.getCount());
                
                // Switch to normal QP if additional variable is small enough
                if (QPMODE == QPMode.QP_AUGMENTED && FastMath.abs(sigma) < getSettings().getEps() && !GRADFAIL) {
                    QPMODE = QPMode.QP_STARDARD;
                    rho = getSettings().getRhoCons();
                }                
            }

            penalty.update(J, JE, JI, x, y, dx, u);

            lagOld= lagrangianGradX(J, JE, JI, x, u);
            double dlan = lagOld.getNorm();
            
            // Convergence check before line search
            KKT = dlan * dlan;
            DHD = dx.dotProduct(H.operate(dx));
            COMPLEMENTARY = complementarySlackness(u, me);
            XNORM = dx.getNorm();
            
            crit0 = KKT <= EPS;
            crit1 = DHD <= EPS2;
            crit2 = XNORM <= EPS * (1.0 + x.getNorm());
            crit3 = VIOLATION <= sqrtEPS;
            crit5 = COMPLEMENTARY <= EPS;
              
            if ((crit0 || crit1) && crit5  && crit3) {
                break;
            }
            
            
            
            if (penalty.getGradient() >= 0) {
                if (QPMODE == QPMode.QP_AUGMENTED)  rho = rhoUp(rho);
                QPMODE = QPMode.QP_AUGMENTED;
                GRADFAIL = true;
            }    
            // If penalty gradient is < 0 proceed with line search
             else  {
                iterations.increment();
                GRADFAIL = false;
                FALLBACK = false;
                
                if (QPMODE == QPMode.QP_AUGMENTED) rho = updateRho(dx, u, sigma);
                
                // Line Search
                alpha = lineSearch.search(penalty);
                FUNEVALCOUNT+=lineSearch.getIteration();
                //UPDATE MULTIPLIER
                if (m > 0) {
                   // y = y.add(u.subtract(y).mapMultiply(alpha));
                    //y.combineToSelf(1.0 - alpha, alpha, u);//this is worse
                    double yk=0;
                    for (int k = 0; k < y.getDimension(); ++k) {
                        yk = y.getEntry(k);
                        y.setEntry(k, yk + (u.getEntry(k) - yk) * alpha);
                    }

                }
                //UPDATE SOLUTION
                //x = x.add(dx.mapMultiply(alpha));
                x.combineToSelf(1.0, alpha, dx);
                // UPDATE FUNCTION AND CONSTRAINT EVAL (Residuals v = eval - lower are stored)
                functionEvalOld = functionEval;
                functionEval = penalty.getObjEval();
                eqEval = penalty.getEqEval();
                ineqEval = penalty.getIqEval();
                bEval = penalty.getBEval();
                
                // Convergence check after line search
                COMPLEMENTARY = complementarySlackness(u, me);
                VIOLATION = constraintViolation();
                lagNew= lagrangianGradX(J, JE, JI, x, u);
                dlan = lagNew.getNorm();
                KKT = dlan * dlan;
                // Step length
                if (alpha > 0) {
                    DHD = DHD * alpha * alpha;
                    XNORM = alpha * dx.getNorm();
                }
                
                FUNDIFF = FastMath.abs(functionEval - functionEvalOld);              
   
                crit0 = (KKT <= EPS);
                crit1 = DHD <= EPS2;
                crit2 = XNORM <= EPS * (1.0 + x.getNorm());
                crit3 = VIOLATION <= sqrtEPS;
                crit4 = FUNDIFF < EPS * (1.0 + FastMath.abs(functionEval));
                crit5 = (COMPLEMENTARY <= EPS);
                
//                if (!lineSearch.isBadStepDetected() && (crit0 || (crit5 && m > 0) || crit4 || crit2 ) && crit1   && crit3) {
//                    break;
//                }
//                if (manualBreak) {
//                    System.out.println("\n[Solver] Manual break triggered by the user.");
//                    break; // Exits the loop exactly at this checkpoint
//                }
                if(customCriteria!=null && customCriteria.converged(alpha, XNORM, DHD, COMPLEMENTARY, KKT, VIOLATION , FUNDIFF,x.getNorm(),FastMath.abs(functionEval))){
                    break;
                }
                
                 if (!lineSearch.isBadStepDetected() && ((crit0 || crit1) && crit5  && crit3)) {
                    break;
                }
//              
                if (lineSearch.isBadStepDetected() && (crit3 && (crit4 || crit2))) {
                    break;
                }
                
                computeGradients();
                RealVector lagnew = lagrangianGradX(J, JE, JI, x, u);
                
                // Hessian update with the logic of line search and internal damping logic
                if (lineSearch.isBadStepFailed()) {
                    FALLBACK = false;
                    bfgs.resetHessian();
                    H = bfgs.getHessian();
                    L = bfgs.getL();
                    lineSearch.resetBadStepCount();
                    BFGSUPDATE = -2;

               } else if (lineSearch.isBadStepDetected()) {
                    BFGSUPDATE = -1;
                    L = bfgs.getL();
                    H = bfgs.getHessian();
                } else {
                    BFGSUPDATE = bfgs.update(dx.mapMultiply(alpha), lagnew.subtract(lagOld));
                    H = bfgs.getHessian();
                    L = bfgs.getL();
                }
                
                formatter.logRow(iterations.getCount(),
                        alpha, lineSearch.getIteration(),
                        XNORM, DHD, COMPLEMENTARY, KKT, VIOLATION, sigma,
                        penalty.getPenaltyEval(), functionEval, FUNDIFF, BFGSUPDATE);
            }
        }
        
        formatter.logRow(iterations.getCount(),
                        alpha, lineSearch.getIteration(),
                        XNORM, DHD, COMPLEMENTARY, KKT, VIOLATION, sigma,
                        penalty.getPenaltyEval(), functionEval, FUNDIFF, -99);
        FUNEVALCOUNT+=iterations.getCount();
        formatter.logRow(FUNEVALCOUNT,crit2, crit1, crit5, crit0, crit3, crit4);

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
     * Compute constraints violations using residuals directly.
     *
     * @return constraints violations
     */
    private double constraintViolation() {
        double crit = 0.0;
        
        if (this.getEqConstraint() != null) {
           
            for (int k = 0; k < eqEval.getDimension(); k++) {
                crit += FastMath.abs(eqEval.getEntry(k));
            }
        }
        
        if (this.getIqConstraint() != null) {
            
            for (int k = 0; k < ineqEval.getDimension(); k++) {
                if (ineqEval.getEntry(k)< 0.0) {
                    crit -= ineqEval.getEntry(k); 
                }
            }
        } 
        
        if (this.BIQ != null) {
            
            for (int k = 0; k < bEval.getDimension(); k++) {
                if (bEval.getEntry(k) < 0.0) {
                    crit -= bEval.getEntry(k);
                }
            }
        } 
        
        return crit;
    }

    /**
     * Compute complementary slackness using residuals directly.
     *
     * @param u  the multiplier vector
     * @param me the number of equality constraints
     * @return complementary slackness
     */
    private double complementarySlackness(RealVector u, int me) {
        if (u == null) {
            return 0.0;
        }
        
        double crit = 0.0;
        int shift = me;
        
         
        if (this.getIqConstraint() != null) {
           
            for (int i = 0; i < ineqEval.getDimension(); i++) {
                crit += FastMath.abs(ineqEval.getEntry(i) * u.getEntry(shift + i));
            }
            shift += ineqEval.getDimension();
        }
        
        if (this.BIQ != null) {
           
            for (int i = 0; i < bEval.getDimension(); i++) {
                crit += FastMath.abs(bEval.getEntry(i) * u.getEntry(shift + i));
            }
        }
        
        return crit;
    }

    /**
     * Increase the rho penalty parameter.
     *
     * @param rho current rho
     * @return increased rho
     */
    private double rhoUp(double rho) {
        return rho * 10.0;
    }

    /**
     * Dynamically update rho based on the current step and multiplier difference.
     *
     * @param dx current step
     * @param dy multiplier step
     * @param additionalVariable slack/augmented variable
     * @return new rho value
     */
    private double updateRho(final RealVector dx, final RealVector dy, final double additionalVariable) {
        int me = JE != null ? JE.getRowDimension() : 0;
        int mi = JI != null ? JI.getRowDimension() : 0;
        int mb = JB != null ? JB.getRowDimension() : 0;
        
        if (me + mi + mb > 0) {
            
            RealVector jacPreDy = new ArrayRealVector(x.getDimension());
            int offset = 0;

            if (me > 0 && JE != null) {
                RealVector dyE = dy.getSubVector(offset, me);
                jacPreDy.combineToSelf(1.0, 1.0, JE.preMultiply(dyE));
                offset += me;
            }
            if (mi > 0 && JI != null) {
                RealVector dyI = dy.getSubVector(offset, mi);
                jacPreDy.combineToSelf(1.0, 1.0, JI.preMultiply(dyI));
                offset += mi;
            }
            if (mb > 0 && JB != null) {
                RealVector dyB = dy.getSubVector(offset, mb);
                jacPreDy.combineToSelf(1.0, 1.0, JB.preMultiply(dyB));
            }

            
            double num = 10.0 * FastMath.pow(dx.dotProduct(jacPreDy), 2);
            double den = (1.0 - additionalVariable) * (1.0 - additionalVariable) * dx.dotProduct(H.operate(dx));
            
            if (den < Precision.SAFE_MIN) {
                den = Precision.SAFE_MIN;
            }

            return FastMath.max(10.0, num / den);
        }
        return 1.0;
    }

    /**
     * Solve augmented problem using residuals.
     *
     * @param y Lagrange multipliers
     * @param rho rho penalty parameter
     * @return problem solution
     */
    /**
 * Solve augmented problem using residuals.
 * <p>
 * Mathematically equivalent to the original version, but optimized to avoid
 * double[][] allocations and implementation-specific casts.
 * </p>
 */
private LagrangeSolution solveAugmentedQP(final RealVector y, final double rho) {
    final int n = x.getDimension();
    final double eps = getSettings().getEps();

    final int me = (getEqConstraint() != null) ? getEqConstraint().dimY() : 0;
    final int mi = (getIqConstraint() != null) ? getIqConstraint().dimY() : 0;
    final int mb = (BIQ != null) ? BIQ.dimY() : 0;
    final int mc = mi + mb;
    
    boolean violated = false;
    if (mi > 0) {
        violated = ineqEval.getMinValue() <= eps || y.getSubVector(me, mi).getMaxValue() > 0;
    }

    final int add = (me > 0 || violated) ? 1 : 0;

    // --- OBJECTIVE ---
    final RealVector g1 = new ArrayRealVector(n + add);
    g1.setSubVector(0, J);

    // L1: [L, 0; 0, sqrt(rho)]
    final RealMatrix L1 = MatrixUtils.createRealMatrix(n + add, n + add);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L1.setEntry(i, j, L.getEntry(i, j));
        }
    }
    if (add == 1) {
        L1.setEntry(n, n, FastMath.sqrt(rho));
    }

    // --- EQUALITIES ---
    LinearEqualityConstraint eqc = null;
    if (me > 0) {
        final RealMatrix Ae = MatrixUtils.createRealMatrix(me, n + add);
        final RealVector eqRhs = new ArrayRealVector(me);
        for (int i = 0; i < me; i++) {
            final double vi = eqEval.getEntry(i);
            eqRhs.setEntry(i, -vi);
            if (add == 1) {
                Ae.setEntry(i, n, -vi);
            }
            for (int j = 0; j < n; j++) {
                Ae.setEntry(i, j, JE.getEntry(i, j));
            }
        }
        eqc = new LinearEqualityConstraint(Ae, eqRhs);
    }
    
    // --- INEQUALITIES & BOUNDS ---
    LinearInequalityConstraint iqc = null;
    if (mc > 0) {
        final RealMatrix A = MatrixUtils.createRealMatrix(mc, n + add);
        final RealVector B = new ArrayRealVector(mc);

        // Part 1: Inequalities
        if (mi > 0) {
            for (int i = 0; i < mi; i++) {
                final double vi = ineqEval.getEntry(i);
                B.setEntry(i, -vi);
                if (add == 1) {
                    if (vi <= eps || y.getEntry(me + i) > 0) {
                        A.setEntry(i, n, -vi);
                    }
                }
                for (int j = 0; j < n; j++) {
                    A.setEntry(i, j, JI.getEntry(i, j));
                }
            }
        }
        
        // Part 2: Bounds (Stacked below inequalities)
        if (mb > 0) {
            for (int i = 0; i < mb; i++) {
                B.setEntry(mi + i, -bEval.getEntry(i));
                for (int j = 0; j < n; j++) {
                    A.setEntry(mi + i, j, JB.getEntry(i, j));
                }
            }
        }
        iqc = new LinearInequalityConstraint(A, B);
    }

    // --- SLACK BOUNDS ---
    LinearBoundedConstraint bc = null;
    if (add == 1) {
        final RealMatrix sigmaA = MatrixUtils.createRealMatrix(1, n + 1);
        sigmaA.setEntry(0, n, 1.0);
        bc = new LinearBoundedConstraint(sigmaA, new ArrayRealVector(1, 0.0), new ArrayRealVector(1, 1.0));
    }
    
    // Solve
    QPDualActiveSolverOption optQP=new QPDualActiveSolverOption();
    optQP.setEps(this.getSettings().getEpsQP());
    optQP.setEpsRelaxation(this.getSettings().getEpsQP());
    optQP.setMatrixMode(QPMatrixMode.CHOLESKY);
    final LagrangeSolution sol = getQPSolver().optimize(optQP,
            new ObjectiveFunction(new QuadraticFunction(L1, g1, 0)), 
            iqc, eqc, bc);
    
    if (sol == null || sol.getX().getDimension() == 0) {
        return sol;
    }
    
    final RealVector fullX = sol.getX();
    final double sigma = (add == 1) ? fullX.getEntry(n) : 0.0;
    
    return new LagrangeSolution(fullX.getSubVector(0, n), 
                                (me + mc == 0) ? null : sol.getLambda().getSubVector(0, me + mc), 
                                sigma);
}
    
   /**
     * Solves the standard Quadratic Programming (QP) subproblem using residuals.
     * * @return a LagrangeSolution containing primal and dual variables, or null if failed.
     */
    private LagrangeSolution solveQP() {
        final int n = x.getDimension();
        final int me = (getEqConstraint() != null) ? getEqConstraint().dimY() : 0;
        final int mi = (getIqConstraint() != null) ? getIqConstraint().dimY() : 0;
        final int mb = (BIQ != null) ? BIQ.dimY() : 0;
        final int mc = mi + mb;

        // --- EQUALITY CONSTRAINTS SETUP ---
        LinearEqualityConstraint eqc = null;
        if (me > 0) {
            // Formulation: JE * d = -eqEval
            eqc = new LinearEqualityConstraint(JE, eqEval.mapMultiply(-1.0));
        }

        // --- INEQUALITY CONSTRAINTS SETUP ---
        LinearInequalityConstraint iqc = null;
        if (mc > 0) {
            // Formulation: [JI; JB] * d >= -[ineqEval; bEval]
            // Pre-allocate global matrix to avoid BlockRealMatrix fragmentation
            final RealMatrix combinedAI = MatrixUtils.createRealMatrix(mc, n);
            final RealVector combinedBI = new ArrayRealVector(mc);

            // Copy Inequality Jacobian (JI) directly into the combined matrix
            if (mi > 0) {
                for (int i = 0; i < mi; i++) {
                    combinedBI.setEntry(i, -ineqEval.getEntry(i));
                    for (int j = 0; j < n; j++) {
                        combinedAI.setEntry(i, j, JI.getEntry(i, j));
                    }
                }
            }

            // Copy Bounds Jacobian (JB) directly into the combined matrix (offset by mi)
            if (mb > 0) {
                for (int i = 0; i < mb; i++) {
                    combinedBI.setEntry(mi + i, -bEval.getEntry(i));
                    for (int j = 0; j < n; j++) {
                        combinedAI.setEntry(mi + i, j, JB.getEntry(i, j));
                    }
                }
            }
            iqc = new LinearInequalityConstraint(combinedAI, combinedBI);
        }

        
        QPDualActiveSolverOption optQP=new QPDualActiveSolverOption();
        optQP.setEps(this.getSettings().getEpsQP());
        optQP.setEpsRelaxation(this.getSettings().getEpsQP());
        optQP.setMatrixMode(QPMatrixMode.CHOLESKY);
        // --- QUADRATIC OBJECTIVE ---
        // Objective: 0.5 * d^T * L * d + J^T * d
        final QuadraticFunction q = new QuadraticFunction(L, J, 0.0);

        // --- SOLVE QP ---
        final LagrangeSolution sol = getQPSolver().optimize(optQP,new ObjectiveFunction(q), iqc, eqc);

        if (sol == null || sol.getX().getDimension() == 0) {
            return null;
        }

        // IMPORTANT: Sigma must be 0.0 in standard QP to avoid triggering slack logic in SQP loop
        return new LagrangeSolution(sol.getX(), 
                                    (me + mc > 0) ? sol.getLambda() : new ArrayRealVector(0), 
                                    0.0);
    }
    

    /**
     * Solves the Quadratic Programming (QP) subproblem in the current SQP
     * iteration if penalty gradient is >= 0
     *
     * @param penaltyGradx gradient of the penalty function
     * @return a {@link LagrangeSolution} representing the QP solution, or
     * {@code null} if the QP failed
     */
    private LagrangeSolution solveQPFallBack(RealVector penaltyGradx) {

        final QuadraticFunction q = new QuadraticFunction(this.H ,penaltyGradx, 0);
        int mb = 0;

        LinearInequalityConstraint iqc = null;

        if (BIQ != null) {
            mb = BIQ.dimY();
            RealMatrix A = JB;
            RealVector B = bEval.mapMultiply(-1.0);
            iqc = new LinearInequalityConstraint(A, B);
        }
        QPDualActiveSolverOption optQP=new QPDualActiveSolverOption();
        optQP.setEps(this.getSettings().getEpsQP());
        optQP.setEpsRelaxation(this.getSettings().getEpsQP());
        optQP.setMatrixMode(QPMatrixMode.FULL);
        LagrangeSolution sol = getQPSolver().optimize(optQP,new ObjectiveFunction(q), iqc);
        
        if (sol == null || sol.getX().getDimension() == 0) {
            return sol;
        }

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
     * L(x, y) = f(x) - yₑᵗ·cₑ(x) - yᵢᵗ·cᵢ(x)
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
     * ∇ₓ L(x, y) = ∇f(x) - JEᵗ·yₑ - JIᵗ·yᵢ
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

        if (otherJE != null) {
            int me = otherJE.getRowDimension();
            RealVector yEq = y.getSubVector(0, me);
            gradL.combineToSelf(1.0, -1.0, otherJE.preMultiply(yEq));
            offset += me;
        }

        if (otherJI != null) {
            int mi = otherJI.getRowDimension();
            RealVector yIq = y.getSubVector(offset, mi);
            gradL.combineToSelf(1.0, -1.0, otherJI.preMultiply(yIq));
            offset += mi;
        }

        if (JB != null) {
            int mb = JB.getRowDimension();
            RealVector yB = y.getSubVector(offset, mb);
            gradL.combineToSelf(1.0, -1.0, JB.preMultiply(yB));
        }

        return gradL;
    }

    /**
     * Computes the optimal finite difference step size and direction for a given variable,
     * strictly respecting its lower and upper bounds.
     * <p>
     * To prevent catastrophic cancellation and numerical noise, this method strictly
     * avoids arbitrary step shrinking. If the remaining room up to the active bound 
     * is strictly less than the required numerical step, the variable is 
     * considered numerically frozen and a zero step is returned.
     * </p>
     * <p>
     * For central differences, if a symmetric step is not feasible, the method attempts 
     * to allocate enough room for an asymmetric 3-point formula step (2 * hBase) 
     * to preserve the <i>O(h^2)</i> truncation error.
     * </p>
     *
     * @param i       the index of the variable
     * @param xi      the current value of the variable
     * @param hBase   the theoretical base step size (e.g., square root or cube root of machine epsilon)
     * @param central if {@code true}, step selection prioritizes central finite differences
     * @return a {@link Pair} where the first element is the signed step size {@code h},
     * and the second element is the direction indicator:
     * <ul>
     * <li>{@code 0} for true central difference</li>
     * <li>{@code +1} for forward one-sided difference</li>
     * <li>{@code -1} for backward one-sided difference</li>
     * </ul>
     * Returns {@code (0.0, 0)} if the variable is strictly bounded and cannot be perturbed safely.
     */
    private Pair<Double, Integer> chooseStep(final int i,
                                             final double xi,
                                             final double hBase,
                                             final boolean central) {

        final double lb = LB.getEntry(i);
        final double ub = UB.getEntry(i);

        final double roomMinus = FastMath.max(0.0,
                Double.isFinite(lb) ? (xi - lb) : Double.POSITIVE_INFINITY);
        final double roomPlus = FastMath.max(0.0,
                Double.isFinite(ub) ? (ub - xi) : Double.POSITIVE_INFINITY);

        if (!central) {
            if (roomPlus >= hBase) {
                return new Pair<>(hBase, +1);
            }
            if (roomMinus >= hBase) {
                return new Pair<>(-hBase, -1);
            }
            return new Pair<>(0.0, 0);
        }

        if (roomPlus >= hBase && roomMinus >= hBase) {
            return new Pair<>(hBase, 0);
        }

        if (roomPlus >= 2.0 * hBase) {
            return new Pair<>(hBase, +1);
        }
        if (roomMinus >= 2.0 * hBase) {
            return new Pair<>(-hBase, -1);
        }

        return new Pair<>(0.0, 0);
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
     * Each variable is perturbed independently by a small step size proportional 
     * to the square root of machine precision. Variables tightly bounded are 
     * gracefully ignored (gradient component set to zero) to prevent out-of-bounds 
     * evaluations and Hessian corruption.
     * </p>
     * <p>
     * Deep copies of the working vectors are instantiated within the loop to 
     * strictly prevent side effects and aliasing issues with external black-box 
     * evaluators that may cache object references.
     * </p>
     */
    private void forwardGradient() {

        int n = x.getDimension();
        double sqrtEps = FastMath.sqrt(getSettings().getGradientEps());
        double xi;
        double h1;
        double fRef = this.functionEval;
        RealVector eqRef = this.eqEval;
        RealVector iqRef = this.ineqEval;

        RealVector gradF = new ArrayRealVector(n);
        RealMatrix gradEq = (getEqConstraint() != null) ? MatrixUtils.createRealMatrix(eqRef.getDimension(), n) : null;
        RealMatrix gradIq = (getIqConstraint() != null) ? MatrixUtils.createRealMatrix(iqRef.getDimension(), n) : null;

        for (int i = 0; i < n; i++) {
            xi = x.getEntry(i);
            h1 = sqrtEps * FastMath.max(1.0, FastMath.abs(xi));
            
            Pair<Double, Integer> stepInfo = chooseStep(i, xi, h1, false);
            double h = stepInfo.getFirst();
            
            if (h == 0.0) {
                gradF.setEntry(i, 0.0);
                if (gradEq != null) {
                    gradEq.setColumnVector(i, new ArrayRealVector(eqEval.getDimension(), 0.0));
                }
                if (gradIq != null) {
                    gradIq.setColumnVector(i, new ArrayRealVector(ineqEval.getDimension(), 0.0));
                }
            } else {
                RealVector xPerturbed = x.copy();
                xPerturbed.setEntry(i, xi + h);
                
                double fPerturbed = getObj().value(xPerturbed);
                gradF.setEntry(i, (fPerturbed - fRef) / h);

                if (gradEq != null) {
                    // Logica residuo: (valPert - lower) - (valRef - lower) = valPert - valRef
                    RealVector eqPerturbed = getEqConstraint().value(xPerturbed).combineToSelf(1.0,-1.0,getEqConstraint().getLowerBound());
                    gradEq.setColumnVector(i, eqPerturbed.combineToSelf(1.0,-1.0,eqRef).mapDivideToSelf(h));
                }

                if (gradIq != null) {
                    RealVector iqPerturbed = getIqConstraint().value(xPerturbed).combineToSelf(1.0,-1.0,getIqConstraint().getLowerBound());
                    gradIq.setColumnVector(i, iqPerturbed.combineToSelf(1.0,-1.0,iqRef).mapDivideToSelf(h));
                }
            }
        }

        this.J = gradF;
        this.JE = gradEq;
        this.JI = gradIq;
    }

    /**
     * Computes the gradient of the objective function and the Jacobians of the
     * constraints using centered finite differences (second-order accurate).
     * <p>
     * The base step size is proportional to the cube root of machine precision. 
     * If a symmetric centered difference is not possible due to variable bounds, 
     * the algorithm falls back to an asymmetric 3-point formula. This industrial-grade 
     * fallback guarantees that the O(h^2) truncation error is preserved, preventing 
     * catastrophic shocks to the BFGS Hessian approximation near active bounds.
     * </p>
     * <p>
     * Vectors are explicitly instantiated to prevent aliasing with stateful objective 
     * functions. Furthermore, evaluations are strictly grouped by the perturbed point 
     * (e.g., all functions evaluated at x+h, then all at x-h) to maximize cache-hit 
     * ratios in complex black-box simulators.
     * </p>
     */
    private void centralGradient() {

        int n = x.getDimension();
        double hBase = FastMath.cbrt(getSettings().getGradientEps());

        RealVector gradF = new ArrayRealVector(n);
        RealMatrix gradEq = (getEqConstraint() != null) ? MatrixUtils.createRealMatrix(eqEval.getDimension(), n) : null;
        RealMatrix gradIq = (getIqConstraint() != null) ? MatrixUtils.createRealMatrix(ineqEval.getDimension(), n) : null;

        for (int i = 0; i < n; i++) {
            double xi = x.getEntry(i);
            double h1 = hBase * FastMath.max(1.0, FastMath.abs(xi));

            Pair<Double, Integer> stepInfo = chooseStep(i, xi, h1, true);
            double h = stepInfo.getFirst();
            int dir = stepInfo.getSecond();
            
            if (h == 0.0) {
                // Variable is completely trapped, numerical derivative is unfeasible
                gradF.setEntry(i, 0.0);
                if (gradEq != null) {
                    gradEq.setColumnVector(i, new ArrayRealVector(eqEval.getDimension(), 0.0));
                }
                if (gradIq != null) {
                    gradIq.setColumnVector(i, new ArrayRealVector(ineqEval.getDimension(), 0.0));
                }
                
            } else if (dir == 0) {
                // -------------------------------------------------
                // CLASSICAL CENTRAL DIFFERENCE
                // -------------------------------------------------
                RealVector xPlus = new ArrayRealVector(x);
                RealVector xMinus = new ArrayRealVector(x);
                xPlus.addToEntry(i, h);
                xMinus.addToEntry(i, -h);

                // --- CACHE-FRIENDLY EVALUATION: Everything at xPlus ---
                double fPlus = getObj().value(xPlus);
                // Calcoliamo i residui v = f(x) - lower per la differenza finita
                RealVector eqPlus = (gradEq != null) ? getEqConstraint().value(xPlus).combineToSelf(1.0,-1.0,getEqConstraint().getLowerBound()) : null;
                RealVector iqPlus = (gradIq != null) ? getIqConstraint().value(xPlus).combineToSelf(1.0,-1.0,getIqConstraint().getLowerBound()) : null;

                // --- CACHE-FRIENDLY EVALUATION: Everything at xMinus ---
                double fMinus = getObj().value(xMinus);
                RealVector eqMinus = (gradEq != null) ? getEqConstraint().value(xMinus).combineToSelf(1.0,-1.0,getEqConstraint().getLowerBound()) : null;
                RealVector iqMinus = (gradIq != null) ? getIqConstraint().value(xMinus).combineToSelf(1.0,-1.0,getIqConstraint().getLowerBound()) : null;

                // Derivatives Calculation (La differenza tra residui è uguale alla differenza tra valori raw)
                gradF.setEntry(i, (fPlus - fMinus) / (2.0 * h));
               
                if (gradEq != null) {
                    // Sfrutto eqPlus gia' allocato per la combinazione: (eqPlus - eqMinus) * inv2h
                    gradEq.setColumnVector(i, eqPlus.combineToSelf(1.0, -1.0, eqMinus).mapDivideToSelf(2.0*h));
                }

                if (gradIq != null) {
                    gradIq.setColumnVector(i, iqPlus.combineToSelf(1.0, -1.0, iqMinus).mapDivideToSelf(2.0*h));
                }

            } else {
                // -------------------------------------------------
                // ASYMMETRIC 3-POINT FORMULA (O(h^2) preservation)
                // -------------------------------------------------
                RealVector x1 = new ArrayRealVector(x);
                RealVector x2 = new ArrayRealVector(x);
                x1.addToEntry(i, h);
                x2.addToEntry(i, 2.0 * h);

                double f0 = functionEval; // Already evaluated at current x (residual if applicable)

                // --- CACHE-FRIENDLY EVALUATION: Everything at x1 ---
                double f1 = getObj().value(x1);
                RealVector eq1 = (gradEq != null) ? getEqConstraint().value(x1).combineToSelf(1.0,-1.0,getEqConstraint().getLowerBound()) : null;
                RealVector iq1 = (gradIq != null) ? getIqConstraint().value(x1).combineToSelf(1.0,-1.0,getIqConstraint().getLowerBound()) : null;

                // --- CACHE-FRIENDLY EVALUATION: Everything at x2 ---
                double f2 = getObj().value(x2);
                RealVector eq2 = (gradEq != null) ? getEqConstraint().value(x2).combineToSelf(1.0,-1.0,getEqConstraint().getLowerBound()) : null;
                RealVector iq2 = (gradIq != null) ? getIqConstraint().value(x2).combineToSelf(1.0,-1.0,getIqConstraint().getLowerBound()) : null;

                // Derivatives Calculation
                gradF.setEntry(i, (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h));

                if (gradEq != null) {
                    // eqEval è già un residuo (f(x) - L) calcolato in precedenza
                    RealVector dEq = eqEval.copy();
                    dEq.combineToSelf(-3.0, 4.0, eq1).combineToSelf(1.0, -1, eq2).mapDivideToSelf(2*h);
                          
                    gradEq.setColumnVector(i, dEq);
                }

                if (gradIq != null) {
                    // ineqEval è già un residuo (f(x) - L) calcolato in precedenza
                   RealVector dIq = ineqEval.copy();
                    dIq.combineToSelf(-3.0, 4.0, iq1).combineToSelf(1.0, -1, iq2).mapDivideToSelf(2*h);
                    gradIq.setColumnVector(i, dIq);
                }
            }
        }

        this.J = gradF;
        this.JE = gradEq;
        this.JI = gradIq;
    }
    /**
     * Set debug printer.
     * @param printer debug printer
     */
    public void setDebugPrinter(final DebugPrinter printer) {
        formatter.setDebugPrinter(printer);
    }

    /**
     * Build bounds constraint in the form of inequality.
     */
    private int buildBoundsAsInequalities() {
        if (getSimpleBounds() == null) {
            int n = getObj().dim();
            this.BIQ = null;
            LB = new ArrayRealVector(n, Double.NEGATIVE_INFINITY);
            UB = new ArrayRealVector(n, Double.POSITIVE_INFINITY);
            this.JB = null;
            return 0;
        }
        int n = getObj().dim();
        LB = (getSimpleBounds().getLower() != null) ? new ArrayRealVector(getSimpleBounds().getLower()) : new ArrayRealVector(n, Double.NEGATIVE_INFINITY);
        UB = (getSimpleBounds().getUpper() != null) ? new ArrayRealVector(getSimpleBounds().getUpper()) : new ArrayRealVector(n, Double.POSITIVE_INFINITY);

        java.util.ArrayList<double[]> rows = new java.util.ArrayList<>();
        java.util.ArrayList<Double> rhs = new java.util.ArrayList<>();

        for (int i = 0; i < n; i++) {
            if (Double.isFinite(LB.getEntry(i))) {
                double[] row = new double[n];
                row[i] = 1.0;
                rows.add(row);
                rhs.add(LB.getEntry(i));
            }
            if (Double.isFinite(UB.getEntry(i))) {
                double[] row = new double[n];
                row[i] = -1.0;
                rows.add(row);
                rhs.add(-UB.getEntry(i));
            }
        }

        if (rows.isEmpty()) {
            this.BIQ = null;
            this.JB = null;
            return 0;
        }

        double[][] amat = rows.toArray(new double[rows.size()][]);
        double[] bvec = new double[rhs.size()];
        for (int k = 0; k < rhs.size(); k++) {
            bvec[k] = rhs.get(k);
        }

        this.BIQ = new LinearInequalityConstraint(new Array2DRowRealMatrix(amat, false), new ArrayRealVector(bvec, false));
        this.JB = this.BIQ.jacobian(null);
        return JB.getRowDimension();
    }

    /**
     * Build initial guess taking account the bounds.
     */
    private RealVector initGuess() {
        int n = this.getObj().dim();        
        SimpleBounds sb = this.getSimpleBounds();
        double[] lbArr = (sb != null) ? sb.getLower() : null;
        double[] ubArr = (sb != null) ? sb.getUpper() : null;

        double[] xArr = (start != null) ? Arrays.copyOf(start, n) : new double[n];

        for (int i = 0; i < n; i++) {
            double li = (lbArr != null && i < lbArr.length && Double.isFinite(lbArr[i])) ? lbArr[i] : Double.NEGATIVE_INFINITY;
            double ui = (ubArr != null && i < ubArr.length && Double.isFinite(ubArr[i])) ? ubArr[i] : Double.POSITIVE_INFINITY;

            if (Double.isFinite(li) && Double.isFinite(ui)) {
                xArr[i] = 0.5 * (li + ui);
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
     * Project direction adjust QP solution direction to be sure respect the bounds.
     */
    private void projectDirectionInPlace(RealVector xVec, RealVector dVec, boolean onDirection) {
        final int n = xVec.getDimension();
        if (onDirection) {
            for (int i = 0; i < n; i++) {
                double xi = xVec.getEntry(i);
                double minStep = LB.getEntry(i) - xi;
                double maxStep = UB.getEntry(i) - xi;
                double di = dVec.getEntry(i);
                if (di < minStep) di = minStep;
                if (di > maxStep) di = maxStep;
                dVec.setEntry(i, di);
            }
        } else {
            for (int i = 0; i < n; i++) {
                if (xVec.getEntry(i) < LB.getEntry(i)) xVec.setEntry(i, LB.getEntry(i));
                if (xVec.getEntry(i) > UB.getEntry(i)) xVec.setEntry(i, UB.getEntry(i));
            }
        }
    }
}
