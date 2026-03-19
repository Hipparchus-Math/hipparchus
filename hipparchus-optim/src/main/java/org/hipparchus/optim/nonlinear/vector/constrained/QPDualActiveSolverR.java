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
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.DecompositionSolver;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Pair;
import org.hipparchus.util.Precision;

/**
 * Implements the dual active-set method by Goldfarb and Idnani (1983) for
 * solving strictly convex quadratic programs of the form:
 *
 * <pre>
 * minimize   (1/2) x^T G x + g0^T x
 * subject to CE^T x  =  ce0     (equality constraints)
 *            CI^T x  &gt;= ci0     (inequality constraints)
 *            lb &lt;= Ax &lt;= ub      (bounded constraints)
 * </pre>
 *
 * <p>
 * Uses an incremental QR factorization updater (QRUpdater) for managing active
 * constraints and maintains multipliers for dual and primal steps.
 * </p>
 *
 * @see <a href="https://doi.org/10.1137/0603006">Goldfarb and Idnani (1983)</a>
 * @since 1.11
 */
public class QPDualActiveSolverR extends QPOptimizer {

    /** No-error code. */
    public static final double ERROR_NONE = 0.0;

    /** Error code reported in lambda[0] when the Hessian decomposition fails. */
    public static final double ERROR_CHOLESKY_DECOMPOSITION = -1.0;

    /** Error code reported in lambda[0] when equality constraints are dependent. */
    public static final double ERROR_DEPENDENT_EQUALITIES = -2.0;

    /** Error code reported in lambda[0] when max iterations are reached. */
    public static final double ERROR_MAX_ITERATIONS = -3.0;

    /** Error code reported when numerical issue. */
    public static final double ERROR_NUMERICAL = -4.0;

    /** Error no Solution Found. */
    public static final double ERROR_INFEASIBLE = -5.0;

    /** Machine epsilon for tolerance checks. */
    private static final double EPS = Math.ulp(1.0);

    /**
     * Maximum number of iterations allowed.Will be adjusted in base of problem
     * dimension.
     */
    private int maxIter;

    /** Quadratic function representing 1/2 x^T G x + g0^T x. */
    private QuadraticFunction function;

    /** Equality constraint data (CE^T x + ce0 = 0). */
    private LinearEqualityConstraint eqConstraints;

    /** Inequality constraint data (CI^T x + ci0 >= 0). */
    private LinearInequalityConstraint iqConstraints;

    /** Bounded constraint data (lower <= Ax <= upper). */
    private LinearBoundedConstraint bConstraints;

    /**
     * Tolerance for symmetric matrix decomposition.
     *
     * @since 4.1
     */
    private MatrixDecompositionTolerance matrixDecompositionTolerance;

    /**
     * If true, function.getP() is interpreted as Cholesky lower factor L
     * (H = L*L^T), not as the Hessian H itself.
     */
    private boolean isCholesky = false;

    /**
     * Parses optimization data to extract the objective function and various
     * constraint sets.
     *
     * @param optData optimization data
     */
    @Override
    protected void parseOptimizationData(final OptimizationData... optData) {
        super.parseOptimizationData(optData);
        // reset QP problem to reuse the same instance of the QP solver;
        this.maxIter = 1000;
        this.function = null;
        this.eqConstraints = null;
        this.iqConstraints = null;
        this.bConstraints = null;
        this.matrixDecompositionTolerance = new MatrixDecompositionTolerance(EPS);
        this.isCholesky = false;

        for (OptimizationData data : optData) {
            if (data instanceof ObjectiveFunction) {
                function = (QuadraticFunction) ((ObjectiveFunction) data).getObjectiveFunction();
            } else if (data instanceof LinearEqualityConstraint) {
                eqConstraints = (LinearEqualityConstraint) data;
            } else if (data instanceof LinearInequalityConstraint) {
                iqConstraints = (LinearInequalityConstraint) data;
            } else if (data instanceof LinearBoundedConstraint) {
                bConstraints = (LinearBoundedConstraint) data;
            } else if (data instanceof IsCholesky) {
                isCholesky = ((IsCholesky) data).isCholesky();
            } else if (data instanceof MatrixDecompositionTolerance) {
                matrixDecompositionTolerance = (MatrixDecompositionTolerance) data;
            }
        }
    }

    /**
     * Finds step for a primal move.
     *
     * @param z        search direction
     * @param ai       constraint row vector
     * @param sv       current violation offset
     * @param equality true if constraint is equality
     * @return the primal step size
     */
    private double findPrimalStep(final RealVector z,
                                  final RealVector ai,
                                  final double sv,
                                  final boolean equality) {
        double norm2 = z.dotProduct(z);
        if (FastMath.abs(norm2) < Precision.EPSILON) {
            return Double.POSITIVE_INFINITY;
        }

        double denom = ai.dotProduct(z);
        if (FastMath.abs(denom) < Precision.SAFE_MIN) {
            return Double.POSITIVE_INFINITY;
        }

        double alpha = -sv / denom;
        // step for inequality should be positive
        if (!equality && alpha < 0) {
            alpha = Double.POSITIVE_INFINITY;
        }

        return alpha;
    }

    /**
     * Finds the blocking step for a dual move.
     *
     * @param u         current multipliers
     * @param r         dual direction
     * @param activeSet map of active indices -> constraint ids
     * @param me        threshold index for equality
     * @return the blocking step and index
     */
    private Pair<Integer, Double> findDualBlockingConstraint(final RealVector u,
                                                             final RealVector r,
                                                             final List<Integer> activeSet,
                                                             final int me) {
        if (activeSet.isEmpty()) {
            return new Pair<>(-1, Double.POSITIVE_INFINITY);
        }

        double alpha = Double.POSITIVE_INFINITY;
        int block = -1;
        int size = activeSet.size();

        for (int i = 0; i < size; i++) {
            double ui = u.getEntry(i);
            double ri = r.getEntry(i);
            // consider only inequality constraint in the active set
            if (ri > Precision.SAFE_MIN && activeSet.get(i) >= me) {
                double cand = ui / ri;
                if (cand < alpha) {
                    alpha = cand;
                    block = i;
                }
            }
        }
        
        return new Pair<>(block, alpha);
    }

    /**
     * Updates multipliers when adding a constraint.
     *
     * @param u       current u or null
     * @param r       dual direction
     * @param alpha   step length
     * @param partial new partial multiplier
     * @return updated multipliers
     */
    private RealVector updateMultipliersOnAddition(final RealVector u,
                                                   final RealVector r,
                                                   final double alpha,
                                                   final double partial) {
        if (u.getDimension() == 0) {
            RealVector v = new ArrayRealVector(1);
            v.set(partial);
            return v;
        }
        return u.add(r.mapMultiply(-alpha)).append(partial);
    }

    /**
     * Updates multipliers when removing a constraint.
     *
     * @param u         current multipliers
     * @param r         dual direction
     * @param alpha     step length
     * @param dropIndex index to remove
     * @return updated multipliers
     */
    private RealVector updateMultipliersOnRemoval(final RealVector u,
                                                  final RealVector r,
                                                  final double alpha,
                                                  final int dropIndex) {
        if (u.getDimension() == 1) {
            return new ArrayRealVector(0, 0);
        }

        RealVector tmp = u.add(r.mapMultiply(-alpha));
        int size = tmp.getDimension();

        if (dropIndex == 0) {
            return tmp.getSubVector(1, size - 1);
        } else if (dropIndex == size - 1) {
            return tmp.getSubVector(0, size - 1);
        }

        RealVector head = tmp.getSubVector(0, dropIndex);
        RealVector tail = tmp.getSubVector(dropIndex + 1, size - dropIndex - 1);
        return head.append(tail);
    }
    
   

    private RealVector removeMultiplierAt(final RealVector u, final int dropIndex) {
        if (u.getDimension() == 1) {
            return new ArrayRealVector(0, 0);
        }
        if (dropIndex == 0) {
            return u.getSubVector(1, u.getDimension() - 1);
        }
        if (dropIndex == u.getDimension() - 1) {
            return u.getSubVector(0, u.getDimension() - 1);
        }
        final RealVector head = u.getSubVector(0, dropIndex);
        final RealVector tail = u.getSubVector(dropIndex + 1, u.getDimension() - dropIndex - 1);
        return head.append(tail);
    }

    private int findEqualitySwapIndex(final RealVector u,
                                      final RealVector r,
                                      final List<Integer> active,
                                      final int p) {
        int swapIndex = -1;
        double minRatio = Double.POSITIVE_INFINITY;

        for (int j = 0; j < active.size(); j++) {
            if (active.get(j) >= p) {
                final double rj = r.getEntry(j);
                if (rj > Precision.SAFE_MIN) {
                    final double ratio = u.getEntry(j) / rj;
                    if (ratio < minRatio) {
                        minRatio = ratio;
                        swapIndex = j;
                    }
                }
            }
        }

        return swapIndex;
    }

    /**
     * Computes constraint weights used for the normalized violation criterion.
     * Weight is 1 / ||a_i|| when possible.
     *
     * @param C full constraint matrix
     * @return weights
     */
    private RealVector computeConstraintWeights(final RealMatrix C) {
        final int m = C.getColumnDimension();
        final RealVector weights = new ArrayRealVector(m);
        for (int i = 0; i < m; i++) {
            final double norm = C.getColumnVector(i).getNorm();
            weights.setEntry(i, norm > Precision.SAFE_MIN ? 1.0 / norm : 1.0);
        }
        return weights;
    }

    /**
     * Find the most violated constraint using normalized violations.
     *
     * Equality constraints use -|c_i(x)| * w_i.
     * Inequality/bound constraints use c_i(x) * w_i.
     *
     * This keeps the original solver convention: a violated constraint must
     * produce a negative value to be selected.
     *
     * @param sv        current constraint evalutation
     * @param weights   constraint weights
     * @param blackList dependent constraint
     * @param activeSet constraint in active set
     * @param me        equality constraint numbers
     * @return the violation value and index
     */
    private Pair<Integer, Double> mostViolatedConstraintNormalized(final RealVector sv,
                                                                  final RealVector weights,
                                                                  final Set<Integer> blackList,
                                                                  final List<Integer> activeSet,
                                                                  final int me) {
        double maxViolation = 0.0;
        int mostViolated = -1;

        for (int i = 0; i < sv.getDimension(); i++) {
            if (blackList.contains(i) || activeSet.contains(i)) {
                continue;
            }

            final double violation;
            if (i < me) {
                violation = -FastMath.abs(sv.getEntry(i)) * weights.getEntry(i);
            } else {
                violation = sv.getEntry(i) * weights.getEntry(i);
            }

            if (violation < maxViolation) {
                maxViolation = violation;
                mostViolated = i;
            }
        }

        return new Pair<>(mostViolated, maxViolation);
    }

    /**
     * Main optimization routine.
     *
     * @return the optimal solution or null if infeasible
     */
    @Override
    public LagrangeSolution doOptimize() {
        RealMatrix G;

        if (!this.isCholesky) {
            G = function.getP();
        } else {
            G = function.getP().multiplyTransposed(function.getP());
        }

        RealVector g0 = function.getQ();
        double g = function.getD();
        int n = G.getColumnDimension();

        RealMatrix CE = null;
        RealVector ce0 = null;
        int p = 0;
        if (eqConstraints != null && eqConstraints.getA().getRowDimension() > 0) {
            CE = eqConstraints.getA().transpose();
            ce0 = eqConstraints.getLowerBound().mapMultiply(-1.0);
            p = CE.getColumnDimension();
        }

        int m1 = (iqConstraints != null) ? iqConstraints.getLowerBound().getDimension() : 0;
        int b1 = (bConstraints != null) ? bConstraints.getLowerBound().getDimension() : 0;
        int m = m1 + 2 * b1;

        final int mc = p + m;
        RealMatrix C = null;
        RealVector c0 = null;
        RealVector weights = null;

        if (mc > 0) {
            C = MatrixUtils.createRealMatrix(n, mc);
            c0 = new ArrayRealVector(mc);

            if (p > 0) {
                C.setSubMatrix(CE.getData(), 0, 0);
                c0.setSubVector(0, ce0);
            }

            if (m1 > 0) {
                RealMatrix Aineq = iqConstraints.jacobian(null);
                RealVector bineq = iqConstraints.getLowerBound();
                C.setSubMatrix(Aineq.transpose().getData(), 0, p);
                c0.setSubVector(p, bineq.mapMultiply(-1.0));
            }

            if (b1 > 0) {
                RealMatrix Abound = bConstraints.jacobian(null);
                RealVector lower = bConstraints.getLowerBound();
                RealVector upper = bConstraints.getUpperBound();
                C.setSubMatrix(Abound.transpose().getData(), 0, p + m1);
                C.setSubMatrix(Abound.scalarMultiply(-1.0).transpose().getData(), 0, p + m1 + b1);
                c0.setSubVector(p + m1, lower.mapMultiply(-1.0));
                c0.setSubVector(p + m1 + b1, upper);
            }

            weights = computeConstraintWeights(C);
        }

        RealVector x;
        RealMatrix L;
        RealMatrix L1;
        QRUpdaterR qrUpdater;
        double tol;
        if (!this.isCholesky) {
            try {
                final double eps = matrixDecompositionTolerance.getEpsMatrixDecomposition();
                final CholeskyDecomposition cholesky = new CholeskyDecomposition(G, eps, eps);
                DecompositionSolver solver = cholesky.getSolver();

                x = solver.solve(g0).mapMultiply(-1.0);
                L = cholesky.getL();
//                if(cholesky.getDiagonalShift()!=0.0)
//                    G=L.multiplyTransposed(L);
//                L=factorizeLQp(G,null);

                L1 = inverseLowerTriangular(L);

//                x = L1.preMultiply(L1.operate(g0)).mapMultiply(-1.0);
//                G=L.multiplyTransposed(L);
//                c1 trace of G matrix
//                double c1 = FastMath.sqrt(G.getTrace());
                double c1 = G.getTrace();
//                c2 trace of inverse of cholesky factorization
//                double c2 = FastMath.sqrt(L1.getTrace());
                double c2 = L1.getTrace();
                tol = mc * c1 * c2 * Precision.EPSILON * 100.0;
                qrUpdater = new QRUpdaterR(L1);
            } catch (MathIllegalArgumentException ex) {
                // matrix is not positive definite return empty solution
                return buildFailureSolution(ERROR_DEPENDENT_EQUALITIES);
            }
        } else {
            L = function.getP();
            L1 = inverseLowerTriangular(L);

            x = L1.preMultiply(L1.operate(g0)).mapMultiply(-1.0);
            double c1 = G.getTrace();
            double c2 = L1.getTrace();
            tol = mc * c1 * c2 * Precision.EPSILON * 100.0;
            qrUpdater = new QRUpdaterR(L1);
        }

        if (mc == 0) {
            return new LagrangeSolution(
                    x,
                    new ArrayRealVector(0, 0),
                    0.5 * x.dotProduct(G.operate(x)) + g0.dotProduct(x) + g
            );
        }

        // max iteration adjusted in base of problem dimension
        this.maxIter = 40 * (n + mc);

        // convergence threshold calculated in base at the matrix conditioning
        // ActiveSet and blackLit(dependent constraints)

        final Set<Integer> blacklist = new HashSet<>();
        List<Integer> active = new ArrayList<>();

        RealVector u = new ArrayRealVector(0, 0);
        RealVector r = new ArrayRealVector(0, 0);
        RealVector d = null;
        RealVector z = null;
        int iteration = 0;

        // Active-set loop for all constraints
        while (mc != 0 && iteration++ < maxIter) {

            RealVector sv;
            // store solution in case constraint can't be added because dependent
            RealVector xOld = x.copy();
            RealVector uOld = u.copy();
            // evaluate constraints
            sv = C.transpose().operate(x).add(c0);

//            // convergence test:
//            // consider only NON-active constraints.
//            // For equalities use absolute normalized violation.
//            // For inequalities/bounds use one-sided normalized violation.
//            double cvMax = 0.0;
//
//            for (int k = 0; k < sv.getDimension(); k++) {
//                if (active.contains(k)) {
//                    continue;
//                }
//
//                final double scaledViolation;
//                if (k < p) {
//                    scaledViolation = FastMath.abs(sv.getEntry(k)) * weights.getEntry(k);
//                } else {
//                    scaledViolation = -sv.getEntry(k) * weights.getEntry(k);
//                }
//
//                if (scaledViolation > cvMax) {
//                    cvMax = scaledViolation;
//                }
//            }
//
//            // Evaluate convergence : no significantly violated NON-active constraint
//            if (cvMax <= Precision.EPSILON) {
//                break; // Optimal solution found
//            }

             // convergence test:
            // consider only NON-active constraints.
            // For equalities use absolute normalized violation.
            // For inequalities/bounds use one-sided normalized violation.
            double cvMax = 0.0;

            for (int k = 0; k < sv.getDimension(); k++) {
                if (active.contains(k)) {
                    continue;
                }

                double residual = sv.getEntry(k);
                double absResidual = FastMath.abs(residual);

                // --- INIZIO FILTRO ANTIRUMORE DI POWELL ---
                // Calcola la somma dei valori assoluti degli operandi (il "TEMP" di Powell)
                // c0 è l'equivalente di B(K), C.getEntry(i, k) è l'equivalente di A(K,I)
                double absSum = FastMath.abs(c0.getEntry(k));
                for (int i = 0; i < x.getDimension(); i++) {
                    // Ricorda: C ha le normali sulle colonne, quindi usiamo C.getEntry(i, k)
                    absSum += FastMath.abs(C.getEntry(i, k) * x.getEntry(i));
                }

                double tempA = absSum + absResidual;
                
                // Test 1: Se sommare il residuo ad absSum non cambia il risultato 
                // in precisione floating-point, il residuo è puro rumore. Ignoralo.
                if (tempA <= absSum) {
                    continue; 
                }
                
                // Test 2: Margine di sicurezza di Powell (1.5 * residuo)
                double tempB = absSum + 1.5 * absResidual;
                if (tempB <= tempA) {
                    continue;
                }
                // --- FINE FILTRO ANTIRUMORE ---

                final double scaledViolation;
                if (k < p) {
                    // Per le uguaglianze, usiamo il valore assoluto già calcolato
                    scaledViolation = absResidual * weights.getEntry(k);
                } else {
                    // Per le disuguaglianze, la violazione è negativa nel tuo framework
                    scaledViolation = -residual * weights.getEntry(k);
                }

                if (scaledViolation > cvMax) {
                    cvMax = scaledViolation;
                }
            }

            // Evaluate convergence : no significantly violated NON-active constraint
            if (cvMax <= Precision.EPSILON) {
                break; // Optimal solution found
            }

////////            sv = C.transpose().operate(x).add(c0);
////////
////////            //calculate norm1-like measure of scaled violations with original flow
////////            double sum = 0;
////////            for (int k = 0; k < sv.getDimension(); k++) {
////////                if (k < p) {
////////                    sum -= FastMath.abs(sv.getEntry(k)) * weights.getEntry(k);
////////                } else {
////////                    sum += FastMath.min(0.0, sv.getEntry(k) * weights.getEntry(k));
////////                }
////////            }
////////
////////            // Evaluate convergence
////////            if (FastMath.abs(sum) <= tol) {
////////                break;// Optimal solution found
////////            }

            // Evaluate most violated constraint, excluding dependent/active loop
            while (iteration++ < maxIter) {

                final Pair<Integer, Double> mostViolated =
                        mostViolatedConstraintNormalized(sv, weights, blacklist, active, p);
                if (mostViolated.getValue() >= 0) {
                    blacklist.clear();
                    break; // reavaluate constraints and optimal condition;
                }

                double t1;
                double t2;
                double t = 0;
                double uPartial = 0;
                int dropIndex;
                RealVector np;
                RealMatrix J2;
                final int constraintIndex = mostViolated.getKey();
                final boolean equality = constraintIndex < p;

                // Dual step loop update multiplier and x (if step is also in primal)
                // until primal step is not done
                while (iteration++ < maxIter) {
                    np = C.getColumnVector(constraintIndex);
                    sv.setEntry(constraintIndex, np.dotProduct(x) + c0.getEntry(constraintIndex));
                    d = qrUpdater.computeD(np);
                    z = qrUpdater.computeZ(d);
                    r = qrUpdater.solveR(d);

                    t1 = findPrimalStep(z, np, sv.getEntry(constraintIndex), equality);
                    Pair<Integer, Double> dualStep = findDualBlockingConstraint(u, r, active, p);
                    t2 = dualStep.getValue();
                    dropIndex = dualStep.getKey();
                    t = FastMath.min(t1, t2);
                    if (!Double.isFinite(t)) {
                        return buildFailureSolution(ERROR_INFEASIBLE);
                    }
                    if (t == t1) {
                        break; // primal full step (exit from dual step loop)
                    } else {
                        // Manage dual step
                        if (t1 < Double.POSITIVE_INFINITY) {
                            // step is also in primal
                            x = x.add(z.mapMultiply(t));
                        }
                        uPartial += t;
                        u = updateMultipliersOnRemoval(u, r, t, dropIndex);
                        
                        qrUpdater.deleteConstraint(dropIndex);

                        active.remove(dropIndex);
                    }
                }

                // Manage full step
                //equality constraint:
                //equality constraint need to be added and if it is no possible some inequality constraint on the active set
                //will be dropped
                if (constraintIndex < p) {
                    // equality case

                    boolean added = false;
                    if (active.size() < n) {
                        added = qrUpdater.addConstraint(d);
                    }

                    if (added) {
                        active.add(constraintIndex);
                        x = x.add(z.mapMultiply(t));
                        uPartial += t;
                        u = updateMultipliersOnAddition(u, r, t, uPartial);
                        
                        blacklist.clear();
                        break; // re-evaluate convergence
                    } else {
                        final int swapIndex = findEqualitySwapIndex(u, r, active, p);
                        
                        if (swapIndex < 0) {
                            return buildFailureSolution(ERROR_DEPENDENT_EQUALITIES);
                        }

                        final int droppedConstraint = active.get(swapIndex);
                        System.out.println("SWAP CONSTRAINT FOUND");
                        qrUpdater.deleteConstraint(swapIndex);
                        active.remove(swapIndex);
                        u = removeMultiplierAt(u, swapIndex);

                        blacklist.add(droppedConstraint);

                        // let the equality be selected again in the next pass
                        break;
                    }
                } else {
                    // inequality / bound case
                    if (active.size() < n && qrUpdater.addConstraint(d)) {
                        active.add(constraintIndex);
                        x = x.add(z.mapMultiply(t));
                        uPartial += t;
                        u = updateMultipliersOnAddition(u, r, t, uPartial);
                        blacklist.clear();
                        break;
                    } else {
                        blacklist.add(constraintIndex);
                    }
                }
                //////MANAGE FULL STEP END
            }
        }

        if (iteration == maxIter) {
            return buildFailureSolution(ERROR_MAX_ITERATIONS);
        }

        return buildSolution(x, u, active, G, g0, g, p, m);
    }

    /**
     * Builds the final solution object with primal x and multipliers lambda.
     *
     * @param x         solution
     * @param u         active set multipliers
     * @param activeSet active set constraints
     * @param G         square matrix of weights for quadratic terms
     * @param g0        vector of weights for linear terms
     * @param g         constant term
     * @param p         number of equalities
     * @param m         number of inequalities
     * @return the optimal solution
     */
    private LagrangeSolution buildSolution(final RealVector x,
                                           final RealVector u,
                                           final List<Integer> activeSet,
                                           final RealMatrix G,
                                           final RealVector g0,
                                           final double g,
                                           final int p,
                                           final int m) {
        final RealVector lambda = new ArrayRealVector(p + m);
        if (!activeSet.isEmpty()) {
            for (int i = 0; i < activeSet.size(); i++) {
                lambda.setEntry(activeSet.get(i), u.getEntry(i));
                if (activeSet.get(i) >= p && u.getEntry(i) < 0) {
                    lambda.setEntry(activeSet.get(i), 0.0);
                }
            }
        }
        final double value = 0.5 * x.dotProduct(G.operate(x)) +
                             g0.dotProduct(x) + g;
        return new LagrangeSolution(x, lambda, value);
    }

    /**
     * Computes the inverse of a lower-triangular matrix via forward
     * substitution.
     *
     * @param L lower triangular matrix
     * @return inverse of Lower Triangular Matrix
     */
    private RealMatrix inverseLowerTriangular(final RealMatrix L) {
        final int n = L.getRowDimension();
        final RealMatrix Linv = MatrixUtils.createRealMatrix(n, n);
        for (int i = 0; i < n; i++) {
            final RealVector e = new ArrayRealVector(n);
            e.setEntry(i, 1.0);
            MatrixUtils.solveLowerTriangularSystem(L, e);
            Linv.setColumnVector(i, e);
        }
        return Linv;
    }

    private LagrangeSolution buildFailureSolution(final double ERROR_CHOLESKY_DECOMPOSITION) {
        return new LagrangeSolution(
                new ArrayRealVector(0, 0),
                new ArrayRealVector(1, ERROR_CHOLESKY_DECOMPOSITION),
                0
        );
    }

}