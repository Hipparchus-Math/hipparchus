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

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.ArrayRealVector;
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
 * CI^T x  &gt;= ci0     (inequality constraints)
 * lb &lt;= Ax &lt;= ub      (bounded constraints)
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
public class QPDualActiveSolver extends QPOptimizer {

    /** No-error code. */
    public static final double ERROR_NONE = 0.0;

    /** Error code reported in lambda[0] when the Hessian decomposition fails. */
    public static final double ERROR_CHOLESKY_DECOMPOSITION = -1.0;

    /** Error code reported in lambda[0] when equality constraints are dependent. */
    public static final double ERROR_DEPENDENT_EQUALITIES = -2.0;

    /** Error code reported in lambda[0] when max iterations are reached. */
    public static final double ERROR_MAX_ITERATIONS = -3.0;

    /** Error code reported when numerical issue occurs. */
    public static final double ERROR_NUMERICAL = -4.0;

    /** Error code reported when no solution is found (infeasible). */
    public static final double ERROR_INFEASIBLE = -5.0;

    /** Machine epsilon for tolerance checks. */
    private static final double EPS = FastMath.ulp(1.0);

    /**
     * Maximum number of iterations allowed. Will be adjusted based on the problem
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

    /** QP Settings */
    private QPDualActiveSolverOption settings;

    /**
     * Tolerance for symmetric matrix decomposition.
     *
     * @since 4.1
     */
    private MatrixDecompositionTolerance matrixDecompositionTolerance;

    /** Simple constructor.
     */
    public QPDualActiveSolver() {
        // nothing to do
    }

    /**
     * Parses optimization data to extract the objective function and various
     * constraint sets.
     *
     * @param optData optimization data
     */
    @Override
    protected void parseOptimizationData(final OptimizationData... optData) {
        super.parseOptimizationData(optData);
        // Reset QP problem to reuse the same instance of the QP solver
        this.maxIter = 1000;
        this.function = null;
        this.eqConstraints = null;
        this.iqConstraints = null;
        this.bConstraints = null;
        this.settings=new QPDualActiveSolverOption();
        this.matrixDecompositionTolerance = new MatrixDecompositionTolerance(EPS);


        for (OptimizationData data : optData) {
            if (data instanceof ObjectiveFunction) {
                function = (QuadraticFunction) ((ObjectiveFunction) data).getObjectiveFunction();
            } else if (data instanceof LinearEqualityConstraint) {
                eqConstraints = (LinearEqualityConstraint) data;
            } else if (data instanceof LinearInequalityConstraint) {
                iqConstraints = (LinearInequalityConstraint) data;
            } else if (data instanceof LinearBoundedConstraint) {
                bConstraints = (LinearBoundedConstraint) data;
            } else if (data instanceof MatrixDecompositionTolerance) {
                matrixDecompositionTolerance = (MatrixDecompositionTolerance) data;
            }
              else if(data instanceof QPDualActiveSolverOption)
            {
               settings = (QPDualActiveSolverOption) data;
            }
        }
    }



    /**
     * Finds the step length for a primal move toward a violated constraint.
     * Modern optimized and numerically stable version.
     *
     * @param z        search direction vector
     * @param ai       constraint row vector
     * @param sv       current violation offset (slack/surplus value)
     * @param equality true if the constraint is an equality constraint
     * @return the primal step size, or Double.POSITIVE_INFINITY if the constraint does not block
     */
    private double findPrimalStep(final RealVector z,
                              final RealVector ai,
                              final double sv,
                              final boolean equality) {

    double sumA = 0.0;
    double sumB = 0.0;
    double aiNormSq = 0.0;
    double zNormSq = 0.0;

    final int dimension = ai.getDimension();

    for (int i = 0; i < dimension; i++) {
        final double aiValue = ai.getEntry(i);
        final double zValue = z.getEntry(i);
        final double product = aiValue * zValue;

        sumA += product;
        sumB += FastMath.abs(product);
        aiNormSq += aiValue * aiValue;
        zNormSq += zValue * zValue;
    }

    if (sumB <= Precision.EPSILON) {
        return Double.POSITIVE_INFINITY;
    }

    final double absSumA = FastMath.abs(sumA);

    double temp = sumB + 0.1 * absSumA;
    double tempA = sumB + 0.2 * absSumA;

    if (temp <= sumB || tempA <= temp) {
        return Double.POSITIVE_INFINITY;
    }

    final double sumC =
            FastMath.sqrt(aiNormSq) * FastMath.sqrt(zNormSq);

    temp = sumC + 0.1 * absSumA;
    tempA = sumC + 0.2 * absSumA;

    if (temp <= sumC || tempA <= temp) {
        return Double.POSITIVE_INFINITY;
    }

    final double alpha = -sv / sumA;

    if (!equality && alpha < 0.0) {
        return Double.POSITIVE_INFINITY;
    }

    return alpha;
}




/**
 * Finds the blocking active inequality for the dual move.
 *
 * If reverseStep is false, the algorithm moves in the usual positive
 * direction and an active inequality can block only if r_i > 0.
 *
 * If reverseStep is true, the algorithm moves in the opposite direction
 * and an active inequality can block only if r_i < 0.
 *
 * The returned step is always a positive length:
 *
 *     alpha = u_i / |r_i|
 *
 * Active equalities are never removed.
 *
 * @param u           current active multipliers
 * @param r           dual direction
 * @param activeSet   active-set global constraint indices
 * @param me          number of equality constraints
 * @param reverseStep true if the primal/dual step is taken in the negative direction
 * @return pair containing active-set position to drop and positive blocking length;
 *         if no blocker exists, returns (-1, +Infinity)
 */
private Pair<Integer, Double> findDualBlockingConstraint(final RealVector u,
                                                          final RealVector r,
                                                          final List<Integer> activeSet,
                                                          final int me,
                                                          final boolean reverseStep) {
    if (activeSet.isEmpty()) {
        return new Pair<>(-1, Double.POSITIVE_INFINITY);
    }

    double alpha = Double.POSITIVE_INFINITY;
    int block = -1;

    for (int i = 0; i < activeSet.size(); i++) {

        /*
         * Only active inequalities can leave.
         * Active equalities must never be dropped.
         */
        if (activeSet.get(i) < me) {
            continue;
        }

        final double ri = r.getEntry(i);
        final double absRi=FastMath.abs(ri);

        if (!reverseStep) {

            /*
             * Forward direction:
             *
             *     u_i(alpha) = u_i - alpha * r_i
             *
             * Only r_i > 0 can drive u_i to zero.
             */
            if (ri <= 0.0) {
                continue;
            }

        } else {

            /*
             * Reverse direction:
             *
             *     signedStep = -alpha
             *     u_i(new) = u_i - signedStep * r_i
             *              = u_i + alpha * r_i
             *
             * Only r_i < 0 can drive u_i to zero.
             */
            if (ri >= 0.0) {
                continue;
            }
        }


          final double cand = u.getEntry(i) / FastMath.abs(ri);





        if (cand < alpha) {
            alpha = cand;
            block = i;
        }
    }

    return new Pair<>(block, alpha);
}

    /**
     * Updates multipliers when adding a constraint.
     *
     * @param u       current multipliers or null
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
            RealVector v = new ArrayRealVector(1,partial);
            return v;
        }

        // Original: return u.add(r.mapMultiply(-alpha)).append(partial);
        // Optimized: copy once, combine in-place, then append

        u.combineToSelf(1.0, -alpha, r);

        return u.append(partial);

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

        // Original: RealVector tmp = u.add(r.mapMultiply(-alpha));
        // Optimized: copy once, combine in-place

        u.combineToSelf(1.0, -alpha, r);

        int size = u.getDimension();

        if (dropIndex == 0) {
            return u.getSubVector(1, size - 1);
        } else if (dropIndex == size - 1) {
            return u.getSubVector(0, size - 1);
        }

        RealVector head = u.getSubVector(0, dropIndex);
        RealVector tail = u.getSubVector(dropIndex + 1, size - dropIndex - 1);
        RealVector result= head.append(tail);
        return result;
    }



    /**
     * Removes a multiplier at a specific index.
     *
     * @param u         the multiplier vector
     * @param dropIndex the index to drop
     * @return the vector without the dropped multiplier
     */
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
            weights.setEntry(i, norm > Precision.EPSILON ? 1.0 / norm : 0.0);
        }
        return weights;
    }

    /**
//     * Finds the most violated constraint using normalized violations and an embedded numerical noise filter.
//     * <p>
//     * Equality constraints use absolute normalized violation.
//     * Inequality/bound constraints use one-sided normalized violation.
//     * The noise filter evaluates the exact numerical floating-point stability of the constraint residual
//     * relative to the magnitude of the terms involved in its calculation.
//     * </p>
//     *
//     * @param sv        current constraint evaluation (residuals)
//     * @param weights   constraint weights
//     * @param blackList dependent constraints to ignore
//     * @param activeSet constraints currently in the active set
//     * @param me        number of equality constraints
//     * @param C         the constraint matrix
//     * @param c0        the constraint constant vector
//     * @param x         the current primal point
//     * @return a Pair containing the index of the most violated constraint and the maximum violation value
//     */
     private Pair<Integer, Double> mostViolatedConstraintNormalized(final RealVector sv,
                                                                   final RealVector weights,
                                                                   final boolean[] active ,
                                                                   final int me,
                                                                   final RealMatrix C,
                                                                   final RealVector c0,
                                                                   final RealVector x) {
        double cvMax = 0.0;
        int kViolated = -1;

        for (int k = 0; k < sv.getDimension(); k++) {
            if (active[k] ) {
                continue;
            }





            final double residual = sv.getEntry(k);
            final double absResidual = FastMath.abs(residual);

            // 1. Scaled Violation
            final double scaledViolation;
            if (k < me) {
                scaledViolation = absResidual * weights.getEntry(k);
            } else {
                scaledViolation = -residual * weights.getEntry(k);
            }

            // 2. Convergence Test: skip if not greater than current maximum
            if (scaledViolation <= cvMax) {
                continue;
            }


            final RealVector rowC=C.getRowVector(k);
            // 3. --- NOISE FILTER ---
            double absSum = FastMath.abs(c0.getEntry(k));
            for (int i = 0; i < x.getDimension(); i++) {
                absSum += FastMath.abs(rowC.getEntry(i) * x.getEntry(i));
            }

            double tempA = absSum + absResidual;
            if (tempA <= absSum) {
                continue; // Noise detected
            }

            double tempB = absSum + 1.5 * absResidual;
            if (tempB <= tempA) {
                continue; // Noise detected
            }

            kViolated = k;
            cvMax = scaledViolation;
        }

        return new Pair<>(kViolated, cvMax);
    }




    /**
     * Main optimization routine.
     *
     * @return the optimal solution or null if infeasible
     */
    @Override
    public LagrangeSolution doOptimize() {



        RealVector g0 = function.getQ();
        double g = function.getD();
        int n = function.getP().getColumnDimension();


        int p =  (eqConstraints != null) ? eqConstraints.dimY():0;
        int m1 = (iqConstraints != null) ? iqConstraints.dimY() : 0;
        int b1 = (bConstraints != null)  ? bConstraints.dimY() : 0;
        int m = m1 + 2 * b1;

        final int mc = p + m;
        RealMatrix C = null;
        RealVector c0 = null;
        RealVector weights = null;

        if (mc > 0) {

            QPConstraintSystem record = buildGlobalConstraints(n);
             C=record.C;
             c0=record.c0;
             weights = record.weights;
        }
        RealMatrix G=null;
        RealVector x=null;
        RealMatrix L=null;
        RealMatrix L1=null;
        QRUpdater qrUpdater=null;

        if (settings.getMatrixMode() == QPMatrixMode.FULL) {
            G = function.getP();
            try {
                final double eps = matrixDecompositionTolerance.getEpsMatrixDecomposition();
                final RegularizedCholeskyDecomposition cholesky = new RegularizedCholeskyDecomposition(G, settings.getEps(), settings.getEps());
                DecompositionSolver solver = cholesky.getSolver();

                x = solver.solve(g0).mapMultiply(-1.0);
                L = cholesky.getL();

                L1 = inverseLowerTriangular(L);
                qrUpdater = new QRUpdater(L1);
            } catch (MathIllegalArgumentException ex) {
                // Matrix is not positive definite, return failure solution
                return buildFailureSolution(ERROR_CHOLESKY_DECOMPOSITION);
            }
        } else if(settings.getMatrixMode() ==QPMatrixMode.CHOLESKY) {

            L = function.getP();
            G = L.multiplyTransposed(L);
            L1 = inverseLowerTriangular(L);

            x = L1.preMultiply(L1.operate(g0)).mapMultiply(-1.0);
            qrUpdater = new QRUpdater(L1);
        }
        else  {

            L1 = function.getP();
            L = inverseLowerTriangular(L1);
            G = L.multiplyTransposed(L);


            x = L1.preMultiply(L1.operate(g0)).mapMultiply(-1.0);
            qrUpdater = new QRUpdater(L1);
        }


        if (mc == 0) {
            double obj = 0.5 * x.dotProduct(G.operate(x)) + g0.dotProduct(x) + g;
            return new LagrangeSolution(x, new ArrayRealVector(0, 0), obj);
        }


        // Maximum iterations adjusted based on problem dimension
        this.maxIter = 40 * (n + mc);

        List<Integer> active = new ArrayList<>();
        final boolean[] activeMask = new boolean[mc];
        RealVector u = new ArrayRealVector(0, 0);
        RealVector r = new ArrayRealVector(0, 0);
        RealVector d = null;
        RealVector z = null;
        int iteration = 0;
        boolean refinementDone = false;
        // Active-set loop for all constraints
        while (mc != 0 && iteration++ < maxIter) {

            // Check Dual Feasibility delete constraint with negative multiplier and recalculate solution
                boolean dualFeasible = false;

                // Continua finché non è dual-fattibile o non svuotiamo il set
                while (!dualFeasible && !active.isEmpty()) {
                    dualFeasible = true;
                    int dropIndex = -1;


//                    // 1. Identifica il primo moltiplicatore negativo puro
                    for (int j = 0; j < active.size(); j++) {
                        if (active.get(j) >= p) { // Se è una disuguaglianza
                            if (u.getEntry(j) < 0.0) {
                                dropIndex = j;
                                dualFeasible = false;
                                break; // Esce al primo negativo incontrato
                            }
                        }
                    }

                    // 2. Se trovato, effettua il drop e riesegui il KKT da capo
                    if (dropIndex >= 0) {
                        qrUpdater.deleteConstraint(dropIndex);
                        final int removedConstraint = active.remove(dropIndex);
                        activeMask[removedConstraint] = false;
                        u=removeMultiplierAt(u, dropIndex);
                        Pair<RealVector, RealVector> xy = recomputeXY(x, G, g0, C, c0, qrUpdater, active);
                        x = xy.getFirst();
                        u = xy.getSecond();

                    }
                }


            RealVector sv;
            // Evaluate constraints sv = C.preMultiply(x).add(c0);
            sv = C.operate(x);
            sv.combineToSelf(1.0, 1.0, c0);
            // Convergence test: find the most violated NON-active constraint
            Pair<Integer, Double> violationInfo = mostViolatedConstraintNormalized(sv, weights, activeMask, p, C, c0, x);
            int kViolated = violationInfo.getKey();
            double cvMax = violationInfo.getValue();


             if (cvMax <= settings.getEps() && refinementDone) {
                break; // Solution found
            }

            // Evaluate convergence: no significantly violated NON-active constraint
            if (cvMax <= settings.getEps()  && !refinementDone) {

                kViolated = -1;
                refinementDone = true;

               Pair<RealVector, RealVector> xy = recomputeXY(x, G, g0, C, c0, qrUpdater, active);
                x = xy.getFirst();
                u = xy.getSecond();



            }



            // Active loop evaluation
            while (iteration++ < maxIter) {
                if (kViolated < 0) {
                    break;
                }

                double t1;
                double t2;
                double t = 0;
                double signedStep=0;
                double uPartial = 0;
                int dropIndex;
                RealVector np;

                final int constraintIndex = kViolated;
                final boolean equality = constraintIndex < p;

                // Dual step loop update multiplier and x (if step is also in primal)
                // until primal step is not done
                while (iteration++ < maxIter) {
                    np = C.getRowVector(constraintIndex);
                    sv.setEntry(constraintIndex, np.dotProduct(x) + c0.getEntry(constraintIndex));
                    d = qrUpdater.computeD(np);
                    z = qrUpdater.computeZ(d);
                    r = qrUpdater.solveR(d);

                    t1 = findPrimalStep(z, np, sv.getEntry(constraintIndex), equality);
                    final boolean reverseStep = equality && sv.getEntry(constraintIndex) > 0.0;
                    t1=FastMath.abs(t1);
                    Pair<Integer, Double> dualStep = findDualBlockingConstraint(u, r, active, p,reverseStep);
                    t2 = dualStep.getValue();
                    dropIndex = dualStep.getKey();
                    t = FastMath.min(t1, t2);
                    signedStep = reverseStep ? -t : t;

                    if (!Double.isFinite(t)) {
                        return cvMax < this.settings.getEpsRelaxation() ?
                               buildSolution(x, u, active, G, g0, g, p, m) :
                               buildFailureSolution(ERROR_INFEASIBLE);
                    }

                    if (t1 <= t2) {
                        break; // Primal full step (exit from dual step loop)
                    } else {
                        // Manage dual step
                        if (t1 < Double.POSITIVE_INFINITY) {
                            // Step is also in primal x = x.add(z.mapMultiply(t));
                            x.combineToSelf(1.0, signedStep, z);
                        }
                        uPartial += signedStep;
                        u = updateMultipliersOnRemoval(u, r, signedStep, dropIndex);

                        qrUpdater.deleteConstraint(dropIndex);
                        final int removedConstraint = active.remove(dropIndex);
                        activeMask[removedConstraint] = false;

                    }
                }


                // Manage full step
                if (active.size() < n && qrUpdater.addConstraint(d)) {
                    active.add(constraintIndex);
                    activeMask[constraintIndex] = true;
                    // Step rimal x = x.add(z.mapMultiply(t));
                    x.combineToSelf(1.0, signedStep, z);
                    uPartial += signedStep;
                    u = updateMultipliersOnAddition(u, r, signedStep, uPartial);

                    break; // Re-evaluate convergence
                } else {

                    return buildFailureSolution(ERROR_DEPENDENT_EQUALITIES);
                }
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
            }
        }
        final double value = 0.5 * x.dotProduct(G.operate(x)) + g0.dotProduct(x) + g;
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

    /**
     * Helper to build a failure solution structure.
     *
     * @param errorCode the failure code
     * @return a LagrangeSolution representing the failure
     */
    private LagrangeSolution buildFailureSolution(final double errorCode) {
        return new LagrangeSolution(
                new ArrayRealVector(0, 0),
                new ArrayRealVector(1, errorCode),
                0
        );
    }

    /**
     * Precomputes the entrywise absolute matrix of G.
     *
     * @param G the symmetric matrix
     * @return the absolute matrix
     */
    private RealMatrix precomputeAbsG(final RealMatrix G) {
        final int n = G.getRowDimension();
        final RealMatrix absG = MatrixUtils.createRealMatrix(n, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                absG.setEntry(i, j, FastMath.abs(G.getEntry(i, j)));
            }
        }

        return absG;
    }

    /**
     * Exact numerical progress check.
     *
     * <p>The test evaluates:
     *
     * <pre>
     * FDIFF  = sum_i ( 2*g0_i + sum_j G_ij*(xOld_j + xNew_j) ) * (xNew_i - xOld_i)
     * FDIFFA = sum_i ( |2*g0_i| + sum_j |G_ij*(xOld_j + xNew_j)| ) * |xNew_i - xOld_i|
     * </pre>
     *
     * and rejects progress if:
     *
     * <pre>
     * FDIFFA + FDIFF     &lt;= FDIFFA
     * FDIFFA + 1.5*FDIFF  &lt;= FDIFFA + FDIFF
     * </pre>
     *
     * <p>To save work:
     * <ul>
     * <li>FDIFF is evaluated exactly using the Cholesky factor {@code L}, with {@code G = L*L^T}</li>
     * <li>FDIFFA is evaluated using a precomputed entrywise absolute matrix {@code absG = |G|}</li>
     * </ul>
     *
     * @param xOld checkpoint point x_old
     * @param xNew current point x_new
     * @param g0 linear term of the quadratic objective
     * @param L lower Cholesky factor of G
     * @param absG precomputed entrywise absolute Hessian |G|
     * @return true if the objective progress is numerically reliable, false otherwise
     */
    private boolean isObjectiveProgressReliableExact(final RealVector xOld,
                                                     final RealVector xNew,
                                                     final RealVector g0,
                                                     final RealMatrix L,
                                                     final RealMatrix absG) {

        // dx = xNew - xOld
        final RealVector dx = xNew.subtract(xOld);

        // Early rejection for exact no-move (stall detected)
        if (dx.getNorm() == 0.0) {
            return false;
        }

        // ------------------------------------------------------------------
        // Exact FDIFF
        //
        //   FDIFF = ( 2*g0 + G*(xOld + xNew) )^T * (xNew - xOld)
        //
        // Using:
        //   gOld = g0 + G*xOld
        // so:
        //   FDIFF = 2*gOld^T*dx + dx^T*G*dx
        //         = 2*gOld^T*dx + ||L^T dx||^2
        // ------------------------------------------------------------------
        final RealVector gOld = L.operate(L.preMultiply(xOld)).add(g0);
        final RealVector ltDx = L.preMultiply(dx);
        final double fDiff = 2.0 * gOld.dotProduct(dx) + ltDx.dotProduct(ltDx);

        // ------------------------------------------------------------------
        // Exact FDIFFA
        //
        //   FDIFFA = ( |2*g0| + |G|*|xOld + xNew| )^T * |dx|
        // ------------------------------------------------------------------
        final RealVector xSumAbs = xOld.add(xNew).map(FastMath::abs);
        final RealVector dxAbs   = dx.map(FastMath::abs);
        final RealVector abs2g0  = g0.mapMultiply(2.0).map(FastMath::abs);

        final RealVector rowAbs = absG.operate(xSumAbs).add(abs2g0);
        final double fDiffA = rowAbs.dotProduct(dxAbs);

        // ------------------------------------------------------------------
        // Rejection logic
        // ------------------------------------------------------------------
        final double s1 = fDiffA + fDiff;
        if (s1 <= fDiffA) {
            return false;
        }

        final double s2 = fDiffA + 1.5 * fDiff;
        if (s2 <= s1) {
            return false;
        }

        return true;
    }

    /**
     * Incremental recomputation of primal x and active multipliers y after a change
     * in the working set (for example after dropping one negative active multiplier).
     *
     * Unlike the absolute KKT reconstruction, this version starts from the current
     * primal point xCurrent and only applies:
     *
     * 1) an active-space correction to reduce the residuals of the current active set,
     * 2) a free-space correction to improve projected stationarity,
     * 3) a recomputation of the active multipliers consistent with the corrected x.
     *
     * Convention:
     * - objective: 0.5 * x^T G x + g0^T x, with G = L*L^T
     * - constraints: C^T x + c0 >= 0
     *
     * @param x  the current primal point
     * @param G         the quadratic objective matrix
     * @param g0        the linear objective vector
     * @param C         the constraint matrix
     * @param c0        the constraint constants
     * @param qrUpdater the factorization updater
     * @param active    the current active set
     * @return a Pair containing the corrected primal x and the active multipliers y
     */
    private Pair<RealVector, RealVector> recomputeXY(final RealVector x,
                                                     final RealMatrix G,
                                                     final RealVector g0,
                                                     final RealMatrix C,
                                                     final RealVector c0,
                                                     final QRUpdater qrUpdater,
                                                     final List<Integer> active) {

        final int nAct = active.size();


        // --------------------------------------------------------------
        // 1) Active-space correction: reduce active residuals starting from xCurrent
        //    rp_i = -(a_i^T x + c0_i)
        // --------------------------------------------------------------
        if (nAct > 0) {
            final double[] rpData = new double[nAct];
            for (int i = 0; i < nAct; i++) {
                final int idx = active.get(i);
                final RealVector ai = C.getRowVector(idx);
                rpData[i] = -(ai.dotProduct(x) + c0.getEntry(idx));
            }
            final RealVector rp = new ArrayRealVector(rpData, false);

            final RealVector tActive = qrUpdater.solveRT(rp);
            if (tActive == null) {
                return null; // singular / numerically unsafe
            }

            // Original: xRec = xRec.add(qrUpdater.applyJActive(tActive));
            // Optimized: in-place addition to avoid object allocation
            x.combineToSelf(1.0, 1.0, qrUpdater.applyJActive(tActive));
        }

        // --------------------------------------------------------------
        // 2) Free-space correction: improve projected stationarity
        //    grad = G*x + g0
        // --------------------------------------------------------------

        // Original: final RealVector gMid = G.operate(xRec).add(g0);
        // Optimized: allocate once for G.operate(xRec), then combine in-place
        final RealVector gMid = G.operate(x);
        gMid.combineToSelf(1.0, 1.0, g0);

        final RealVector dMid = qrUpdater.computeD(gMid);

        // Original: xRec = xRec.subtract(qrUpdater.computeZ(dMid));
        // Optimized: in-place subtraction (xRec = 1.0 * xRec + (-1.0) * Z)
        x.combineToSelf(1.0, -1.0, qrUpdater.computeZ(dMid));

        // --------------------------------------------------------------
        // 3) Recompute active multipliers from corrected xRec
        // --------------------------------------------------------------
        if (nAct == 0) {
            return new Pair<>(x, new ArrayRealVector(0));
        }

        RealVector yActive = qrUpdater.solveR(dMid);

        if (yActive == null) {
            return null;
        }

        if (yActive.getDimension() > nAct) {
            yActive = yActive.getSubVector(0, nAct);
        }

        return new Pair<>(x, yActive);
    }

    /**
     * Private container for the integrated QP system.
     * Designed to be stack-allocated.
     */
    private static final class QPConstraintSystem {
        /** Global constraint matrix C (n x mc) where columns are constraint gradients. */
        private final RealMatrix C;
        /** Right-hand side vector c0 for the QP problem. */
        private final RealVector c0;
        /** Normalization weights based on the 2-norm of each constraint gradient. */
        private final RealVector weights;

        private QPConstraintSystem(RealMatrix C, RealVector c0, RealVector weights) {
            this.C = C;
            this.c0 = c0;
            this.weights = weights;
        }
    }

        /**
     * Builds the virtual row-oriented QP constraint system.
     *
     * <p>The global constraint matrix is represented as:</p>
     *
     * <pre>
     * A = [ Ae  ]
     *     [ Ai  ]
     *     [ Ab  ]
     *     [-Ab  ]
     * </pre>
     *
     * <p>No global physical matrix and no transpose are created.</p>
     *
     * @param n number of primal variables
     * @return integrated constraint system
     */
    private QPConstraintSystem buildGlobalConstraints(final int n) {

        final int p = eqConstraints != null ? eqConstraints.dimY() : 0;

        final int m1 = iqConstraints != null ? iqConstraints.dimY() : 0;

        final int b1 = bConstraints != null ? bConstraints.dimY() : 0;

        final int mc = p + m1 + 2 * b1;

        /*
         * Original row-oriented matrices.
         *
         * No transpose and no coefficient copy are performed.
         */
        final RealMatrix Ae = eqConstraints != null ? eqConstraints.getA() : null;

        final RealMatrix Ai = iqConstraints != null ? iqConstraints.jacobian(null) : null;

        final RealMatrix Ab = bConstraints != null ? bConstraints.jacobian(null) : null;

        /*
         * Virtual global matrix:
         *
         *     [ Ae  ]
         * A = [ Ai  ]
         *     [ Ab  ]
         *     [-Ab  ]
         */
        final RealMatrix A = new QPConstraintMatrix(n,Ae, Ai, Ab);

        final double[] c0Data = new double[mc];
        final double[] wData = new double[mc];

        final double safeMin = Precision.EPSILON;

        /*
         * Equality constraints:
         *
         *     Ae_j x - be_j = 0
         */
        if (p > 0) {

            final RealVector be = eqConstraints.getLowerBound();

            for (int j = 0; j < p; j++) {

                c0Data[j] = -be.getEntry(j);

                double sumSq = 0.0;

                for (int i = 0; i < n; i++) {
                    final double value = Ae.getEntry(j, i);
                    sumSq += value * value;
                }

                final double norm = FastMath.sqrt(sumSq);

                wData[j] = norm > safeMin ?
                        1.0 / norm : 0.0;
            }
        }

        /*
         * Inequality constraints:
         *
         *     Ai_j x - bi_j >= 0
         */
        if (m1 > 0) {

            final RealVector bi = iqConstraints.getLowerBound();

            for (int j = 0; j < m1; j++) {

                final int globalIndex = p + j;

                c0Data[globalIndex] = -bi.getEntry(j);

                double sumSq = 0.0;

                for (int i = 0; i < n; i++) {
                    final double value = Ai.getEntry(j, i);
                    sumSq += value * value;
                }

                final double norm = FastMath.sqrt(sumSq);

                wData[globalIndex] = norm > safeMin ?
                        1.0 / norm : 0.0;
            }
        }

        /*
         * Bounded constraints:
         *
         * Lower:
         *
         *     Ab_j x - lb_j >= 0
         *
         * Upper:
         *
         *     -Ab_j x + ub_j >= 0
         *
         * Both rows have the same Euclidean norm.
         */
        if (b1 > 0) {

            final RealVector lb = bConstraints.getLowerBound();
            final RealVector ub = bConstraints.getUpperBound();

            final int lowerOffset = p + m1;
            final int upperOffset = lowerOffset + b1;

            for (int j = 0; j < b1; j++) {

                final int lowerIndex = lowerOffset + j;
                final int upperIndex = upperOffset + j;

                c0Data[lowerIndex] = -lb.getEntry(j);
                c0Data[upperIndex] = ub.getEntry(j);

                double sumSq = 0.0;

                for (int i = 0; i < n; i++) {
                    final double value = Ab.getEntry(j, i);
                    sumSq += value * value;
                }

                final double norm = FastMath.sqrt(sumSq);

                final double inverseNorm = norm > safeMin ?
                        1.0 / norm : 0.0;

                wData[lowerIndex] = inverseNorm;
                wData[upperIndex] = inverseNorm;
            }
        }

        return new QPConstraintSystem(
                A,
                new ArrayRealVector(c0Data, false),
                new ArrayRealVector(wData, false));
    }

}
