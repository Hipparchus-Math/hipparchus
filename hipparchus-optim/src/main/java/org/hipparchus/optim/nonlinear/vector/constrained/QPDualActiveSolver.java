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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


/**
 * Implements the dual active-set method by Goldfarb and Idnani (1983) for
 * solving strictly convex quadratic programs of the form:
 * <pre>
 * minimize   (1/2) x^T G x + g0^T x
 * subject to CE^T x  =  ce0     (equality constraints)
 *            CI^T x  &gt;= ci0     (inequality constraints)
 *            lb &lt;= Ax &lt;= ub      (bounded constraints)
 * </pre> <p>
 * Uses an incremental QR factorization updater (QRUpdater) for managing active
 * constraints and maintains multipliers for dual and primal steps.
 * </p>
 *
 * @see <a href="https://doi.org/10.1137/0603006">Goldfarb and Idnani (1983)</a>
 * @since 1.11
 */
public class QPDualActiveSolver extends QPOptimizer {

    /**
     * Machine epsilon for tolerance checks.
     */
    private static final double EPS = Math.ulp(1.0);

    /**
     * Maximum number of iterations allowed.Will be adjusted in base of problem
     * dimension
     */
    private int maxIter;

    /**
     * Quadratic function representing 1/2 x^T G x + g0^T x.
     */
    private QuadraticFunction function;

    /**
     * Equality constraint data (CE^T x + ce0 = 0).
     */
    private LinearEqualityConstraint eqConstraints;

    /**
     * Inequality constraint data (CI^T x + ci0 >= 0).
     */
    private LinearInequalityConstraint iqConstraints;

    /**
     * Bounded constraint data (lower <= Ax <= upper).
     */
    private LinearBoundedConstraint bConstraints;

    /** Tolerance for symmetric matrix decomposition.
     * @since 4.1
     */
    private MatrixDecompositionTolerance matrixDecompositionTolerance;

     /**
     * Inverse of Cholesky factorization if passed from external.
     */
    private RealMatrix inverseL;

    /**
     * Parses optimization data to extract the objective function and various
     * constraint sets.
     * @param optData optimization data
     */
    @Override
    protected void parseOptimizationData(OptimizationData... optData) {
        super.parseOptimizationData(optData);
        //reset QP problem to reuse the same instance of the QP solver;
        this.maxIter = 1000;
        this.function = null;
        this.eqConstraints = null;
        this.iqConstraints = null;
        this.bConstraints = null;
        this.matrixDecompositionTolerance = new MatrixDecompositionTolerance(EPS);
        this.inverseL = null;
        for (OptimizationData data : optData) {
            if (data instanceof ObjectiveFunction) {
                function = (QuadraticFunction) ((ObjectiveFunction) data).getObjectiveFunction();
            } else if (data instanceof LinearEqualityConstraint) {
                eqConstraints = (LinearEqualityConstraint) data;
            } else if (data instanceof LinearInequalityConstraint) {
                iqConstraints = (LinearInequalityConstraint) data;
            } else if (data instanceof LinearBoundedConstraint) {
                bConstraints = (LinearBoundedConstraint) data;
            } else if (data instanceof InverseCholesky) {
                inverseL = ((InverseCholesky) data).getInverseL();
            } else if (data instanceof MatrixDecompositionTolerance) {
                matrixDecompositionTolerance = (MatrixDecompositionTolerance) data;
            }
        }
    }

    /**
     * Finds step for a primal move.
     *
     * @param z search direction
     * @param ai constraint row vector
     * @param sv current violation offset
     * @param equality true if constraint is equality
     * @return the primal step size
     */
    private double findPrimalStep(final RealVector z, final RealVector ai, final double sv, final boolean equality) {
        double norm2 = z.dotProduct(z);
        if (FastMath.abs(norm2) < Precision.EPSILON) {
            return Double.POSITIVE_INFINITY;
        }

        double denom = ai.dotProduct(z);
        double alpha = -sv / denom;
        //step for inequality should be positive
        if (!equality && alpha < 0) {
            alpha = Double.POSITIVE_INFINITY;
        }

        return alpha;
    }

    /**
     * Finds the blocking step for a dual move.
     *
     * @param u current multipliers
     * @param r dual direction
     * @param activeSet map of active indices -> constraint ids
     * @param me threshold index for equality
     * @return the blocking step and index
     */
    private Pair<Integer, Double> findDualBlockingConstraint(final RealVector u, final RealVector r,
                                                             final List<Integer> activeSet, final int me) {
        if (activeSet.isEmpty()) {
            return new Pair<>(-1, Double.POSITIVE_INFINITY);
        }
        double alpha = Double.POSITIVE_INFINITY;
        int block = -1;
        int size = activeSet.size();
        for (int i = 0; i < size; i++) {
            double ui = u.getEntry(i);
            double ri = r.getEntry(i);
            //consider only inequality constraint in the active set
            if (ri > 0 && activeSet.get(i) >= me) {
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
     * @param u current u or null
     * @param r dual direction
     * @param alpha step length
     * @param partial new partial multiplier
     * @return updated multipliers
     */
    private RealVector updateMultipliersOnAddition(final RealVector u, final RealVector r,
                                                   final double alpha, final double partial) {
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
     * @param u current multipliers
     * @param r dual direction
     * @param alpha step length
     * @param dropIndex index to remove
     * @return updated multipliers
     */
    private RealVector updateMultipliersOnRemoval(final RealVector u, final RealVector r,
                                                  final double alpha, final int dropIndex) {
        if (u.getDimension() == 1) {
            return new ArrayRealVector(0,0);
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
    /**
     * Find the most violated constraint.
     *
     * @param sv current constraint evalutation
     * @param blackList dependent constraint
     * @param activeSet constraint in active set
     * @param me equality constraint numbers
     * @return the violation value and index
     */
    private Pair<Integer, Double> mostViolatedConstraint(RealVector sv, Set<Integer> blackList, List<Integer> activeSet, int me) {
        double maxViolation = 0;
        int mostViolated = -1;
        for (int i = 0; i < sv.getDimension(); i++) {
            if (blackList.contains(me + i) || activeSet.contains(me + i)) {
                continue;
            }
            double violation = sv.getEntry(i);
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
        RealMatrix G = function.getP();
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
        RealMatrix CI = null;
        RealVector ci0 = null;
        if (m > 0) {
            CI  = MatrixUtils.createRealMatrix(n, m);
            ci0 = new ArrayRealVector(m);
            if (m1 > 0) {
                RealMatrix Aineq = iqConstraints.jacobian(null);
                RealVector bineq = iqConstraints.getLowerBound();
                CI.setSubMatrix(Aineq.transpose().getData(), 0, 0);
                ci0.setSubVector(0, bineq.mapMultiply(-1.0));
            }
            if (b1 > 0) {
                RealMatrix Abound = bConstraints.jacobian(null);
                RealVector lower  = bConstraints.getLowerBound();
                RealVector upper  = bConstraints.getUpperBound();
                CI.setSubMatrix(Abound.transpose().getData(), 0, m1);
                CI.setSubMatrix(Abound.scalarMultiply(-1.0).transpose().getData(), 0, m1 + b1);
                ci0.setSubVector(m1, lower.mapMultiply(-1.0));
                ci0.setSubVector(m1 + b1, upper);
            }
        }
        RealVector x;
        RealMatrix L;
        RealMatrix L1;
        QRUpdater qrUpdater;
        double tol;
        if (this.inverseL == null) {
            try {
                final double eps = matrixDecompositionTolerance.getEpsMatrixDecomposition();
                final CholeskyDecomposition cholesky = new CholeskyDecomposition(G, eps, eps);
                DecompositionSolver solver = cholesky.getSolver();
                x = solver.solve(g0).mapMultiply(-1.0);
                L = cholesky.getL();
                L1 = inverseLowerTriangular(L);
                //c1 trace of G matrix
               double c1 = FastMath.sqrt(G.getTrace());
               //c2 trace of inverse of cholesky factorization
               double c2 = FastMath.sqrt(L1.getTrace());
               tol = m * c1 * c2 * Precision.EPSILON * 100.0;
                qrUpdater = new QRUpdater(L1);
            } catch (MathIllegalArgumentException ex) {
                // matrix is not positive definite return empty solution
                return new LagrangeSolution(new ArrayRealVector(0,0), new ArrayRealVector(0,0), 0.0);
            }
        } else {
            L = this.inverseL;
            L1 = inverseLowerTriangular(L);
            RealMatrix G1 = L1.multiplyTransposed(L1);
            x = L.preMultiply(L.operate(g0)).mapMultiply(-1.0);
            double c1 = FastMath.sqrt(G1.getTrace());
            double c2 = FastMath.sqrt(L.getTrace());
            tol = m * c1 * c2 * Precision.EPSILON * 100.0;
            qrUpdater = new QRUpdater(L);
        }
        if (m + p == 0) {
            return new LagrangeSolution(x,
                                        new ArrayRealVector(0,0),
                                        0.5 * x.dotProduct(G.operate(x)) + g0.dotProduct(x) + g);
        }
        //max iteration adjusted in base of problem dimension
        this.maxIter = 40 * (n + m + p);

        //convergence threshold calculated in base at the matrix conditioning
        //ActiveSet and blackLit(dependent constraints)

        final Set<Integer> blacklist = new HashSet<>();
        List<Integer>      active    = new ArrayList<>();

        RealVector u = new ArrayRealVector(0,0);
        RealVector r = new ArrayRealVector(0,0);
        RealVector d = null;
        RealVector z = null;
        // Add equality constraints in the active set updating x solution and multipliers
        for (int i = 0; i < p; i++) {
            RealVector ai = CE.getColumnVector(i);
            double sve = ai.dotProduct(x) + ce0.getEntry(i);
            RealMatrix Q = qrUpdater.getJ();
            d = Q.transpose().operate(ai);
            RealMatrix J2 = qrUpdater.getJ2();
            z = (n - active.size() > 0) ?
                J2.operate(d.getSubVector(active.size(), n - active.size())) :
                new ArrayRealVector(n);
            if (!active.isEmpty()) {
                r = qrUpdater.getRInv().operate(d.getSubVector(0, active.size()));
            }
            double alpha = findPrimalStep(z, ai, sve, true);
            x = x.add(z.mapMultiply(alpha));
            u = updateMultipliersOnAddition(u, r, alpha, alpha);
            if (!qrUpdater.addConstraint(d)) {
                return null;//equality constraint are linearly dependent
            }
            active.add(i);
        }
        int iteration = 0;

        // Active-set loop for inequalities
        while (m != 0 && iteration++ < maxIter) {

            RealVector sv;
            //store solution in case constraint can't be added because dependent
            RealVector xOld = x;
            RealVector uOld = u;
            //evaluate inequality constraints
            sv = CI.transpose().operate(x).add(ci0);

            //calculate norm1 of the constraints
            double sum = 0;
            for (int k = 0; k < sv.getDimension(); k++) {
                sum += FastMath.min(0.0, sv.getEntry(k));
            }

            // Evaluate convergence
            if (FastMath.abs(sum) <= tol) {
                break;// Optimal solution found
            }

            // Evaluate most violated constraint, excluding dependent/active loop
            while (iteration++ < maxIter) {
                final Pair<Integer, Double> mostViolated = mostViolatedConstraint(sv, blacklist, active, p);
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
                // Dual step loop update multiplier and x (if step is also in primal) until primal step is not done
                while (iteration++ < maxIter) {
                    np = CI.getColumnVector(mostViolated.getKey());
                    sv.setEntry(mostViolated.getKey(), np.dotProduct(x) + ci0.getEntry(mostViolated.getKey()));
                    d = qrUpdater.getJ().transpose().operate(np);
                    J2 = qrUpdater.getJ2();
                    z = (n - active.size() > 0) ?
                        J2.operate(d.getSubVector(active.size(), n - active.size())) :
                        new ArrayRealVector(n);

                    if (!active.isEmpty()) {
                        r = qrUpdater.getRInv().operate(d.getSubVector(0, active.size()));
                    }

                    t1 = findPrimalStep(z, np, sv.getEntry(mostViolated.getKey()), false);
                    Pair<Integer, Double> dualStep = findDualBlockingConstraint(u, r, active, p);
                    t2 = dualStep.getValue();
                    dropIndex = dualStep.getKey();
                    t = FastMath.min(t1, t2);
                    if (t == t1) {
                        break; // primal full step (exit from dual step loop)
                    } else {
                        //Manage dual step
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
                if (qrUpdater.addConstraint(d)&& active.size()<n) {

                    active.add(p + mostViolated.getKey());
                    x = x.add(z.mapMultiply(t));
                    uPartial += t;
                    u = updateMultipliersOnAddition(u, r, t, uPartial);
                    blacklist.clear();
                    break; //revaluate convergence (exit from most violated constraint loop)
                } else {
                    // dependent constraint -> add in blacklist and revert state
                    // revaluate only violated constraint without recalculate them;
                    blacklist.add(p + mostViolated.getKey());
                    x = xOld;
                    u = uOld;
                }
            }
        }
        if (iteration == maxIter) {
            return new LagrangeSolution(new ArrayRealVector(0,0),
                                        new ArrayRealVector(0,0),
                                        0.0); // no optimal solution is found
        }
        return buildSolution(x, u, active, G, g0, g, p, m);
    }

    /**
     * Builds the final solution object with primal x and multipliers lambda.
     * @param x solution
     * @param u active set multipliers
     * @param activeSet active set constraints
     * @param G square matrix of weights for quadratic terms
     * @param g0 vector of weights for linear terms
     * @param g constant term
     * @param p number of equalities
     * @param m number of inequalities
     * @return the optimal solution
     */
    private LagrangeSolution buildSolution(final RealVector x, final RealVector u, final List<Integer> activeSet,
                                           final RealMatrix G, final RealVector g0, final double g, final int p, final int m) {
        final RealVector lambda = new ArrayRealVector(p + m);
        if (!activeSet.isEmpty()) {
            for (int i = 0; i < activeSet.size(); i++) {
                lambda.setEntry(activeSet.get(i), u.getEntry(i));
            }
        }
        final double value = 0.5 * x.dotProduct(G.operate(x)) +
                             g0.dotProduct(x) + g;
        return new LagrangeSolution(x, lambda, value);
    }

    /**
     * Computes the inverse of a lower-triangular matrix via forward
     * substitution.
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

}
