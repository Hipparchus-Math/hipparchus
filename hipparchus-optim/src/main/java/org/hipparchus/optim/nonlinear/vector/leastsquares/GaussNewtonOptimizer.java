/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.optim.nonlinear.vector.leastsquares;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.CholeskyDecomposer;
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.LUDecomposer;
import org.hipparchus.linear.LUDecomposition;
import org.hipparchus.linear.MatrixDecomposer;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.QRDecomposer;
import org.hipparchus.linear.QRDecomposition;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.SingularValueDecomposer;
import org.hipparchus.linear.SingularValueDecomposition;
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.LocalizedOptimFormats;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem.Evaluation;
import org.hipparchus.util.Incrementor;
import org.hipparchus.util.Pair;

/**
 * Gauss-Newton least-squares solver.
 * <p>
 * This class solve a least-square problem by solving the normal equations
 * of the linearized problem at each iteration. Either LU decomposition or
 * Cholesky decomposition can be used to solve the normal equations, or QR
 * decomposition or SVD decomposition can be used to solve the linear system.
 * Cholesky/LU decomposition is faster but QR decomposition is more robust for difficult
 * problems, and SVD can compute a solution for rank-deficient problems.
 */
public class GaussNewtonOptimizer implements LeastSquaresOptimizer {

    /**
     * The decomposition algorithm to use to solve the normal equations.
     *
     * @deprecated Use {@link MatrixDecomposer} instead.
     */
    @Deprecated
    public enum Decomposition {
        /**
         * Solve by forming the normal equations (J<sup>T</sup>Jx=J<sup>T</sup>r) and
         * using the {@link LUDecomposition}.
         *
         * <p> Theoretically this method takes mn<sup>2</sup>/2 operations to compute the
         * normal matrix and n<sup>3</sup>/3 operations (m &gt; n) to solve the system using
         * the LU decomposition. </p>
         */
        LU {
            @Override
            protected RealVector solve(final RealMatrix jacobian,
                                       final RealVector residuals) {
                try {
                    final Pair<RealMatrix, RealVector> normalEquation =
                            computeNormalMatrix(jacobian, residuals);
                    final RealMatrix normal = normalEquation.getFirst();
                    final RealVector jTr = normalEquation.getSecond();
                    return new LUDecomposition(normal, SINGULARITY_THRESHOLD)
                            .getSolver()
                            .solve(jTr);
                } catch (MathIllegalArgumentException e) {
                    throw new MathIllegalStateException(LocalizedOptimFormats.UNABLE_TO_SOLVE_SINGULAR_PROBLEM, e);
                }
            }

            @Override
            protected MatrixDecomposer getDecomposer() {
                return new LUDecomposer(SINGULARITY_THRESHOLD);
            }

            @Override
            protected boolean isFormNormalEquations() {
                return true;
            }
        },
        /**
         * Solve the linear least squares problem (Jx=r) using the {@link
         * QRDecomposition}.
         *
         * <p> Theoretically this method takes mn<sup>2</sup> - n<sup>3</sup>/3 operations
         * (m &gt; n) and has better numerical accuracy than any method that forms the normal
         * equations. </p>
         */
        QR {
            @Override
            protected RealVector solve(final RealMatrix jacobian,
                                       final RealVector residuals) {
                try {
                    return new QRDecomposition(jacobian, SINGULARITY_THRESHOLD)
                            .getSolver()
                            .solve(residuals);
                } catch (MathIllegalArgumentException e) {
                    throw new MathIllegalStateException(LocalizedOptimFormats.UNABLE_TO_SOLVE_SINGULAR_PROBLEM, e);
                }
            }

            @Override
            protected MatrixDecomposer getDecomposer() {
                return new QRDecomposer(SINGULARITY_THRESHOLD);
            }

            @Override
            protected boolean isFormNormalEquations() {
                return false;
            }
        },
        /**
         * Solve by forming the normal equations (J<sup>T</sup>Jx=J<sup>T</sup>r) and
         * using the {@link CholeskyDecomposition}.
         *
         * <p> Theoretically this method takes mn<sup>2</sup>/2 operations to compute the
         * normal matrix and n<sup>3</sup>/6 operations (m &gt; n) to solve the system using
         * the Cholesky decomposition. </p>
         */
        CHOLESKY {
            @Override
            protected RealVector solve(final RealMatrix jacobian,
                                       final RealVector residuals) {
                try {
                    final Pair<RealMatrix, RealVector> normalEquation =
                            computeNormalMatrix(jacobian, residuals);
                    final RealMatrix normal = normalEquation.getFirst();
                    final RealVector jTr = normalEquation.getSecond();
                    return new CholeskyDecomposition(
                            normal, SINGULARITY_THRESHOLD, SINGULARITY_THRESHOLD)
                            .getSolver()
                            .solve(jTr);
                } catch (MathIllegalArgumentException e) {
                    if (e.getSpecifier() == LocalizedCoreFormats.NOT_POSITIVE_DEFINITE_MATRIX) {
                        throw new MathIllegalStateException(LocalizedOptimFormats.UNABLE_TO_SOLVE_SINGULAR_PROBLEM, e);
                    } else {
                        throw e;
                    }
                }
            }

            @Override
            protected MatrixDecomposer getDecomposer() {
                return new CholeskyDecomposer(SINGULARITY_THRESHOLD,
                        SINGULARITY_THRESHOLD);
            }

            @Override
            protected boolean isFormNormalEquations() {
                return true;
            }
        },
        /**
         * Solve the linear least squares problem using the {@link
         * SingularValueDecomposition}.
         *
         * <p> This method is slower, but can provide a solution for rank deficient and
         * nearly singular systems.
         */
        SVD {
            @Override
            protected RealVector solve(final RealMatrix jacobian,
                                       final RealVector residuals) {
                return new SingularValueDecomposition(jacobian)
                        .getSolver()
                        .solve(residuals);
            }

            @Override
            protected MatrixDecomposer getDecomposer() {
                return new SingularValueDecomposer();
            }

            @Override
            protected boolean isFormNormalEquations() {
                return false;
            }
        };

        /**
         * Solve the linear least squares problem Jx=r.
         *
         * @param jacobian  the Jacobian matrix, J. the number of rows &gt;= the number or
         *                  columns.
         * @param residuals the computed residuals, r.
         * @return the solution x, to the linear least squares problem Jx=r.
         * @throws MathIllegalStateException if the matrix properties (e.g. singular) do not
         *                              permit a solution.
         */
        protected abstract RealVector solve(RealMatrix jacobian,
                                            RealVector residuals);

        /**
         * Get the equivalent matrix decomposer.
         *
         * @return the decomposer.
         */
        protected abstract MatrixDecomposer getDecomposer();

        /**
         * Get if this decomposition forms the normal equations explicitly.
         *
         * @return {@code true} if the normal equations are formed explicitly, {@code
         * false} otherwise.
         */
        protected abstract boolean isFormNormalEquations();

    }

    /**
     * The singularity threshold for matrix decompositions. Determines when a {@link
     * MathIllegalStateException} is thrown. The current value was the default value for {@link
     * LUDecomposition}.
     */
    private static final double SINGULARITY_THRESHOLD = 1e-11;

    /** Decomposition. */
    @Deprecated
    private final Decomposition decomposition;
    /** Decomposer */
    private final MatrixDecomposer decomposer;
    /** Indicates if normal equations should be formed explicitly. */
    private final boolean formNormalEquations;

    /**
     * Creates a Gauss Newton optimizer.
     * <p/>
     * The default for the algorithm is to solve the normal equations using QR
     * decomposition.
     */
    public GaussNewtonOptimizer() {
        this(Decomposition.QR);
    }

    /**
     * Create a Gauss Newton optimizer that uses the given decomposition algorithm to
     * solve the normal equations.
     *
     * @param decomposition the {@link Decomposition} algorithm.
     * @deprecated Use {@link #GaussNewtonOptimizer(MatrixDecomposer, boolean)} instead.
     * The new constructor provides control of the numerical tolerances in the
     * decomposition.
     */
    @Deprecated
    public GaussNewtonOptimizer(final Decomposition decomposition) {
        this(decomposition, decomposition.getDecomposer(),
                decomposition.isFormNormalEquations());
    }

    /**
     * Create a Gauss Newton optimizer that uses the given matrix decomposition algorithm
     * to solve the normal equations.
     *
     * @param decomposer          the decomposition algorithm to use.
     * @param formNormalEquations whether the normal equations should be explicitly
     *                            formed. If {@code true} then {@code decomposer} is used
     *                            to solve J<sup>T</sup>Jx=J<sup>T</sup>r, otherwise
     *                            {@code decomposer} is used to solve Jx=r. If {@code
     *                            decomposer} can only solve square systems then this
     *                            parameter should be {@code true}.
     */
    public GaussNewtonOptimizer(final MatrixDecomposer decomposer,
                                final boolean formNormalEquations) {
        this(null, decomposer, formNormalEquations);
    }

    /**
     * Bridge constructor until the deprecated field can be removed.
     *
     * @param decomposition       the value of the deprecated field.
     * @param decomposer          to use.
     * @param formNormalEquations explicitly or not.
     */
    @Deprecated
    private GaussNewtonOptimizer(final Decomposition decomposition,
                                 final MatrixDecomposer decomposer,
                                 final boolean formNormalEquations) {
        this.decomposition = decomposition;
        this.decomposer = decomposer;
        this.formNormalEquations = formNormalEquations;
    }

    /**
     * Get the matrix decomposition algorithm used to solve the normal equations.
     *
     * @return the matrix {@link Decomposition} algorithm. May be {@code null}.
     * @deprecated Use {@link #getDecomposer()} and {@link #isFormNormalEquations()}
     * instead.
     */
    @Deprecated
    public Decomposition getDecomposition() {
        return this.decomposition;
    }

    /**
     * Configure the decomposition algorithm.
     *
     * @param newDecomposition the {@link Decomposition} algorithm to use.
     * @return a new instance.
     * @deprecated Use {@link #withDecomposer(MatrixDecomposer)} and {@link
     * #withFormNormalEquations(boolean)} instead. the new methods allow the numerical
     * tolerance of the decomposition to be set.
     */
    @Deprecated
    public GaussNewtonOptimizer withDecomposition(final Decomposition newDecomposition) {
        return new GaussNewtonOptimizer(newDecomposition);
    }

    /**
     * Get the matrix decomposition algorithm.
     *
     * @return the decomposition algorithm.
     */
    public MatrixDecomposer getDecomposer() {
        return decomposer;
    }

    /**
     * Configure the matrix decomposition algorithm.
     *
     * @param newDecomposer the decomposition algorithm to use.
     * @return a new instance.
     */
    public GaussNewtonOptimizer withDecomposer(final MatrixDecomposer newDecomposer) {
        return new GaussNewtonOptimizer(newDecomposer, this.isFormNormalEquations());
    }

    /**
     * Get if the normal equations are explicitly formed.
     *
     * @return if the normal equations should be explicitly formed. If {@code true} then
     * {@code decomposer} is used to solve J<sup>T</sup>Jx=J<sup>T</sup>r, otherwise
     * {@code decomposer} is used to solve Jx=r.
     */
    public boolean isFormNormalEquations() {
        return formNormalEquations;
    }

    /**
     * Configure if the normal equations should be explicitly formed.
     *
     * @param newFormNormalEquations whether the normal equations should be explicitly
     *                               formed. If {@code true} then {@code decomposer} is used
     *                               to solve J<sup>T</sup>Jx=J<sup>T</sup>r, otherwise
     *                               {@code decomposer} is used to solve Jx=r. If {@code
     *                               decomposer} can only solve square systems then this
     *                               parameter should be {@code true}.
     * @return a new instance.
     */
    public GaussNewtonOptimizer withFormNormalEquations(final boolean newFormNormalEquations) {
        return new GaussNewtonOptimizer(this.getDecomposer(), newFormNormalEquations);
    }

    /** {@inheritDoc} */
    @Override
    public Optimum optimize(final LeastSquaresProblem lsp) {
        //create local evaluation and iteration counts
        final Incrementor evaluationCounter = lsp.getEvaluationCounter();
        final Incrementor iterationCounter = lsp.getIterationCounter();
        final ConvergenceChecker<Evaluation> checker
                = lsp.getConvergenceChecker();

        // Computation will be useless without a checker (see "for-loop").
        if (checker == null) {
            throw new NullArgumentException();
        }

        RealVector currentPoint = lsp.getStart();

        // iterate until convergence is reached
        Evaluation current = null;
        while (true) {
            iterationCounter.increment();

            // evaluate the objective function and its jacobian
            Evaluation previous = current;
            // Value of the objective function at "currentPoint".
            evaluationCounter.increment();
            current = lsp.evaluate(currentPoint);
            final RealVector currentResiduals = current.getResiduals();
            final RealMatrix weightedJacobian = current.getJacobian();
            currentPoint = current.getPoint();

            // Check convergence.
            if (previous != null &&
                checker.converged(iterationCounter.getCount(), previous, current)) {
                return Optimum.of(current,
                                  evaluationCounter.getCount(),
                                  iterationCounter.getCount());
            }

            // solve the linearized least squares problem
            final RealMatrix lhs; // left hand side
            final RealVector rhs; // right hand side
            if (this.formNormalEquations) {
                final Pair<RealMatrix, RealVector> normalEquation =
                        computeNormalMatrix(weightedJacobian, currentResiduals);
                lhs = normalEquation.getFirst();
                rhs = normalEquation.getSecond();
            } else {
                lhs = weightedJacobian;
                rhs = currentResiduals;
            }
            final RealVector dX;
            try {
                dX = this.decomposer.decompose(lhs).solve(rhs);
            } catch (MathIllegalArgumentException e) {
                // change exception message
                throw new MathIllegalStateException(
                        LocalizedOptimFormats.UNABLE_TO_SOLVE_SINGULAR_PROBLEM, e);
            }
            // update the estimated parameters
            currentPoint = currentPoint.add(dX);
        }
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {
        return "GaussNewtonOptimizer{" +
                "decomposer=" + decomposer +
                ", formNormalEquations=" + formNormalEquations +
                '}';
    }

    /**
     * Compute the normal matrix, J<sup>T</sup>J.
     *
     * @param jacobian  the m by n jacobian matrix, J. Input.
     * @param residuals the m by 1 residual vector, r. Input.
     * @return  the n by n normal matrix and  the n by 1 J<sup>Tr vector.
     */
    private static Pair<RealMatrix, RealVector> computeNormalMatrix(final RealMatrix jacobian,
                                                                    final RealVector residuals) {
        //since the normal matrix is symmetric, we only need to compute half of it.
        final int nR = jacobian.getRowDimension();
        final int nC = jacobian.getColumnDimension();
        //allocate space for return values
        final RealMatrix normal = MatrixUtils.createRealMatrix(nC, nC);
        final RealVector jTr = new ArrayRealVector(nC);
        //for each measurement
        for (int i = 0; i < nR; ++i) {
            //compute JTr for measurement i
            for (int j = 0; j < nC; j++) {
                jTr.setEntry(j, jTr.getEntry(j) +
                        residuals.getEntry(i) * jacobian.getEntry(i, j));
            }

            // add the the contribution to the normal matrix for measurement i
            for (int k = 0; k < nC; ++k) {
                //only compute the upper triangular part
                for (int l = k; l < nC; ++l) {
                    normal.setEntry(k, l, normal.getEntry(k, l) +
                            jacobian.getEntry(i, k) * jacobian.getEntry(i, l));
                }
            }
        }
        //copy the upper triangular part to the lower triangular part.
        for (int i = 0; i < nC; i++) {
            for (int j = 0; j < i; j++) {
                normal.setEntry(i, j, normal.getEntry(j, i));
            }
        }
        return new Pair<RealMatrix, RealVector>(normal, jTr);
    }

}
