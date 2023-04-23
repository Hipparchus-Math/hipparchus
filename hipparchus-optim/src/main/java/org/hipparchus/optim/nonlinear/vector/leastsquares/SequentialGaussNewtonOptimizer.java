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
package org.hipparchus.optim.nonlinear.vector.leastsquares;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.MatrixDecomposer;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.QRDecomposer;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.LocalizedOptimFormats;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem.Evaluation;
import org.hipparchus.util.Incrementor;
import org.hipparchus.util.Pair;

/**
 * Sequential Gauss-Newton least-squares solver.
 * <p>
 * This class solve a least-square problem by solving the normal equations of
 * the linearized problem at each iteration.
 * </p>
 *
 */
public class SequentialGaussNewtonOptimizer implements LeastSquaresOptimizer {

    /**
     * The singularity threshold for matrix decompositions. Determines when a
     * {@link MathIllegalStateException} is thrown. The current value was the
     * default value for {@link LUDecomposition}.
     */
    private static final double SINGULARITY_THRESHOLD = 1e-11;

    /** Decomposer. */
    private final MatrixDecomposer decomposer;

    /** Indicates if normal equations should be formed explicitly. */
    private final boolean formNormalEquations;

    /** Old evaluation previously computed. */
    private final Evaluation oldEvaluation;

    /** Old jacobian previously computed. */
    private final RealMatrix oldLhs;

    /** Old residuals previously computed. */
    private final RealVector oldRhs;

    /**
     * Create a sequential Gauss Newton optimizer.
     * <p>
     * The default for the algorithm is to use QR decomposition, not
     * form normal equations and have no previous evaluation
     * </p>
     *
     */
    public SequentialGaussNewtonOptimizer() {
        this(new QRDecomposer(SINGULARITY_THRESHOLD), false, null);
    }

    /**
     * Create a sequential Gauss Newton optimizer that uses the given matrix
     * decomposition algorithm to solve the normal equations.
     * <p>
     * The {@code decomposer} is used to solve J<sup>T</sup>Jx=J<sup>T</sup>r.
     * </p>
     *
     * @param decomposer the decomposition algorithm to use.
     * @param formNormalEquations whether the normal equations should be explicitly
     *                            formed. If {@code true} then {@code decomposer} is used
     *                            to solve J<sup>T</sup>Jx=J<sup>T</sup>r, otherwise
     *                            {@code decomposer} is used to solve Jx=r. If {@code
     *                            decomposer} can only solve square systems then this
     *                            parameter should be {@code true}.
     * @param evaluation old evaluation previously computed, null if there are no previous evaluations.
     */
    public SequentialGaussNewtonOptimizer(final MatrixDecomposer decomposer,
                                          final boolean formNormalEquations,
                                          final Evaluation evaluation) {
        this.decomposer          = decomposer;
        this.formNormalEquations = formNormalEquations;
        this.oldEvaluation       = evaluation;
        if (evaluation == null) {
            this.oldLhs = null;
            this.oldRhs = null;
        } else {
            if (formNormalEquations) {
                final Pair<RealMatrix, RealVector> normalEquation =
                                computeNormalMatrix(evaluation.getJacobian(), evaluation.getResiduals());
                // solve the linearized least squares problem
                this.oldLhs = normalEquation.getFirst();
                this.oldRhs = normalEquation.getSecond();
            } else {
                this.oldLhs = evaluation.getJacobian();
                this.oldRhs = evaluation.getResiduals();
            }
        }
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
    public SequentialGaussNewtonOptimizer withDecomposer(final MatrixDecomposer newDecomposer) {
        return new SequentialGaussNewtonOptimizer(newDecomposer,
                                                  this.isFormNormalEquations(),
                                                  this.getOldEvaluation());
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
    public SequentialGaussNewtonOptimizer withFormNormalEquations(final boolean newFormNormalEquations) {
        return new SequentialGaussNewtonOptimizer(this.getDecomposer(),
                                                  newFormNormalEquations,
                                                  this.getOldEvaluation());
    }

    /**
     * Get the previous evaluation used by the optimizer.
     *
     * @return the previous evaluation.
     */
    public Evaluation getOldEvaluation() {
        return oldEvaluation;
    }

    /**
     * Configure the previous evaluation used by the optimizer.
     * <p>
     * This building method uses a complete evaluation to retrieve
     * a priori data. Note that as {@link #withAPrioriData(RealVector, RealMatrix)}
     * generates a fake evaluation and calls this method, either
     * {@link #withAPrioriData(RealVector, RealMatrix)} or {@link #withEvaluation(LeastSquaresProblem.Evaluation)}
     * should be called, but not both as the last one called will override the previous one.
     * </p>
     * @param previousEvaluation the previous evaluation used by the optimizer.
     * @return a new instance.
     */
    public SequentialGaussNewtonOptimizer withEvaluation(final Evaluation previousEvaluation) {
        return new SequentialGaussNewtonOptimizer(this.getDecomposer(),
                                                  this.isFormNormalEquations(),
                                                  previousEvaluation);
    }

    /**
     * Configure from a priori state and covariance.
     * <p>
     * This building method generates a fake evaluation and calls
     * {@link #withEvaluation(LeastSquaresProblem.Evaluation)}, so either
     * {@link #withAPrioriData(RealVector, RealMatrix)} or {@link #withEvaluation(LeastSquaresProblem.Evaluation)}
     * should be called, but not both as the last one called will override the previous one.
     * <p>
     * A Cholesky decomposition is used to compute the weighted jacobian from the
     * a priori covariance. This method uses the default thresholds of the decomposition.
     * </p>
     * @param aPrioriState a priori state to use
     * @param aPrioriCovariance a priori covariance to use
     * @return a new instance.
     * @see #withAPrioriData(RealVector, RealMatrix, double, double)
     */
    public SequentialGaussNewtonOptimizer withAPrioriData(final RealVector aPrioriState,
                                                          final RealMatrix aPrioriCovariance) {
        return withAPrioriData(aPrioriState, aPrioriCovariance,
                               CholeskyDecomposition.DEFAULT_RELATIVE_SYMMETRY_THRESHOLD,
                               CholeskyDecomposition.DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD);
    }

    /**
     * Configure from a priori state and covariance.
     * <p>
     * This building method generates a fake evaluation and calls
     * {@link #withEvaluation(LeastSquaresProblem.Evaluation)}, so either
     * {@link #withAPrioriData(RealVector, RealMatrix)} or {@link #withEvaluation(LeastSquaresProblem.Evaluation)}
     * should be called, but not both as the last one called will override the previous one.
     * <p>
     * A Cholesky decomposition is used to compute the weighted jacobian from the
     * a priori covariance.
     * </p>
     * @param aPrioriState a priori state to use
     * @param aPrioriCovariance a priori covariance to use
     * @param relativeSymmetryThreshold Cholesky decomposition threshold above which off-diagonal
     * elements are considered too different and matrix not symmetric
     * @param absolutePositivityThreshold Cholesky decomposition threshold below which diagonal
     * elements are considered null and matrix not positive definite
     * @return a new instance.
     * @since 2.3
     */
    public SequentialGaussNewtonOptimizer withAPrioriData(final RealVector aPrioriState,
                                                          final RealMatrix aPrioriCovariance,
                                                          final double relativeSymmetryThreshold,
                                                          final double absolutePositivityThreshold) {

        // we consider the a priori state and covariance come from a
        // previous estimation with exactly one observation of each state
        // component, so partials are the identity matrix, weight is the
        // square root of inverse of covariance, and residuals are zero

        // create a fake weighted Jacobian
        final RealMatrix jTj              = getDecomposer().decompose(aPrioriCovariance).getInverse();
        final RealMatrix weightedJacobian = new CholeskyDecomposition(jTj,
                                                                      relativeSymmetryThreshold,
                                                                      absolutePositivityThreshold).getLT();

        // create fake zero residuals
        final RealVector residuals        = MatrixUtils.createRealVector(aPrioriState.getDimension());

        // combine everything as an evaluation
        final Evaluation fakeEvaluation   = new AbstractEvaluation(aPrioriState.getDimension()) {

            /** {@inheritDoc} */
            @Override
            public RealVector getResiduals() {
                return residuals;
            }

            /** {@inheritDoc} */
            @Override
            public RealVector getPoint() {
                return aPrioriState;
            }

            /** {@inheritDoc} */
            @Override
            public RealMatrix getJacobian() {
                return weightedJacobian;
            }
        };

        return withEvaluation(fakeEvaluation);

    }

    /** {@inheritDoc} */
    @Override
    public Optimum optimize(final LeastSquaresProblem lsp) {
        // create local evaluation and iteration counts
        final Incrementor evaluationCounter = lsp.getEvaluationCounter();
        final Incrementor iterationCounter = lsp.getIterationCounter();
        final ConvergenceChecker<Evaluation> checker =
            lsp.getConvergenceChecker();

        // Computation will be useless without a checker (see "for-loop").
        if (checker == null) {
            throw new NullArgumentException();
        }

        RealVector currentPoint = lsp.getStart();

        if (oldEvaluation != null &&
            currentPoint.getDimension() != oldEvaluation.getPoint().getDimension()) {
            throw new MathIllegalStateException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                      currentPoint.getDimension(), oldEvaluation.getPoint().getDimension());
        }

        // iterate until convergence is reached
        Evaluation current = null;
        while (true) {
            iterationCounter.increment();

            // evaluate the objective function and its jacobian
            final Evaluation previous = current;

            // Value of the objective function at "currentPoint".
            evaluationCounter.increment();
            current = lsp.evaluate(currentPoint);
            final RealVector currentResiduals = current.getResiduals();
            final RealMatrix weightedJacobian = current.getJacobian();

            currentPoint = current.getPoint();

            // Check convergence.
            if (previous != null &&
                checker.converged(iterationCounter.getCount(), previous,
                                  current)) {
                // combine old and new evaluations
                final Evaluation combinedEvaluation = oldEvaluation == null ?
                                                      current :
                                                      new CombinedEvaluation(oldEvaluation, current);
                return Optimum.of(combinedEvaluation, evaluationCounter.getCount(),
                                  iterationCounter.getCount());
            }

           // solve the linearized least squares problem
            final RealMatrix lhs; // left hand side
            final RealVector rhs; // right hand side
            if (this.formNormalEquations) {
                final Pair<RealMatrix, RealVector> normalEquation =
                                computeNormalMatrix(weightedJacobian, currentResiduals);

                lhs = oldLhs == null ?
                      normalEquation.getFirst() :
                      normalEquation.getFirst().add(oldLhs); // left hand side
                rhs = oldRhs == null ?
                      normalEquation.getSecond() :
                      normalEquation.getSecond().add(oldRhs); // right hand side
            } else {
                lhs = oldLhs == null ?
                      weightedJacobian :
                      combineJacobians(oldLhs, weightedJacobian);
                rhs = oldRhs == null ?
                      currentResiduals :
                      combineResiduals(oldRhs, currentResiduals);
            }

            final RealVector dX;
            try {
                dX = this.decomposer.decompose(lhs).solve(rhs);
            } catch (MathIllegalArgumentException e) {
                // change exception message
                throw new MathIllegalStateException(LocalizedOptimFormats.UNABLE_TO_SOLVE_SINGULAR_PROBLEM,
                                                    e);
            }
            // update the estimated parameters
            currentPoint = currentPoint.add(dX);

        }
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {
        return "SequentialGaussNewtonOptimizer{" +
               "decomposer=" + decomposer + '}';
    }

    /**
     * Compute the normal matrix, J<sup>T</sup>J.
     *
     * @param jacobian  the m by n jacobian matrix, J. Input.
     * @param residuals the m by 1 residual vector, r. Input.
     * @return  the n by n normal matrix and the n by 1 J<sup>Tr</sup> vector.
     */
    private static Pair<RealMatrix, RealVector>
        computeNormalMatrix(final RealMatrix jacobian,
                            final RealVector residuals) {
        // since the normal matrix is symmetric, we only need to compute half of
        // it.
        final int nR = jacobian.getRowDimension();
        final int nC = jacobian.getColumnDimension();
        // allocate space for return values
        final RealMatrix normal = MatrixUtils.createRealMatrix(nC, nC);
        final RealVector jTr = new ArrayRealVector(nC);
        // for each measurement
        for (int i = 0; i < nR; ++i) {
            // compute JTr for measurement i
            for (int j = 0; j < nC; j++) {
                jTr.setEntry(j,
                             jTr.getEntry(j) +
                                residuals.getEntry(i) *
                                               jacobian.getEntry(i, j));
            }

            // add the the contribution to the normal matrix for measurement i
            for (int k = 0; k < nC; ++k) {
                // only compute the upper triangular part
                for (int l = k; l < nC; ++l) {
                    normal
                        .setEntry(k, l,
                                  normal.getEntry(k,
                                                  l) +
                                        jacobian.getEntry(i, k) *
                                                       jacobian.getEntry(i, l));
                }
            }
        }
        // copy the upper triangular part to the lower triangular part.
        for (int i = 0; i < nC; i++) {
            for (int j = 0; j < i; j++) {
                normal.setEntry(i, j, normal.getEntry(j, i));
            }
        }
        return new Pair<RealMatrix, RealVector>(normal, jTr);
    }

    /** Combine Jacobian matrices
     * @param oldJacobian old Jacobian matrix
     * @param newJacobian new Jacobian matrix
     * @return combined Jacobian matrix
     */
    private static RealMatrix combineJacobians(final RealMatrix oldJacobian,
                                               final RealMatrix newJacobian) {
        final int oldRowDimension    = oldJacobian.getRowDimension();
        final int oldColumnDimension = oldJacobian.getColumnDimension();
        final RealMatrix jacobian =
                        MatrixUtils.createRealMatrix(oldRowDimension + newJacobian.getRowDimension(),
                                                     oldColumnDimension);
        jacobian.setSubMatrix(oldJacobian.getData(), 0,               0);
        jacobian.setSubMatrix(newJacobian.getData(), oldRowDimension, 0);
        return jacobian;
    }

    /** Combine residuals vectors
     * @param oldResiduals old residuals vector
     * @param newResiduals new residuals vector
     * @return combined residuals vector
     */
    private static RealVector combineResiduals(final RealVector oldResiduals,
                                               final RealVector newResiduals) {
        return oldResiduals.append(newResiduals);
    }

    /**
     * Container with an old and a new evaluation and combine both of them
     */
    private static class CombinedEvaluation extends AbstractEvaluation {

        /** Point of evaluation. */
        private final RealVector point;

        /** Derivative at point. */
        private final RealMatrix jacobian;

        /** Computed residuals. */
        private final RealVector residuals;

        /**
         * Create an {@link Evaluation} with no weights.
         *
         * @param oldEvaluation the old evaluation.
         * @param newEvaluation the new evaluation
         */
        private CombinedEvaluation(final Evaluation oldEvaluation,
                                   final Evaluation newEvaluation) {

            super(oldEvaluation.getResiduals().getDimension() +
                  newEvaluation.getResiduals().getDimension());

            this.point    = newEvaluation.getPoint();
            this.jacobian = combineJacobians(oldEvaluation.getJacobian(),
                                             newEvaluation.getJacobian());
            this.residuals = combineResiduals(oldEvaluation.getResiduals(),
                                              newEvaluation.getResiduals());
        }

        /** {@inheritDoc} */
        @Override
        public RealMatrix getJacobian() {
            return jacobian;
        }

        /** {@inheritDoc} */
        @Override
        public RealVector getPoint() {
            return point;
        }

        /** {@inheritDoc} */
        @Override
        public RealVector getResiduals() {
            return residuals;
        }

    }

}
