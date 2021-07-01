/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.optim.nonlinear.vector.leastsquares;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.linear.ArrayRealVector;
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

    /** Old evaluation previously computed. */
    private final Evaluation oldEvaluation;

    /** Old jacobian previously computed. */
    private final RealMatrix oldLhs;

    /** Old residuals previously computed. */
    private final RealVector oldRhs;

    /**
     * Create a sequential Gauss Newton optimizer.
     * <p>
     * The default for the algorithm is to solve the normal equations
     * J<sup>T</sup>Jx=J<sup>T</sup>r using QR decomposition.
     * </p>
     *
     * @param evaluation old evaluation previously computed, null if there are no previous evaluations.
     */
    public SequentialGaussNewtonOptimizer(final Evaluation evaluation) {
        this(new QRDecomposer(SINGULARITY_THRESHOLD), evaluation);
    }

    /**
     * Create a sequential Gauss Newton optimizer that uses the given matrix
     * decomposition algorithm to solve the normal equations.
     * <p>
     * The {@code decomposer} is used to solve J<sup>T</sup>Jx=J<sup>T</sup>r.
     * </p>
     *
     * @param decomposer the decomposition algorithm to use.
     * @param evaluation old evaluation previously computed, null if there are no previous evaluations.
     */
    public SequentialGaussNewtonOptimizer(final MatrixDecomposer decomposer,
                                          final Evaluation evaluation) {
        this.decomposer = decomposer;
        this.oldEvaluation = evaluation;
        if (evaluation == null) {
            this.oldLhs = null;
            this.oldRhs = null;
        } else {
            final Pair<RealMatrix, RealVector> normalEquation =
                            computeNormalMatrix(evaluation.getJacobian(), evaluation.getResiduals());
            // solve the linearized least squares problem
            this.oldLhs = normalEquation.getFirst();
            this.oldRhs = normalEquation.getSecond();
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
     * Get the previous evaluation used by the optimizer.
     *
     * @return the previous evaluation.
     */
    public Evaluation getOldEvaluation() {
        return oldEvaluation;
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

            final Pair<RealMatrix, RealVector> normalEquation =
                computeNormalMatrix(weightedJacobian, currentResiduals);
            // solve the linearized least squares problem

            final RealMatrix lhs = oldLhs == null ?
                                   normalEquation.getFirst() :
                                   normalEquation.getFirst().add(oldLhs); // left hand side
            final RealVector rhs = oldRhs == null ?
                                   normalEquation.getSecond() :
                                   normalEquation.getSecond().add(oldRhs); // right hand side

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

            final RealMatrix oldJacobian = oldEvaluation.getJacobian();
            final RealMatrix newJacobian = newEvaluation.getJacobian();

            final int oldRowDimension    = oldJacobian.getRowDimension();
            final int oldColumnDimension = oldJacobian.getColumnDimension();

            this.jacobian = MatrixUtils.createRealMatrix(oldRowDimension + newJacobian.getRowDimension(),
                                                         oldColumnDimension);
            jacobian.setSubMatrix(oldJacobian.getData(), 0,               0);
            jacobian.setSubMatrix(newJacobian.getData(), oldRowDimension, 0);

            this.point     = newEvaluation.getPoint();
            this.residuals = oldEvaluation.getResiduals().append(newEvaluation.getResiduals());
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
