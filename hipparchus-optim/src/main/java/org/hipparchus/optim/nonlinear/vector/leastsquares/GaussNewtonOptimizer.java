/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.optim.nonlinear.vector.leastsquares;

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
     * The singularity threshold for matrix decompositions. Determines when a {@link
     * MathIllegalStateException} is thrown. The current value was the default value for {@link
     * org.hipparchus.linear.LUDecomposition}.
     */
    private static final double SINGULARITY_THRESHOLD = 1e-11;

    /** Decomposer */
    private final MatrixDecomposer decomposer;

    /** Indicates if normal equations should be formed explicitly. */
    private final boolean formNormalEquations;

    /**
     * Creates a Gauss Newton optimizer.
     * <p>
     * The default for the algorithm is to use QR decomposition and not
     * form normal equations.
     * </p>
     */
    public GaussNewtonOptimizer() {
        this(new QRDecomposer(SINGULARITY_THRESHOLD), false);
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
        this.decomposer          = decomposer;
        this.formNormalEquations = formNormalEquations;
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
