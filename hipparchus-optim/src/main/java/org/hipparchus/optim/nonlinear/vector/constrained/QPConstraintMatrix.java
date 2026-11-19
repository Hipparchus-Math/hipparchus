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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.AbstractRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Read-only, zero-copy row-oriented view of all linear constraints of a QP.
 *
 * <p>The represented matrix has the form:</p>
 *
 * <pre>
 *     A = [ Ae          ]
 *         [ Ai          ]
 *         [ Ab          ]
 *         [ -Ab         ]
 * </pre>
 *
 * <p>where each constraint gradient is stored as a row. The class keeps only
 * references to the original equality, inequality and bounded-constraint
 * matrices. In particular, the upper-bound block {@code -Ab} is generated on
 * demand and is never stored.</p>
 *
 * <p>The global row ordering is:</p>
 * <ol>
 *   <li>equality constraints,</li>
 *   <li>inequality constraints,</li>
 *   <li>lower bounded constraints,</li>
 *   <li>upper bounded constraints.</li>
 * </ol>
 *
 * <p>The matrix is intentionally read-only. Calling {@link #copy()} creates a
 * physical independent matrix, but normal solver operations such as
 * {@link #getEntry(int, int)}, {@link #getRowVector(int)} and
 * {@link #operate(RealVector)} do not materialize the global matrix.</p>
 */
public final class QPConstraintMatrix extends AbstractRealMatrix {

    /** Number of primal variables, i.e. matrix columns. */
    private final int variableCount;

    /** Equality matrix, one equality per row. */
    private final RealMatrix equalityMatrix;

    /** Inequality matrix, one inequality per row. */
    private final RealMatrix inequalityMatrix;

    /** Bounded-constraint matrix, one bounded expression per row. */
    private final RealMatrix boundedMatrix;

    /** Number of equality rows. */
    private final int equalityCount;

    /** Number of inequality rows. */
    private final int inequalityCount;

    /** Number of bounded-expression rows. */
    private final int boundedCount;

    /** First global row of the inequality block. */
    private final int inequalityOffset;

    /** First global row of the lower-bound block. */
    private final int lowerBoundOffset;

    /** First global row of the upper-bound block. */
    private final int upperBoundOffset;

    /** Total number of represented constraints. */
    private final int constraintCount;

    /**
     * Build a row-oriented virtual QP constraint matrix.
     *
     * @param variableCount number of primal variables
     * @param equalityMatrix equality matrix, or {@code null}
     * @param inequalityMatrix inequality matrix, or {@code null}
     * @param boundedMatrix bounded-expression matrix, or {@code null}
     */
    public QPConstraintMatrix(final int variableCount,
                              final RealMatrix equalityMatrix,
                              final RealMatrix inequalityMatrix,
                              final RealMatrix boundedMatrix) {

        if (variableCount < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NEGATIVE_VALUE, variableCount);
        }

        this.variableCount = variableCount;
        this.equalityMatrix = equalityMatrix;
        this.inequalityMatrix = inequalityMatrix;
        this.boundedMatrix = boundedMatrix;

        checkColumnDimension(equalityMatrix,   variableCount);
        checkColumnDimension(inequalityMatrix, variableCount);
        checkColumnDimension(boundedMatrix,    variableCount);

        this.equalityCount = equalityMatrix == null ? 0 : equalityMatrix.getRowDimension();
        this.inequalityCount = inequalityMatrix == null ? 0 : inequalityMatrix.getRowDimension();
        this.boundedCount = boundedMatrix == null ? 0 : boundedMatrix.getRowDimension();

        this.inequalityOffset = equalityCount;
        this.lowerBoundOffset = equalityCount + inequalityCount;
        this.upperBoundOffset = lowerBoundOffset + boundedCount;
        this.constraintCount = upperBoundOffset + boundedCount;
    }

    /** Validate the number of columns of a backing block.
     * @param matrix constraint matrix
     * @param expected expected number of columns
     */
    private static void checkColumnDimension(final RealMatrix matrix,
                                             final int expected) {
        if (matrix != null && matrix.getColumnDimension() != expected) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   matrix.getColumnDimension(), expected);
        }
    }

    /** {@inheritDoc} */
    @Override
    public int getRowDimension() {
        return constraintCount;
    }

    /** {@inheritDoc} */
    @Override
    public int getColumnDimension() {
        return variableCount;
    }

    /**
     * Get a coefficient without constructing the global matrix.
     *
     * @param row global constraint index
     * @param column primal-variable index
     * @return coefficient A[row, column]
     */
    @Override
    public double getEntry(final int row, final int column) {
        MatrixUtils.checkMatrixIndex(this, row, column);

        if (row < inequalityOffset) {
            return equalityMatrix.getEntry(row, column);
        }
        if (row < lowerBoundOffset) {
            return inequalityMatrix.getEntry(row - inequalityOffset, column);
        }
        if (row < upperBoundOffset) {
            return boundedMatrix.getEntry(row - lowerBoundOffset, column);
        }
        return -boundedMatrix.getEntry(row - upperBoundOffset, column);
    }

    /**
     * Get a constraint gradient as a row vector.
     *
     * <p>The backing Hipparchus matrix already returns an independent vector.
     * For an upper bound that same newly allocated vector is negated in place,
     * avoiding an additional vector allocation.</p>
     *
     * @param row global constraint index
     * @return constraint gradient
     */
    @Override
    public RealVector getRowVector(final int row) {
        MatrixUtils.checkRowIndex(this, row);

        if (row < inequalityOffset) {
            return equalityMatrix.getRowVector(row);
        }
        if (row < lowerBoundOffset) {
            return inequalityMatrix.getRowVector(row - inequalityOffset);
        }
        if (row < upperBoundOffset) {
            return boundedMatrix.getRowVector(row - lowerBoundOffset);
        }

        return boundedMatrix.getRowVector(row - upperBoundOffset)
                            .mapMultiplyToSelf(-1.0);
    }

    /**
     * Get a constraint gradient as an array.
     *
     * @param row global constraint index
     * @return independent coefficient array
     */
    @Override
    public double[] getRow(final int row) {
        MatrixUtils.checkRowIndex(this, row);

        if (row < inequalityOffset) {
            return equalityMatrix.getRow(row);
        }
        if (row < lowerBoundOffset) {
            return inequalityMatrix.getRow(row - inequalityOffset);
        }
        if (row < upperBoundOffset) {
            return boundedMatrix.getRow(row - lowerBoundOffset);
        }

        final double[] data = boundedMatrix.getRow(row - upperBoundOffset);
        for (int i = 0; i < data.length; i++) {
            data[i] = -data[i];
        }
        return data;
    }

    /**
     * Evaluate all linear constraint expressions A x.
     *
     * <p>Each physical backing block is multiplied only once. The product of
     * the bounded matrix is reused for both the lower and upper blocks, with
     * opposite signs.</p>
     *
     * @param x primal vector
     * @return A x
     */
    @Override
    public RealVector operate(final RealVector x) {
        final double[] result = new double[constraintCount];

        if (equalityCount > 0) {
            copy(equalityMatrix.operate(x), result, 0, 1.0);
        }

        if (inequalityCount > 0) {
            copy(inequalityMatrix.operate(x), result, inequalityOffset, 1.0);
        }

        if (boundedCount > 0) {
            final RealVector boundedValue = boundedMatrix.operate(x);
            copy(boundedValue, result, lowerBoundOffset, 1.0);
            copy(boundedValue, result, upperBoundOffset, -1.0);
        }

        return new ArrayRealVector(result, false);
    }

    /** Copy a vector into an array, optionally changing its sign.
     * @param source      source vector
     * @param destination destination array
     * @param offset      offset in destination array
     * @param scale       scaling factor applied to all components
     */
    private static void copy(final RealVector source,
                             final double[] destination,
                             final int offset,
                             final double scale) {
        for (int i = 0; i < source.getDimension(); i++) {
            destination[offset + i] = scale * source.getEntry(i);
        }
    }


    /**
     * This view is read-only.
     */
    @Override
    public void setEntry(final int row, final int column, final double value) {
        throw new UnsupportedOperationException("QPConstraintMatrix is read-only");
    }

    /**
     * Create an independent physical matrix for results of generic operations.
     */
    @Override
    public RealMatrix createMatrix(final int rowDimension, final int columnDimension) {
        return MatrixUtils.createRealMatrix(rowDimension, columnDimension);
    }

    /**
     * Materialize an independent deep copy.
     */
    @Override
    public RealMatrix copy() {
        final RealMatrix result = MatrixUtils.createRealMatrix(constraintCount, variableCount);
        for (int row = 0; row < constraintCount; row++) {
            result.setRow(row, getRow(row));
        }
        return result;
    }
}
