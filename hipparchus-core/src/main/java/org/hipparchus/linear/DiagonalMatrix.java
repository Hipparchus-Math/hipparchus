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
package org.hipparchus.linear;

import java.io.Serializable;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;

/**
 * Implementation of a diagonal matrix.
 *
 */
public class DiagonalMatrix extends AbstractRealMatrix
    implements Serializable {
    /** Serializable version identifier. */
    private static final long serialVersionUID = 20121229L;
    /** Entries of the diagonal. */
    private final double[] data;

    /**
     * Creates a matrix with the supplied dimension.
     *
     * @param dimension Number of rows and columns in the new matrix.
     * @throws MathIllegalArgumentException if the dimension is
     * not positive.
     */
    public DiagonalMatrix(final int dimension)
        throws MathIllegalArgumentException {
        super(dimension, dimension);
        data = new double[dimension];
    }

    /**
     * Creates a matrix using the input array as the underlying data.
     * <br>
     * The input array is copied, not referenced.
     *
     * @param d Data for the new matrix.
     */
    public DiagonalMatrix(final double[] d) {
        this(d, true);
    }

    /**
     * Creates a matrix using the input array as the underlying data.
     * <br>
     * If an array is created specially in order to be embedded in a
     * this instance and not used directly, the {@code copyArray} may be
     * set to {@code false}.
     * This will prevent the copying and improve performance as no new
     * array will be built and no data will be copied.
     *
     * @param d Data for new matrix.
     * @param copyArray if {@code true}, the input array will be copied,
     * otherwise it will be referenced.
     * @exception NullArgumentException if d is null
     */
    public DiagonalMatrix(final double[] d, final boolean copyArray)
        throws NullArgumentException {
        MathUtils.checkNotNull(d);
        data = copyArray ? d.clone() : d;
    }

    /**
     * {@inheritDoc}
     *
     * @throws MathIllegalArgumentException if the requested dimensions are not equal.
     */
    @Override
    public RealMatrix createMatrix(final int rowDimension,
                                   final int columnDimension)
        throws MathIllegalArgumentException {
        if (rowDimension != columnDimension) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   rowDimension, columnDimension);
        }

        return new DiagonalMatrix(rowDimension);
    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix copy() {
        return new DiagonalMatrix(data);
    }

    /**
     * Compute the sum of {@code this} and {@code m}.
     *
     * @param m Matrix to be added.
     * @return {@code this + m}.
     * @throws MathIllegalArgumentException if {@code m} is not the same
     * size as {@code this}.
     */
    public DiagonalMatrix add(final DiagonalMatrix m)
        throws MathIllegalArgumentException {
        // Safety check.
        MatrixUtils.checkAdditionCompatible(this, m);

        final int dim = getRowDimension();
        final double[] outData = new double[dim];
        for (int i = 0; i < dim; i++) {
            outData[i] = data[i] + m.data[i];
        }

        return new DiagonalMatrix(outData, false);
    }

    /**
     * Returns {@code this} minus {@code m}.
     *
     * @param m Matrix to be subtracted.
     * @return {@code this - m}
     * @throws MathIllegalArgumentException if {@code m} is not the same
     * size as {@code this}.
     */
    public DiagonalMatrix subtract(final DiagonalMatrix m)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSubtractionCompatible(this, m);

        final int dim = getRowDimension();
        final double[] outData = new double[dim];
        for (int i = 0; i < dim; i++) {
            outData[i] = data[i] - m.data[i];
        }

        return new DiagonalMatrix(outData, false);
    }

    /**
     * Returns the result of postmultiplying {@code this} by {@code m}.
     *
     * @param m matrix to postmultiply by
     * @return {@code this * m}
     * @throws MathIllegalArgumentException if
     * {@code columnDimension(this) != rowDimension(m)}
     */
    public DiagonalMatrix multiply(final DiagonalMatrix m)
        throws MathIllegalArgumentException {
        MatrixUtils.checkMultiplicationCompatible(this, m);

        final int dim = getRowDimension();
        final double[] outData = new double[dim];
        for (int i = 0; i < dim; i++) {
            outData[i] = data[i] * m.data[i];
        }

        return new DiagonalMatrix(outData, false);
    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix multiply(final RealMatrix m)
        throws MathIllegalArgumentException {
        if (m instanceof DiagonalMatrix) {
            return multiply((DiagonalMatrix) m);
        } else {
            MatrixUtils.checkMultiplicationCompatible(this, m);
            final RealMatrix product = m.createMatrix(m.getRowDimension(), m.getColumnDimension());
            product.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                /** {@inheritDoc} */
                @Override
                public double visit(int row, int column, double value) {
                    return data[row] * m.getEntry(row, column);
                }
            });
            return product;
        }
    }

    /**
     * Returns the result of postmultiplying {@code this} by {@code m^T}.
     * @param m matrix to first transpose and second postmultiply by
     * @return {@code this * m}
     * @throws MathIllegalArgumentException if
     * {@code columnDimension(this) != columnDimension(m)}
     * @since 1.3
     */
    public DiagonalMatrix multiplyTransposed(final DiagonalMatrix m)
        throws MathIllegalArgumentException {
        // transposition is no-op for diagonal matrices
        return multiply(m);
    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix multiplyTransposed(final RealMatrix m)
        throws MathIllegalArgumentException {
        if (m instanceof DiagonalMatrix) {
            return multiplyTransposed((DiagonalMatrix) m);
        } else {
            MatrixUtils.checkSameColumnDimension(this, m);
            final RealMatrix product = m.createMatrix(m.getColumnDimension(), m.getRowDimension());
            product.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                /** {@inheritDoc} */
                @Override
                public double visit(int row, int column, double value) {
                    return data[row] * m.getEntry(column, row);
                }
            });
            return product;
        }
    }

    /**
     * Returns the result of postmultiplying {@code this^T} by {@code m}.
     * @param m matrix to first transpose and second postmultiply by
     * @return {@code this^T * m}
     * @throws MathIllegalArgumentException if
     * {@code columnDimension(this) != columnDimension(m)}
     * @since 1.3
     */
    public DiagonalMatrix transposeMultiply(final DiagonalMatrix m)
        throws MathIllegalArgumentException {
        // transposition is no-op for diagonal matrices
        return multiply(m);
    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix transposeMultiply(final RealMatrix m) {
        if (m instanceof DiagonalMatrix) {
            return transposeMultiply((DiagonalMatrix) m);
        } else {
            // transposition is no-op for diagonal matrices
            return multiply(m);
        }
    }

    /** {@inheritDoc} */
    @Override
    public double[][] getData() {
        final int dim = getRowDimension();
        final double[][] out = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            out[i][i] = data[i];
        }

        return out;
    }

    /**
     * Gets a reference to the underlying data array.
     *
     * @return 1-dimensional array of entries.
     */
    public double[] getDataRef() {
        return data; // NOPMD - returning an internal array is intentional and documented here
    }

    /** {@inheritDoc} */
    @Override
    public double getEntry(final int row, final int column)
        throws MathIllegalArgumentException {
        MatrixUtils.checkMatrixIndex(this, row, column);
        return row == column ? data[row] : 0;
    }

    /** {@inheritDoc}
     * @throws MathIllegalArgumentException if {@code row != column} and value is non-zero.
     */
    @Override
    public void setEntry(final int row, final int column, final double value)
        throws MathIllegalArgumentException {
        if (row == column) {
            MatrixUtils.checkRowIndex(this, row);
            data[row] = value;
        } else {
            ensureZero(value);
        }
    }

    /** {@inheritDoc}
     * @throws MathIllegalArgumentException if {@code row != column} and increment is non-zero.
     */
    @Override
    public void addToEntry(final int row,
                           final int column,
                           final double increment)
        throws MathIllegalArgumentException {
        if (row == column) {
            MatrixUtils.checkRowIndex(this, row);
            data[row] += increment;
        } else {
            ensureZero(increment);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void multiplyEntry(final int row,
                              final int column,
                              final double factor)
        throws MathIllegalArgumentException {
        // we don't care about non-diagonal elements for multiplication
        if (row == column) {
            MatrixUtils.checkRowIndex(this, row);
            data[row] *= factor;
        }
    }

    /** {@inheritDoc} */
    @Override
    public int getRowDimension() {
        return data.length;
    }

    /** {@inheritDoc} */
    @Override
    public int getColumnDimension() {
        return data.length;
    }

    /** {@inheritDoc} */
    @Override
    public double[] operate(final double[] v)
        throws MathIllegalArgumentException {
        return multiply(new DiagonalMatrix(v, false)).getDataRef();
    }

    /** {@inheritDoc} */
    @Override
    public double[] preMultiply(final double[] v)
        throws MathIllegalArgumentException {
        return operate(v);
    }

    /** {@inheritDoc} */
    @Override
    public RealVector preMultiply(final RealVector v) throws MathIllegalArgumentException {
        final double[] vectorData;
        if (v instanceof ArrayRealVector) {
            vectorData = ((ArrayRealVector) v).getDataRef();
        } else {
            vectorData = v.toArray();
        }
        return MatrixUtils.createRealVector(preMultiply(vectorData));
    }

    /** Ensure a value is zero.
     * @param value value to check
     * @exception MathIllegalArgumentException if value is not zero
     */
    private void ensureZero(final double value) throws MathIllegalArgumentException {
        if (!Precision.equals(0.0, value, 1)) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE,
                                                   FastMath.abs(value), 0);
        }
    }

    /**
     * Computes the inverse of this diagonal matrix.
     * <p>
     * Note: this method will use a singularity threshold of 0,
     * use {@link #inverse(double)} if a different threshold is needed.
     *
     * @return the inverse of {@code m}
     * @throws MathIllegalArgumentException if the matrix is singular
     */
    public DiagonalMatrix inverse() throws MathIllegalArgumentException {
        return inverse(0);
    }

    /**
     * Computes the inverse of this diagonal matrix.
     *
     * @param threshold Singularity threshold.
     * @return the inverse of {@code m}
     * @throws MathIllegalArgumentException if the matrix is singular
     */
    public DiagonalMatrix inverse(double threshold) throws MathIllegalArgumentException {
        if (isSingular(threshold)) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SINGULAR_MATRIX);
        }

        final double[] result = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            result[i] = 1.0 / data[i];
        }
        return new DiagonalMatrix(result, false);
    }

    /** Returns whether this diagonal matrix is singular, i.e. any diagonal entry
     * is equal to {@code 0} within the given threshold.
     *
     * @param threshold Singularity threshold.
     * @return {@code true} if the matrix is singular, {@code false} otherwise
     */
    public boolean isSingular(double threshold) {
        for (double datum : data) {
            if (Precision.equals(datum, 0.0, threshold)) {
                return true;
            }
        }
        return false;
    }
}
