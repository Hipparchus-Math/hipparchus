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

package org.hipparchus.linear;

import java.io.Serializable;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Implementation of {@link RealMatrix} using a {@code double[][]} array to
 * store entries.
 *
 */
public class Array2DRowRealMatrix extends AbstractRealMatrix implements Serializable {
    /** Serializable version identifier. */
    private static final long serialVersionUID = -1067294169172445528L;

    /** Entries of the matrix. */
    private double data[][];

    /**
     * Creates a matrix with no data
     */
    public Array2DRowRealMatrix() {
        // This constructor is intentionally empty. Nothing special is needed here.
    }

    /**
     * Create a new RealMatrix with the supplied row and column dimensions.
     *
     * @param rowDimension Number of rows in the new matrix.
     * @param columnDimension Number of columns in the new matrix.
     * @throws MathIllegalArgumentException if the row or column dimension is
     * not positive.
     */
    public Array2DRowRealMatrix(final int rowDimension,
                                final int columnDimension)
        throws MathIllegalArgumentException {
        super(rowDimension, columnDimension);
        data = new double[rowDimension][columnDimension];
    }

    /**
     * Create a new {@code RealMatrix} using the input array as the underlying
     * data array.
     * <p>The input array is copied, not referenced. This constructor has
     * the same effect as calling {@link #Array2DRowRealMatrix(double[][], boolean)}
     * with the second argument set to {@code true}.</p>
     *
     * @param d Data for the new matrix.
     * @throws MathIllegalArgumentException if {@code d} is not rectangular.
     * @throws MathIllegalArgumentException if {@code d} row or column dimension is zero.
     * @throws NullArgumentException if {@code d} is {@code null}.
     * @see #Array2DRowRealMatrix(double[][], boolean)
     */
    public Array2DRowRealMatrix(final double[][] d)
        throws MathIllegalArgumentException, NullArgumentException {
        copyIn(d);
    }

    /**
     * Create a new RealMatrix using the input array as the underlying
     * data array.
     * If an array is built specially in order to be embedded in a
     * RealMatrix and not used directly, the {@code copyArray} may be
     * set to {@code false}. This will prevent the copying and improve
     * performance as no new array will be built and no data will be copied.
     *
     * @param d Data for new matrix.
     * @param copyArray if {@code true}, the input array will be copied,
     * otherwise it will be referenced.
     * @throws MathIllegalArgumentException if {@code d} is not rectangular.
     * @throws MathIllegalArgumentException if {@code d} row or column dimension is zero.
     * @throws NullArgumentException if {@code d} is {@code null}.
     * @see #Array2DRowRealMatrix(double[][])
     */
    public Array2DRowRealMatrix(final double[][] d, final boolean copyArray) // NOPMD - array copy is taken care of by parameter
        throws MathIllegalArgumentException,
        NullArgumentException {
        if (copyArray) {
            copyIn(d);
        } else {
            if (d == null) {
                throw new NullArgumentException();
            }
            final int nRows = d.length;
            if (nRows == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.AT_LEAST_ONE_ROW);
            }
            final int nCols = d[0].length;
            if (nCols == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.AT_LEAST_ONE_COLUMN);
            }
            for (int r = 1; r < nRows; r++) {
                if (d[r].length != nCols) {
                    throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                           d[r].length, nCols);
                }
            }
            data = d;
        }
    }

    /**
     * Create a new (column) RealMatrix using {@code v} as the
     * data for the unique column of the created matrix.
     * The input array is copied.
     *
     * @param v Column vector holding data for new matrix.
     */
    public Array2DRowRealMatrix(final double[] v) {
        final int nRows = v.length;
        data = new double[nRows][1];
        for (int row = 0; row < nRows; row++) {
            data[row][0] = v[row];
        }
    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix createMatrix(final int rowDimension,
                                   final int columnDimension)
        throws MathIllegalArgumentException {
        return new Array2DRowRealMatrix(rowDimension, columnDimension);
    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix copy() {
        return new Array2DRowRealMatrix(copyOut(), false);
    }

    /**
     * Compute the sum of {@code this} and {@code m}.
     *
     * @param m Matrix to be added.
     * @return {@code this + m}.
     * @throws MathIllegalArgumentException if {@code m} is not the same
     * size as {@code this}.
     */
    public Array2DRowRealMatrix add(final Array2DRowRealMatrix m)
        throws MathIllegalArgumentException {
        // Safety check.
        MatrixUtils.checkAdditionCompatible(this, m);

        final int rowCount    = getRowDimension();
        final int columnCount = getColumnDimension();
        final double[][] outData = new double[rowCount][columnCount];
        for (int row = 0; row < rowCount; row++) {
            final double[] dataRow    = data[row];
            final double[] mRow       = m.data[row];
            final double[] outDataRow = outData[row];
            for (int col = 0; col < columnCount; col++) {
                outDataRow[col] = dataRow[col] + mRow[col];
            }
        }

        return new Array2DRowRealMatrix(outData, false);
    }

    /**
     * Returns {@code this} minus {@code m}.
     *
     * @param m Matrix to be subtracted.
     * @return {@code this - m}
     * @throws MathIllegalArgumentException if {@code m} is not the same
     * size as {@code this}.
     */
    public Array2DRowRealMatrix subtract(final Array2DRowRealMatrix m)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSubtractionCompatible(this, m);

        final int rowCount    = getRowDimension();
        final int columnCount = getColumnDimension();
        final double[][] outData = new double[rowCount][columnCount];
        for (int row = 0; row < rowCount; row++) {
            final double[] dataRow    = data[row];
            final double[] mRow       = m.data[row];
            final double[] outDataRow = outData[row];
            for (int col = 0; col < columnCount; col++) {
                outDataRow[col] = dataRow[col] - mRow[col];
            }
        }

        return new Array2DRowRealMatrix(outData, false);
    }

    /**
     * Returns the result of postmultiplying {@code this} by {@code m}.
     *
     * @param m matrix to postmultiply by
     * @return {@code this * m}
     * @throws MathIllegalArgumentException if
     * {@code columnDimension(this) != rowDimension(m)}
     */
    public Array2DRowRealMatrix multiply(final Array2DRowRealMatrix m)
        throws MathIllegalArgumentException {
        MatrixUtils.checkMultiplicationCompatible(this, m);

        final int nRows = this.getRowDimension();
        final int nCols = m.getColumnDimension();
        final int nSum = this.getColumnDimension();

        final double[][] outData = new double[nRows][nCols];
        // Will hold a column of "m".
        final double[] mCol = new double[nSum];
        final double[][] mData = m.data;

        // Multiply.
        for (int col = 0; col < nCols; col++) {
            // Copy all elements of column "col" of "m" so that
            // will be in contiguous memory.
            for (int mRow = 0; mRow < nSum; mRow++) {
                mCol[mRow] = mData[mRow][col];
            }

            for (int row = 0; row < nRows; row++) {
                final double[] dataRow = data[row];
                double sum = 0;
                for (int i = 0; i < nSum; i++) {
                    sum += dataRow[i] * mCol[i];
                }
                outData[row][col] = sum;
            }
        }

        return new Array2DRowRealMatrix(outData, false);
    }

    /**
     * Returns the result of postmultiplying {@code this} by {@code m^T}.
     * @param m matrix to first transpose and second postmultiply by
     * @return {@code this * m^T}
     * @throws MathIllegalArgumentException if
     * {@code columnDimension(this) != columnDimension(m)}
     * @since 1.3
     */
    public RealMatrix multiplyTransposed(final Array2DRowRealMatrix m)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSameColumnDimension(this, m);

        final int nRows = this.getRowDimension();
        final int nCols = m.getRowDimension();
        final int nSum  = this.getColumnDimension();

        final RealMatrix out = MatrixUtils.createRealMatrix(nRows, nCols);
        final double[][] mData   = m.data;

        // Multiply.
        for (int col = 0; col < nCols; col++) {
            for (int row = 0; row < nRows; row++) {
                final double[] dataRow = data[row];
                final double[] mRow    = mData[col];
                double sum = 0;
                for (int i = 0; i < nSum; i++) {
                    sum += dataRow[i] * mRow[i];
                }
                out.setEntry(row, col, sum);
            }
        }

        return out;

    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix multiplyTransposed(final RealMatrix m) {
        if (m instanceof Array2DRowRealMatrix) {
            return multiplyTransposed((Array2DRowRealMatrix) m);
        } else {
            MatrixUtils.checkSameColumnDimension(this, m);

            final int nRows = this.getRowDimension();
            final int nCols = m.getRowDimension();
            final int nSum  = this.getColumnDimension();

            final RealMatrix out = MatrixUtils.createRealMatrix(nRows, nCols);

            // Multiply.
            for (int col = 0; col < nCols; col++) {
                for (int row = 0; row < nRows; row++) {
                    final double[] dataRow = data[row];
                    double sum = 0;
                    for (int i = 0; i < nSum; i++) {
                        sum += dataRow[i] * m.getEntry(col, i);
                    }
                    out.setEntry(row, col, sum);
                }
            }

            return out;

        }
    }

    /**
     * Returns the result of postmultiplying {@code this^T} by {@code m}.
     * @param m matrix to postmultiply by
     * @return {@code this^T * m}
     * @throws MathIllegalArgumentException if
     * {@code columnDimension(this) != columnDimension(m)}
     * @since 1.3
     */
    public RealMatrix transposeMultiply(final Array2DRowRealMatrix m)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSameRowDimension(this, m);

        final int nRows = this.getColumnDimension();
        final int nCols = m.getColumnDimension();
        final int nSum  = this.getRowDimension();

        final RealMatrix out = MatrixUtils.createRealMatrix(nRows, nCols);
        final double[][] mData   = m.data;

        // Multiply.
        for (int k = 0; k < nSum; k++) {
            final double[] dataK = data[k];
            final double[] mK    = mData[k];
            for (int row = 0; row < nRows; row++) {
                final double dataIRow = dataK[row];
                for (int col = 0; col < nCols; col++) {
                    out.addToEntry(row, col, dataIRow * mK[col]);
                }
            }
        }

        return out;

    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix transposeMultiply(final RealMatrix m) {
        if (m instanceof Array2DRowRealMatrix) {
            return transposeMultiply((Array2DRowRealMatrix) m);
        } else {
            MatrixUtils.checkSameRowDimension(this, m);

            final int nRows = this.getColumnDimension();
            final int nCols = m.getColumnDimension();
            final int nSum  = this.getRowDimension();

            final RealMatrix out = MatrixUtils.createRealMatrix(nRows, nCols);

            // Multiply.
            for (int k = 0; k < nSum; k++) {
                final double[] dataK = data[k];
                for (int row = 0; row < nRows; row++) {
                    final double dataIRow = dataK[row];
                    for (int col = 0; col < nCols; col++) {
                        out.addToEntry(row, col, dataIRow * m.getEntry(k, col));
                    }
                }
            }

            return out;

        }
    }

    /** {@inheritDoc} */
    @Override
    public double[][] getData() {
        return copyOut();
    }

    /**
     * Get a reference to the underlying data array.
     *
     * @return 2-dimensional array of entries.
     */
    public double[][] getDataRef() {
        return data; // NOPMD - returning an internal array is intentional and documented here
    }

    /** {@inheritDoc} */
    @Override
    public void setSubMatrix(final double[][] subMatrix, final int row,
                             final int column)
        throws MathIllegalArgumentException, NullArgumentException {
        if (data == null) {
            if (row > 0) {
                throw new MathIllegalStateException(LocalizedCoreFormats.FIRST_ROWS_NOT_INITIALIZED_YET, row);
            }
            if (column > 0) {
                throw new MathIllegalStateException(LocalizedCoreFormats.FIRST_COLUMNS_NOT_INITIALIZED_YET, column);
            }
            MathUtils.checkNotNull(subMatrix);
            final int nRows = subMatrix.length;
            if (nRows == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.AT_LEAST_ONE_ROW);
            }

            final int nCols = subMatrix[0].length;
            if (nCols == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.AT_LEAST_ONE_COLUMN);
            }
            data = new double[subMatrix.length][nCols];
            for (int i = 0; i < data.length; ++i) {
                if (subMatrix[i].length != nCols) {
                    throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                           subMatrix[i].length, nCols);
                }
                System.arraycopy(subMatrix[i], 0, data[i + row], column, nCols);
            }
        } else {
            super.setSubMatrix(subMatrix, row, column);
        }

    }

    /** {@inheritDoc} */
    @Override
    public double getEntry(final int row, final int column)
        throws MathIllegalArgumentException {
        MatrixUtils.checkMatrixIndex(this, row, column);
        return data[row][column];
    }

    /** {@inheritDoc} */
    @Override
    public void setEntry(final int row, final int column, final double value)
        throws MathIllegalArgumentException {
        MatrixUtils.checkMatrixIndex(this, row, column);
        data[row][column] = value;
    }

    /** {@inheritDoc} */
    @Override
    public void addToEntry(final int row, final int column,
                           final double increment)
        throws MathIllegalArgumentException {
        MatrixUtils.checkMatrixIndex(this, row, column);
        data[row][column] += increment;
    }

    /** {@inheritDoc} */
    @Override
    public void multiplyEntry(final int row, final int column,
                              final double factor)
        throws MathIllegalArgumentException {
        MatrixUtils.checkMatrixIndex(this, row, column);
        data[row][column] *= factor;
    }

    /** {@inheritDoc} */
    @Override
    public int getRowDimension() {
        return (data == null) ? 0 : data.length;
    }

    /** {@inheritDoc} */
    @Override
    public int getColumnDimension() {
        return ((data == null) || (data[0] == null)) ? 0 : data[0].length;
    }

    /** {@inheritDoc} */
    @Override
    public double[] operate(final double[] v)
        throws MathIllegalArgumentException {
        final int nRows = this.getRowDimension();
        final int nCols = this.getColumnDimension();
        if (v.length != nCols) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   v.length, nCols);
        }
        final double[] out = new double[nRows];
        for (int row = 0; row < nRows; row++) {
            final double[] dataRow = data[row];
            double sum = 0;
            for (int i = 0; i < nCols; i++) {
                sum += dataRow[i] * v[i];
            }
            out[row] = sum;
        }
        return out;
    }

    /** {@inheritDoc} */
    @Override
    public double[] preMultiply(final double[] v)
        throws MathIllegalArgumentException {
        final int nRows = getRowDimension();
        final int nCols = getColumnDimension();
        if (v.length != nRows) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   v.length, nRows);
        }

        final double[] out = new double[nCols];
        for (int col = 0; col < nCols; ++col) {
            double sum = 0;
            for (int i = 0; i < nRows; ++i) {
                sum += data[i][col] * v[i];
            }
            out[col] = sum;
        }

        return out;

    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix getSubMatrix(final int startRow, final int endRow,
                                   final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
        final int rowCount = endRow - startRow + 1;
        final int columnCount = endColumn - startColumn + 1;
        final double[][] outData = new double[rowCount][columnCount];
        for (int i = 0; i < rowCount; ++i) {
            System.arraycopy(data[startRow + i], startColumn, outData[i], 0, columnCount);
        }

        Array2DRowRealMatrix subMatrix = new Array2DRowRealMatrix();
        subMatrix.data = outData;
        return subMatrix;
    }

    /** {@inheritDoc} */
    @Override
    public double walkInRowOrder(final RealMatrixChangingVisitor visitor) {
        final int rows    = getRowDimension();
        final int columns = getColumnDimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int i = 0; i < rows; ++i) {
            final double[] rowI = data[i];
            for (int j = 0; j < columns; ++j) {
                rowI[j] = visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public double walkInRowOrder(final RealMatrixPreservingVisitor visitor) {
        final int rows    = getRowDimension();
        final int columns = getColumnDimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int i = 0; i < rows; ++i) {
            final double[] rowI = data[i];
            for (int j = 0; j < columns; ++j) {
                visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public double walkInRowOrder(final RealMatrixChangingVisitor visitor,
                                 final int startRow, final int endRow,
                                 final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
        visitor.start(getRowDimension(), getColumnDimension(),
                      startRow, endRow, startColumn, endColumn);
        for (int i = startRow; i <= endRow; ++i) {
            final double[] rowI = data[i];
            for (int j = startColumn; j <= endColumn; ++j) {
                rowI[j] = visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public double walkInRowOrder(final RealMatrixPreservingVisitor visitor,
                                 final int startRow, final int endRow,
                                 final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
        visitor.start(getRowDimension(), getColumnDimension(),
                      startRow, endRow, startColumn, endColumn);
        for (int i = startRow; i <= endRow; ++i) {
            final double[] rowI = data[i];
            for (int j = startColumn; j <= endColumn; ++j) {
                visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public double walkInColumnOrder(final RealMatrixChangingVisitor visitor) {
        final int rows    = getRowDimension();
        final int columns = getColumnDimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int j = 0; j < columns; ++j) {
            for (int i = 0; i < rows; ++i) {
                final double[] rowI = data[i];
                rowI[j] = visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public double walkInColumnOrder(final RealMatrixPreservingVisitor visitor) {
        final int rows    = getRowDimension();
        final int columns = getColumnDimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int j = 0; j < columns; ++j) {
            for (int i = 0; i < rows; ++i) {
                visitor.visit(i, j, data[i][j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public double walkInColumnOrder(final RealMatrixChangingVisitor visitor,
                                    final int startRow, final int endRow,
                                    final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
        visitor.start(getRowDimension(), getColumnDimension(),
                      startRow, endRow, startColumn, endColumn);
        for (int j = startColumn; j <= endColumn; ++j) {
            for (int i = startRow; i <= endRow; ++i) {
                final double[] rowI = data[i];
                rowI[j] = visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public double walkInColumnOrder(final RealMatrixPreservingVisitor visitor,
                                    final int startRow, final int endRow,
                                    final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
        visitor.start(getRowDimension(), getColumnDimension(),
                      startRow, endRow, startColumn, endColumn);
        for (int j = startColumn; j <= endColumn; ++j) {
            for (int i = startRow; i <= endRow; ++i) {
                visitor.visit(i, j, data[i][j]);
            }
        }
        return visitor.end();
    }

    /**
     * Get a fresh copy of the underlying data array.
     *
     * @return a copy of the underlying data array.
     */
    private double[][] copyOut() {
        final int nRows = this.getRowDimension();
        final double[][] out = new double[nRows][this.getColumnDimension()];
        // can't copy 2-d array in one shot, otherwise get row references
        for (int i = 0; i < nRows; i++) {
            System.arraycopy(data[i], 0, out[i], 0, data[i].length);
        }
        return out;
    }

    /**
     * Replace data with a fresh copy of the input array.
     *
     * @param in Data to copy.
     * @throws MathIllegalArgumentException if the input array is empty.
     * @throws MathIllegalArgumentException if the input array is not rectangular.
     * @throws NullArgumentException if the input array is {@code null}.
     */
    private void copyIn(final double[][] in)
        throws MathIllegalArgumentException, NullArgumentException {
        setSubMatrix(in, 0, 0);
    }

    /** {@inheritDoc} */
    @Override
    public double[] getRow(final int row) throws MathIllegalArgumentException {
        MatrixUtils.checkRowIndex(this, row);
        final int nCols = getColumnDimension();
        final double[] out = new double[nCols];
        System.arraycopy(data[row], 0, out, 0, nCols);
        return out;
    }

    /** {@inheritDoc} */
    @Override
    public void setRow(final int row, final double[] array)
        throws MathIllegalArgumentException {
        MatrixUtils.checkRowIndex(this, row);
        final int nCols = getColumnDimension();
        if (array.length != nCols) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH_2x2,
                                                   1, array.length, 1, nCols);
        }
        System.arraycopy(array, 0, data[row], 0, nCols);
    }

    /**
     * Kronecker product of the current matrix and the parameter matrix.
     *
     * @param b matrix to post Kronecker-multiply by
     * @return this ⨂ b
     */
    public RealMatrix kroneckerProduct(final RealMatrix b) {
        final int m = getRowDimension();
        final int n = getColumnDimension();

        final int p = b.getRowDimension();
        final int q = b.getColumnDimension();

        final RealMatrix kroneckerProduct = MatrixUtils.createRealMatrix(m * p, n * q);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                kroneckerProduct.setSubMatrix(b.scalarMultiply(getEntry(i, j)) .getData(), i * p, j * q);
            }
        }

        return kroneckerProduct;
    }

    /**
     * Transforms a matrix in a vector (Vectorization).
     * @return a one column matrix
     */
    public RealMatrix stack() {
        final int m = getRowDimension();
        final int n = getColumnDimension();

        final RealMatrix stacked = MatrixUtils.createRealMatrix(m * n, 1);

        for (int i = 0; i < m; i++) {
            stacked.setSubMatrix(getColumnMatrix(i).getData(), i * n, 0);
        }

        return stacked;
    }

    /**
     * Transforms a one-column stacked matrix into a squared matrix (devectorization).
     * @return square matrix
     */
    public RealMatrix unstackSquare() {
        final int m = getRowDimension();
        final int n = getColumnDimension();
        final int s = (int) FastMath.round(FastMath.sqrt(m));

        if (n != 1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH, n, 1);
        }
        if (s * s != m) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX, s, ((double) m) / s);
        }

        final RealMatrix unstacked = MatrixUtils.createRealMatrix(s, s);

        for (int i = 0; i < s; i++) {
            unstacked.setColumnMatrix(i, getSubMatrix(i * s, i * s + s - 1, 0, 0));
        }

        return unstacked;
    }

}
