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
import org.hipparchus.Field;
import org.hipparchus.FieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/**
 * Implementation of {@link FieldMatrix} using a {@link FieldElement}[][] array to store entries.
 * <p>
 * As specified in the {@link FieldMatrix} interface, matrix element indexing
 * is 0-based -- e.g., {@code getEntry(0, 0)}
 * returns the element in the first row, first column of the matrix
 * </p>
 *
 * @param <T> the type of the field elements
 */
public class Array2DRowFieldMatrix<T extends FieldElement<T>>
    extends AbstractFieldMatrix<T>
    implements Serializable {
    /** Serializable version identifier */
    private static final long serialVersionUID = 7260756672015356458L;
    /** Entries of the matrix */
    private T[][] data;

    /**
     * Creates a matrix with no data
     * @param field field to which the elements belong
     */
    public Array2DRowFieldMatrix(final Field<T> field) {
        super(field);
    }

    /**
     * Create a new {@code FieldMatrix<T>} with the supplied row and column dimensions.
     *
     * @param field Field to which the elements belong.
     * @param rowDimension Number of rows in the new matrix.
     * @param columnDimension Number of columns in the new matrix.
     * @throws MathIllegalArgumentException if row or column dimension is not positive.
     */
    public Array2DRowFieldMatrix(final Field<T> field, final int rowDimension,
                                 final int columnDimension)
        throws MathIllegalArgumentException {
        super(field, rowDimension, columnDimension);
        data = MathArrays.buildArray(field, rowDimension, columnDimension);
    }

    /**
     * Create a new {@code FieldMatrix<T>} using the input array as the underlying
     * data array.
     * <p>The input array is copied, not referenced. This constructor has
     * the same effect as calling {@link #Array2DRowFieldMatrix(FieldElement[][], boolean)}
     * with the second argument set to {@code true}.</p>
     *
     * @param d Data for the new matrix.
     * @throws MathIllegalArgumentException if {@code d} is not rectangular.
     * @throws NullArgumentException if {@code d} is {@code null}.
     * @throws MathIllegalArgumentException if there are not at least one row and one column.
     * @see #Array2DRowFieldMatrix(FieldElement[][], boolean)
     */
    public Array2DRowFieldMatrix(final T[][] d)
        throws MathIllegalArgumentException, NullArgumentException {
        this(extractField(d), d);
    }

    /**
     * Create a new {@code FieldMatrix<T>} using the input array as the underlying
     * data array.
     * <p>The input array is copied, not referenced. This constructor has
     * the same effect as calling {@link #Array2DRowFieldMatrix(FieldElement[][], boolean)}
     * with the second argument set to {@code true}.</p>
     *
     * @param field Field to which the elements belong.
     * @param d Data for the new matrix.
     * @throws MathIllegalArgumentException if {@code d} is not rectangular.
     * @throws NullArgumentException if {@code d} is {@code null}.
     * @throws MathIllegalArgumentException if there are not at least one row and one column.
     * @see #Array2DRowFieldMatrix(FieldElement[][], boolean)
     */
    public Array2DRowFieldMatrix(final Field<T> field, final T[][] d)
        throws MathIllegalArgumentException, NullArgumentException {
        super(field);
        copyIn(d);
    }

    /**
     * Create a new {@code FieldMatrix<T>} using the input array as the underlying
     * data array.
     * <p>If an array is built specially in order to be embedded in a
     * {@code FieldMatrix<T>} and not used directly, the {@code copyArray} may be
     * set to {@code false}. This will prevent the copying and improve
     * performance as no new array will be built and no data will be copied.</p>
     *
     * @param d Data for the new matrix.
     * @param copyArray Whether to copy or reference the input array.
     * @throws MathIllegalArgumentException if {@code d} is not rectangular.
     * @throws MathIllegalArgumentException if there are not at least one row and one column.
     * @throws NullArgumentException if {@code d} is {@code null}.
     * @see #Array2DRowFieldMatrix(FieldElement[][])
     */
    public Array2DRowFieldMatrix(final T[][] d, final boolean copyArray)
        throws MathIllegalArgumentException,
        NullArgumentException {
        this(extractField(d), d, copyArray);
    }

    /**
     * Create a new {@code FieldMatrix<T>} using the input array as the underlying
     * data array.
     * <p>If an array is built specially in order to be embedded in a
     * {@code FieldMatrix<T>} and not used directly, the {@code copyArray} may be
     * set to {@code false}. This will prevent the copying and improve
     * performance as no new array will be built and no data will be copied.</p>
     *
     * @param field Field to which the elements belong.
     * @param d Data for the new matrix.
     * @param copyArray Whether to copy or reference the input array.
     * @throws MathIllegalArgumentException if {@code d} is not rectangular.
     * @throws MathIllegalArgumentException if there are not at least one row and one column.
     * @throws NullArgumentException if {@code d} is {@code null}.
     * @see #Array2DRowFieldMatrix(FieldElement[][])
     */
    public Array2DRowFieldMatrix(final Field<T> field, final T[][] d, final boolean copyArray)
        throws MathIllegalArgumentException, NullArgumentException {
        super(field);
        if (copyArray) {
            copyIn(d);
        } else {
            MathUtils.checkNotNull(d);
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
                                                           nCols, d[r].length);
                }
            }
            data = d; // NOPMD - array copy is taken care of by parameter
        }
    }

    /**
     * Create a new (column) {@code FieldMatrix<T>} using {@code v} as the
     * data for the unique column of the created matrix.
     * The input array is copied.
     *
     * @param v Column vector holding data for new matrix.
     * @throws MathIllegalArgumentException if v is empty
     */
    public Array2DRowFieldMatrix(final T[] v) throws MathIllegalArgumentException {
        this(extractField(v), v);
    }

    /**
     * Create a new (column) {@code FieldMatrix<T>} using {@code v} as the
     * data for the unique column of the created matrix.
     * The input array is copied.
     *
     * @param field Field to which the elements belong.
     * @param v Column vector holding data for new matrix.
     */
    public Array2DRowFieldMatrix(final Field<T> field, final T[] v) {
        super(field);
        final int nRows = v.length;
        data = MathArrays.buildArray(getField(), nRows, 1);
        for (int row = 0; row < nRows; row++) {
            data[row][0] = v[row];
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldMatrix<T> createMatrix(final int rowDimension,
                                       final int columnDimension)
        throws MathIllegalArgumentException {
        return new Array2DRowFieldMatrix<>(getField(), rowDimension, columnDimension);
    }

    /** {@inheritDoc} */
    @Override
    public FieldMatrix<T> copy() {
        return new Array2DRowFieldMatrix<>(getField(), copyOut(), false);
    }

    /**
     * Add {@code m} to this matrix.
     *
     * @param m Matrix to be added.
     * @return {@code this} + m.
     * @throws MathIllegalArgumentException if {@code m} is not the same
     * size as this matrix.
     */
    public Array2DRowFieldMatrix<T> add(final Array2DRowFieldMatrix<T> m)
        throws MathIllegalArgumentException {
        // safety check
        checkAdditionCompatible(m);

        final int rowCount    = getRowDimension();
        final int columnCount = getColumnDimension();
        final T[][] outData = MathArrays.buildArray(getField(), rowCount, columnCount);
        for (int row = 0; row < rowCount; row++) {
            final T[] dataRow    = data[row];
            final T[] mRow       = m.data[row];
            final T[] outDataRow = outData[row];
            for (int col = 0; col < columnCount; col++) {
                outDataRow[col] = dataRow[col].add(mRow[col]);
            }
        }

        return new Array2DRowFieldMatrix<>(getField(), outData, false);
    }

    /**
     * Subtract {@code m} from this matrix.
     *
     * @param m Matrix to be subtracted.
     * @return {@code this} + m.
     * @throws MathIllegalArgumentException if {@code m} is not the same
     * size as this matrix.
     */
    public Array2DRowFieldMatrix<T> subtract(final Array2DRowFieldMatrix<T> m)
        throws MathIllegalArgumentException {
        // safety check
        checkSubtractionCompatible(m);

        final int rowCount    = getRowDimension();
        final int columnCount = getColumnDimension();
        final T[][] outData = MathArrays.buildArray(getField(), rowCount, columnCount);
        for (int row = 0; row < rowCount; row++) {
            final T[] dataRow    = data[row];
            final T[] mRow       = m.data[row];
            final T[] outDataRow = outData[row];
            for (int col = 0; col < columnCount; col++) {
                outDataRow[col] = dataRow[col].subtract(mRow[col]);
            }
        }

        return new Array2DRowFieldMatrix<>(getField(), outData, false);

    }

    /**
     * Postmultiplying this matrix by {@code m}.
     *
     * @param m Matrix to postmultiply by.
     * @return {@code this} * m.
     * @throws MathIllegalArgumentException if the number of columns of this
     * matrix is not equal to the number of rows of {@code m}.
     */
    public Array2DRowFieldMatrix<T> multiply(final Array2DRowFieldMatrix<T> m)
        throws MathIllegalArgumentException {
        // safety check
        checkMultiplicationCompatible(m);

        final int nRows = this.getRowDimension();
        final int nCols = m.getColumnDimension();
        final int nSum = this.getColumnDimension();
        final T[][] outData = MathArrays.buildArray(getField(), nRows, nCols);
        for (int row = 0; row < nRows; row++) {
            final T[] dataRow    = data[row];
            final T[] outDataRow = outData[row];
            for (int col = 0; col < nCols; col++) {
                T sum = getField().getZero();
                for (int i = 0; i < nSum; i++) {
                    sum = sum.add(dataRow[i].multiply(m.data[i][col]));
                }
                outDataRow[col] = sum;
            }
        }

        return new Array2DRowFieldMatrix<>(getField(), outData, false);

    }

    /**
     * Returns the result of postmultiplying {@code this} by {@code m^T}.
     * @param m matrix to first transpose and second postmultiply by
     * @return {@code this * m^T}
     * @throws MathIllegalArgumentException if
     * {@code columnDimension(this) != columnDimension(m)}
     * @since 1.3
     */
    public FieldMatrix<T> multiplyTransposed(final Array2DRowFieldMatrix<T> m)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSameColumnDimension(this, m);

        final int nRows = this.getRowDimension();
        final int nCols = m.getRowDimension();
        final int nSum  = this.getColumnDimension();

        final FieldMatrix<T> out   = MatrixUtils.createFieldMatrix(getField(), nRows, nCols);
        final T[][]          mData = m.data;

        // Multiply.
        for (int col = 0; col < nCols; col++) {
            for (int row = 0; row < nRows; row++) {
                final T[] dataRow = data[row];
                final T[] mRow    = mData[col];
                T sum = getField().getZero();
                for (int i = 0; i < nSum; i++) {
                    sum = sum.add(dataRow[i].multiply(mRow[i]));
                }
                out.setEntry(row, col, sum);
            }
        }

        return out;

    }

    /** {@inheritDoc} */
    @Override
    public FieldMatrix<T> multiplyTransposed(final FieldMatrix<T> m) {
        if (m instanceof Array2DRowFieldMatrix) {
            return multiplyTransposed((Array2DRowFieldMatrix<T>) m);
        } else {
            MatrixUtils.checkSameColumnDimension(this, m);

            final int nRows = this.getRowDimension();
            final int nCols = m.getRowDimension();
            final int nSum  = this.getColumnDimension();

            final FieldMatrix<T> out = MatrixUtils.createFieldMatrix(getField(), nRows, nCols);

            // Multiply.
            for (int col = 0; col < nCols; col++) {
                for (int row = 0; row < nRows; row++) {
                    final T[] dataRow = data[row];
                    T sum = getField().getZero();
                    for (int i = 0; i < nSum; i++) {
                        sum = sum.add(dataRow[i].multiply(m.getEntry(col, i)));
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
    public FieldMatrix<T> transposeMultiply(final Array2DRowFieldMatrix<T> m)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSameRowDimension(this, m);

        final int nRows = this.getColumnDimension();
        final int nCols = m.getColumnDimension();
        final int nSum  = this.getRowDimension();

        final FieldMatrix<T> out   = MatrixUtils.createFieldMatrix(getField(), nRows, nCols);
        final T[][]          mData = m.data;

        // Multiply.
        for (int k = 0; k < nSum; k++) {
            final T[] dataK = data[k];
            final T[] mK    = mData[k];
            for (int row = 0; row < nRows; row++) {
                final T dataIRow = dataK[row];
                for (int col = 0; col < nCols; col++) {
                    out.addToEntry(row, col, dataIRow.multiply(mK[col]));
                }
            }
        }

        return out;
    }

    /** {@inheritDoc} */
    @Override
    public FieldMatrix<T> transposeMultiply(final FieldMatrix<T> m) {
        if (m instanceof Array2DRowFieldMatrix) {
            return transposeMultiply((Array2DRowFieldMatrix<T>) m);
        } else {
            MatrixUtils.checkSameRowDimension(this, m);

            final int nRows = this.getColumnDimension();
            final int nCols = m.getColumnDimension();
            final int nSum  = this.getRowDimension();

            final FieldMatrix<T> out = MatrixUtils.createFieldMatrix(getField(), nRows, nCols);

            // Multiply.
            for (int k = 0; k < nSum; k++) {
                final T[] dataK = data[k];
                for (int row = 0; row < nRows; row++) {
                    final T dataIRow = dataK[row];
                    for (int col = 0; col < nCols; col++) {
                        out.addToEntry(row, col, dataIRow.multiply(m.getEntry(k, col)));
                    }
                }
            }

            return out;

        }
    }

    /** {@inheritDoc} */
    @Override
    public T[][] getData() {
        return copyOut();
    }

    /**
     * Get a reference to the underlying data array.
     * This methods returns internal data, <strong>not</strong> fresh copy of it.
     *
     * @return the 2-dimensional array of entries.
     */
    public T[][] getDataRef() {
        return data; // NOPMD - returning an internal array is intentional and documented here
    }

    /** {@inheritDoc} */
    @Override
    public void setSubMatrix(final T[][] subMatrix, final int row,
                             final int column)
        throws MathIllegalArgumentException, NullArgumentException {
        if (data == null) {
            if (row > 0) {
                throw new MathIllegalStateException(LocalizedCoreFormats.FIRST_ROWS_NOT_INITIALIZED_YET, row);
            }
            if (column > 0) {
                throw new MathIllegalStateException(LocalizedCoreFormats.FIRST_COLUMNS_NOT_INITIALIZED_YET, column);
            }
            final int nRows = subMatrix.length;
            if (nRows == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.AT_LEAST_ONE_ROW);
            }

            final int nCols = subMatrix[0].length;
            if (nCols == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.AT_LEAST_ONE_COLUMN);
            }
            data = MathArrays.buildArray(getField(), subMatrix.length, nCols);
            for (int i = 0; i < data.length; ++i) {
                if (subMatrix[i].length != nCols) {
                    throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                           nCols, subMatrix[i].length);
                }
                System.arraycopy(subMatrix[i], 0, data[i + row], column, nCols);
            }
        } else {
            super.setSubMatrix(subMatrix, row, column);
        }

    }

    /** {@inheritDoc} */
    @Override
    public T getEntry(final int row, final int column)
        throws MathIllegalArgumentException {
      try {
        return data[row][column];
      } catch (IndexOutOfBoundsException e) {
        // throw the exact cause of the exception
        checkRowIndex(row);
        checkColumnIndex(column);
        // should never happen
        throw e;
      }
    }

    /** {@inheritDoc} */
    @Override
    public void setEntry(final int row, final int column, final T value)
        throws MathIllegalArgumentException {
        checkRowIndex(row);
        checkColumnIndex(column);

        data[row][column] = value;
    }

    /** {@inheritDoc} */
    @Override
    public void addToEntry(final int row, final int column, final T increment)
        throws MathIllegalArgumentException {
        checkRowIndex(row);
        checkColumnIndex(column);

        data[row][column] = data[row][column].add(increment);
    }

    /** {@inheritDoc} */
    @Override
    public void multiplyEntry(final int row, final int column, final T factor)
        throws MathIllegalArgumentException {
        checkRowIndex(row);
        checkColumnIndex(column);

        data[row][column] = data[row][column].multiply(factor);
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
    public T[] operate(final T[] v) throws MathIllegalArgumentException {
        final int nRows = this.getRowDimension();
        final int nCols = this.getColumnDimension();
        if (v.length != nCols) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   v.length, nCols);
        }
        final T[] out = MathArrays.buildArray(getField(), nRows);
        for (int row = 0; row < nRows; row++) {
            final T[] dataRow = data[row];
            T sum = getField().getZero();
            for (int i = 0; i < nCols; i++) {
                sum = sum.add(dataRow[i].multiply(v[i]));
            }
            out[row] = sum;
        }
        return out;
    }

    /** {@inheritDoc} */
    @Override
    public T[] preMultiply(final T[] v) throws MathIllegalArgumentException {
        final int nRows = getRowDimension();
        final int nCols = getColumnDimension();
        if (v.length != nRows) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   v.length, nRows);
        }

        final T[] out = MathArrays.buildArray(getField(), nCols);
        for (int col = 0; col < nCols; ++col) {
            T sum = getField().getZero();
            for (int i = 0; i < nRows; ++i) {
                sum = sum.add(data[i][col].multiply(v[i]));
            }
            out[col] = sum;
        }

        return out;
    }

    /** {@inheritDoc} */
    @Override
    public FieldMatrix<T> getSubMatrix(final int startRow, final int endRow,
                                       final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn);
        final int rowCount = endRow - startRow + 1;
        final int columnCount = endColumn - startColumn + 1;
        final T[][] outData = MathArrays.buildArray(getField(), rowCount, columnCount);
        for (int i = 0; i < rowCount; ++i) {
            System.arraycopy(data[startRow + i], startColumn, outData[i], 0, columnCount);
        }

        Array2DRowFieldMatrix<T> subMatrix = new Array2DRowFieldMatrix<>(getField());
        subMatrix.data = outData;
        return subMatrix;
    }

    /** {@inheritDoc} */
    @Override
    public T walkInRowOrder(final FieldMatrixChangingVisitor<T> visitor) {
        final int rows    = getRowDimension();
        final int columns = getColumnDimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int i = 0; i < rows; ++i) {
            final T[] rowI = data[i];
            for (int j = 0; j < columns; ++j) {
                rowI[j] = visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public T walkInRowOrder(final FieldMatrixPreservingVisitor<T> visitor) {
        final int rows    = getRowDimension();
        final int columns = getColumnDimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int i = 0; i < rows; ++i) {
            final T[] rowI = data[i];
            for (int j = 0; j < columns; ++j) {
                visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public T walkInRowOrder(final FieldMatrixChangingVisitor<T> visitor,
                            final int startRow, final int endRow,
                            final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
        visitor.start(getRowDimension(), getColumnDimension(),
                      startRow, endRow, startColumn, endColumn);
        for (int i = startRow; i <= endRow; ++i) {
            final T[] rowI = data[i];
            for (int j = startColumn; j <= endColumn; ++j) {
                rowI[j] = visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public T walkInRowOrder(final FieldMatrixPreservingVisitor<T> visitor,
                            final int startRow, final int endRow,
                            final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
        visitor.start(getRowDimension(), getColumnDimension(),
                      startRow, endRow, startColumn, endColumn);
        for (int i = startRow; i <= endRow; ++i) {
            final T[] rowI = data[i];
            for (int j = startColumn; j <= endColumn; ++j) {
                visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public T walkInColumnOrder(final FieldMatrixChangingVisitor<T> visitor) {
        final int rows    = getRowDimension();
        final int columns = getColumnDimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int j = 0; j < columns; ++j) {
            for (int i = 0; i < rows; ++i) {
                final T[] rowI = data[i];
                rowI[j] = visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public T walkInColumnOrder(final FieldMatrixPreservingVisitor<T> visitor) {
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
    public T walkInColumnOrder(final FieldMatrixChangingVisitor<T> visitor,
                               final int startRow, final int endRow,
                               final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
    checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
        visitor.start(getRowDimension(), getColumnDimension(),
                      startRow, endRow, startColumn, endColumn);
        for (int j = startColumn; j <= endColumn; ++j) {
            for (int i = startRow; i <= endRow; ++i) {
                final T[] rowI = data[i];
                rowI[j] = visitor.visit(i, j, rowI[j]);
            }
        }
        return visitor.end();
    }

    /** {@inheritDoc} */
    @Override
    public T walkInColumnOrder(final FieldMatrixPreservingVisitor<T> visitor,
                               final int startRow, final int endRow,
                               final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        checkSubMatrixIndex(startRow, endRow, startColumn, endColumn);
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
    private T[][] copyOut() {
        final int nRows = this.getRowDimension();
        final T[][] out = MathArrays.buildArray(getField(), nRows, getColumnDimension());
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
    private void copyIn(final T[][] in)
        throws MathIllegalArgumentException, NullArgumentException {
        setSubMatrix(in, 0, 0);
    }

    /** {@inheritDoc} */
    @Override
    public T[] getRow(final int row) throws MathIllegalArgumentException {
        MatrixUtils.checkRowIndex(this, row);
        final int nCols = getColumnDimension();
        final T[] out = MathArrays.buildArray(getField(), nCols);
        System.arraycopy(data[row], 0, out, 0, nCols);
        return out;
    }

    /** {@inheritDoc} */
    @Override
    public void setRow(final int row, final T[] array)
        throws MathIllegalArgumentException {
        MatrixUtils.checkRowIndex(this, row);
        final int nCols = getColumnDimension();
        if (array.length != nCols) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH_2x2,
                                                   1, array.length, 1, nCols);
        }
        System.arraycopy(array, 0, data[row], 0, nCols);
    }

}
