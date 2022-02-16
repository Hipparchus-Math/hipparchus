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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.FieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.fraction.BigFraction;
import org.hipparchus.fraction.Fraction;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;

/**
 * A collection of static methods that operate on or return matrices.
 *
 */
public class MatrixUtils {

    /**
     * The default format for {@link RealMatrix} objects.
     */
    public static final RealMatrixFormat DEFAULT_FORMAT = RealMatrixFormat.getRealMatrixFormat();

    /**
     * A format for {@link RealMatrix} objects compatible with octave.
     */
    public static final RealMatrixFormat OCTAVE_FORMAT = new RealMatrixFormat("[", "]", "", "", "; ", ", ");

    /** Pade coefficients required for the matrix exponential calculation. */
    private static final double[] PADE_COEFFICIENTS_3 = {
            120.0,
            60.0,
            12.0,
            1.0
    };

    /** Pade coefficients required for the matrix exponential calculation. */
    private static final double[] PADE_COEFFICIENTS_5 = {
            30240.0,
            15120.0,
            3360.0,
            420.0,
            30.0,
            1
    };

    /** Pade coefficients required for the matrix exponential calculation. */
    private static final double[] PADE_COEFFICIENTS_7 = {
            17297280.0,
            8648640.0,
            1995840.0,
            277200.0,
            25200.0,
            1512.0,
            56.0,
            1.0
    };

    /** Pade coefficients required for the matrix exponential calculation. */
    private static final double[] PADE_COEFFICIENTS_9 = {
            17643225600.0,
            8821612800.0,
            2075673600.0,
            302702400.0,
            30270240.0,
            2162160.0,
            110880.0,
            3960.0,
            90.0,
            1.0
    };

    /** Pade coefficients required for the matrix exponential calculation. */
    private static final double[] PADE_COEFFICIENTS_13 = {
            6.476475253248e+16,
            3.238237626624e+16,
            7.7717703038976e+15,
            1.1873537964288e+15,
            129060195264000.0,
            10559470521600.0,
            670442572800.0,
            33522128640.0,
            1323241920.0,
            40840800.0,
            960960.0,
            16380.0,
            182.0,
            1.0
    };

    /**
     * Private constructor.
     */
    private MatrixUtils() {
        super();
    }

    /**
     * Returns a {@link RealMatrix} with specified dimensions.
     * <p>The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 64&times;64 for a
     * square matrix) which can be stored in a 32kB array, a {@link
     * Array2DRowRealMatrix} instance is built. Above this threshold a {@link
     * BlockRealMatrix} instance is built.</p>
     * <p>The matrix elements are all set to 0.0.</p>
     * @param rows number of rows of the matrix
     * @param columns number of columns of the matrix
     * @return  RealMatrix with specified dimensions
     * @see #createRealMatrix(double[][])
     */
    public static RealMatrix createRealMatrix(final int rows, final int columns) {
        return (rows * columns <= 4096) ?
                new Array2DRowRealMatrix(rows, columns) : new BlockRealMatrix(rows, columns);
    }

    /**
     * Returns a {@link FieldMatrix} with specified dimensions.
     * <p>The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 64&times;64 for a
     * square matrix), a {@link FieldMatrix} instance is built. Above
     * this threshold a {@link BlockFieldMatrix} instance is built.</p>
     * <p>The matrix elements are all set to field.getZero().</p>
     * @param <T> the type of the field elements
     * @param field field to which the matrix elements belong
     * @param rows number of rows of the matrix
     * @param columns number of columns of the matrix
     * @return  FieldMatrix with specified dimensions
     * @see #createFieldMatrix(FieldElement[][])
     */
    public static <T extends FieldElement<T>> FieldMatrix<T> createFieldMatrix(final Field<T> field,
                                                                               final int rows,
                                                                               final int columns) {
        return (rows * columns <= 4096) ?
                new Array2DRowFieldMatrix<T>(field, rows, columns) : new BlockFieldMatrix<T>(field, rows, columns);
    }

    /**
     * Returns a {@link RealMatrix} whose entries are the the values in the
     * the input array.
     * <p>The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 64&times;64 for a
     * square matrix) which can be stored in a 32kB array, a {@link
     * Array2DRowRealMatrix} instance is built. Above this threshold a {@link
     * BlockRealMatrix} instance is built.</p>
     * <p>The input array is copied, not referenced.</p>
     *
     * @param data input array
     * @return  RealMatrix containing the values of the array
     * @throws org.hipparchus.exception.MathIllegalArgumentException
     * if {@code data} is not rectangular (not all rows have the same length).
     * @throws MathIllegalArgumentException if a row or column is empty.
     * @throws NullArgumentException if either {@code data} or {@code data[0]}
     * is {@code null}.
     * @throws MathIllegalArgumentException if {@code data} is not rectangular.
     * @see #createRealMatrix(int, int)
     */
    public static RealMatrix createRealMatrix(double[][] data)
        throws MathIllegalArgumentException, NullArgumentException {
        if (data == null ||
            data[0] == null) {
            throw new NullArgumentException();
        }
        return (data.length * data[0].length <= 4096) ?
                new Array2DRowRealMatrix(data) : new BlockRealMatrix(data);
    }

    /**
     * Returns a {@link FieldMatrix} whose entries are the the values in the
     * the input array.
     * <p>The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 64&times;64 for a
     * square matrix), a {@link FieldMatrix} instance is built. Above
     * this threshold a {@link BlockFieldMatrix} instance is built.</p>
     * <p>The input array is copied, not referenced.</p>
     * @param <T> the type of the field elements
     * @param data input array
     * @return a matrix containing the values of the array.
     * @throws org.hipparchus.exception.MathIllegalArgumentException
     * if {@code data} is not rectangular (not all rows have the same length).
     * @throws MathIllegalArgumentException if a row or column is empty.
     * @throws NullArgumentException if either {@code data} or {@code data[0]}
     * is {@code null}.
     * @see #createFieldMatrix(Field, int, int)
     */
    public static <T extends FieldElement<T>> FieldMatrix<T> createFieldMatrix(T[][] data)
        throws MathIllegalArgumentException, NullArgumentException {
        if (data == null ||
            data[0] == null) {
            throw new NullArgumentException();
        }
        return (data.length * data[0].length <= 4096) ?
                new Array2DRowFieldMatrix<T>(data) : new BlockFieldMatrix<T>(data);
    }

    /**
     * Returns <code>dimension x dimension</code> identity matrix.
     *
     * @param dimension dimension of identity matrix to generate
     * @return identity matrix
     * @throws IllegalArgumentException if dimension is not positive
     */
    public static RealMatrix createRealIdentityMatrix(int dimension) {
        final RealMatrix m = createRealMatrix(dimension, dimension);
        for (int i = 0; i < dimension; ++i) {
            m.setEntry(i, i, 1.0);
        }
        return m;
    }

    /**
     * Returns <code>dimension x dimension</code> identity matrix.
     *
     * @param <T> the type of the field elements
     * @param field field to which the elements belong
     * @param dimension dimension of identity matrix to generate
     * @return identity matrix
     * @throws IllegalArgumentException if dimension is not positive
     */
    public static <T extends FieldElement<T>> FieldMatrix<T>
        createFieldIdentityMatrix(final Field<T> field, final int dimension) {
        final T zero = field.getZero();
        final T one  = field.getOne();
        final T[][] d = MathArrays.buildArray(field, dimension, dimension);
        for (int row = 0; row < dimension; row++) {
            final T[] dRow = d[row];
            Arrays.fill(dRow, zero);
            dRow[row] = one;
        }
        return new Array2DRowFieldMatrix<T>(field, d, false);
    }

    /**
     * Returns a diagonal matrix with specified elements.
     *
     * @param diagonal diagonal elements of the matrix (the array elements
     * will be copied)
     * @return diagonal matrix
     */
    public static RealMatrix createRealDiagonalMatrix(final double[] diagonal) {
        final RealMatrix m = createRealMatrix(diagonal.length, diagonal.length);
        for (int i = 0; i < diagonal.length; ++i) {
            m.setEntry(i, i, diagonal[i]);
        }
        return m;
    }

    /**
     * Returns a diagonal matrix with specified elements.
     *
     * @param <T> the type of the field elements
     * @param diagonal diagonal elements of the matrix (the array elements
     * will be copied)
     * @return diagonal matrix
     */
    public static <T extends FieldElement<T>> FieldMatrix<T>
        createFieldDiagonalMatrix(final T[] diagonal) {
        final FieldMatrix<T> m =
            createFieldMatrix(diagonal[0].getField(), diagonal.length, diagonal.length);
        for (int i = 0; i < diagonal.length; ++i) {
            m.setEntry(i, i, diagonal[i]);
        }
        return m;
    }

    /**
     * Creates a {@link RealVector} using the data from the input array.
     *
     * @param data the input data
     * @return a data.length RealVector
     * @throws MathIllegalArgumentException if {@code data} is empty.
     * @throws NullArgumentException if {@code data} is {@code null}.
     */
    public static RealVector createRealVector(double[] data)
        throws MathIllegalArgumentException, NullArgumentException {
        if (data == null) {
            throw new NullArgumentException();
        }
        return new ArrayRealVector(data, true);
    }

    /**
     * Creates a {@link RealVector} with specified dimensions.
     *
     * @param dimension dimension of the vector
     * @return a new vector
     * @since 1.3
     */
    public static RealVector createRealVector(final int dimension) {
        return new ArrayRealVector(new double[dimension]);
    }

    /**
     * Creates a {@link FieldVector} using the data from the input array.
     *
     * @param <T> the type of the field elements
     * @param data the input data
     * @return a data.length FieldVector
     * @throws MathIllegalArgumentException if {@code data} is empty.
     * @throws NullArgumentException if {@code data} is {@code null}.
     * @throws MathIllegalArgumentException if {@code data} has 0 elements
     */
    public static <T extends FieldElement<T>> FieldVector<T> createFieldVector(final T[] data)
        throws MathIllegalArgumentException, NullArgumentException {
        if (data == null) {
            throw new NullArgumentException();
        }
        if (data.length == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT);
        }
        return new ArrayFieldVector<T>(data[0].getField(), data, true);
    }

    /**
     * Creates a {@link FieldVector} with specified dimensions.
     *
     * @param <T> the type of the field elements
     * @param field field to which array elements belong
     * @param dimension dimension of the vector
     * @return a new vector
     * @since 1.3
     */
    public static <T extends FieldElement<T>> FieldVector<T> createFieldVector(final Field<T> field, final int dimension) {
        return new ArrayFieldVector<>(MathArrays.buildArray(field, dimension));
    }

    /**
     * Create a row {@link RealMatrix} using the data from the input
     * array.
     *
     * @param rowData the input row data
     * @return a 1 x rowData.length RealMatrix
     * @throws MathIllegalArgumentException if {@code rowData} is empty.
     * @throws NullArgumentException if {@code rowData} is {@code null}.
     */
    public static RealMatrix createRowRealMatrix(double[] rowData)
        throws MathIllegalArgumentException, NullArgumentException {
        if (rowData == null) {
            throw new NullArgumentException();
        }
        final int nCols = rowData.length;
        final RealMatrix m = createRealMatrix(1, nCols);
        for (int i = 0; i < nCols; ++i) {
            m.setEntry(0, i, rowData[i]);
        }
        return m;
    }

    /**
     * Create a row {@link FieldMatrix} using the data from the input
     * array.
     *
     * @param <T> the type of the field elements
     * @param rowData the input row data
     * @return a 1 x rowData.length FieldMatrix
     * @throws MathIllegalArgumentException if {@code rowData} is empty.
     * @throws NullArgumentException if {@code rowData} is {@code null}.
     */
    public static <T extends FieldElement<T>> FieldMatrix<T>
        createRowFieldMatrix(final T[] rowData)
        throws MathIllegalArgumentException, NullArgumentException {
        if (rowData == null) {
            throw new NullArgumentException();
        }
        final int nCols = rowData.length;
        if (nCols == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.AT_LEAST_ONE_COLUMN);
        }
        final FieldMatrix<T> m = createFieldMatrix(rowData[0].getField(), 1, nCols);
        for (int i = 0; i < nCols; ++i) {
            m.setEntry(0, i, rowData[i]);
        }
        return m;
    }

    /**
     * Creates a column {@link RealMatrix} using the data from the input
     * array.
     *
     * @param columnData  the input column data
     * @return a columnData x 1 RealMatrix
     * @throws MathIllegalArgumentException if {@code columnData} is empty.
     * @throws NullArgumentException if {@code columnData} is {@code null}.
     */
    public static RealMatrix createColumnRealMatrix(double[] columnData)
        throws MathIllegalArgumentException, NullArgumentException {
        if (columnData == null) {
            throw new NullArgumentException();
        }
        final int nRows = columnData.length;
        final RealMatrix m = createRealMatrix(nRows, 1);
        for (int i = 0; i < nRows; ++i) {
            m.setEntry(i, 0, columnData[i]);
        }
        return m;
    }

    /**
     * Creates a column {@link FieldMatrix} using the data from the input
     * array.
     *
     * @param <T> the type of the field elements
     * @param columnData  the input column data
     * @return a columnData x 1 FieldMatrix
     * @throws MathIllegalArgumentException if {@code data} is empty.
     * @throws NullArgumentException if {@code columnData} is {@code null}.
     */
    public static <T extends FieldElement<T>> FieldMatrix<T>
        createColumnFieldMatrix(final T[] columnData)
        throws MathIllegalArgumentException, NullArgumentException {
        if (columnData == null) {
            throw new NullArgumentException();
        }
        final int nRows = columnData.length;
        if (nRows == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.AT_LEAST_ONE_ROW);
        }
        final FieldMatrix<T> m = createFieldMatrix(columnData[0].getField(), nRows, 1);
        for (int i = 0; i < nRows; ++i) {
            m.setEntry(i, 0, columnData[i]);
        }
        return m;
    }

    /**
     * Checks whether a matrix is symmetric, within a given relative tolerance.
     *
     * @param matrix Matrix to check.
     * @param relativeTolerance Tolerance of the symmetry check.
     * @param raiseException If {@code true}, an exception will be raised if
     * the matrix is not symmetric.
     * @return {@code true} if {@code matrix} is symmetric.
     * @throws MathIllegalArgumentException if the matrix is not square.
     * @throws MathIllegalArgumentException if the matrix is not symmetric.
     */
    private static boolean isSymmetricInternal(RealMatrix matrix,
                                               double relativeTolerance,
                                               boolean raiseException) {
        final int rows = matrix.getRowDimension();
        if (rows != matrix.getColumnDimension()) {
            if (raiseException) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                       rows, matrix.getColumnDimension());
            } else {
                return false;
            }
        }
        for (int i = 0; i < rows; i++) {
            for (int j = i + 1; j < rows; j++) {
                final double mij = matrix.getEntry(i, j);
                final double mji = matrix.getEntry(j, i);
                if (FastMath.abs(mij - mji) >
                    FastMath.max(FastMath.abs(mij), FastMath.abs(mji)) * relativeTolerance) {
                    if (raiseException) {
                        throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SYMMETRIC_MATRIX,
                                                               i, j, relativeTolerance);
                    } else {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     * Checks whether a matrix is symmetric.
     *
     * @param matrix Matrix to check.
     * @param eps Relative tolerance.
     * @throws MathIllegalArgumentException if the matrix is not square.
     * @throws MathIllegalArgumentException if the matrix is not symmetric.
     */
    public static void checkSymmetric(RealMatrix matrix,
                                      double eps) {
        isSymmetricInternal(matrix, eps, true);
    }

    /**
     * Checks whether a matrix is symmetric.
     *
     * @param matrix Matrix to check.
     * @param eps Relative tolerance.
     * @return {@code true} if {@code matrix} is symmetric.
     */
    public static boolean isSymmetric(RealMatrix matrix,
                                      double eps) {
        return isSymmetricInternal(matrix, eps, false);
    }

    /**
     * Check if matrix indices are valid.
     *
     * @param m Matrix.
     * @param row Row index to check.
     * @param column Column index to check.
     * @throws MathIllegalArgumentException if {@code row} or {@code column} is not
     * a valid index.
     */
    public static void checkMatrixIndex(final AnyMatrix m,
                                        final int row, final int column)
        throws MathIllegalArgumentException {
        checkRowIndex(m, row);
        checkColumnIndex(m, column);
    }

    /**
     * Check if a row index is valid.
     *
     * @param m Matrix.
     * @param row Row index to check.
     * @throws MathIllegalArgumentException if {@code row} is not a valid index.
     */
    public static void checkRowIndex(final AnyMatrix m, final int row)
        throws MathIllegalArgumentException {
        if (row < 0 ||
            row >= m.getRowDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ROW_INDEX,
                                          row, 0, m.getRowDimension() - 1);
        }
    }

    /**
     * Check if a column index is valid.
     *
     * @param m Matrix.
     * @param column Column index to check.
     * @throws MathIllegalArgumentException if {@code column} is not a valid index.
     */
    public static void checkColumnIndex(final AnyMatrix m, final int column)
        throws MathIllegalArgumentException {
        if (column < 0 || column >= m.getColumnDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.COLUMN_INDEX,
                                           column, 0, m.getColumnDimension() - 1);
        }
    }

    /**
     * Check if submatrix ranges indices are valid.
     * Rows and columns are indicated counting from 0 to {@code n - 1}.
     *
     * @param m Matrix.
     * @param startRow Initial row index.
     * @param endRow Final row index.
     * @param startColumn Initial column index.
     * @param endColumn Final column index.
     * @throws MathIllegalArgumentException if the indices are invalid.
     * @throws MathIllegalArgumentException if {@code endRow < startRow} or
     * {@code endColumn < startColumn}.
     */
    public static void checkSubMatrixIndex(final AnyMatrix m,
                                           final int startRow, final int endRow,
                                           final int startColumn, final int endColumn)
        throws MathIllegalArgumentException {
        checkRowIndex(m, startRow);
        checkRowIndex(m, endRow);
        if (endRow < startRow) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INITIAL_ROW_AFTER_FINAL_ROW,
                                                endRow, startRow, false);
        }

        checkColumnIndex(m, startColumn);
        checkColumnIndex(m, endColumn);
        if (endColumn < startColumn) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INITIAL_COLUMN_AFTER_FINAL_COLUMN,
                                                endColumn, startColumn, false);
        }


    }

    /**
     * Check if submatrix ranges indices are valid.
     * Rows and columns are indicated counting from 0 to n-1.
     *
     * @param m Matrix.
     * @param selectedRows Array of row indices.
     * @param selectedColumns Array of column indices.
     * @throws NullArgumentException if {@code selectedRows} or
     * {@code selectedColumns} are {@code null}.
     * @throws MathIllegalArgumentException if the row or column selections are empty (zero
     * length).
     * @throws MathIllegalArgumentException if row or column selections are not valid.
     */
    public static void checkSubMatrixIndex(final AnyMatrix m,
                                           final int[] selectedRows,
                                           final int[] selectedColumns)
        throws MathIllegalArgumentException, NullArgumentException {
        if (selectedRows == null) {
            throw new NullArgumentException();
        }
        if (selectedColumns == null) {
            throw new NullArgumentException();
        }
        if (selectedRows.length == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.EMPTY_SELECTED_ROW_INDEX_ARRAY);
        }
        if (selectedColumns.length == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.EMPTY_SELECTED_COLUMN_INDEX_ARRAY);
        }

        for (final int row : selectedRows) {
            checkRowIndex(m, row);
        }
        for (final int column : selectedColumns) {
            checkColumnIndex(m, column);
        }
    }

    /**
     * Check if matrices are addition compatible.
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @throws MathIllegalArgumentException if the matrices are not addition
     * compatible.
     */
    public static void checkAdditionCompatible(final AnyMatrix left, final AnyMatrix right)
        throws MathIllegalArgumentException {
        if ((left.getRowDimension()    != right.getRowDimension()) ||
            (left.getColumnDimension() != right.getColumnDimension())) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH_2x2,
                                                   left.getRowDimension(), left.getColumnDimension(),
                                                   right.getRowDimension(), right.getColumnDimension());
        }
    }

    /**
     * Check if matrices are subtraction compatible
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @throws MathIllegalArgumentException if the matrices are not addition
     * compatible.
     */
    public static void checkSubtractionCompatible(final AnyMatrix left, final AnyMatrix right)
        throws MathIllegalArgumentException {
        if ((left.getRowDimension()    != right.getRowDimension()) ||
            (left.getColumnDimension() != right.getColumnDimension())) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH_2x2,
                                                   left.getRowDimension(), left.getColumnDimension(),
                                                   right.getRowDimension(), right.getColumnDimension());
        }
    }

    /**
     * Check if matrices are multiplication compatible
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @throws MathIllegalArgumentException if matrices are not multiplication
     * compatible.
     */
    public static void checkMultiplicationCompatible(final AnyMatrix left, final AnyMatrix right)
        throws MathIllegalArgumentException {
        if (left.getColumnDimension() != right.getRowDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   left.getColumnDimension(), right.getRowDimension());
        }
    }

    /**
     * Check if matrices have the same number of columns.
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @throws MathIllegalArgumentException if matrices don't have the same number of columns.
     * @since 1.3
     */
    public static void checkSameColumnDimension(final AnyMatrix left, final AnyMatrix right)
        throws MathIllegalArgumentException {
        if (left.getColumnDimension() != right.getColumnDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   left.getColumnDimension(), right.getColumnDimension());
        }
    }

    /**
     * Check if matrices have the same number of rows.
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @throws MathIllegalArgumentException if matrices don't have the same number of rows.
     * @since 1.3
     */
    public static void checkSameRowDimension(final AnyMatrix left, final AnyMatrix right)
        throws MathIllegalArgumentException {
        if (left.getRowDimension() != right.getRowDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   left.getRowDimension(), right.getRowDimension());
        }
    }

    /**
     * Convert a {@link FieldMatrix}/{@link Fraction} matrix to a {@link RealMatrix}.
     * @param m Matrix to convert.
     * @return the converted matrix.
     */
    public static Array2DRowRealMatrix fractionMatrixToRealMatrix(final FieldMatrix<Fraction> m) {
        final FractionMatrixConverter converter = new FractionMatrixConverter();
        m.walkInOptimizedOrder(converter);
        return converter.getConvertedMatrix();
    }

    /** Converter for {@link FieldMatrix}/{@link Fraction}. */
    private static class FractionMatrixConverter extends DefaultFieldMatrixPreservingVisitor<Fraction> {
        /** Converted array. */
        private double[][] data;
        /** Simple constructor. */
        FractionMatrixConverter() {
            super(Fraction.ZERO);
        }

        /** {@inheritDoc} */
        @Override
        public void start(int rows, int columns,
                          int startRow, int endRow, int startColumn, int endColumn) {
            data = new double[rows][columns];
        }

        /** {@inheritDoc} */
        @Override
        public void visit(int row, int column, Fraction value) {
            data[row][column] = value.doubleValue();
        }

        /**
         * Get the converted matrix.
         *
         * @return the converted matrix.
         */
        Array2DRowRealMatrix getConvertedMatrix() {
            return new Array2DRowRealMatrix(data, false);
        }

    }

    /**
     * Convert a {@link FieldMatrix}/{@link BigFraction} matrix to a {@link RealMatrix}.
     *
     * @param m Matrix to convert.
     * @return the converted matrix.
     */
    public static Array2DRowRealMatrix bigFractionMatrixToRealMatrix(final FieldMatrix<BigFraction> m) {
        final BigFractionMatrixConverter converter = new BigFractionMatrixConverter();
        m.walkInOptimizedOrder(converter);
        return converter.getConvertedMatrix();
    }

    /** Converter for {@link FieldMatrix}/{@link BigFraction}. */
    private static class BigFractionMatrixConverter extends DefaultFieldMatrixPreservingVisitor<BigFraction> {
        /** Converted array. */
        private double[][] data;
        /** Simple constructor. */
        BigFractionMatrixConverter() {
            super(BigFraction.ZERO);
        }

        /** {@inheritDoc} */
        @Override
        public void start(int rows, int columns,
                          int startRow, int endRow, int startColumn, int endColumn) {
            data = new double[rows][columns];
        }

        /** {@inheritDoc} */
        @Override
        public void visit(int row, int column, BigFraction value) {
            data[row][column] = value.doubleValue();
        }

        /**
         * Get the converted matrix.
         *
         * @return the converted matrix.
         */
        Array2DRowRealMatrix getConvertedMatrix() {
            return new Array2DRowRealMatrix(data, false);
        }
    }

    /** Serialize a {@link RealVector}.
     * <p>
     * This method is intended to be called from within a private
     * <code>writeObject</code> method (after a call to
     * <code>oos.defaultWriteObject()</code>) in a class that has a
     * {@link RealVector} field, which should be declared <code>transient</code>.
     * This way, the default handling does not serialize the vector (the {@link
     * RealVector} interface is not serializable by default) but this method does
     * serialize it specifically.
     * </p>
     * <p>
     * The following example shows how a simple class with a name and a real vector
     * should be written:
     * <pre><code>
     * public class NamedVector implements Serializable {
     *
     *     private final String name;
     *     private final transient RealVector coefficients;
     *
     *     // omitted constructors, getters ...
     *
     *     private void writeObject(ObjectOutputStream oos) throws IOException {
     *         oos.defaultWriteObject();  // takes care of name field
     *         MatrixUtils.serializeRealVector(coefficients, oos);
     *     }
     *
     *     private void readObject(ObjectInputStream ois) throws ClassNotFoundException, IOException {
     *         ois.defaultReadObject();  // takes care of name field
     *         MatrixUtils.deserializeRealVector(this, "coefficients", ois);
     *     }
     *
     * }
     * </code></pre>
     * </p>
     *
     * @param vector real vector to serialize
     * @param oos stream where the real vector should be written
     * @exception IOException if object cannot be written to stream
     * @see #deserializeRealVector(Object, String, ObjectInputStream)
     */
    public static void serializeRealVector(final RealVector vector,
                                           final ObjectOutputStream oos)
        throws IOException {
        final int n = vector.getDimension();
        oos.writeInt(n);
        for (int i = 0; i < n; ++i) {
            oos.writeDouble(vector.getEntry(i));
        }
    }

    /** Deserialize  a {@link RealVector} field in a class.
     * <p>
     * This method is intended to be called from within a private
     * <code>readObject</code> method (after a call to
     * <code>ois.defaultReadObject()</code>) in a class that has a
     * {@link RealVector} field, which should be declared <code>transient</code>.
     * This way, the default handling does not deserialize the vector (the {@link
     * RealVector} interface is not serializable by default) but this method does
     * deserialize it specifically.
     * </p>
     * @param instance instance in which the field must be set up
     * @param fieldName name of the field within the class (may be private and final)
     * @param ois stream from which the real vector should be read
     * @exception ClassNotFoundException if a class in the stream cannot be found
     * @exception IOException if object cannot be read from the stream
     * @see #serializeRealVector(RealVector, ObjectOutputStream)
     */
    public static void deserializeRealVector(final Object instance,
                                             final String fieldName,
                                             final ObjectInputStream ois)
      throws ClassNotFoundException, IOException {
        try {

            // read the vector data
            final int n = ois.readInt();
            final double[] data = new double[n];
            for (int i = 0; i < n; ++i) {
                data[i] = ois.readDouble();
            }

            // create the instance
            final RealVector vector = new ArrayRealVector(data, false);

            // set up the field
            final java.lang.reflect.Field f =
                instance.getClass().getDeclaredField(fieldName);
            f.setAccessible(true);
            f.set(instance, vector);

        } catch (NoSuchFieldException | IllegalAccessException e) {
            IOException ioe = new IOException();
            ioe.initCause(e);
            throw ioe;
        }

    }

    /** Serialize a {@link RealMatrix}.
     * <p>
     * This method is intended to be called from within a private
     * <code>writeObject</code> method (after a call to
     * <code>oos.defaultWriteObject()</code>) in a class that has a
     * {@link RealMatrix} field, which should be declared <code>transient</code>.
     * This way, the default handling does not serialize the matrix (the {@link
     * RealMatrix} interface is not serializable by default) but this method does
     * serialize it specifically.
     * </p>
     * <p>
     * The following example shows how a simple class with a name and a real matrix
     * should be written:
     * <pre><code>
     * public class NamedMatrix implements Serializable {
     *
     *     private final String name;
     *     private final transient RealMatrix coefficients;
     *
     *     // omitted constructors, getters ...
     *
     *     private void writeObject(ObjectOutputStream oos) throws IOException {
     *         oos.defaultWriteObject();  // takes care of name field
     *         MatrixUtils.serializeRealMatrix(coefficients, oos);
     *     }
     *
     *     private void readObject(ObjectInputStream ois) throws ClassNotFoundException, IOException {
     *         ois.defaultReadObject();  // takes care of name field
     *         MatrixUtils.deserializeRealMatrix(this, "coefficients", ois);
     *     }
     *
     * }
     * </code></pre>
     * </p>
     *
     * @param matrix real matrix to serialize
     * @param oos stream where the real matrix should be written
     * @exception IOException if object cannot be written to stream
     * @see #deserializeRealMatrix(Object, String, ObjectInputStream)
     */
    public static void serializeRealMatrix(final RealMatrix matrix,
                                           final ObjectOutputStream oos)
        throws IOException {
        final int n = matrix.getRowDimension();
        final int m = matrix.getColumnDimension();
        oos.writeInt(n);
        oos.writeInt(m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                oos.writeDouble(matrix.getEntry(i, j));
            }
        }
    }

    /** Deserialize  a {@link RealMatrix} field in a class.
     * <p>
     * This method is intended to be called from within a private
     * <code>readObject</code> method (after a call to
     * <code>ois.defaultReadObject()</code>) in a class that has a
     * {@link RealMatrix} field, which should be declared <code>transient</code>.
     * This way, the default handling does not deserialize the matrix (the {@link
     * RealMatrix} interface is not serializable by default) but this method does
     * deserialize it specifically.
     * </p>
     * @param instance instance in which the field must be set up
     * @param fieldName name of the field within the class (may be private and final)
     * @param ois stream from which the real matrix should be read
     * @exception ClassNotFoundException if a class in the stream cannot be found
     * @exception IOException if object cannot be read from the stream
     * @see #serializeRealMatrix(RealMatrix, ObjectOutputStream)
     */
    public static void deserializeRealMatrix(final Object instance,
                                             final String fieldName,
                                             final ObjectInputStream ois)
      throws ClassNotFoundException, IOException {
        try {

            // read the matrix data
            final int n = ois.readInt();
            final int m = ois.readInt();
            final double[][] data = new double[n][m];
            for (int i = 0; i < n; ++i) {
                final double[] dataI = data[i];
                for (int j = 0; j < m; ++j) {
                    dataI[j] = ois.readDouble();
                }
            }

            // create the instance
            final RealMatrix matrix = new Array2DRowRealMatrix(data, false);

            // set up the field
            final java.lang.reflect.Field f =
                instance.getClass().getDeclaredField(fieldName);
            f.setAccessible(true);
            f.set(instance, matrix);

        } catch (NoSuchFieldException | IllegalAccessException e) {
            IOException ioe = new IOException();
            ioe.initCause(e);
            throw ioe;
        }
    }

    /**Solve  a  system of composed of a Lower Triangular Matrix
     * {@link RealMatrix}.
     * <p>
     * This method is called to solve systems of equations which are
     * of the lower triangular form. The matrix {@link RealMatrix}
     * is assumed, though not checked, to be in lower triangular form.
     * The vector {@link RealVector} is overwritten with the solution.
     * The matrix is checked that it is square and its dimensions match
     * the length of the vector.
     * </p>
     * @param rm RealMatrix which is lower triangular
     * @param b  RealVector this is overwritten
     * @throws MathIllegalArgumentException if the matrix and vector are not
     * conformable
     * @throws MathIllegalArgumentException if the matrix {@code rm} is not square
     * @throws MathRuntimeException if the absolute value of one of the diagonal
     * coefficient of {@code rm} is lower than {@link Precision#SAFE_MIN}
     */
    public static void solveLowerTriangularSystem(RealMatrix rm, RealVector b)
        throws MathIllegalArgumentException, MathRuntimeException {
        if ((rm == null) || (b == null) || ( rm.getRowDimension() != b.getDimension())) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   (rm == null) ? 0 : rm.getRowDimension(),
                                                   (b  == null) ? 0 : b.getDimension());
        }
        if( rm.getColumnDimension() != rm.getRowDimension() ){
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   rm.getRowDimension(), rm.getColumnDimension());
        }
        int rows = rm.getRowDimension();
        for( int i = 0 ; i < rows ; i++ ){
            double diag = rm.getEntry(i, i);
            if( FastMath.abs(diag) < Precision.SAFE_MIN ){
                throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
            }
            double bi = b.getEntry(i)/diag;
            b.setEntry(i,  bi );
            for( int j = i+1; j< rows; j++ ){
                b.setEntry(j, b.getEntry(j)-bi*rm.getEntry(j,i)  );
            }
        }
    }

    /** Solver a  system composed  of an Upper Triangular Matrix
     * {@link RealMatrix}.
     * <p>
     * This method is called to solve systems of equations which are
     * of the lower triangular form. The matrix {@link RealMatrix}
     * is assumed, though not checked, to be in upper triangular form.
     * The vector {@link RealVector} is overwritten with the solution.
     * The matrix is checked that it is square and its dimensions match
     * the length of the vector.
     * </p>
     * @param rm RealMatrix which is upper triangular
     * @param b  RealVector this is overwritten
     * @throws MathIllegalArgumentException if the matrix and vector are not
     * conformable
     * @throws MathIllegalArgumentException if the matrix {@code rm} is not
     * square
     * @throws MathRuntimeException if the absolute value of one of the diagonal
     * coefficient of {@code rm} is lower than {@link Precision#SAFE_MIN}
     */
    public static void solveUpperTriangularSystem(RealMatrix rm, RealVector b)
        throws MathIllegalArgumentException, MathRuntimeException {
        if ((rm == null) || (b == null) || ( rm.getRowDimension() != b.getDimension())) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   (rm == null) ? 0 : rm.getRowDimension(),
                                                   (b  == null) ? 0 : b.getDimension());
        }
        if( rm.getColumnDimension() != rm.getRowDimension() ){
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   rm.getRowDimension(), rm.getColumnDimension());
        }
        int rows = rm.getRowDimension();
        for( int i = rows-1 ; i >-1 ; i-- ){
            double diag = rm.getEntry(i, i);
            if( FastMath.abs(diag) < Precision.SAFE_MIN ){
                throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
            }
            double bi = b.getEntry(i)/diag;
            b.setEntry(i,  bi );
            for( int j = i-1; j>-1; j-- ){
                b.setEntry(j, b.getEntry(j)-bi*rm.getEntry(j,i)  );
            }
        }
    }

    /**
     * Computes the inverse of the given matrix by splitting it into
     * 4 sub-matrices.
     *
     * @param m Matrix whose inverse must be computed.
     * @param splitIndex Index that determines the "split" line and
     * column.
     * The element corresponding to this index will part of the
     * upper-left sub-matrix.
     * @return the inverse of {@code m}.
     * @throws MathIllegalArgumentException if {@code m} is not square.
     */
    public static RealMatrix blockInverse(RealMatrix m,
                                          int splitIndex) {
        final int n = m.getRowDimension();
        if (m.getColumnDimension() != n) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   m.getRowDimension(), m.getColumnDimension());
        }

        final int splitIndex1 = splitIndex + 1;

        final RealMatrix a = m.getSubMatrix(0, splitIndex, 0, splitIndex);
        final RealMatrix b = m.getSubMatrix(0, splitIndex, splitIndex1, n - 1);
        final RealMatrix c = m.getSubMatrix(splitIndex1, n - 1, 0, splitIndex);
        final RealMatrix d = m.getSubMatrix(splitIndex1, n - 1, splitIndex1, n - 1);

        final SingularValueDecomposition aDec = new SingularValueDecomposition(a);
        final DecompositionSolver aSolver = aDec.getSolver();
        if (!aSolver.isNonSingular()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SINGULAR_MATRIX);
        }
        final RealMatrix aInv = aSolver.getInverse();

        final SingularValueDecomposition dDec = new SingularValueDecomposition(d);
        final DecompositionSolver dSolver = dDec.getSolver();
        if (!dSolver.isNonSingular()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SINGULAR_MATRIX);
        }
        final RealMatrix dInv = dSolver.getInverse();

        final RealMatrix tmp1 = a.subtract(b.multiply(dInv).multiply(c));
        final SingularValueDecomposition tmp1Dec = new SingularValueDecomposition(tmp1);
        final DecompositionSolver tmp1Solver = tmp1Dec.getSolver();
        if (!tmp1Solver.isNonSingular()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SINGULAR_MATRIX);
        }
        final RealMatrix result00 = tmp1Solver.getInverse();

        final RealMatrix tmp2 = d.subtract(c.multiply(aInv).multiply(b));
        final SingularValueDecomposition tmp2Dec = new SingularValueDecomposition(tmp2);
        final DecompositionSolver tmp2Solver = tmp2Dec.getSolver();
        if (!tmp2Solver.isNonSingular()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SINGULAR_MATRIX);
        }
        final RealMatrix result11 = tmp2Solver.getInverse();

        final RealMatrix result01 = aInv.multiply(b).multiply(result11).scalarMultiply(-1);
        final RealMatrix result10 = dInv.multiply(c).multiply(result00).scalarMultiply(-1);

        final RealMatrix result = new Array2DRowRealMatrix(n, n);
        result.setSubMatrix(result00.getData(), 0, 0);
        result.setSubMatrix(result01.getData(), 0, splitIndex1);
        result.setSubMatrix(result10.getData(), splitIndex1, 0);
        result.setSubMatrix(result11.getData(), splitIndex1, splitIndex1);

        return result;
    }

    /**
     * Computes the inverse of the given matrix.
     * <p>
     * By default, the inverse of the matrix is computed using the QR-decomposition,
     * unless a more efficient method can be determined for the input matrix.
     * <p>
     * Note: this method will use a singularity threshold of 0,
     * use {@link #inverse(RealMatrix, double)} if a different threshold is needed.
     *
     * @param matrix Matrix whose inverse shall be computed
     * @return the inverse of {@code matrix}
     * @throws NullArgumentException if {@code matrix} is {@code null}
     * @throws MathIllegalArgumentException if m is singular
     * @throws MathIllegalArgumentException if matrix is not square
     */
    public static RealMatrix inverse(RealMatrix matrix)
            throws MathIllegalArgumentException, NullArgumentException {
        return inverse(matrix, 0);
    }

    /**
     * Computes the inverse of the given matrix.
     * <p>
     * By default, the inverse of the matrix is computed using the QR-decomposition,
     * unless a more efficient method can be determined for the input matrix.
     *
     * @param matrix Matrix whose inverse shall be computed
     * @param threshold Singularity threshold
     * @return the inverse of {@code m}
     * @throws NullArgumentException if {@code matrix} is {@code null}
     * @throws MathIllegalArgumentException if matrix is singular
     * @throws MathIllegalArgumentException if matrix is not square
     */
    public static RealMatrix inverse(RealMatrix matrix, double threshold)
            throws MathIllegalArgumentException, NullArgumentException {

        MathUtils.checkNotNull(matrix);

        if (!matrix.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   matrix.getRowDimension(), matrix.getColumnDimension());
        }

        if (matrix instanceof DiagonalMatrix) {
            return ((DiagonalMatrix) matrix).inverse(threshold);
        } else {
            QRDecomposition decomposition = new QRDecomposition(matrix, threshold);
            return decomposition.getSolver().getInverse();
        }
    }

    /**
     * Computes the <a href="https://mathworld.wolfram.com/MatrixExponential.html">
     * matrix exponential</a> of the given matrix.
     *
     * The algorithm implementation follows the Pade approximant method of
     * <p>Higham, Nicholas J. “The Scaling and Squaring Method for the Matrix Exponential
     * Revisited.” SIAM Journal on Matrix Analysis and Applications 26, no. 4 (January 2005): 1179–93.</p>
     *
     * @param rm RealMatrix whose inverse shall be computed
     * @return The inverse of {@code rm}
     * @throws MathIllegalArgumentException if matrix is not square
     */
    public static RealMatrix matrixExponential(final RealMatrix rm) {

        // Check that the input matrix is square
        if (!rm.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                    rm.getRowDimension(), rm.getColumnDimension());
        }

        // Preprocessing to reduce the norm
        int dim = rm.getRowDimension();
        final RealMatrix identity = MatrixUtils.createRealIdentityMatrix(dim);
        final double preprocessScale = rm.getTrace() / dim;
        RealMatrix scaledMatrix = rm.copy();
        scaledMatrix = scaledMatrix.subtract(identity.scalarMultiply(preprocessScale));

        // Select pade degree required
        final double l1Norm = rm.getNorm1();
        double[] padeCoefficients;
        int squaringCount = 0;

        if (l1Norm < 1.495585217958292e-2) {
            padeCoefficients = PADE_COEFFICIENTS_3;
        } else if (l1Norm < 2.539398330063230e-1) {
            padeCoefficients = PADE_COEFFICIENTS_5;
        } else if (l1Norm < 9.504178996162932e-1) {
            padeCoefficients = PADE_COEFFICIENTS_7;
        } else if (l1Norm < 2.097847961257068) {
            padeCoefficients = PADE_COEFFICIENTS_9;
        } else {
            padeCoefficients = PADE_COEFFICIENTS_13;

            // Calculate scaling factor
            final double normScale = 5.371920351148152;
            squaringCount = Math.max(0, Math.getExponent(l1Norm / normScale));

            // Scale matrix by power of 2
            final int finalSquaringCount = squaringCount;
            scaledMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int row, int column, double value) {
                    return Math.scalb(value, -finalSquaringCount);
                }
            });
        }

        // Calculate U and V using Horner
        // See Golub, Gene H., and Charles F. van Loan. Matrix Computations. 4th ed.
        // John Hopkins University Press, 2013.  pages 530/531
        final RealMatrix scaledMatrix2 = scaledMatrix.multiply(scaledMatrix);
        final int coeffLength = padeCoefficients.length;

        // Calculate V
        RealMatrix padeV = MatrixUtils.createRealMatrix(dim, dim);
        for (int i = coeffLength - 1; i > 1; i -= 2) {
            padeV = scaledMatrix2.multiply(padeV.add(identity.scalarMultiply(padeCoefficients[i])));
        }
        padeV = scaledMatrix.multiply(padeV.add(identity.scalarMultiply(padeCoefficients[1])));

        // Calculate U
        RealMatrix padeU = MatrixUtils.createRealMatrix(dim, dim);
        for (int i = coeffLength - 2; i > 1; i -= 2) {
            padeU = scaledMatrix2.multiply(padeU.add(identity.scalarMultiply(padeCoefficients[i])));
        }
        padeU = padeU.add(identity.scalarMultiply(padeCoefficients[0]));

        // Calculate pade approximate by solving (U-V) F = (U+V) for F
        RealMatrix padeNumer = padeU.add(padeV);
        RealMatrix padeDenom = padeU.subtract(padeV);

        // Calculate the matrix ratio
        QRDecomposition decomposition = new QRDecomposition(padeDenom);
        RealMatrix result = decomposition.getSolver().solve(padeNumer);

        // Repeated squaring if matrix was scaled
        for (int i = 0; i < squaringCount; i++) {
            result = result.multiply(result);
        }

        // Undo preprocessing
        result = result.scalarMultiply(Math.exp(preprocessScale));

        return result;
    }

    /** Orthonormalize a list of vectors.
     * <p>
     * Orthonormalization is performed by using the Modified Gram-Schmidt process.
     * </p>
     * @param independent list of independent vectors
     * @param threshold projected vectors with a norm less than or equal to this threshold
     * are considered to have zero norm, hence the vectors they come from are not independent from
     * previous vectors
     * @param dependentVectorsHandler handler for dependent vectors
     * @return orthonormal basis having the same span as {@code independent}
     * @since 2.1
     */
    public static List<RealVector> orthonormalize(final List<RealVector> independent,
                                                  final double threshold, final DependentVectorsHandler handler) {

        // create separate list
        final List<RealVector> basis = new ArrayList<>(independent);

        // loop over basis vectors
        int index = 0;
        while (index < basis.size()) {

            // check dependency
            final RealVector vi = basis.get(index);
            final double norm = vi.getNorm();
            if (norm <= threshold) {
                // the current vector is dependent from the previous ones
                index = handler.manageDependent(index, basis);
            } else {

                // normalize basis vector in place
                vi.mapDivideToSelf(vi.getNorm());

                // project remaining vectors in place
                for (int j = index + 1; j < basis.size(); ++j) {
                    final RealVector vj  = basis.get(j);
                    final double     dot = vi.dotProduct(vj);
                    for (int k = 0; k < vj.getDimension(); ++k) {
                        vj.setEntry(k, vj.getEntry(k) - dot * vi.getEntry(k));
                    }
                }

                ++index;

            }

        }

        return basis;

    }

    /** Orthonormalize a list of vectors.
     * <p>
     * Orthonormalization is performed by using the Modified Gram-Schmidt process.
     * </p>
     * @param <T> type of the field elements
     * @param independent list of independent vectors
     * @param threshold projected vectors with a norm less than or equal to this threshold
     * are considered to have zero norm, hence the vectors they come from are not independent from
     * previous vectors
     * @param dependentVectorsHandler handler for dependent vectors
     * @return orthonormal basis having the same span as {@code independent}
     * @since 2.1
     */
    public static <T extends CalculusFieldElement<T>> List<FieldVector<T>> orthonormalize(final Field<T> field,
                                                                                          final List<FieldVector<T>> independent,
                                                                                          final T threshold,
                                                                                          final DependentVectorsHandler handler) {

        // create separate list
        final List<FieldVector<T>> basis = new ArrayList<>(independent);

        // loop over basis vectors
        int index = 0;
        while (index < basis.size()) {

            // check dependency
            final FieldVector<T> vi = basis.get(index);
            final T norm = vi.dotProduct(vi).sqrt();
            if (norm.subtract(threshold).getReal() <= 0) {
                // the current vector is dependent from the previous ones
                index = handler.manageDependent(field, index, basis);
            } else {

                // normalize basis vector in place
                vi.mapDivideToSelf(norm);

                // project remaining vectors in place
                for (int j = index + 1; j < basis.size(); ++j) {
                    final FieldVector<T> vj  = basis.get(j);
                    final T              dot = vi.dotProduct(vj);
                    for (int k = 0; k < vj.getDimension(); ++k) {
                        vj.setEntry(k, vj.getEntry(k).subtract(dot.multiply(vi.getEntry(k))));
                    }
                }

                ++index;

            }
        }

        return basis;

    }

}
