/*
 * Licensed to the Hipparchus project under one or more
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
package org.hipparchus.linear;

import org.hipparchus.Field;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.junit.Assert;
import org.junit.Test;

public class FieldMatrixTest {

    @Test
    public void testDefaultMultiplyTransposed() {
        FieldMatrix<Binary64> a = createMatrix(new double[][] { {1d,2d,3d}, {2d,5d,3d}, {1d,0d,8d} });
        FieldMatrix<Binary64> b = createMatrix(new double[][] { {4d, -5d, 6d} });
        FieldMatrix<Binary64> abTRef = a.multiplyTransposed(b);
        DefaultMatrix dma = new DefaultMatrix(a);
        DefaultMatrix dmb = new DefaultMatrix(b);
        FieldMatrix<Binary64> abT = dma.multiplyTransposed(dmb);
        FieldMatrix<Binary64> diff = abT.subtract(abTRef);
        for (int i = 0; i < diff.getRowDimension(); ++i) {
            for (int j = 0; j < diff.getColumnDimension(); ++j) {
                Assert.assertEquals(0.0, diff.getEntry(i, j).doubleValue(), 1.0e-10);
            }
        }
    }

    @Test
    public void testDefaultTransposeMultiply() {
        FieldMatrix<Binary64> a = createMatrix(new double[][] { {1d,2d,3d}, {2d,5d,3d}, {1d,0d,8d} });
        FieldMatrix<Binary64> b = createMatrix(new double[][] { {4d}, {-5d}, {6d} });
        FieldMatrix<Binary64> aTbRef = a.transposeMultiply(b);
        DefaultMatrix dma = new DefaultMatrix(a);
        DefaultMatrix dmb = new DefaultMatrix(b);
        FieldMatrix<Binary64> aTb = dma.transposeMultiply(dmb);
        FieldMatrix<Binary64> diff = aTb.subtract(aTbRef);
        for (int i = 0; i < diff.getRowDimension(); ++i) {
            for (int j = 0; j < diff.getColumnDimension(); ++j) {
                Assert.assertEquals(0.0, diff.getEntry(i, j).doubleValue(), 1.0e-10);
            }
        }
    }

    @Test
    public void testDefaultMap() {
        FieldMatrix<Binary64> a = createMatrix(new double[][] { {1d,2d,3d}, {2d,5d,3d}, {1d,0d,8d} });
        FieldMatrix<Binary64> result = a.add(a.map(x -> x.negate()));
        result.walkInOptimizedOrder(new FieldMatrixPreservingVisitor<Binary64>() {
            
            @Override
            public void visit(int row, int column, Binary64 value) {
                Assert.assertEquals(0.0, value.getReal(), 1.0e-10);
            }
            
            @Override
            public void start(int rows, int columns, int startRow, int endRow,
                              int startColumn, int endColumn) {
            }
            
            @Override
            public Binary64 end() {
                return Binary64Field.getInstance().getZero();
            }
        });
    }

    @Test
    public void testArithmeticalBlending() {
        // Given
        final Field<Binary64> field = Binary64Field.getInstance();

        final FieldMatrix<Binary64> matrix1 = MatrixUtils.createFieldMatrix(field, 2, 2);
        matrix1.setEntry(0, 0, new Binary64(1));
        matrix1.setEntry(0, 1, new Binary64(2));
        matrix1.setEntry(1, 0, new Binary64(3));
        matrix1.setEntry(1, 1, new Binary64(4));

        final FieldMatrix<Binary64> matrix2 = MatrixUtils.createFieldMatrix(field, 2, 2);
        matrix2.setEntry(0, 0, new Binary64(2));
        matrix2.setEntry(0, 1, new Binary64(4));
        matrix2.setEntry(1, 0, new Binary64(9));
        matrix2.setEntry(1, 1, new Binary64(16));

        final Binary64 blendingValue = new Binary64(0.65);

        // When
        final FieldMatrix<Binary64> blendedMatrix = matrix1.blendArithmeticallyWith(matrix2, blendingValue);

        // Then
        Assert.assertEquals(1.65 , blendedMatrix.getEntry(0,0).getReal(), 1.0e-15);
        Assert.assertEquals(3.3  , blendedMatrix.getEntry(0,1).getReal(), 1.0e-15);
        Assert.assertEquals(6.9  , blendedMatrix.getEntry(1,0).getReal(), 1.0e-15);
        Assert.assertEquals(11.8 , blendedMatrix.getEntry(1,1).getReal(), 1.0e-15);
    }

    // local class that does NOT override multiplyTransposed nor transposeMultiply nor map nor mapToSelf
    // so the default methods are called
    private class DefaultMatrix extends AbstractFieldMatrix<Binary64> {

        FieldMatrix<Binary64> m;
        public DefaultMatrix(FieldMatrix<Binary64> m) {
            super(Binary64Field.getInstance());
            this.m = m;
        }

        @Override
        public int getRowDimension() {
            return m.getRowDimension();
        }

        @Override
        public int getColumnDimension() {
            return m.getColumnDimension();
        }

        @Override
        public FieldMatrix<Binary64> createMatrix(int rowDimension, int columnDimension) {
            return m.createMatrix(rowDimension, columnDimension);
        }

        @Override
        public FieldMatrix<Binary64> copy() {
            return m.copy();
        }

        @Override
        public void setEntry(int row, int column, Binary64 value) {
            m.setEntry(row, column, value);
        }

        @Override
        public void addToEntry(int row, int column, Binary64 increment) {
            m.addToEntry(row, column, increment);
        }

        @Override
        public void multiplyEntry(int row, int column, Binary64 factor) {
            m.multiplyEntry(row, column, factor);
        }

        @Override
        public Binary64 getEntry(int row, int column) {
            return m.getEntry(row, column);
        }

    }

    private FieldMatrix<Binary64> createMatrix(double[][] a) {
        FieldMatrix<Binary64> m = MatrixUtils.createFieldMatrix(Binary64Field.getInstance(),
                                                                 a.length, a[0].length);
        for (int i = 0; i < m.getRowDimension(); ++i) {
            for (int j = 0; j < m.getColumnDimension(); ++j) {
                m.setEntry(i, j, new Binary64(a[i][j]));
            }
        }
        return m;
    }

}
