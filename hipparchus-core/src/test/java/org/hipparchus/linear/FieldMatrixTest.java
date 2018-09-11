/*
 * Licensed to the Hipparchus project under one or more
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
package org.hipparchus.linear;

import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.junit.Assert;
import org.junit.Test;

public class FieldMatrixTest {

    @Test
    public void testDefaultMultiplyTransposed() {
        FieldMatrix<Decimal64> a = createMatrix(new double[][] { {1d,2d,3d}, {2d,5d,3d}, {1d,0d,8d} });
        FieldMatrix<Decimal64> b = createMatrix(new double[][] { {4d, -5d, 6d} });
        FieldMatrix<Decimal64> abTRef = a.multiplyTransposed(b);
        DefaultMatrix dma = new DefaultMatrix(a);
        DefaultMatrix dmb = new DefaultMatrix(b);
        FieldMatrix<Decimal64> abT = dma.multiplyTransposed(dmb);
        FieldMatrix<Decimal64> diff = abT.subtract(abTRef);
        for (int i = 0; i < diff.getRowDimension(); ++i) {
            for (int j = 0; j < diff.getColumnDimension(); ++j) {
                Assert.assertEquals(0.0, diff.getEntry(i, j).doubleValue(), 1.0e-10);
            }
        }
    }

    @Test
    public void testDefaultTransposeMultiply() {
        FieldMatrix<Decimal64> a = createMatrix(new double[][] { {1d,2d,3d}, {2d,5d,3d}, {1d,0d,8d} });
        FieldMatrix<Decimal64> b = createMatrix(new double[][] { {4d}, {-5d}, {6d} });
        FieldMatrix<Decimal64> aTbRef = a.transposeMultiply(b);
        DefaultMatrix dma = new DefaultMatrix(a);
        DefaultMatrix dmb = new DefaultMatrix(b);
        FieldMatrix<Decimal64> aTb = dma.transposeMultiply(dmb);
        FieldMatrix<Decimal64> diff = aTb.subtract(aTbRef);
        for (int i = 0; i < diff.getRowDimension(); ++i) {
            for (int j = 0; j < diff.getColumnDimension(); ++j) {
                Assert.assertEquals(0.0, diff.getEntry(i, j).doubleValue(), 1.0e-10);
            }
        }
    }

    // local class that does NOT override multiplyTransposed nor transposeMultiply
    // so the default methods are called
    private class DefaultMatrix extends AbstractFieldMatrix<Decimal64> {

        FieldMatrix<Decimal64> m;
        public DefaultMatrix(FieldMatrix<Decimal64> m) {
            super(Decimal64Field.getInstance());
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
        public FieldMatrix<Decimal64> createMatrix(int rowDimension, int columnDimension) {
            return m.createMatrix(rowDimension, columnDimension);
        }

        @Override
        public FieldMatrix<Decimal64> copy() {
            return m.copy();
        }

        @Override
        public void setEntry(int row, int column, Decimal64 value) {
            m.setEntry(row, column, value);
        }

        @Override
        public void addToEntry(int row, int column, Decimal64 increment) {
            m.addToEntry(row, column, increment);
        }

        @Override
        public void multiplyEntry(int row, int column, Decimal64 factor) {
            m.multiplyEntry(row, column, factor);
        }

        @Override
        public Decimal64 getEntry(int row, int column) {
            return m.getEntry(row, column);
        }

    }

    private FieldMatrix<Decimal64> createMatrix(double[][] a) {
        FieldMatrix<Decimal64> m = MatrixUtils.createFieldMatrix(Decimal64Field.getInstance(),
                                                                 a.length, a[0].length);
        for (int i = 0; i < m.getRowDimension(); ++i) {
            for (int j = 0; j < m.getColumnDimension(); ++j) {
                m.setEntry(i, j, new Decimal64(a[i][j]));
            }
        }
        return m;
    }

}
