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

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link DiagonalMatrix} class.
 */
class DiagonalMatrixTest {
    @Test
    void testConstructor1() {
        final int dim = 3;
        final DiagonalMatrix m = new DiagonalMatrix(dim);
        assertEquals(dim, m.getRowDimension());
        assertEquals(dim, m.getColumnDimension());
    }

    @Test
    void testConstructor2() {
        final double[] d = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(d);
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j = 0; j < m.getRowDimension(); j++) {
                if (i == j) {
                    assertEquals(d[i], m.getEntry(i, j), 0d);
                } else {
                    assertEquals(0d, m.getEntry(i, j), 0d);
                }
            }
        }

        // Check that the underlying was copied.
        d[0] = 0;
        assertFalse(d[0] == m.getEntry(0, 0));
    }

    @Test
    void testConstructor3() {
        final double[] d = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(d, false);
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j = 0; j < m.getRowDimension(); j++) {
                if (i == j) {
                    assertEquals(d[i], m.getEntry(i, j), 0d);
                } else {
                    assertEquals(0d, m.getEntry(i, j), 0d);
                }
            }
        }

        // Check that the underlying is referenced.
        d[0] = 0;
        assertEquals(d[0], m.getEntry(0, 0));

    }

    @Test
    void testCreateError() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final double[] d = {-1.2, 3.4, 5};
            final DiagonalMatrix m = new DiagonalMatrix(d, false);
            m.createMatrix(5, 3);
        });
    }

    @Test
    void testCreate() {
        final double[] d = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(d, false);
        final RealMatrix p = m.createMatrix(5, 5);
        assertTrue(p instanceof DiagonalMatrix);
        assertEquals(5, p.getRowDimension());
        assertEquals(5, p.getColumnDimension());
    }

    @Test
    void testCopy() {
        final double[] d = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(d, false);
        final DiagonalMatrix p = (DiagonalMatrix) m.copy();
        for (int i = 0; i < m.getRowDimension(); ++i) {
            assertEquals(m.getEntry(i, i), p.getEntry(i, i), 1.0e-20);
        }
    }

    @Test
    void testGetData() {
        final double[] data = { -1.2, 3.4, 5 };
        final int dim = 3;
        final DiagonalMatrix m = new DiagonalMatrix(dim);
        for (int i = 0; i < dim; i++) {
            m.setEntry(i, i, data[i]);
        }

        final double[][] out = m.getData();
        assertEquals(dim, out.length);
        for (int i = 0; i < m.getRowDimension(); i++) {
            assertEquals(dim, out[i].length);
            for (int j = 0; j < m.getRowDimension(); j++) {
                if (i == j) {
                    assertEquals(data[i], out[i][j], 0d);
                } else {
                    assertEquals(0d, out[i][j], 0d);
                }
            }
        }
    }

    @Test
    void testAdd() {
        final double[] data1 = { -1.2, 3.4, 5 };
        final DiagonalMatrix m1 = new DiagonalMatrix(data1);

        final double[] data2 = { 10.1, 2.3, 45 };
        final DiagonalMatrix m2 = new DiagonalMatrix(data2);

        final DiagonalMatrix result = m1.add(m2);
        assertEquals(m1.getRowDimension(), result.getRowDimension());
        for (int i = 0; i < result.getRowDimension(); i++) {
            for (int j = 0; j < result.getRowDimension(); j++) {
                if (i == j) {
                    assertEquals(data1[i] + data2[i], result.getEntry(i, j), 0d);
                } else {
                    assertEquals(0d, result.getEntry(i, j), 0d);
                }
            }
        }
    }

    @Test
    void testSubtract() {
        final double[] data1 = { -1.2, 3.4, 5 };
        final DiagonalMatrix m1 = new DiagonalMatrix(data1);

        final double[] data2 = { 10.1, 2.3, 45 };
        final DiagonalMatrix m2 = new DiagonalMatrix(data2);

        final DiagonalMatrix result = m1.subtract(m2);
        assertEquals(m1.getRowDimension(), result.getRowDimension());
        for (int i = 0; i < result.getRowDimension(); i++) {
            for (int j = 0; j < result.getRowDimension(); j++) {
                if (i == j) {
                    assertEquals(data1[i] - data2[i], result.getEntry(i, j), 0d);
                } else {
                    assertEquals(0d, result.getEntry(i, j), 0d);
                }
            }
        }
    }

    @Test
    void testAddToEntry() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(data);

        for (int i = 0; i < m.getRowDimension(); i++) {
            m.addToEntry(i, i, i);
            assertEquals(data[i] + i, m.getEntry(i, i), 0d);
        }
    }

    @Test
    void testMultiplyEntry() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(data);

        for (int i = 0; i < m.getRowDimension(); i++) {
            m.multiplyEntry(i, i, i);
            assertEquals(data[i] * i, m.getEntry(i, i), 0d);
        }
    }

    @Test
    void testMultiply1() {
        final double[] data1 = { -1.2, 3.4, 5 };
        final DiagonalMatrix m1 = new DiagonalMatrix(data1);
        final double[] data2 = { 10.1, 2.3, 45 };
        final DiagonalMatrix m2 = new DiagonalMatrix(data2);

        final DiagonalMatrix result = (DiagonalMatrix) m1.multiply((RealMatrix) m2);
        assertEquals(m1.getRowDimension(), result.getRowDimension());
        for (int i = 0; i < result.getRowDimension(); i++) {
            for (int j = 0; j < result.getRowDimension(); j++) {
                if (i == j) {
                    assertEquals(data1[i] * data2[i], result.getEntry(i, j), 0d);
                } else {
                    assertEquals(0d, result.getEntry(i, j), 0d);
                }
            }
        }
    }

    @Test
    void testMultiply2() {
        final double[] data1 = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag1 = new DiagonalMatrix(data1);

        final double[][] data2 = { { -1.2, 3.4 },
                                   { -5.6, 7.8 },
                                   {  9.1, 2.3 } };
        final RealMatrix dense2 = new Array2DRowRealMatrix(data2);
        final RealMatrix dense1 = new Array2DRowRealMatrix(diag1.getData());

        final RealMatrix diagResult = diag1.multiply(dense2);
        final RealMatrix denseResult = dense1.multiply(dense2);

        for (int i = 0; i < dense1.getRowDimension(); i++) {
            for (int j = 0; j < dense2.getColumnDimension(); j++) {
                assertEquals(denseResult.getEntry(i, j),
                                    diagResult.getEntry(i, j), 0d);
            }
        }
    }

    @Test
    void testMultiplyTransposedDiagonalMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x4b20cb5a0440c929l);
        for (int rows = 1; rows <= 64; rows += 7) {
            final DiagonalMatrix a = new DiagonalMatrix(rows);
            for (int i = 0; i < rows; ++i) {
                a.setEntry(i, i, randomGenerator.nextDouble());
            }
            final DiagonalMatrix b = new DiagonalMatrix(rows);
            for (int i = 0; i < rows; ++i) {
                b.setEntry(i, i, randomGenerator.nextDouble());
            }
            assertEquals(0.0,
                                a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).getNorm1(),
                                1.0e-15);
        }
    }

    @Test
    void testMultiplyTransposedArray2DRowRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x0fa7b97d4826cd43l);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            final DiagonalMatrix a = new DiagonalMatrix(rows);
            for (int i = 0; i < rows; ++i) {
                a.setEntry(i, i, randomGenerator.nextDouble());
            }
            for (int interm = 1; interm <= 64; interm += 7) {
                final Array2DRowRealMatrix b = new Array2DRowRealMatrix(interm, rows);
                b.walkInOptimizedOrder(randomSetter);
                assertEquals(0.0,
                                    a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).getNorm1(),
                                    1.0e-15);
            }
        }
    }

    @Test
    void testMultiplyTransposedWrongDimensions() {
        try {
            new DiagonalMatrix(3).multiplyTransposed(new DiagonalMatrix(2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testTransposeMultiplyDiagonalMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x4b20cb5a0440c929l);
        for (int rows = 1; rows <= 64; rows += 7) {
            final DiagonalMatrix a = new DiagonalMatrix(rows);
            for (int i = 0; i < rows; ++i) {
                a.setEntry(i, i, randomGenerator.nextDouble());
            }
            final DiagonalMatrix b = new DiagonalMatrix(rows);
            for (int i = 0; i < rows; ++i) {
                b.setEntry(i, i, randomGenerator.nextDouble());
            }
            assertEquals(0.0,
                                a.transposeMultiply(b).subtract(a.transpose().multiply(b)).getNorm1(),
                                1.0e-15);
        }
    }

    @Test
    void testTransposeMultiplyArray2DRowRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x0fa7b97d4826cd43l);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            final DiagonalMatrix a = new DiagonalMatrix(rows);
            for (int i = 0; i < rows; ++i) {
                a.setEntry(i, i, randomGenerator.nextDouble());
            }
            for (int interm = 1; interm <= 64; interm += 7) {
                final Array2DRowRealMatrix b = new Array2DRowRealMatrix(rows, interm);
                b.walkInOptimizedOrder(randomSetter);
                assertEquals(0.0,
                                    a.transposeMultiply(b).subtract(a.transpose().multiply(b)).getNorm1(),
                                    1.0e-15);
            }
        }
    }

    @Test
    void testTransposeMultiplyWrongDimensions() {
        try {
            new DiagonalMatrix(3).transposeMultiply(new DiagonalMatrix(2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testOperate() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        final RealMatrix dense = new Array2DRowRealMatrix(diag.getData());

        final double[] v = { 6.7, 890.1, 23.4 };
        final double[] diagResult = diag.operate(v);
        final double[] denseResult = dense.operate(v);

        UnitTestUtils.customAssertEquals(diagResult, denseResult, 0d);
    }

    @Test
    void testPreMultiply() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        final RealMatrix dense = new Array2DRowRealMatrix(diag.getData());

        final double[] v = { 6.7, 890.1, 23.4 };
        final double[] diagResult = diag.preMultiply(v);
        final double[] denseResult = dense.preMultiply(v);

        UnitTestUtils.customAssertEquals(diagResult, denseResult, 0d);
    }

    @Test
    void testPreMultiplyVector() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        final RealMatrix dense = new Array2DRowRealMatrix(diag.getData());

        final double[] v = { 6.7, 890.1, 23.4 };
        final RealVector vector = MatrixUtils.createRealVector(v);
        final RealVector diagResult = diag.preMultiply(vector);
        final RealVector denseResult = dense.preMultiply(vector);

        UnitTestUtils.customAssertEquals("preMultiply(Vector) returns wrong result", diagResult, denseResult, 0d);
    }

    @Test
    void testSetNonDiagonalEntry() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final DiagonalMatrix diag = new DiagonalMatrix(3);
            diag.setEntry(1, 2, 3.4);
        });
    }

    @Test
    void testSetNonDiagonalZero() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.setEntry(1, 2, 0.0);
        assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test
    void testAddNonDiagonalEntry() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final DiagonalMatrix diag = new DiagonalMatrix(3);
            diag.addToEntry(1, 2, 3.4);
        });
    }

    @Test
    void testAddNonDiagonalZero() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.addToEntry(1, 2, 0.0);
        assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test
    void testMultiplyNonDiagonalEntry() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.multiplyEntry(1, 2, 3.4);
        assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test
    void testMultiplyNonDiagonalZero() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.multiplyEntry(1, 2, 0.0);
        assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test
    void testSetEntryOutOfRange() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final DiagonalMatrix diag = new DiagonalMatrix(3);
            diag.setEntry(3, 3, 3.4);
        });
    }

    @Test
    void testNull() {
        assertThrows(NullArgumentException.class, () -> {
            new DiagonalMatrix(null, false);
        });
    }

    @Test
    void testSetSubMatrixError() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final double[] data = {-1.2, 3.4, 5};
            final DiagonalMatrix diag = new DiagonalMatrix(data);
            diag.setSubMatrix(new double[][]{{1.0, 1.0}, {1.0, 1.0}}, 1, 1);
        });
    }

    @Test
    void testSetSubMatrix() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        diag.setSubMatrix(new double[][] { {0.0, 5.0, 0.0}, {0.0, 0.0, 6.0}}, 1, 0);
        assertEquals(-1.2, diag.getEntry(0, 0), 1.0e-20);
        assertEquals( 5.0, diag.getEntry(1, 1), 1.0e-20);
        assertEquals( 6.0, diag.getEntry(2, 2), 1.0e-20);
    }

    @Test
    void testInverseError() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final double[] data = {1, 2, 0};
            final DiagonalMatrix diag = new DiagonalMatrix(data);
            diag.inverse();
        });
    }

    @Test
    void testInverseError2() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final double[] data = {1, 2, 1e-6};
            final DiagonalMatrix diag = new DiagonalMatrix(data);
            diag.inverse(1e-5);
        });
    }

    @Test
    void testInverse() {
        final double[] data = { 1, 2, 3 };
        final DiagonalMatrix m = new DiagonalMatrix(data);
        final DiagonalMatrix inverse = m.inverse();

        final DiagonalMatrix result = m.multiply(inverse);
        UnitTestUtils.customAssertEquals("DiagonalMatrix.inverse() returns wrong result",
                                         MatrixUtils.createRealIdentityMatrix(data.length), result, Math.ulp(1d));
    }

}
