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

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test cases for the {@link DiagonalMatrix} class.
 */
public class DiagonalMatrixTest {
    @Test
    public void testConstructor1() {
        final int dim = 3;
        final DiagonalMatrix m = new DiagonalMatrix(dim);
        Assert.assertEquals(dim, m.getRowDimension());
        Assert.assertEquals(dim, m.getColumnDimension());
    }

    @Test
    public void testConstructor2() {
        final double[] d = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(d);
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j = 0; j < m.getRowDimension(); j++) {
                if (i == j) {
                    Assert.assertEquals(d[i], m.getEntry(i, j), 0d);
                } else {
                    Assert.assertEquals(0d, m.getEntry(i, j), 0d);
                }
            }
        }

        // Check that the underlying was copied.
        d[0] = 0;
        Assert.assertFalse(d[0] == m.getEntry(0, 0));
    }

    @Test
    public void testConstructor3() {
        final double[] d = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(d, false);
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j = 0; j < m.getRowDimension(); j++) {
                if (i == j) {
                    Assert.assertEquals(d[i], m.getEntry(i, j), 0d);
                } else {
                    Assert.assertEquals(0d, m.getEntry(i, j), 0d);
                }
            }
        }

        // Check that the underlying is referenced.
        d[0] = 0;
        Assert.assertTrue(d[0] == m.getEntry(0, 0));

    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testCreateError() {
        final double[] d = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(d, false);
        m.createMatrix(5, 3);
    }

    @Test
    public void testCreate() {
        final double[] d = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(d, false);
        final RealMatrix p = m.createMatrix(5, 5);
        Assert.assertTrue(p instanceof DiagonalMatrix);
        Assert.assertEquals(5, p.getRowDimension());
        Assert.assertEquals(5, p.getColumnDimension());
    }

    @Test
    public void testCopy() {
        final double[] d = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(d, false);
        final DiagonalMatrix p = (DiagonalMatrix) m.copy();
        for (int i = 0; i < m.getRowDimension(); ++i) {
            Assert.assertEquals(m.getEntry(i, i), p.getEntry(i, i), 1.0e-20);
        }
    }

    @Test
    public void testGetData() {
        final double[] data = { -1.2, 3.4, 5 };
        final int dim = 3;
        final DiagonalMatrix m = new DiagonalMatrix(dim);
        for (int i = 0; i < dim; i++) {
            m.setEntry(i, i, data[i]);
        }

        final double[][] out = m.getData();
        Assert.assertEquals(dim, out.length);
        for (int i = 0; i < m.getRowDimension(); i++) {
            Assert.assertEquals(dim, out[i].length);
            for (int j = 0; j < m.getRowDimension(); j++) {
                if (i == j) {
                    Assert.assertEquals(data[i], out[i][j], 0d);
                } else {
                    Assert.assertEquals(0d, out[i][j], 0d);
                }
            }
        }
    }

    @Test
    public void testAdd() {
        final double[] data1 = { -1.2, 3.4, 5 };
        final DiagonalMatrix m1 = new DiagonalMatrix(data1);

        final double[] data2 = { 10.1, 2.3, 45 };
        final DiagonalMatrix m2 = new DiagonalMatrix(data2);

        final DiagonalMatrix result = m1.add(m2);
        Assert.assertEquals(m1.getRowDimension(), result.getRowDimension());
        for (int i = 0; i < result.getRowDimension(); i++) {
            for (int j = 0; j < result.getRowDimension(); j++) {
                if (i == j) {
                    Assert.assertEquals(data1[i] + data2[i], result.getEntry(i, j), 0d);
                } else {
                    Assert.assertEquals(0d, result.getEntry(i, j), 0d);
                }
            }
        }
    }

    @Test
    public void testSubtract() {
        final double[] data1 = { -1.2, 3.4, 5 };
        final DiagonalMatrix m1 = new DiagonalMatrix(data1);

        final double[] data2 = { 10.1, 2.3, 45 };
        final DiagonalMatrix m2 = new DiagonalMatrix(data2);

        final DiagonalMatrix result = m1.subtract(m2);
        Assert.assertEquals(m1.getRowDimension(), result.getRowDimension());
        for (int i = 0; i < result.getRowDimension(); i++) {
            for (int j = 0; j < result.getRowDimension(); j++) {
                if (i == j) {
                    Assert.assertEquals(data1[i] - data2[i], result.getEntry(i, j), 0d);
                } else {
                    Assert.assertEquals(0d, result.getEntry(i, j), 0d);
                }
            }
        }
    }

    @Test
    public void testAddToEntry() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(data);

        for (int i = 0; i < m.getRowDimension(); i++) {
            m.addToEntry(i, i, i);
            Assert.assertEquals(data[i] + i, m.getEntry(i, i), 0d);
        }
    }

    @Test
    public void testMultiplyEntry() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix m = new DiagonalMatrix(data);

        for (int i = 0; i < m.getRowDimension(); i++) {
            m.multiplyEntry(i, i, i);
            Assert.assertEquals(data[i] * i, m.getEntry(i, i), 0d);
        }
    }

    @Test
    public void testMultiply1() {
        final double[] data1 = { -1.2, 3.4, 5 };
        final DiagonalMatrix m1 = new DiagonalMatrix(data1);
        final double[] data2 = { 10.1, 2.3, 45 };
        final DiagonalMatrix m2 = new DiagonalMatrix(data2);

        final DiagonalMatrix result = (DiagonalMatrix) m1.multiply((RealMatrix) m2);
        Assert.assertEquals(m1.getRowDimension(), result.getRowDimension());
        for (int i = 0; i < result.getRowDimension(); i++) {
            for (int j = 0; j < result.getRowDimension(); j++) {
                if (i == j) {
                    Assert.assertEquals(data1[i] * data2[i], result.getEntry(i, j), 0d);
                } else {
                    Assert.assertEquals(0d, result.getEntry(i, j), 0d);
                }
            }
        }
    }

    @Test
    public void testMultiply2() {
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
                Assert.assertEquals(denseResult.getEntry(i, j),
                                    diagResult.getEntry(i, j), 0d);
            }
        }
    }

    @Test
    public void testMultiplyTransposedDiagonalMatrix() {
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
            Assert.assertEquals(0.0,
                                a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).getNorm(),
                                1.0e-15);
        }
    }

    @Test
    public void testMultiplyTransposedArray2DRowRealMatrix() {
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
                Assert.assertEquals(0.0,
                                    a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).getNorm(),
                                    1.0e-15);
            }
        }
    }

    @Test
    public void testMultiplyTransposedWrongDimensions() {
        try {
            new DiagonalMatrix(3).multiplyTransposed(new DiagonalMatrix(2));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            Assert.assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    public void testTransposeMultiplyDiagonalMatrix() {
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
            Assert.assertEquals(0.0,
                                a.transposeMultiply(b).subtract(a.transpose().multiply(b)).getNorm(),
                                1.0e-15);
        }
    }

    @Test
    public void testTransposeMultiplyArray2DRowRealMatrix() {
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
                Assert.assertEquals(0.0,
                                    a.transposeMultiply(b).subtract(a.transpose().multiply(b)).getNorm(),
                                    1.0e-15);
            }
        }
    }

    @Test
    public void testTransposeMultiplyWrongDimensions() {
        try {
            new DiagonalMatrix(3).transposeMultiply(new DiagonalMatrix(2));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            Assert.assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    public void testOperate() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        final RealMatrix dense = new Array2DRowRealMatrix(diag.getData());

        final double[] v = { 6.7, 890.1, 23.4 };
        final double[] diagResult = diag.operate(v);
        final double[] denseResult = dense.operate(v);

        UnitTestUtils.assertEquals(diagResult, denseResult, 0d);
    }

    @Test
    public void testPreMultiply() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        final RealMatrix dense = new Array2DRowRealMatrix(diag.getData());

        final double[] v = { 6.7, 890.1, 23.4 };
        final double[] diagResult = diag.preMultiply(v);
        final double[] denseResult = dense.preMultiply(v);

        UnitTestUtils.assertEquals(diagResult, denseResult, 0d);
    }

    @Test
    public void testPreMultiplyVector() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        final RealMatrix dense = new Array2DRowRealMatrix(diag.getData());

        final double[] v = { 6.7, 890.1, 23.4 };
        final RealVector vector = MatrixUtils.createRealVector(v);
        final RealVector diagResult = diag.preMultiply(vector);
        final RealVector denseResult = dense.preMultiply(vector);

        UnitTestUtils.assertEquals("preMultiply(Vector) returns wrong result", diagResult, denseResult, 0d);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testSetNonDiagonalEntry() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.setEntry(1, 2, 3.4);
    }

    @Test
    public void testSetNonDiagonalZero() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.setEntry(1, 2, 0.0);
        Assert.assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testAddNonDiagonalEntry() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.addToEntry(1, 2, 3.4);
    }

    @Test
    public void testAddNonDiagonalZero() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.addToEntry(1, 2, 0.0);
        Assert.assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test
    public void testMultiplyNonDiagonalEntry() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.multiplyEntry(1, 2, 3.4);
        Assert.assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test
    public void testMultiplyNonDiagonalZero() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.multiplyEntry(1, 2, 0.0);
        Assert.assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testSetEntryOutOfRange() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.setEntry(3, 3, 3.4);
    }

    @Test(expected=NullArgumentException.class)
    public void testNull() {
        new DiagonalMatrix(null, false);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testSetSubMatrixError() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        diag.setSubMatrix(new double[][] { {1.0, 1.0}, {1.0, 1.0}}, 1, 1);
    }

    @Test
    public void testSetSubMatrix() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        diag.setSubMatrix(new double[][] { {0.0, 5.0, 0.0}, {0.0, 0.0, 6.0}}, 1, 0);
        Assert.assertEquals(-1.2, diag.getEntry(0, 0), 1.0e-20);
        Assert.assertEquals( 5.0, diag.getEntry(1, 1), 1.0e-20);
        Assert.assertEquals( 6.0, diag.getEntry(2, 2), 1.0e-20);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testInverseError() {
        final double[] data = { 1, 2, 0 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        diag.inverse();
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testInverseError2() {
        final double[] data = { 1, 2, 1e-6 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        diag.inverse(1e-5);
    }

    @Test
    public void testInverse() {
        final double[] data = { 1, 2, 3 };
        final DiagonalMatrix m = new DiagonalMatrix(data);
        final DiagonalMatrix inverse = m.inverse();

        final DiagonalMatrix result = m.multiply(inverse);
        UnitTestUtils.assertEquals("DiagonalMatrix.inverse() returns wrong result",
                MatrixUtils.createRealIdentityMatrix(data.length), result, Math.ulp(1d));
    }

}
