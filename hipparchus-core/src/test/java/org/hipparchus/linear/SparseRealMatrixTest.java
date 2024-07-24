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
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link OpenMapRealMatrix} class.
 *
 */
final class SparseRealMatrixTest {

    // 3 x 3 identity matrix
    protected double[][] id = { { 1d, 0d, 0d }, { 0d, 1d, 0d }, { 0d, 0d, 1d } };
    // Test data for group operations
    protected double[][] testData = { { 1d, 2d, 3d }, { 2d, 5d, 3d },
            { 1d, 0d, 8d } };
    protected double[][] testDataLU = { { 2d, 5d, 3d }, { .5d, -2.5d, 6.5d },
            { 0.5d, 0.2d, .2d } };
    protected double[][] testDataPlus2 = { { 3d, 4d, 5d }, { 4d, 7d, 5d },
            { 3d, 2d, 10d } };
    protected double[][] testDataMinus = { { -1d, -2d, -3d },
            { -2d, -5d, -3d }, { -1d, 0d, -8d } };
    protected double[] testDataRow1 = { 1d, 2d, 3d };
    protected double[] testDataCol3 = { 3d, 3d, 8d };
    protected double[][] testDataInv = { { -40d, 16d, 9d }, { 13d, -5d, -3d },
            { 5d, -2d, -1d } };
    protected double[] preMultTest = { 8, 12, 33 };
    protected double[][] testData2 = { { 1d, 2d, 3d }, { 2d, 5d, 3d } };
    protected double[][] testData2T = { { 1d, 2d }, { 2d, 5d }, { 3d, 3d } };
    protected double[][] testDataPlusInv = { { -39d, 18d, 12d },
            { 15d, 0d, 0d }, { 6d, -2d, 7d } };

    // lu decomposition tests
    protected double[][] luData = { { 2d, 3d, 3d }, { 0d, 5d, 7d }, { 6d, 9d, 8d } };
    protected double[][] luDataLUDecomposition = { { 6d, 9d, 8d },
            { 0d, 5d, 7d }, { 0.33333333333333, 0d, 0.33333333333333 } };

    // singular matrices
    protected double[][] singular = { { 2d, 3d }, { 2d, 3d } };
    protected double[][] bigSingular = { { 1d, 2d, 3d, 4d },
            { 2d, 5d, 3d, 4d }, { 7d, 3d, 256d, 1930d }, { 3d, 7d, 6d, 8d } }; // 4th

    // row
    // =
    // 1st
    // +
    // 2nd
    protected double[][] detData = { { 1d, 2d, 3d }, { 4d, 5d, 6d },
            { 7d, 8d, 10d } };
    protected double[][] detData2 = { { 1d, 3d }, { 2d, 4d } };

    // vectors
    protected double[] testVector = { 1, 2, 3 };
    protected double[] testVector2 = { 1, 2, 3, 4 };

    // submatrix accessor tests
    protected double[][] subTestData = { { 1, 2, 3, 4 },
            { 1.5, 2.5, 3.5, 4.5 }, { 2, 4, 6, 8 }, { 4, 5, 6, 7 } };

    // array selections
    protected double[][] subRows02Cols13 = { { 2, 4 }, { 4, 8 } };
    protected double[][] subRows03Cols12 = { { 2, 3 }, { 5, 6 } };
    protected double[][] subRows03Cols123 = { { 2, 3, 4 }, { 5, 6, 7 } };

    // effective permutations
    protected double[][] subRows20Cols123 = { { 4, 6, 8 }, { 2, 3, 4 } };
    protected double[][] subRows31Cols31 = { { 7, 5 }, { 4.5, 2.5 } };

    // contiguous ranges
    protected double[][] subRows01Cols23 = { { 3, 4 }, { 3.5, 4.5 } };
    protected double[][] subRows23Cols00 = { { 2 }, { 4 } };
    protected double[][] subRows00Cols33 = { { 4 } };

    // row matrices
    protected double[][] subRow0 = { { 1, 2, 3, 4 } };
    protected double[][] subRow3 = { { 4, 5, 6, 7 } };

    // column matrices
    protected double[][] subColumn1 = { { 2 }, { 2.5 }, { 4 }, { 5 } };
    protected double[][] subColumn3 = { { 4 }, { 4.5 }, { 8 }, { 7 } };

    // tolerances
    protected double entryTolerance = 10E-16;
    protected double normTolerance = 10E-14;

    /** test dimensions */
    @Test
    void testDimensions() {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        OpenMapRealMatrix m2 = createSparseMatrix(testData2);
        assertEquals(3, m.getRowDimension(), "testData row dimension");
        assertEquals(3, m.getColumnDimension(), "testData column dimension");
        assertTrue(m.isSquare(), "testData is square");
        assertEquals(2, m2.getRowDimension(), "testData2 row dimension");
        assertEquals(3, m2.getColumnDimension(), "testData2 column dimension");
        assertFalse(m2.isSquare(), "testData2 is not square");
    }

    /** test copy functions */
    @Test
    void testCopyFunctions() {
        OpenMapRealMatrix m1 = createSparseMatrix(testData);
        RealMatrix m2 = m1.copy();
        assertEquals(m1.getClass(), m2.getClass());
        assertEquals((m2), m1);
        OpenMapRealMatrix m3 = createSparseMatrix(testData);
        RealMatrix m4 = m3.copy();
        assertEquals(m3.getClass(), m4.getClass());
        assertEquals((m4), m3);
    }

    /** test add */
    @Test
    void testAdd() {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        OpenMapRealMatrix mInv = createSparseMatrix(testDataInv);
        OpenMapRealMatrix mDataPlusInv = createSparseMatrix(testDataPlusInv);
        RealMatrix mPlusMInv = m.add(mInv);
        for (int row = 0; row < m.getRowDimension(); row++) {
            for (int col = 0; col < m.getColumnDimension(); col++) {
                assertEquals(mDataPlusInv.getEntry(row, col), mPlusMInv.getEntry(row, col),
                    entryTolerance,
                    "sum entry entry");
            }
        }
    }

    /** test add failure */
    @Test
    void testAddFail() {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        OpenMapRealMatrix m2 = createSparseMatrix(testData2);
        try {
            m.add(m2);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test norm */
    @Test
    void testNorm() {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        OpenMapRealMatrix m2 = createSparseMatrix(testData2);
        assertEquals(14d, m.getNorm1(), entryTolerance, "testData norm");
        assertEquals(7d, m2.getNorm1(), entryTolerance, "testData2 norm");
        assertEquals(10d, m.getNormInfty(), entryTolerance, "testData norm");
        assertEquals(10d, m2.getNormInfty(), entryTolerance, "testData2 norm");
    }

    /** test m-n = m + -n */
    @Test
    void testPlusMinus() {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        OpenMapRealMatrix n = createSparseMatrix(testDataInv);
        customAssertClose("m-n = m + -n", m.subtract(n),
                          n.scalarMultiply(-1d).add(m), entryTolerance);
        try {
            m.subtract(createSparseMatrix(testData2));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test multiply */
    @Test
    void testMultiply() {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        OpenMapRealMatrix mInv = createSparseMatrix(testDataInv);
        OpenMapRealMatrix identity = createSparseMatrix(id);
        OpenMapRealMatrix m2 = createSparseMatrix(testData2);
        customAssertClose("inverse multiply", m.multiply(mInv), identity,
                          entryTolerance);
        customAssertClose("inverse multiply", m.multiply(new BlockRealMatrix(testDataInv)), identity,
                          entryTolerance);
        customAssertClose("inverse multiply", mInv.multiply(m), identity,
                          entryTolerance);
        customAssertClose("identity multiply", m.multiply(identity), m,
                          entryTolerance);
        customAssertClose("identity multiply", identity.multiply(mInv), mInv,
                          entryTolerance);
        customAssertClose("identity multiply", m2.multiply(identity), m2,
                          entryTolerance);
        try {
            m.multiply(createSparseMatrix(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    // Additional Test for Array2DRowRealMatrixTest.testMultiply

    private double[][] d3 = new double[][] { { 1, 2, 3, 4 }, { 5, 6, 7, 8 } };
    private double[][] d4 = new double[][] { { 1 }, { 2 }, { 3 }, { 4 } };
    private double[][] d5 = new double[][] { { 30 }, { 70 } };

    @Test
    void testMultiply2() {
        RealMatrix m3 = createSparseMatrix(d3);
        RealMatrix m4 = createSparseMatrix(d4);
        RealMatrix m5 = createSparseMatrix(d5);
        customAssertClose("m3*m4=m5", m3.multiply(m4), m5, entryTolerance);
    }

    @Test
    void testMultiplyTransposedSparseRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x5f31d5645cf821efl);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final OpenMapRealMatrix a = new OpenMapRealMatrix(rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final OpenMapRealMatrix b = new OpenMapRealMatrix(interm, cols);
                    b.walkInOptimizedOrder(randomSetter);
                    assertEquals(0.0,
                                        a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).getNorm1(),
                                        1.0e-15);
                }
            }
        }
    }

    @Test
    void testMultiplyTransposedWrongDimensions() {
        try {
            new OpenMapRealMatrix(2, 3).multiplyTransposed(new OpenMapRealMatrix(3, 2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testTransposeMultiplySparseRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x5f31d5645cf821efl);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final OpenMapRealMatrix a = new OpenMapRealMatrix(rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final OpenMapRealMatrix b = new OpenMapRealMatrix(rows, interm);
                    b.walkInOptimizedOrder(randomSetter);
                    assertEquals(0.0,
                                        a.transposeMultiply(b).subtract(a.transpose().multiply(b)).getNorm1(),
                                        4.0e-13);
                }
            }
        }
    }

    @Test
    void testTransposeMultiplyWrongDimensions() {
        try {
            new OpenMapRealMatrix(2, 3).transposeMultiply(new OpenMapRealMatrix(3, 2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(2, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(3, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    /** test trace */
    @Test
    void testTrace() {
        RealMatrix m = createSparseMatrix(id);
        assertEquals(3d, m.getTrace(), entryTolerance, "identity trace");
        m = createSparseMatrix(testData2);
        try {
            m.getTrace();
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test sclarAdd */
    @Test
    void testScalarAdd() {
        RealMatrix m = createSparseMatrix(testData);
        customAssertClose("scalar add", createSparseMatrix(testDataPlus2),
                          m.scalarAdd(2d), entryTolerance);
    }

    /** test operate */
    @Test
    void testOperate() {
        RealMatrix m = createSparseMatrix(id);
        customAssertClose("identity operate", testVector, m.operate(testVector),
                          entryTolerance);
        customAssertClose("identity operate", testVector, m.operate(
                new ArrayRealVector(testVector)).toArray(), entryTolerance);
        m = createSparseMatrix(bigSingular);
        try {
            m.operate(testVector);
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test issue MATH-209 */
    @Test
    void testMath209() {
        RealMatrix a = createSparseMatrix(new double[][] {
                { 1, 2 }, { 3, 4 }, { 5, 6 } });
        double[] b = a.operate(new double[] { 1, 1 });
        assertEquals(a.getRowDimension(), b.length);
        assertEquals(3.0, b[0], 1.0e-12);
        assertEquals(7.0, b[1], 1.0e-12);
        assertEquals(11.0, b[2], 1.0e-12);
    }

    /** test transpose */
    @Test
    void testTranspose() {
        RealMatrix m = createSparseMatrix(testData);
        RealMatrix mIT = new LUDecomposition(m).getSolver().getInverse().transpose();
        RealMatrix mTI = new LUDecomposition(m.transpose()).getSolver().getInverse();
        customAssertClose("inverse-transpose", mIT, mTI, normTolerance);
        m = createSparseMatrix(testData2);
        RealMatrix mt = createSparseMatrix(testData2T);
        customAssertClose("transpose", mt, m.transpose(), normTolerance);
    }

    /** test preMultiply by vector */
    @Test
    void testPremultiplyVector() {
        RealMatrix m = createSparseMatrix(testData);
        customAssertClose("premultiply", m.preMultiply(testVector), preMultTest,
                          normTolerance);
        customAssertClose("premultiply", m.preMultiply(
            new ArrayRealVector(testVector).toArray()), preMultTest, normTolerance);
        m = createSparseMatrix(bigSingular);
        try {
            m.preMultiply(testVector);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testPremultiply() {
        RealMatrix m3 = createSparseMatrix(d3);
        RealMatrix m4 = createSparseMatrix(d4);
        RealMatrix m5 = createSparseMatrix(d5);
        customAssertClose("m3*m4=m5", m4.preMultiply(m3), m5, entryTolerance);

        OpenMapRealMatrix m = createSparseMatrix(testData);
        OpenMapRealMatrix mInv = createSparseMatrix(testDataInv);
        OpenMapRealMatrix identity = createSparseMatrix(id);
        customAssertClose("inverse multiply", m.preMultiply(mInv), identity,
                          entryTolerance);
        customAssertClose("inverse multiply", mInv.preMultiply(m), identity,
                          entryTolerance);
        customAssertClose("identity multiply", m.preMultiply(identity), m,
                          entryTolerance);
        customAssertClose("identity multiply", identity.preMultiply(mInv), mInv,
                          entryTolerance);
        try {
            m.preMultiply(createSparseMatrix(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testGetVectors() {
        RealMatrix m = createSparseMatrix(testData);
        customAssertClose("get row", m.getRow(0), testDataRow1, entryTolerance);
        customAssertClose("get col", m.getColumn(2), testDataCol3, entryTolerance);
        try {
            m.getRow(10);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            m.getColumn(-1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testGetEntry() {
        RealMatrix m = createSparseMatrix(testData);
        assertEquals(2d, m.getEntry(0, 1), entryTolerance, "get entry");
        try {
            m.getEntry(10, 4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    /** test examples in user guide */
    @Test
    void testExamples() {
        // Create a real matrix with two rows and three columns
        double[][] matrixData = { { 1d, 2d, 3d }, { 2d, 5d, 3d } };
        RealMatrix m = createSparseMatrix(matrixData);
        // One more with three rows, two columns
        double[][] matrixData2 = { { 1d, 2d }, { 2d, 5d }, { 1d, 7d } };
        RealMatrix n = createSparseMatrix(matrixData2);
        // Now multiply m by n
        RealMatrix p = m.multiply(n);
        assertEquals(2, p.getRowDimension());
        assertEquals(2, p.getColumnDimension());
        // Invert p
        RealMatrix pInverse = new LUDecomposition(p).getSolver().getInverse();
        assertEquals(2, pInverse.getRowDimension());
        assertEquals(2, pInverse.getColumnDimension());

        // Solve example
        double[][] coefficientsData = { { 2, 3, -2 }, { -1, 7, 6 },
                { 4, -3, -5 } };
        RealMatrix coefficients = createSparseMatrix(coefficientsData);
        RealVector constants = new ArrayRealVector(new double[]{ 1, -2, 1 }, false);
        RealVector solution = new LUDecomposition(coefficients).getSolver().solve(constants);
        final double cst0 = constants.getEntry(0);
        final double cst1 = constants.getEntry(1);
        final double cst2 = constants.getEntry(2);
        final double sol0 = solution.getEntry(0);
        final double sol1 = solution.getEntry(1);
        final double sol2 = solution.getEntry(2);
        assertEquals(2 * sol0 + 3 * sol1 - 2 * sol2, cst0, 1E-12);
        assertEquals(-1 * sol0 + 7 * sol1 + 6 * sol2, cst1, 1E-12);
        assertEquals(4 * sol0 - 3 * sol1 - 5 * sol2, cst2, 1E-12);

    }

    // test submatrix accessors
    @Test
    void testSubMatrix() {
        RealMatrix m = createSparseMatrix(subTestData);
        RealMatrix mRows23Cols00 = createSparseMatrix(subRows23Cols00);
        RealMatrix mRows00Cols33 = createSparseMatrix(subRows00Cols33);
        RealMatrix mRows01Cols23 = createSparseMatrix(subRows01Cols23);
        RealMatrix mRows02Cols13 = createSparseMatrix(subRows02Cols13);
        RealMatrix mRows03Cols12 = createSparseMatrix(subRows03Cols12);
        RealMatrix mRows03Cols123 = createSparseMatrix(subRows03Cols123);
        RealMatrix mRows20Cols123 = createSparseMatrix(subRows20Cols123);
        RealMatrix mRows31Cols31 = createSparseMatrix(subRows31Cols31);
        assertEquals(mRows23Cols00, m.getSubMatrix(2, 3, 0, 0), "Rows23Cols00");
        assertEquals(mRows00Cols33, m.getSubMatrix(0, 0, 3, 3), "Rows00Cols33");
        assertEquals(mRows01Cols23, m.getSubMatrix(0, 1, 2, 3), "Rows01Cols23");
        assertEquals(mRows02Cols13,
            m.getSubMatrix(new int[] { 0, 2 }, new int[] { 1, 3 }),
            "Rows02Cols13");
        assertEquals(mRows03Cols12,
            m.getSubMatrix(new int[] { 0, 3 }, new int[] { 1, 2 }),
            "Rows03Cols12");
        assertEquals(mRows03Cols123,
            m.getSubMatrix(new int[] { 0, 3 }, new int[] { 1, 2, 3 }),
            "Rows03Cols123");
        assertEquals(mRows20Cols123,
            m.getSubMatrix(new int[] { 2, 0 }, new int[] { 1, 2, 3 }),
            "Rows20Cols123");
        assertEquals(mRows31Cols31,
            m.getSubMatrix(new int[] { 3, 1 }, new int[] { 3, 1 }),
            "Rows31Cols31");
        assertEquals(mRows31Cols31,
            m.getSubMatrix(new int[] { 3, 1 }, new int[] { 3, 1 }),
            "Rows31Cols31");

        try {
            m.getSubMatrix(1, 0, 2, 4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getSubMatrix(-1, 1, 2, 2);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getSubMatrix(1, 0, 2, 2);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getSubMatrix(1, 0, 2, 4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getSubMatrix(new int[] {}, new int[] { 0 });
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getSubMatrix(new int[] { 0 }, new int[] { 4 });
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetRowMatrix() {
        RealMatrix m = createSparseMatrix(subTestData);
        RealMatrix mRow0 = createSparseMatrix(subRow0);
        RealMatrix mRow3 = createSparseMatrix(subRow3);
        assertEquals(mRow0, m.getRowMatrix(0), "Row0");
        assertEquals(mRow3, m.getRowMatrix(3), "Row3");
        try {
            m.getRowMatrix(-1);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getRowMatrix(4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetColumnMatrix() {
        RealMatrix m = createSparseMatrix(subTestData);
        RealMatrix mColumn1 = createSparseMatrix(subColumn1);
        RealMatrix mColumn3 = createSparseMatrix(subColumn3);
        assertEquals(mColumn1, m.getColumnMatrix(1), "Column1");
        assertEquals(mColumn3, m.getColumnMatrix(3), "Column3");
        try {
            m.getColumnMatrix(-1);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getColumnMatrix(4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetRowVector() {
        RealMatrix m = createSparseMatrix(subTestData);
        RealVector mRow0 = new ArrayRealVector(subRow0[0]);
        RealVector mRow3 = new ArrayRealVector(subRow3[0]);
        assertEquals(mRow0, m.getRowVector(0), "Row0");
        assertEquals(mRow3, m.getRowVector(3), "Row3");
        try {
            m.getRowVector(-1);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getRowVector(4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetColumnVector() {
        RealMatrix m = createSparseMatrix(subTestData);
        RealVector mColumn1 = columnToVector(subColumn1);
        RealVector mColumn3 = columnToVector(subColumn3);
        assertEquals(mColumn1, m.getColumnVector(1), "Column1");
        assertEquals(mColumn3, m.getColumnVector(3), "Column3");
        try {
            m.getColumnVector(-1);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getColumnVector(4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    private RealVector columnToVector(double[][] column) {
        double[] data = new double[column.length];
        for (int i = 0; i < data.length; ++i) {
            data[i] = column[i][0];
        }
        return new ArrayRealVector(data, false);
    }

    @Test
    void testEqualsAndHashCode() {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        OpenMapRealMatrix m1 = m.copy();
        OpenMapRealMatrix mt = (OpenMapRealMatrix) m.transpose();
        assertTrue(m.hashCode() != mt.hashCode());
        assertEquals(m.hashCode(), m1.hashCode());
        assertEquals(m, m);
        assertEquals(m, m1);
        assertNotEquals(null, m);
        assertNotEquals(m, mt);
        assertNotEquals(m, createSparseMatrix(bigSingular));
    }

    @Test
    void testToString() {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        assertEquals("OpenMapRealMatrix{{1.0,2.0,3.0},{2.0,5.0,3.0},{1.0,0.0,8.0}}",
            m.toString());
        m = new OpenMapRealMatrix(1, 1);
        assertEquals("OpenMapRealMatrix{{0.0}}", m.toString());
    }

    @Test
    void testSetSubMatrix() {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        m.setSubMatrix(detData2, 1, 1);
        RealMatrix expected = createSparseMatrix(new double[][] {
                { 1.0, 2.0, 3.0 }, { 2.0, 1.0, 3.0 }, { 1.0, 2.0, 4.0 } });
        assertEquals(expected, m);

        m.setSubMatrix(detData2, 0, 0);
        expected = createSparseMatrix(new double[][] {
                { 1.0, 3.0, 3.0 }, { 2.0, 4.0, 3.0 }, { 1.0, 2.0, 4.0 } });
        assertEquals(expected, m);

        m.setSubMatrix(testDataPlus2, 0, 0);
        expected = createSparseMatrix(new double[][] {
                { 3.0, 4.0, 5.0 }, { 4.0, 7.0, 5.0 }, { 3.0, 2.0, 10.0 } });
        assertEquals(expected, m);

        // javadoc example
        OpenMapRealMatrix matrix =
            createSparseMatrix(new double[][] {
        { 1, 2, 3, 4 }, { 5, 6, 7, 8 }, { 9, 0, 1, 2 } });
        matrix.setSubMatrix(new double[][] { { 3, 4 }, { 5, 6 } }, 1, 1);
        expected = createSparseMatrix(new double[][] {
                { 1, 2, 3, 4 }, { 5, 3, 4, 8 }, { 9, 5, 6, 2 } });
        assertEquals(expected, matrix);

        // dimension overflow
        try {
            m.setSubMatrix(testData, 1, 1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
        // dimension underflow
        try {
            m.setSubMatrix(testData, -1, 1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
        try {
            m.setSubMatrix(testData, 1, -1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        // null
        try {
            m.setSubMatrix(null, 1, 1);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException e) {
            // expected
        }
        try {
            new OpenMapRealMatrix(0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        // ragged
        try {
            m.setSubMatrix(new double[][] { { 1 }, { 2, 3 } }, 0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        // empty
        try {
            m.setSubMatrix(new double[][] { {} }, 0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

    }

    @Test
    void testSerial()  {
        OpenMapRealMatrix m = createSparseMatrix(testData);
        assertEquals(m,UnitTestUtils.serializeAndRecover(m));
    }

    // --------------- -----------------Protected methods

    /** verifies that two matrices are close (1-norm) */
    protected void customAssertClose(String msg, RealMatrix m, RealMatrix n,
                                     double tolerance) {
        assertTrue(m.subtract(n).getNorm1() < tolerance, msg);
    }

    /** verifies that two vectors are close (sup norm) */
    protected void customAssertClose(String msg, double[] m, double[] n,
                                     double tolerance) {
        if (m.length != n.length) {
            fail("vectors not same length");
        }
        for (int i = 0; i < m.length; i++) {
            assertEquals(m[i], n[i],
                    tolerance,
                    msg + " " + i + " elements differ");
        }
    }

    private OpenMapRealMatrix createSparseMatrix(double[][] data) {
        OpenMapRealMatrix matrix = new OpenMapRealMatrix(data.length, data[0].length);
        for (int row = 0; row < data.length; row++) {
            for (int col = 0; col < data[row].length; col++) {
                matrix.setEntry(row, col, data[row][col]);
            }
        }
        return matrix;
    }
}
