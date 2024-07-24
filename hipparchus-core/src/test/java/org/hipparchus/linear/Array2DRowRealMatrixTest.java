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
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import java.text.DecimalFormat;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link Array2DRowRealMatrix} class.
 *
 */

final class Array2DRowRealMatrixTest {

    // 3 x 3 identity matrix
    protected double[][] id = { {1d,0d,0d}, {0d,1d,0d}, {0d,0d,1d} };

    // Test data for group operations
    protected double[][] testData = { {1d,2d,3d}, {2d,5d,3d}, {1d,0d,8d} };
    protected double[][] testDataLU = {{2d, 5d, 3d}, {.5d, -2.5d, 6.5d}, {0.5d, 0.2d, .2d}};
    protected double[][] testDataPlus2 = { {3d,4d,5d}, {4d,7d,5d}, {3d,2d,10d} };
    protected double[][] testDataMinus = { {-1d,-2d,-3d}, {-2d,-5d,-3d},
       {-1d,0d,-8d} };
    protected double[] testDataRow1 = {1d,2d,3d};
    protected double[] testDataCol3 = {3d,3d,8d};
    protected double[][] testDataInv =
        { {-40d,16d,9d}, {13d,-5d,-3d}, {5d,-2d,-1d} };
    protected double[] preMultTest = {8,12,33};
    protected double[][] testData2 ={ {1d,2d,3d}, {2d,5d,3d}};
    protected double[][] testData2T = { {1d,2d}, {2d,5d}, {3d,3d}};
    protected double[][] testDataPlusInv =
        { {-39d,18d,12d}, {15d,0d,0d}, {6d,-2d,7d} };

    // lu decomposition tests
    protected double[][] luData = { {2d,3d,3d}, {0d,5d,7d}, {6d,9d,8d} };
    protected double[][] luDataLUDecomposition = { {6d,9d,8d}, {0d,5d,7d},
            {0.33333333333333,0d,0.33333333333333} };

    // singular matrices
    protected double[][] singular = { {2d,3d}, {2d,3d} };
    protected double[][] bigSingular = {{1d,2d,3d,4d}, {2d,5d,3d,4d},
        {7d,3d,256d,1930d}, {3d,7d,6d,8d}}; // 4th row = 1st + 2nd
    protected double[][] detData = { {1d,2d,3d}, {4d,5d,6d}, {7d,8d,10d} };
    protected double[][] detData2 = { {1d, 3d}, {2d, 4d}};

    // vectors
    protected double[] testVector = {1,2,3};
    protected double[] testVector2 = {1,2,3,4};

    // submatrix accessor tests
    protected double[][] subTestData = {{1, 2, 3, 4}, {1.5, 2.5, 3.5, 4.5},
            {2, 4, 6, 8}, {4, 5, 6, 7}};
    // array selections
    protected double[][] subRows02Cols13 = { {2, 4}, {4, 8}};
    protected double[][] subRows03Cols12 = { {2, 3}, {5, 6}};
    protected double[][] subRows03Cols123 = { {2, 3, 4} , {5, 6, 7}};
    // effective permutations
    protected double[][] subRows20Cols123 = { {4, 6, 8} , {2, 3, 4}};
    protected double[][] subRows31Cols31 = {{7, 5}, {4.5, 2.5}};
    // contiguous ranges
    protected double[][] subRows01Cols23 = {{3,4} , {3.5, 4.5}};
    protected double[][] subRows23Cols00 = {{2} , {4}};
    protected double[][] subRows00Cols33 = {{4}};
    // row matrices
    protected double[][] subRow0 = {{1,2,3,4}};
    protected double[][] subRow3 = {{4,5,6,7}};
    // column matrices
    protected double[][] subColumn1 = {{2}, {2.5}, {4}, {5}};
    protected double[][] subColumn3 = {{4}, {4.5}, {8}, {7}};

    // tolerances
    protected double entryTolerance = 10E-16;
    protected double normTolerance = 10E-14;
    protected double powerTolerance = 10E-16;

    /** test dimensions */
    @Test
    void testDimensions() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testData2);
        assertEquals(3,m.getRowDimension(),"testData row dimension");
        assertEquals(3,m.getColumnDimension(),"testData column dimension");
        assertTrue(m.isSquare(),"testData is square");
        assertEquals(2, m2.getRowDimension(), "testData2 row dimension");
        assertEquals(3, m2.getColumnDimension(), "testData2 column dimension");
        assertFalse(m2.isSquare(), "testData2 is not square");
    }

    /** test copy functions */
    @Test
    void testCopyFunctions() {
        Array2DRowRealMatrix m1 = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(m1.getData());
        assertEquals(m2,m1);
        Array2DRowRealMatrix m3 = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m4 = new Array2DRowRealMatrix(m3.getData(), false);
        assertEquals(m4,m3);
    }

    /** test add */
    @Test
    void testAdd() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix mInv = new Array2DRowRealMatrix(testDataInv);
        RealMatrix mPlusMInv = m.add(mInv);
        double[][] sumEntries = mPlusMInv.getData();
        for (int row = 0; row < m.getRowDimension(); row++) {
            for (int col = 0; col < m.getColumnDimension(); col++) {
                assertEquals(testDataPlusInv[row][col],sumEntries[row][col],
                        entryTolerance,
                        "sum entry entry");
            }
        }
    }

    /** test add failure */
    @Test
    void testAddFail() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testData2);
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
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testData2);
        assertEquals(10d,m.getNormInfty(),entryTolerance,"testData norm");
        assertEquals(10d,m2.getNormInfty(),entryTolerance,"testData2 norm");
        assertEquals(14d,m.getNorm1(),entryTolerance,"testData norm");
        assertEquals(7d,m2.getNorm1(),entryTolerance,"testData2 norm");
    }

    /** test Frobenius norm */
    @Test
    void testFrobeniusNorm() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testData2);
        assertEquals(FastMath.sqrt(117.0), m.getFrobeniusNorm(), entryTolerance, "testData Frobenius norm");
        assertEquals(FastMath.sqrt(52.0), m2.getFrobeniusNorm(), entryTolerance, "testData2 Frobenius norm");
    }

    /** test m-n = m + -n */
    @Test
    void testPlusMinus() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testDataInv);
        UnitTestUtils.customAssertEquals("m-n = m + -n", m.subtract(m2),
                                         m2.scalarMultiply(-1d).add(m), entryTolerance);
        try {
            m.subtract(new Array2DRowRealMatrix(testData2));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test multiply */
    @Test
    void testMultiply() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix mInv = new Array2DRowRealMatrix(testDataInv);
        Array2DRowRealMatrix identity = new Array2DRowRealMatrix(id);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testData2);
        UnitTestUtils.customAssertEquals("inverse multiply", m.multiply(mInv),
                                         identity, entryTolerance);
        UnitTestUtils.customAssertEquals("inverse multiply", mInv.multiply(m),
                                         identity, entryTolerance);
        UnitTestUtils.customAssertEquals("identity multiply", m.multiply(identity),
                                         m, entryTolerance);
        UnitTestUtils.customAssertEquals("identity multiply", identity.multiply(mInv),
                                         mInv, entryTolerance);
        UnitTestUtils.customAssertEquals("identity multiply", m2.multiply(identity),
                                         m2, entryTolerance);
        try {
            m.multiply(new Array2DRowRealMatrix(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    //Additional Test for Array2DRowRealMatrixTest.testMultiply

    private final double[][] d3 = new double[][] {{1,2,3,4},{5,6,7,8}};
    private final double[][] d4 = new double[][] {{1},{2},{3},{4}};
    private final double[][] d5 = new double[][] {{30},{70}};

    @Test
    void testMultiply2() {
       RealMatrix m3 = new Array2DRowRealMatrix(d3);
       RealMatrix m4 = new Array2DRowRealMatrix(d4);
       RealMatrix m5 = new Array2DRowRealMatrix(d5);
       UnitTestUtils.customAssertEquals("m3*m4=m5", m3.multiply(m4), m5, entryTolerance);
    }

    @Test
    void testMultiplyTransposedArray2DRowRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xdeff3d383a112763l);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final Array2DRowRealMatrix a = new Array2DRowRealMatrix(rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final Array2DRowRealMatrix b = new Array2DRowRealMatrix(interm, cols);
                    b.walkInOptimizedOrder(randomSetter);
                    assertEquals(0.0,
                                        a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).getNorm1(),
                                        1.0e-15);
                }
            }
        }
    }

    @Test
    void testMultiplyTransposedBlockRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x463e54fb50b900fel);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final Array2DRowRealMatrix a = new Array2DRowRealMatrix(rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final BlockRealMatrix b = new BlockRealMatrix(interm, cols);
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
            new Array2DRowRealMatrix(2, 3).multiplyTransposed(new Array2DRowRealMatrix(3, 2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testTransposeMultiplyArray2DRowRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xdeff3d383a112763l);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final Array2DRowRealMatrix a = new Array2DRowRealMatrix(rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final Array2DRowRealMatrix b = new Array2DRowRealMatrix(rows, interm);
                    b.walkInOptimizedOrder(randomSetter);
                    assertEquals(0.0,
                                        a.transposeMultiply(b).subtract(a.transpose().multiply(b)).getNorm1(),
                                        1.0e-15);
                }
            }
        }
    }

    @Test
    void testTransposeMultiplyBlockRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x463e54fb50b900fel);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final Array2DRowRealMatrix a = new Array2DRowRealMatrix(rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final BlockRealMatrix b = new BlockRealMatrix(rows, interm);
                    b.walkInOptimizedOrder(randomSetter);
                    assertEquals(0.0,
                                        a.transposeMultiply(b).subtract(a.transpose().multiply(b)).getNorm1(),
                                        1.0e-15);
                }
            }
        }
    }

    @Test
    void testTransposeMultiplyWrongDimensions() {
        try {
            new Array2DRowRealMatrix(2, 3).transposeMultiply(new Array2DRowRealMatrix(3, 2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(2, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(3, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testPower() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix mInv = new Array2DRowRealMatrix(testDataInv);
        Array2DRowRealMatrix mPlusInv = new Array2DRowRealMatrix(testDataPlusInv);
        Array2DRowRealMatrix identity = new Array2DRowRealMatrix(id);

        UnitTestUtils.customAssertEquals("m^0", m.power(0),
                                         identity, entryTolerance);
        UnitTestUtils.customAssertEquals("mInv^0", mInv.power(0),
                                         identity, entryTolerance);
        UnitTestUtils.customAssertEquals("mPlusInv^0", mPlusInv.power(0),
                                         identity, entryTolerance);

        UnitTestUtils.customAssertEquals("m^1", m.power(1),
                                         m, entryTolerance);
        UnitTestUtils.customAssertEquals("mInv^1", mInv.power(1),
                                         mInv, entryTolerance);
        UnitTestUtils.customAssertEquals("mPlusInv^1", mPlusInv.power(1),
                                         mPlusInv, entryTolerance);

        RealMatrix C1 = m.copy();
        RealMatrix C2 = mInv.copy();
        RealMatrix C3 = mPlusInv.copy();

        for (int i = 2; i <= 10; ++i) {
            C1 = C1.multiply(m);
            C2 = C2.multiply(mInv);
            C3 = C3.multiply(mPlusInv);

            UnitTestUtils.customAssertEquals("m^" + i, m.power(i),
                                             C1, entryTolerance);
            UnitTestUtils.customAssertEquals("mInv^" + i, mInv.power(i),
                                             C2, entryTolerance);
            UnitTestUtils.customAssertEquals("mPlusInv^" + i, mPlusInv.power(i),
                                             C3, entryTolerance);
        }

        try {
            Array2DRowRealMatrix mNotSquare = new Array2DRowRealMatrix(testData2T);
            mNotSquare.power(2);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }

        try {
            m.power(-1);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test trace */
    @Test
    void testTrace() {
        RealMatrix m = new Array2DRowRealMatrix(id);
        assertEquals(3d,m.getTrace(),entryTolerance,"identity trace");
        m = new Array2DRowRealMatrix(testData2);
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
        RealMatrix m = new Array2DRowRealMatrix(testData);
        UnitTestUtils.customAssertEquals("scalar add", new Array2DRowRealMatrix(testDataPlus2),
                                         m.scalarAdd(2d), entryTolerance);
    }

    /** test operate */
    @Test
    void testOperate() {
        RealMatrix m = new Array2DRowRealMatrix(id);
        UnitTestUtils.customAssertEquals("identity operate", testVector,
                                         m.operate(testVector), entryTolerance);
        UnitTestUtils.customAssertEquals("identity operate", testVector,
                                         m.operate(new ArrayRealVector(testVector)).toArray(), entryTolerance);
        m = new Array2DRowRealMatrix(bigSingular);
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
        RealMatrix a = new Array2DRowRealMatrix(new double[][] {
                { 1, 2 }, { 3, 4 }, { 5, 6 }
        }, false);
        double[] b = a.operate(new double[] { 1, 1 });
        assertEquals(a.getRowDimension(), b.length);
        assertEquals( 3.0, b[0], 1.0e-12);
        assertEquals( 7.0, b[1], 1.0e-12);
        assertEquals(11.0, b[2], 1.0e-12);
    }

    /** test transpose */
    @Test
    void testTranspose() {
        RealMatrix m = new Array2DRowRealMatrix(testData);
        RealMatrix mIT = new LUDecomposition(m).getSolver().getInverse().transpose();
        RealMatrix mTI = new LUDecomposition(m.transpose()).getSolver().getInverse();
        UnitTestUtils.customAssertEquals("inverse-transpose", mIT, mTI, normTolerance);
        m = new Array2DRowRealMatrix(testData2);
        RealMatrix mt = new Array2DRowRealMatrix(testData2T);
        UnitTestUtils.customAssertEquals("transpose", mt, m.transpose(), normTolerance);
    }

    /** test preMultiply by vector */
    @Test
    void testPremultiplyVector() {
        RealMatrix m = new Array2DRowRealMatrix(testData);
        UnitTestUtils.customAssertEquals("premultiply", m.preMultiply(testVector),
                                         preMultTest, normTolerance);
        UnitTestUtils.customAssertEquals("premultiply", m.preMultiply(new ArrayRealVector(testVector).toArray()),
                                         preMultTest, normTolerance);
        m = new Array2DRowRealMatrix(bigSingular);
        try {
            m.preMultiply(testVector);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testPremultiply() {
        RealMatrix m3 = new Array2DRowRealMatrix(d3);
        RealMatrix m4 = new Array2DRowRealMatrix(d4);
        RealMatrix m5 = new Array2DRowRealMatrix(d5);
        UnitTestUtils.customAssertEquals("m3*m4=m5", m4.preMultiply(m3), m5, entryTolerance);

        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix mInv = new Array2DRowRealMatrix(testDataInv);
        Array2DRowRealMatrix identity = new Array2DRowRealMatrix(id);
        UnitTestUtils.customAssertEquals("inverse multiply", m.preMultiply(mInv),
                                         identity, entryTolerance);
        UnitTestUtils.customAssertEquals("inverse multiply", mInv.preMultiply(m),
                                         identity, entryTolerance);
        UnitTestUtils.customAssertEquals("identity multiply", m.preMultiply(identity),
                                         m, entryTolerance);
        UnitTestUtils.customAssertEquals("identity multiply", identity.preMultiply(mInv),
                                         mInv, entryTolerance);
        try {
            m.preMultiply(new Array2DRowRealMatrix(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testGetVectors() {
        RealMatrix m = new Array2DRowRealMatrix(testData);
        UnitTestUtils.customAssertEquals("get row", m.getRow(0), testDataRow1, entryTolerance);
        UnitTestUtils.customAssertEquals("get col", m.getColumn(2), testDataCol3, entryTolerance);
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
        RealMatrix m = new Array2DRowRealMatrix(testData);
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
        double[][] matrixData = { {1d,2d,3d}, {2d,5d,3d}};
        RealMatrix m = new Array2DRowRealMatrix(matrixData);
        // One more with three rows, two columns
        double[][] matrixData2 = { {1d,2d}, {2d,5d}, {1d, 7d}};
        RealMatrix n = new Array2DRowRealMatrix(matrixData2);
        // Now multiply m by n
        RealMatrix p = m.multiply(n);
        assertEquals(2, p.getRowDimension());
        assertEquals(2, p.getColumnDimension());
        // Invert p
        RealMatrix pInverse = new LUDecomposition(p).getSolver().getInverse();
        assertEquals(2, pInverse.getRowDimension());
        assertEquals(2, pInverse.getColumnDimension());

        // Solve example
        double[][] coefficientsData = {{2, 3, -2}, {-1, 7, 6}, {4, -3, -5}};
        RealMatrix coefficients = new Array2DRowRealMatrix(coefficientsData);
        RealVector constants = new ArrayRealVector(new double[]{1, -2, 1}, false);
        RealVector solution = new LUDecomposition(coefficients).getSolver().solve(constants);
        final double cst0 = constants.getEntry(0);
        final double cst1 = constants.getEntry(1);
        final double cst2 = constants.getEntry(2);
        final double sol0 = solution.getEntry(0);
        final double sol1 = solution.getEntry(1);
        final double sol2 = solution.getEntry(2);
        assertEquals(2 * sol0 + 3 * sol1 -2 * sol2, cst0, 1E-12);
        assertEquals(-1 * sol0 + 7 * sol1 + 6 * sol2, cst1, 1E-12);
        assertEquals(4 * sol0 - 3 * sol1 -5 * sol2, cst2, 1E-12);
    }

    // test submatrix accessors
    @Test
    void testGetSubMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        checkGetSubMatrix(m, subRows23Cols00,  2 , 3 , 0, 0, false);
        checkGetSubMatrix(m, subRows00Cols33,  0 , 0 , 3, 3, false);
        checkGetSubMatrix(m, subRows01Cols23,  0 , 1 , 2, 3, false);
        checkGetSubMatrix(m, subRows02Cols13,  new int[] { 0, 2 }, new int[] { 1, 3 },    false);
        checkGetSubMatrix(m, subRows03Cols12,  new int[] { 0, 3 }, new int[] { 1, 2 },    false);
        checkGetSubMatrix(m, subRows03Cols123, new int[] { 0, 3 }, new int[] { 1, 2, 3 }, false);
        checkGetSubMatrix(m, subRows20Cols123, new int[] { 2, 0 }, new int[] { 1, 2, 3 }, false);
        checkGetSubMatrix(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 },    false);
        checkGetSubMatrix(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 },    false);
        checkGetSubMatrix(m, null,  1, 0, 2, 4, true);
        checkGetSubMatrix(m, null, -1, 1, 2, 2, true);
        checkGetSubMatrix(m, null,  1, 0, 2, 2, true);
        checkGetSubMatrix(m, null,  1, 0, 2, 4, true);
        checkGetSubMatrix(m, null, new int[] {},    new int[] { 0 }, true);
        checkGetSubMatrix(m, null, new int[] { 0 }, new int[] { 4 }, true);
    }

    private void checkGetSubMatrix(RealMatrix m, double[][] reference,
                                   int startRow, int endRow, int startColumn, int endColumn,
                                   boolean mustFail) {
        try {
            RealMatrix sub = m.getSubMatrix(startRow, endRow, startColumn, endColumn);
            assertEquals(new Array2DRowRealMatrix(reference), sub);
            if (mustFail) {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentException or MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (!mustFail) {
                throw e;
            }
        }
    }

    private void checkGetSubMatrix(RealMatrix m, double[][] reference,
                                   int[] selectedRows, int[] selectedColumns,
                                   boolean mustFail) {
        try {
            RealMatrix sub = m.getSubMatrix(selectedRows, selectedColumns);
            assertEquals(new Array2DRowRealMatrix(reference), sub);
            if (mustFail) {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentException or MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (!mustFail) {
                throw e;
            }
        }
    }

    @Test
    void testCopySubMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        checkCopy(m, subRows23Cols00,  2 , 3 , 0, 0, false);
        checkCopy(m, subRows00Cols33,  0 , 0 , 3, 3, false);
        checkCopy(m, subRows01Cols23,  0 , 1 , 2, 3, false);
        checkCopy(m, subRows02Cols13,  new int[] { 0, 2 }, new int[] { 1, 3 },    false);
        checkCopy(m, subRows03Cols12,  new int[] { 0, 3 }, new int[] { 1, 2 },    false);
        checkCopy(m, subRows03Cols123, new int[] { 0, 3 }, new int[] { 1, 2, 3 }, false);
        checkCopy(m, subRows20Cols123, new int[] { 2, 0 }, new int[] { 1, 2, 3 }, false);
        checkCopy(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 },    false);
        checkCopy(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 },    false);

        checkCopy(m, null,  1, 0, 2, 4, true);
        checkCopy(m, null, -1, 1, 2, 2, true);
        checkCopy(m, null,  1, 0, 2, 2, true);
        checkCopy(m, null,  1, 0, 2, 4, true);
        checkCopy(m, null, new int[] {},    new int[] { 0 }, true);
        checkCopy(m, null, new int[] { 0 }, new int[] { 4 }, true);

        // rectangular check
        double[][] copy = new double[][] { { 0, 0, 0 }, { 0, 0 } };
        checkCopy(m, copy, 0, 1, 0, 2, true);
        checkCopy(m, copy, new int[] { 0, 1 }, new int[] { 0, 1, 2 }, true);
    }

    private void checkCopy(RealMatrix m, double[][] reference,
                           int startRow, int endRow, int startColumn, int endColumn,
                           boolean mustFail) {
        try {
            double[][] sub = (reference == null) ?
                             new double[1][1] : createIdenticalCopy(reference);
            m.copySubMatrix(startRow, endRow, startColumn, endColumn, sub);
            assertEquals(new Array2DRowRealMatrix(reference), new Array2DRowRealMatrix(sub));
            if (mustFail) {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentException or MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (!mustFail) {
                throw e;
            }
        }
    }

    private void checkCopy(RealMatrix m, double[][] reference,
                           int[] selectedRows, int[] selectedColumns,
                           boolean mustFail) {
        try {
            double[][] sub = (reference == null) ?
                    new double[1][1] : createIdenticalCopy(reference);
            m.copySubMatrix(selectedRows, selectedColumns, sub);
            assertEquals(new Array2DRowRealMatrix(reference), new Array2DRowRealMatrix(sub));
            if (mustFail) {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentException or MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (!mustFail) {
                throw e;
            }
        }
    }

    private double[][] createIdenticalCopy(final double[][] matrix) {
        final double[][] matrixCopy = new double[matrix.length][];
        for (int i = 0; i < matrixCopy.length; i++) {
            matrixCopy[i] = new double[matrix[i].length];
        }
        return matrixCopy;
    }

    @Test
    void testGetRowMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealMatrix mRow0 = new Array2DRowRealMatrix(subRow0);
        RealMatrix mRow3 = new Array2DRowRealMatrix(subRow3);
        assertEquals(mRow0,
                m.getRowMatrix(0),
                "Row0");
        assertEquals(mRow3,
                m.getRowMatrix(3),
                "Row3");
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
    void testSetRowMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealMatrix mRow3 = new Array2DRowRealMatrix(subRow3);
        assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowMatrix(0, mRow3);
        assertEquals(mRow3, m.getRowMatrix(0));
        try {
            m.setRowMatrix(-1, mRow3);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.setRowMatrix(0, m);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetColumnMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealMatrix mColumn1 = new Array2DRowRealMatrix(subColumn1);
        RealMatrix mColumn3 = new Array2DRowRealMatrix(subColumn3);
        assertEquals(mColumn1,
                m.getColumnMatrix(1),
                "Column1");
        assertEquals(mColumn3,
                m.getColumnMatrix(3),
                "Column3");
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
    void testSetColumnMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealMatrix mColumn3 = new Array2DRowRealMatrix(subColumn3);
        assertNotSame(mColumn3, m.getColumnMatrix(1));
        m.setColumnMatrix(1, mColumn3);
        assertEquals(mColumn3, m.getColumnMatrix(1));
        try {
            m.setColumnMatrix(-1, mColumn3);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.setColumnMatrix(0, m);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetRowVector() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
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
    void testSetRowVector() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealVector mRow3 = new ArrayRealVector(subRow3[0]);
        assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowVector(0, mRow3);
        assertEquals(mRow3, m.getRowVector(0));
        try {
            m.setRowVector(-1, mRow3);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.setRowVector(0, new ArrayRealVector(5));
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetColumnVector() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
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

    @Test
    void testSetColumnVector() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealVector mColumn3 = columnToVector(subColumn3);
        assertNotSame(mColumn3, m.getColumnVector(1));
        m.setColumnVector(1, mColumn3);
        assertEquals(mColumn3, m.getColumnVector(1));
        try {
            m.setColumnVector(-1, mColumn3);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.setColumnVector(0, new ArrayRealVector(5));
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
    void testGetRow() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        checkArrays(subRow0[0], m.getRow(0));
        checkArrays(subRow3[0], m.getRow(3));
        try {
            m.getRow(-1);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getRow(4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testSetRow() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        assertTrue(subRow3[0][0] != m.getRow(0)[0]);
        m.setRow(0, subRow3[0]);
        checkArrays(subRow3[0], m.getRow(0));
        try {
            m.setRow(-1, subRow3[0]);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.setRow(0, new double[5]);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetColumn() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        double[] mColumn1 = columnToArray(subColumn1);
        double[] mColumn3 = columnToArray(subColumn3);
        checkArrays(mColumn1, m.getColumn(1));
        checkArrays(mColumn3, m.getColumn(3));
        try {
            m.getColumn(-1);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.getColumn(4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testSetColumn() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        double[] mColumn3 = columnToArray(subColumn3);
        assertTrue(mColumn3[0] != m.getColumn(1)[0]);
        m.setColumn(1, mColumn3);
        checkArrays(mColumn3, m.getColumn(1));
        try {
            m.setColumn(-1, mColumn3);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            m.setColumn(0, new double[5]);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    private double[] columnToArray(double[][] column) {
        double[] data = new double[column.length];
        for (int i = 0; i < data.length; ++i) {
            data[i] = column[i][0];
        }
        return data;
    }

    private void checkArrays(double[] expected, double[] actual) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(expected[i], actual[i], 0);
        }
    }

    @Test
    void testEqualsAndHashCode() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m1 = (Array2DRowRealMatrix) m.copy();
        Array2DRowRealMatrix mt = (Array2DRowRealMatrix) m.transpose();
        assertTrue(m.hashCode() != mt.hashCode());
        assertEquals(m.hashCode(), m1.hashCode());
        assertEquals(m, m);
        assertEquals(m, m1);
        assertNotEquals(null, m);
        assertNotEquals(m, mt);
        assertNotEquals(m, new Array2DRowRealMatrix(bigSingular));
    }

    @Test
    void testToString() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        assertEquals("Array2DRowRealMatrix{{1.0,2.0,3.0},{2.0,5.0,3.0},{1.0,0.0,8.0}}",
                m.toString());
        m = new Array2DRowRealMatrix();
        assertEquals("Array2DRowRealMatrix{}",
                m.toString());
    }

    @Test
    void testSetSubMatrix() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        m.setSubMatrix(detData2,1,1);
        RealMatrix expected = MatrixUtils.createRealMatrix
            (new double[][] {{1.0,2.0,3.0},{2.0,1.0,3.0},{1.0,2.0,4.0}});
        assertEquals(expected, m);

        m.setSubMatrix(detData2,0,0);
        expected = MatrixUtils.createRealMatrix
            (new double[][] {{1.0,3.0,3.0},{2.0,4.0,3.0},{1.0,2.0,4.0}});
        assertEquals(expected, m);

        m.setSubMatrix(testDataPlus2,0,0);
        expected = MatrixUtils.createRealMatrix
            (new double[][] {{3.0,4.0,5.0},{4.0,7.0,5.0},{3.0,2.0,10.0}});
        assertEquals(expected, m);

        // dimension overflow
        try {
            m.setSubMatrix(testData,1,1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
        // dimension underflow
        try {
            m.setSubMatrix(testData,-1,1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
        try {
            m.setSubMatrix(testData,1,-1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        // null
        try {
            m.setSubMatrix(null,1,1);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException e) {
            // expected
        }
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix();
        try {
            m2.setSubMatrix(testData,0,1);
            fail("expecting MathIllegalStateException");
        } catch (MathIllegalStateException e) {
            // expected
        }
        try {
            m2.setSubMatrix(testData,1,0);
            fail("expecting MathIllegalStateException");
        } catch (MathIllegalStateException e) {
            // expected
        }

        // ragged
        try {
            m.setSubMatrix(new double[][] {{1}, {2, 3}}, 0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        // empty
        try {
            m.setSubMatrix(new double[][] {{}}, 0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
    }

    @Test
    void testWalk() {
        int rows    = 150;
        int columns = 75;

        RealMatrix m = new Array2DRowRealMatrix(rows, columns);
        m.walkInRowOrder(new SetVisitor());
        GetVisitor getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInRowOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            assertEquals(0.0, m.getEntry(i, 0), 0);
            assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            assertEquals(0.0, m.getEntry(0, j), 0);
            assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInColumnOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInColumnOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            assertEquals(0.0, m.getEntry(i, 0), 0);
            assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            assertEquals(0.0, m.getEntry(0, j), 0);
            assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            assertEquals(0.0, m.getEntry(i, 0), 0);
            assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            assertEquals(0.0, m.getEntry(0, j), 0);
            assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            assertEquals(0.0, m.getEntry(i, 0), 0);
            assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            assertEquals(0.0, m.getEntry(0, j), 0);
            assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }
    }

    @Test
    void testSerial()  {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        assertEquals(m,UnitTestUtils.serializeAndRecover(m));
    }


    private static class SetVisitor extends DefaultRealMatrixChangingVisitor {
        @Override
        public double visit(int i, int j, double value) {
            return i + j / 1024.0;
        }
    }

    private static class GetVisitor extends DefaultRealMatrixPreservingVisitor {
        private int count = 0;
        @Override
        public void visit(int i, int j, double value) {
            ++count;
            assertEquals(i + j / 1024.0, value, 0.0);
        }
        public int getCount() {
            return count;
        }
    }

    //--------------- -----------------Protected methods

    /** extracts the l  and u matrices from compact lu representation */
    protected void splitLU(RealMatrix lu, double[][] lowerData, double[][] upperData) {
        if (!lu.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   lu.getRowDimension(), lu.getColumnDimension());
        }
        if (lowerData.length != lowerData[0].length) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   lowerData.length, lowerData[0].length);
        }
        if (upperData.length != upperData[0].length) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   upperData.length, upperData[0].length);
        }
        if (lowerData.length != upperData.length) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   lowerData.length, upperData.length);
        }
        if (lowerData.length != lu.getRowDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   lowerData.length, lu.getRowDimension());
        }

        int n = lu.getRowDimension();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j < i) {
                    lowerData[i][j] = lu.getEntry(i, j);
                    upperData[i][j] = 0d;
                } else if (i == j) {
                    lowerData[i][j] = 1d;
                    upperData[i][j] = lu.getEntry(i, j);
                } else {
                    lowerData[i][j] = 0d;
                    upperData[i][j] = lu.getEntry(i, j);
                }
            }
        }
    }

    /** Returns the result of applying the given row permutation to the matrix */
    protected RealMatrix permuteRows(RealMatrix matrix, int[] permutation) {
        if (!matrix.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   matrix.getRowDimension(), matrix.getColumnDimension());
        }
        if (matrix.getRowDimension() != permutation.length) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   matrix.getRowDimension(), permutation.length);
        }

        int n = matrix.getRowDimension();
        int m = matrix.getColumnDimension();
        double[][] out = new double[m][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out[i][j] = matrix.getEntry(permutation[i], j);
            }
        }
        return new Array2DRowRealMatrix(out);
    }

//    /** Useful for debugging */
//    private void dumpMatrix(RealMatrix m) {
//          for (int i = 0; i < m.getRowDimension(); i++) {
//              String os = "";
//              for (int j = 0; j < m.getColumnDimension(); j++) {
//                  os += m.getEntry(i, j) + " ";
//              }
//              System.out.println(os);
//          }
//    }

    static private final RealMatrixFormat f = new RealMatrixFormat("", "",
                                                                   "\n", "", "", "\t\t", new DecimalFormat(
                                                                                   " ##############0.0000;-##############0.0000"));

    @Test
    void testKroneckerProduct() {
        // AA = [-3 2;1 1];
        // BB = [1 0;0 1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { -3, 2 },
            { 1, 1 } });
        Array2DRowRealMatrix A_ = (Array2DRowRealMatrix) A;

        RealMatrix B = MatrixUtils.createRealMatrix(new double[][] { { 1, 0 },
            { 0, 1 } });

        RealMatrix C = A_.kroneckerProduct(B);

        assertNotNull(C);

        RealMatrix C_expected = MatrixUtils.createRealMatrix(new double[][] {
            { -3, 0, 2, 0 }, { -0, -3, 0, 2 }, { 1, 0, 1, 0 },
            { 0, 1, 0, 1 } });
        assertEquals(f.format(C_expected), f.format(C.scalarAdd(0)));

    }

    @Test
    void testStack() {
        // AA = [-3 2;1 1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { -3, 2 },
            { 1, 1 } });
        Array2DRowRealMatrix A_ = (Array2DRowRealMatrix) A;

        RealMatrix C = A_.stack();

        assertNotNull(C);

        RealMatrix C_expected = MatrixUtils.createRealMatrix(new double[][] {
            { -3 }, { 1 }, { 2 }, { 1 } });
        assertEquals(f.format(C_expected), f.format(C));

    }

    @Test
    void testUnstackSquare() {
        // AA = [-3 ;2;1; 1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { -3 },
            { 1 }, { 2 }, { 1 } });
        Array2DRowRealMatrix A_ = (Array2DRowRealMatrix) A;

        RealMatrix C = A_.unstackSquare();

        assertNotNull(C);

        RealMatrix C_expected = MatrixUtils.createRealMatrix(new double[][] {
            { -3, 2 }, { 1, 1 } });
        assertEquals(f.format(C_expected), f.format(C));

    }

    @Test
    void testUnstackNotsquare() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            // AA = [-3 ;2;1; 1;1];

            RealMatrix A = MatrixUtils.createRealMatrix(new double[][]{{-3},
                {1}, {2}, {1}, {1}});
            Array2DRowRealMatrix A_ = (Array2DRowRealMatrix) A;

            A_.unstackSquare();
        });
    }
}

