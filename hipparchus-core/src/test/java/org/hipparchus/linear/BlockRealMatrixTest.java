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
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link BlockRealMatrix} class.
 *
 */

final class BlockRealMatrixTest {

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

    /** test dimensions */
    @Test
    void testDimensions() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m2 = new BlockRealMatrix(testData2);
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
        Random r = new Random(66636328996002l);
        BlockRealMatrix m1 = createRandomMatrix(r, 47, 83);
        BlockRealMatrix m2 = new BlockRealMatrix(m1.getData());
        assertEquals(m1, m2);
        BlockRealMatrix m3 = new BlockRealMatrix(testData);
        BlockRealMatrix m4 = new BlockRealMatrix(m3.getData());
        assertEquals(m3, m4);
    }

    /** test add */
    @Test
    void testAdd() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix mInv = new BlockRealMatrix(testDataInv);
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
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m2 = new BlockRealMatrix(testData2);
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
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m2 = new BlockRealMatrix(testData2);
        assertEquals(14d,m.getNorm1(),entryTolerance,"testData norm");
        assertEquals(7d,m2.getNorm1(),entryTolerance,"testData2 norm");
        assertEquals(10d,m.getNormInfty(),entryTolerance,"testData norm");
        assertEquals(10d,m2.getNormInfty(),entryTolerance,"testData2 norm");
    }

    /** test Frobenius norm */
    @Test
    void testFrobeniusNorm() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m2 = new BlockRealMatrix(testData2);
        assertEquals(FastMath.sqrt(117.0), m.getFrobeniusNorm(), entryTolerance, "testData Frobenius norm");
        assertEquals(FastMath.sqrt(52.0), m2.getFrobeniusNorm(), entryTolerance, "testData2 Frobenius norm");
    }

    /** test m-n = m + -n */
    @Test
    void testPlusMinus() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m2 = new BlockRealMatrix(testDataInv);
        customAssertClose(m.subtract(m2), m2.scalarMultiply(-1d).add(m), entryTolerance);
        try {
            m.subtract(new BlockRealMatrix(testData2));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test multiply */
    @Test
    void testMultiply() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix mInv = new BlockRealMatrix(testDataInv);
        BlockRealMatrix identity = new BlockRealMatrix(id);
        BlockRealMatrix m2 = new BlockRealMatrix(testData2);
        customAssertClose(m.multiply(mInv), identity, entryTolerance);
        customAssertClose(mInv.multiply(m), identity, entryTolerance);
        customAssertClose(m.multiply(identity), m, entryTolerance);
        customAssertClose(identity.multiply(mInv), mInv, entryTolerance);
        customAssertClose(m2.multiply(identity), m2, entryTolerance);
        try {
            m.multiply(new BlockRealMatrix(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testSeveralBlocks() {
        RealMatrix m = new BlockRealMatrix(35, 71);
        for (int i = 0; i < m.getRowDimension(); ++i) {
            for (int j = 0; j < m.getColumnDimension(); ++j) {
                m.setEntry(i, j, i + j / 1024.0);
            }
        }

        RealMatrix mT = m.transpose();
        assertEquals(m.getRowDimension(), mT.getColumnDimension());
        assertEquals(m.getColumnDimension(), mT.getRowDimension());
        for (int i = 0; i < mT.getRowDimension(); ++i) {
            for (int j = 0; j < mT.getColumnDimension(); ++j) {
                assertEquals(m.getEntry(j, i), mT.getEntry(i, j), 0);
            }
        }

        RealMatrix mPm = m.add(m);
        for (int i = 0; i < mPm.getRowDimension(); ++i) {
            for (int j = 0; j < mPm.getColumnDimension(); ++j) {
                assertEquals(2 * m.getEntry(i, j), mPm.getEntry(i, j), 0);
            }
        }

        RealMatrix mPmMm = mPm.subtract(m);
        for (int i = 0; i < mPmMm.getRowDimension(); ++i) {
            for (int j = 0; j < mPmMm.getColumnDimension(); ++j) {
                assertEquals(m.getEntry(i, j), mPmMm.getEntry(i, j), 0);
            }
        }

        RealMatrix mTm = mT.multiply(m);
        for (int i = 0; i < mTm.getRowDimension(); ++i) {
            for (int j = 0; j < mTm.getColumnDimension(); ++j) {
                double sum = 0;
                for (int k = 0; k < mT.getColumnDimension(); ++k) {
                    sum += (k + i / 1024.0) * (k + j / 1024.0);
                }
                assertEquals(sum, mTm.getEntry(i, j), 0);
            }
        }

        RealMatrix mmT = m.multiply(mT);
        for (int i = 0; i < mmT.getRowDimension(); ++i) {
            for (int j = 0; j < mmT.getColumnDimension(); ++j) {
                double sum = 0;
                for (int k = 0; k < m.getColumnDimension(); ++k) {
                    sum += (i + k / 1024.0) * (j + k / 1024.0);
                }
                assertEquals(sum, mmT.getEntry(i, j), 0);
            }
        }

        RealMatrix sub1 = m.getSubMatrix(2, 9, 5, 20);
        for (int i = 0; i < sub1.getRowDimension(); ++i) {
            for (int j = 0; j < sub1.getColumnDimension(); ++j) {
                assertEquals((i + 2) + (j + 5) / 1024.0, sub1.getEntry(i, j), 0);
            }
        }

        RealMatrix sub2 = m.getSubMatrix(10, 12, 3, 70);
        for (int i = 0; i < sub2.getRowDimension(); ++i) {
            for (int j = 0; j < sub2.getColumnDimension(); ++j) {
                assertEquals((i + 10) + (j + 3) / 1024.0, sub2.getEntry(i, j), 0);
            }
        }

        RealMatrix sub3 = m.getSubMatrix(30, 34, 0, 5);
        for (int i = 0; i < sub3.getRowDimension(); ++i) {
            for (int j = 0; j < sub3.getColumnDimension(); ++j) {
                assertEquals((i + 30) + (j + 0) / 1024.0, sub3.getEntry(i, j), 0);
            }
        }

        RealMatrix sub4 = m.getSubMatrix(30, 32, 62, 65);
        for (int i = 0; i < sub4.getRowDimension(); ++i) {
            for (int j = 0; j < sub4.getColumnDimension(); ++j) {
                assertEquals((i + 30) + (j + 62) / 1024.0, sub4.getEntry(i, j), 0);
            }
        }

    }

    //Additional Test for BlockRealMatrixTest.testMultiply

    private double[][] d3 = new double[][] {{1,2,3,4},{5,6,7,8}};
    private double[][] d4 = new double[][] {{1},{2},{3},{4}};
    private double[][] d5 = new double[][] {{30},{70}};

    @Test
    void testMultiply2() {
        RealMatrix m3 = new BlockRealMatrix(d3);
        RealMatrix m4 = new BlockRealMatrix(d4);
        RealMatrix m5 = new BlockRealMatrix(d5);
        customAssertClose(m3.multiply(m4), m5, entryTolerance);
    }

    @Test
    void testMultiplyTransposedBlockRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xfaa1594a49a1359el);
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final BlockRealMatrix a = new BlockRealMatrix(rows, cols);
                a.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                    public double visit(final int row, final int column, final double value) {
                        return randomGenerator.nextDouble();
                    }
                });
                for (int interm = 1; interm <= 64; interm += 7) {
                    final BlockRealMatrix b = new BlockRealMatrix(interm, cols);
                    b.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                        public double visit(final int row, final int column, final double value) {
                            return randomGenerator.nextDouble();
                        }
                    });
                    assertEquals(0.0,
                                        a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).getNorm1(),
                                        1.0e-15);
                }
            }
        }
    }

    @Test
    void testMultiplyTransposedArray2DRowRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xac2d0185fc69670bl);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final BlockRealMatrix a = new BlockRealMatrix(rows, cols);
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
    void testMultiplyTransposedWrongDimensions() {
        try {
            new BlockRealMatrix(2, 3).multiplyTransposed(new BlockRealMatrix(3, 2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testTransposeMultiplyBlockRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xfaa1594a49a1359el);
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final BlockRealMatrix a = new BlockRealMatrix(rows, cols);
                a.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                    public double visit(final int row, final int column, final double value) {
                        return randomGenerator.nextDouble();
                    }
                });
                for (int interm = 1; interm <= 64; interm += 7) {
                    final BlockRealMatrix b = new BlockRealMatrix(rows, interm);
                    b.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                        public double visit(final int row, final int column, final double value) {
                            return randomGenerator.nextDouble();
                        }
                    });
                    assertEquals(0.0,
                                        a.transposeMultiply(b).subtract(a.transpose().multiply(b)).getNorm1(),
                                        1.0e-15);
                }
            }
        }
    }

    @Test
    void testTransposeMultiplyArray2DRowRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xac2d0185fc69670bl);
        final RealMatrixChangingVisitor randomSetter = new DefaultRealMatrixChangingVisitor() {
            public double visit(final int row, final int column, final double value) {
                return randomGenerator.nextDouble();
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final BlockRealMatrix a = new BlockRealMatrix(rows, cols);
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
    void testTransposeMultiplyWrongDimensions() {
        try {
            new BlockRealMatrix(2, 3).transposeMultiply(new BlockRealMatrix(3, 2));
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
        RealMatrix m = new BlockRealMatrix(id);
        assertEquals(3d,m.getTrace(),entryTolerance,"identity trace");
        m = new BlockRealMatrix(testData2);
        try {
            m.getTrace();
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test scalarAdd */
    @Test
    void testScalarAdd() {
        RealMatrix m = new BlockRealMatrix(testData);
        customAssertClose(new BlockRealMatrix(testDataPlus2), m.scalarAdd(2d), entryTolerance);
    }

    /** test operate */
    @Test
    void testOperate() {
        RealMatrix m = new BlockRealMatrix(id);
        customAssertClose(testVector, m.operate(testVector), entryTolerance);
        customAssertClose(testVector, m.operate(new ArrayRealVector(testVector)).toArray(), entryTolerance);
        m = new BlockRealMatrix(bigSingular);
        try {
            m.operate(testVector);
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testOperateLarge() {
        int p = (7 * BlockRealMatrix.BLOCK_SIZE) / 2;
        int q = (5 * BlockRealMatrix.BLOCK_SIZE) / 2;
        int r =  3 * BlockRealMatrix.BLOCK_SIZE;
        Random random = new Random(111007463902334l);
        RealMatrix m1 = createRandomMatrix(random, p, q);
        RealMatrix m2 = createRandomMatrix(random, q, r);
        RealMatrix m1m2 = m1.multiply(m2);
        for (int i = 0; i < r; ++i) {
            checkArrays(m1m2.getColumn(i), m1.operate(m2.getColumn(i)));
        }
    }

    @Test
    void testOperatePremultiplyLarge() {
        int p = (7 * BlockRealMatrix.BLOCK_SIZE) / 2;
        int q = (5 * BlockRealMatrix.BLOCK_SIZE) / 2;
        int r =  3 * BlockRealMatrix.BLOCK_SIZE;
        Random random = new Random(111007463902334l);
        RealMatrix m1 = createRandomMatrix(random, p, q);
        RealMatrix m2 = createRandomMatrix(random, q, r);
        RealMatrix m1m2 = m1.multiply(m2);
        for (int i = 0; i < p; ++i) {
            checkArrays(m1m2.getRow(i), m2.preMultiply(m1.getRow(i)));
        }
    }

    /** test issue MATH-209 */
    @Test
    void testMath209() {
        RealMatrix a = new BlockRealMatrix(new double[][] {
                { 1, 2 }, { 3, 4 }, { 5, 6 }
        });
        double[] b = a.operate(new double[] { 1, 1 });
        assertEquals(a.getRowDimension(), b.length);
        assertEquals( 3.0, b[0], 1.0e-12);
        assertEquals( 7.0, b[1], 1.0e-12);
        assertEquals(11.0, b[2], 1.0e-12);
    }

    /** test transpose */
    @Test
    void testTranspose() {
        RealMatrix m = new BlockRealMatrix(testData);
        RealMatrix mIT = new LUDecomposition(m).getSolver().getInverse().transpose();
        RealMatrix mTI = new LUDecomposition(m.transpose()).getSolver().getInverse();
        customAssertClose(mIT, mTI, normTolerance);
        m = new BlockRealMatrix(testData2);
        RealMatrix mt = new BlockRealMatrix(testData2T);
        customAssertClose(mt, m.transpose(), normTolerance);
    }

    /** test preMultiply by vector */
    @Test
    void testPremultiplyVector() {
        RealMatrix m = new BlockRealMatrix(testData);
        customAssertClose(m.preMultiply(testVector), preMultTest, normTolerance);
        customAssertClose(m.preMultiply(new ArrayRealVector(testVector).toArray()),
                          preMultTest, normTolerance);
        m = new BlockRealMatrix(bigSingular);
        try {
            m.preMultiply(testVector);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testArithmeticBlending() {

        // Given
        final RealMatrix m1 = new BlockRealMatrix(new double[][] {
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 },
        });
        final RealMatrix m2 = new BlockRealMatrix(new double[][] {
                { 10, 11, 12 },
                { 13, 14, 15 },
                { 16, 17, 18 },
        });

        final double blendingValue = 0.2;

        // When
        final RealMatrix blendedMatrix = m1.blendArithmeticallyWith(m2, blendingValue);

        // Then
        final RealMatrix expectedMatrix = new BlockRealMatrix(new double[][] {
                { 2.8, 3.8, 4.8 },
                { 5.8, 6.8, 7.8 },
                { 8.8, 9.8, 10.8 }
        });

        customAssertClose(expectedMatrix, blendedMatrix, normTolerance);
    }

    @Test
    void testPremultiply() {
        RealMatrix m3 = new BlockRealMatrix(d3);
        RealMatrix m4 = new BlockRealMatrix(d4);
        RealMatrix m5 = new BlockRealMatrix(d5);
        customAssertClose(m4.preMultiply(m3), m5, entryTolerance);

        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix mInv = new BlockRealMatrix(testDataInv);
        BlockRealMatrix identity = new BlockRealMatrix(id);
        customAssertClose(m.preMultiply(mInv), identity, entryTolerance);
        customAssertClose(mInv.preMultiply(m), identity, entryTolerance);
        customAssertClose(m.preMultiply(identity), m, entryTolerance);
        customAssertClose(identity.preMultiply(mInv), mInv, entryTolerance);
        try {
            m.preMultiply(new BlockRealMatrix(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testGetVectors() {
        RealMatrix m = new BlockRealMatrix(testData);
        customAssertClose(m.getRow(0), testDataRow1, entryTolerance);
        customAssertClose(m.getColumn(2), testDataCol3, entryTolerance);
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
        RealMatrix m = new BlockRealMatrix(testData);
        assertEquals(2d, m.getEntry(0, 1), entryTolerance, "get entry");
        try {
            m.getEntry(10, 4);
            fail ("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    /** test examples in user guide */
    @Test
    void testExamples() {
        // Create a real matrix with two rows and three columns
        double[][] matrixData = { {1d,2d,3d}, {2d,5d,3d}};
        RealMatrix m = new BlockRealMatrix(matrixData);
        // One more with three rows, two columns
        double[][] matrixData2 = { {1d,2d}, {2d,5d}, {1d, 7d}};
        RealMatrix n = new BlockRealMatrix(matrixData2);
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
        RealMatrix coefficients = new BlockRealMatrix(coefficientsData);
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
        RealMatrix m = new BlockRealMatrix(subTestData);
        checkGetSubMatrix(m, subRows23Cols00,  2 , 3 , 0, 0);
        checkGetSubMatrix(m, subRows00Cols33,  0 , 0 , 3, 3);
        checkGetSubMatrix(m, subRows01Cols23,  0 , 1 , 2, 3);
        checkGetSubMatrix(m, subRows02Cols13,  new int[] { 0, 2 }, new int[] { 1, 3 });
        checkGetSubMatrix(m, subRows03Cols12,  new int[] { 0, 3 }, new int[] { 1, 2 });
        checkGetSubMatrix(m, subRows03Cols123, new int[] { 0, 3 }, new int[] { 1, 2, 3 });
        checkGetSubMatrix(m, subRows20Cols123, new int[] { 2, 0 }, new int[] { 1, 2, 3 });
        checkGetSubMatrix(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 });
        checkGetSubMatrix(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 });
        checkGetSubMatrix(m, null,  1, 0, 2, 4);
        checkGetSubMatrix(m, null, -1, 1, 2, 2);
        checkGetSubMatrix(m, null,  1, 0, 2, 2);
        checkGetSubMatrix(m, null,  1, 0, 2, 4);
        checkGetSubMatrix(m, null, new int[] {},    new int[] { 0 });
        checkGetSubMatrix(m, null, new int[] { 0 }, new int[] { 4 });
    }

    private void checkGetSubMatrix(RealMatrix m, double[][] reference,
                                   int startRow, int endRow, int startColumn, int endColumn) {
        try {
            RealMatrix sub = m.getSubMatrix(startRow, endRow, startColumn, endColumn);
            if (reference != null) {
                assertEquals(new BlockRealMatrix(reference), sub);
            } else {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentException or MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (reference != null) {
                throw e;
            }
        }
    }

    private void checkGetSubMatrix(RealMatrix m, double[][] reference,
                                   int[] selectedRows, int[] selectedColumns) {
        try {
            RealMatrix sub = m.getSubMatrix(selectedRows, selectedColumns);
            if (reference != null) {
                assertEquals(new BlockRealMatrix(reference), sub);
            } else {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentExceptiono r MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (reference != null) {
                throw e;
            }
        }
    }

    @Test
    void testGetSetMatrixLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        RealMatrix sub = new BlockRealMatrix(n - 4, n - 4).scalarAdd(1);

        m.setSubMatrix(sub.getData(), 2, 2);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if ((i < 2) || (i > n - 3) || (j < 2) || (j > n - 3)) {
                    assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        assertEquals(sub, m.getSubMatrix(2, n - 3, 2, n - 3));

    }

    @Test
    void testCopySubMatrix() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        checkCopy(m, subRows23Cols00,  2 , 3 , 0, 0);
        checkCopy(m, subRows00Cols33,  0 , 0 , 3, 3);
        checkCopy(m, subRows01Cols23,  0 , 1 , 2, 3);
        checkCopy(m, subRows02Cols13,  new int[] { 0, 2 }, new int[] { 1, 3 });
        checkCopy(m, subRows03Cols12,  new int[] { 0, 3 }, new int[] { 1, 2 });
        checkCopy(m, subRows03Cols123, new int[] { 0, 3 }, new int[] { 1, 2, 3 });
        checkCopy(m, subRows20Cols123, new int[] { 2, 0 }, new int[] { 1, 2, 3 });
        checkCopy(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 });
        checkCopy(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 });

        checkCopy(m, null,  1, 0, 2, 4);
        checkCopy(m, null, -1, 1, 2, 2);
        checkCopy(m, null,  1, 0, 2, 2);
        checkCopy(m, null,  1, 0, 2, 4);
        checkCopy(m, null, new int[] {},    new int[] { 0 });
        checkCopy(m, null, new int[] { 0 }, new int[] { 4 });
    }

    private void checkCopy(RealMatrix m, double[][] reference,
                           int startRow, int endRow, int startColumn, int endColumn) {
        try {
            double[][] sub = (reference == null) ?
                             new double[1][1] :
                             new double[reference.length][reference[0].length];
            m.copySubMatrix(startRow, endRow, startColumn, endColumn, sub);
            if (reference != null) {
                assertEquals(new BlockRealMatrix(reference), new BlockRealMatrix(sub));
            } else {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentException or MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (reference != null) {
                throw e;
            }
        }
    }

    private void checkCopy(RealMatrix m, double[][] reference,
                           int[] selectedRows, int[] selectedColumns) {
        try {
            double[][] sub = (reference == null) ?
                    new double[1][1] :
                    new double[reference.length][reference[0].length];
            m.copySubMatrix(selectedRows, selectedColumns, sub);
            if (reference != null) {
                assertEquals(new BlockRealMatrix(reference), new BlockRealMatrix(sub));
            } else {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentException or MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (reference != null) {
                throw e;
            }
        }
    }

    @Test
    void testGetRowMatrix() {
        RealMatrix m     = new BlockRealMatrix(subTestData);
        RealMatrix mRow0 = new BlockRealMatrix(subRow0);
        RealMatrix mRow3 = new BlockRealMatrix(subRow3);
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
    void testSetRowMatrix() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        RealMatrix mRow3 = new BlockRealMatrix(subRow3);
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
    void testGetSetRowMatrixLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        RealMatrix sub = new BlockRealMatrix(1, n).scalarAdd(1);

        m.setRowMatrix(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        assertEquals(sub, m.getRowMatrix(2));
    }

    @Test
    void testGetColumnMatrix() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        RealMatrix mColumn1 = new BlockRealMatrix(subColumn1);
        RealMatrix mColumn3 = new BlockRealMatrix(subColumn3);
        assertEquals(mColumn1, m.getColumnMatrix(1));
        assertEquals(mColumn3, m.getColumnMatrix(3));
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
        RealMatrix m = new BlockRealMatrix(subTestData);
        RealMatrix mColumn3 = new BlockRealMatrix(subColumn3);
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
    void testGetSetColumnMatrixLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        RealMatrix sub = new BlockRealMatrix(n, 1).scalarAdd(1);

        m.setColumnMatrix(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        assertEquals(sub, m.getColumnMatrix(2));

    }

    @Test
    void testGetRowVector() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        RealVector mRow0 = new ArrayRealVector(subRow0[0]);
        RealVector mRow3 = new ArrayRealVector(subRow3[0]);
        assertEquals(mRow0, m.getRowVector(0));
        assertEquals(mRow3, m.getRowVector(3));
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
        RealMatrix m = new BlockRealMatrix(subTestData);
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
    void testGetSetRowVectorLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        RealVector sub = new ArrayRealVector(n, 1.0);

        m.setRowVector(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        assertEquals(sub, m.getRowVector(2));
    }

    @Test
    void testGetColumnVector() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        RealVector mColumn1 = columnToVector(subColumn1);
        RealVector mColumn3 = columnToVector(subColumn3);
        assertEquals(mColumn1, m.getColumnVector(1));
        assertEquals(mColumn3, m.getColumnVector(3));
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
        RealMatrix m = new BlockRealMatrix(subTestData);
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

    @Test
    void testGetSetColumnVectorLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        RealVector sub = new ArrayRealVector(n, 1.0);

        m.setColumnVector(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        assertEquals(sub, m.getColumnVector(2));
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
        RealMatrix m = new BlockRealMatrix(subTestData);
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
        RealMatrix m = new BlockRealMatrix(subTestData);
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
    void testGetSetRowLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        double[] sub = new double[n];
        Arrays.fill(sub, 1.0);

        m.setRow(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        checkArrays(sub, m.getRow(2));
    }

    @Test
    void testGetColumn() {
        RealMatrix m = new BlockRealMatrix(subTestData);
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
        RealMatrix m = new BlockRealMatrix(subTestData);
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

    @Test
    void testGetSetColumnLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        double[] sub = new double[n];
        Arrays.fill(sub, 1.0);

        m.setColumn(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        checkArrays(sub, m.getColumn(2));
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
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m1 = m.copy();
        BlockRealMatrix mt = m.transpose();
        assertTrue(m.hashCode() != mt.hashCode());
        assertEquals(m.hashCode(), m1.hashCode());
        assertEquals(m, m);
        assertEquals(m, m1);
        assertNotEquals(null, m);
        assertNotEquals(m, mt);
        assertNotEquals(m, new BlockRealMatrix(bigSingular));
    }

    @Test
    void testToString() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        assertEquals("BlockRealMatrix{{1.0,2.0,3.0},{2.0,5.0,3.0},{1.0,0.0,8.0}}",
                m.toString());
    }

    @Test
    void testSetSubMatrix() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        m.setSubMatrix(detData2,1,1);
        RealMatrix expected = new BlockRealMatrix
            (new double[][] {{1.0,2.0,3.0},{2.0,1.0,3.0},{1.0,2.0,4.0}});
        assertEquals(expected, m);

        m.setSubMatrix(detData2,0,0);
        expected = new BlockRealMatrix
            (new double[][] {{1.0,3.0,3.0},{2.0,4.0,3.0},{1.0,2.0,4.0}});
        assertEquals(expected, m);

        m.setSubMatrix(testDataPlus2,0,0);
        expected = new BlockRealMatrix
            (new double[][] {{3.0,4.0,5.0},{4.0,7.0,5.0},{3.0,2.0,10.0}});
        assertEquals(expected, m);

        // javadoc example
        BlockRealMatrix matrix = new BlockRealMatrix
            (new double[][] {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 0, 1 , 2}});
        matrix.setSubMatrix(new double[][] {{3, 4}, {5, 6}}, 1, 1);
        expected = new BlockRealMatrix
            (new double[][] {{1, 2, 3, 4}, {5, 3, 4, 8}, {9, 5 ,6, 2}});
        assertEquals(expected, matrix);

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

        RealMatrix m = new BlockRealMatrix(rows, columns);
        m.walkInRowOrder(new SetVisitor());
        GetVisitor getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockRealMatrix(rows, columns);
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

        m = new BlockRealMatrix(rows, columns);
        m.walkInColumnOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockRealMatrix(rows, columns);
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

        m = new BlockRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockRealMatrix(rows, columns);
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

        m = new BlockRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockRealMatrix(rows, columns);
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
        BlockRealMatrix m = new BlockRealMatrix(testData);
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

    /** verifies that two matrices are close (1-norm) */
    protected void customAssertClose(RealMatrix m, RealMatrix n, double tolerance) {
        assertTrue(m.subtract(n).getNorm1() < tolerance);
    }

    /** verifies that two vectors are close (sup norm) */
    protected void customAssertClose(double[] m, double[] n, double tolerance) {
        if (m.length != n.length) {
            fail("vectors not same length");
        }
        for (int i = 0; i < m.length; i++) {
            assertEquals(m[i], n[i], tolerance);
        }
    }

    private BlockRealMatrix createRandomMatrix(Random r, int rows, int columns) {
        BlockRealMatrix m = new BlockRealMatrix(rows, columns);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                m.setEntry(i, j, 200 * r.nextDouble() - 100);
            }
        }
        return m;
    }
}
