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
import org.hipparchus.fraction.Fraction;
import org.hipparchus.fraction.FractionField;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
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
 * Test cases for the {@link BlockFieldMatrix} class.
 *
 */

final class BlockFieldMatrixTest {

    // 3 x 3 identity matrix
    protected Fraction[][] id = {
            {new Fraction(1),new Fraction(0),new Fraction(0)},
            {new Fraction(0),new Fraction(1),new Fraction(0)},
            {new Fraction(0),new Fraction(0),new Fraction(1)}
    };

    // Test data for group operations
    protected Fraction[][] testData = {
            {new Fraction(1),new Fraction(2),new Fraction(3)},
            {new Fraction(2),new Fraction(5),new Fraction(3)},
            {new Fraction(1),new Fraction(0),new Fraction(8)}
    };
    protected Fraction[][] testDataLU = {
            {new Fraction(2), new Fraction(5), new Fraction(3)},
            {new Fraction(1, 2), new Fraction(-5, 2), new Fraction(13, 2)},
            {new Fraction(1, 2), new Fraction(1, 5), new Fraction(1, 5)}
    };
    protected Fraction[][] testDataPlus2 = {
            {new Fraction(3),new Fraction(4),new Fraction(5)},
            {new Fraction(4),new Fraction(7),new Fraction(5)},
            {new Fraction(3),new Fraction(2),new Fraction(10)}
    };
    protected Fraction[][] testDataMinus = {
            {new Fraction(-1),new Fraction(-2),new Fraction(-3)},
            {new Fraction(-2),new Fraction(-5),new Fraction(-3)},
            {new Fraction(-1),new Fraction(0),new Fraction(-8)}
    };
    protected Fraction[] testDataRow1 = {new Fraction(1),new Fraction(2),new Fraction(3)};
    protected Fraction[] testDataCol3 = {new Fraction(3),new Fraction(3),new Fraction(8)};
    protected Fraction[][] testDataInv = {
            {new Fraction(-40),new Fraction(16),new Fraction(9)},
            {new Fraction(13),new Fraction(-5),new Fraction(-3)},
            {new Fraction(5),new Fraction(-2),new Fraction(-1)}
    };
    protected Fraction[] preMultTest = {new Fraction(8), new Fraction(12), new Fraction(33)};
    protected Fraction[][] testData2 = {
            {new Fraction(1),new Fraction(2),new Fraction(3)},
            {new Fraction(2),new Fraction(5),new Fraction(3)}
    };
    protected Fraction[][] testData2T = {
            {new Fraction(1),new Fraction(2)},
            {new Fraction(2),new Fraction(5)},
            {new Fraction(3),new Fraction(3)}
    };
    protected Fraction[][] testDataPlusInv = {
            {new Fraction(-39),new Fraction(18),new Fraction(12)},
            {new Fraction(15),new Fraction(0),new Fraction(0)},
            {new Fraction(6),new Fraction(-2),new Fraction(7)}
    };

    // lu decomposition tests
    protected Fraction[][] luData = {
            {new Fraction(2),new Fraction(3),new Fraction(3)},
            {new Fraction(0),new Fraction(5),new Fraction(7)},
            {new Fraction(6),new Fraction(9),new Fraction(8)}
    };
    protected Fraction[][] luDataLUDecomposition = {
            {new Fraction(6),new Fraction(9),new Fraction(8)},
            {new Fraction(0),new Fraction(5),new Fraction(7)},
            {new Fraction(1, 3),new Fraction(0),new Fraction(1, 3)}
    };

    // singular matrices
    protected Fraction[][] singular = { {new Fraction(2),new Fraction(3)}, {new Fraction(2),new Fraction(3)} };
    protected Fraction[][] bigSingular = {
            {new Fraction(1),new Fraction(2),new Fraction(3),new Fraction(4)},
            {new Fraction(2),new Fraction(5),new Fraction(3),new Fraction(4)},
            {new Fraction(7),new Fraction(3),new Fraction(256),new Fraction(1930)},
            {new Fraction(3),new Fraction(7),new Fraction(6),new Fraction(8)}
    }; // 4th row = 1st + 2nd
    protected Fraction[][] detData = {
            {new Fraction(1),new Fraction(2),new Fraction(3)},
            {new Fraction(4),new Fraction(5),new Fraction(6)},
            {new Fraction(7),new Fraction(8),new Fraction(10)}
    };
    protected Fraction[][] detData2 = { {new Fraction(1), new Fraction(3)}, {new Fraction(2), new Fraction(4)}};

    // vectors
    protected Fraction[] testVector = {new Fraction(1),new Fraction(2),new Fraction(3)};
    protected Fraction[] testVector2 = {new Fraction(1),new Fraction(2),new Fraction(3),new Fraction(4)};

    // submatrix accessor tests
    protected Fraction[][] subTestData = {
            {new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4)},
            {new Fraction(3, 2), new Fraction(5, 2), new Fraction(7, 2), new Fraction(9, 2)},
            {new Fraction(2), new Fraction(4), new Fraction(6), new Fraction(8)},
            {new Fraction(4), new Fraction(5), new Fraction(6), new Fraction(7)}
    };
    // array selections
    protected Fraction[][] subRows02Cols13 = { {new Fraction(2), new Fraction(4)}, {new Fraction(4), new Fraction(8)}};
    protected Fraction[][] subRows03Cols12 = { {new Fraction(2), new Fraction(3)}, {new Fraction(5), new Fraction(6)}};
    protected Fraction[][] subRows03Cols123 = {
            {new Fraction(2), new Fraction(3), new Fraction(4)},
            {new Fraction(5), new Fraction(6), new Fraction(7)}
    };
    // effective permutations
    protected Fraction[][] subRows20Cols123 = {
            {new Fraction(4), new Fraction(6), new Fraction(8)},
            {new Fraction(2), new Fraction(3), new Fraction(4)}
    };
    protected Fraction[][] subRows31Cols31 = {{new Fraction(7), new Fraction(5)}, {new Fraction(9, 2), new Fraction(5, 2)}};
    // contiguous ranges
    protected Fraction[][] subRows01Cols23 = {{new Fraction(3),new Fraction(4)} , {new Fraction(7, 2), new Fraction(9, 2)}};
    protected Fraction[][] subRows23Cols00 = {{new Fraction(2)} , {new Fraction(4)}};
    protected Fraction[][] subRows00Cols33 = {{new Fraction(4)}};
    // row matrices
    protected Fraction[][] subRow0 = {{new Fraction(1),new Fraction(2),new Fraction(3),new Fraction(4)}};
    protected Fraction[][] subRow3 = {{new Fraction(4),new Fraction(5),new Fraction(6),new Fraction(7)}};
    // column matrices
    protected Fraction[][] subColumn1 = {{new Fraction(2)}, {new Fraction(5, 2)}, {new Fraction(4)}, {new Fraction(5)}};
    protected Fraction[][] subColumn3 = {{new Fraction(4)}, {new Fraction(9, 2)}, {new Fraction(8)}, {new Fraction(7)}};

    // tolerances
    protected double entryTolerance = 10E-16;
    protected double normTolerance = 10E-14;

    /** test dimensions */
    @Test
    void testDimensions() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> m2 = new BlockFieldMatrix<Fraction>(testData2);
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
        BlockFieldMatrix<Fraction> m1 = createRandomMatrix(r, 47, 83);
        BlockFieldMatrix<Fraction> m2 = new BlockFieldMatrix<Fraction>(m1.getData());
        assertEquals(m1, m2);
        BlockFieldMatrix<Fraction> m3 = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> m4 = new BlockFieldMatrix<Fraction>(m3.getData());
        assertEquals(m3, m4);
    }

    /** test add */
    @Test
    void testAdd() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> mInv = new BlockFieldMatrix<Fraction>(testDataInv);
        FieldMatrix<Fraction> mPlusMInv = m.add(mInv);
        Fraction[][] sumEntries = mPlusMInv.getData();
        for (int row = 0; row < m.getRowDimension(); row++) {
            for (int col = 0; col < m.getColumnDimension(); col++) {
                assertEquals(testDataPlusInv[row][col],sumEntries[row][col]);
            }
        }
    }

    /** test add failure */
    @Test
    void testAddFail() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> m2 = new BlockFieldMatrix<Fraction>(testData2);
        try {
            m.add(m2);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test m-n = m + -n */
    @Test
    void testPlusMinus() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> m2 = new BlockFieldMatrix<Fraction>(testDataInv);
        UnitTestUtils.customAssertEquals(m.subtract(m2), m2.scalarMultiply(new Fraction(-1)).add(m));
        try {
            m.subtract(new BlockFieldMatrix<Fraction>(testData2));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test multiply */
    @Test
    void testMultiply() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> mInv = new BlockFieldMatrix<Fraction>(testDataInv);
        BlockFieldMatrix<Fraction> identity = new BlockFieldMatrix<Fraction>(id);
        BlockFieldMatrix<Fraction> m2 = new BlockFieldMatrix<Fraction>(testData2);
        UnitTestUtils.customAssertEquals(m.multiply(mInv), identity);
        UnitTestUtils.customAssertEquals(mInv.multiply(m), identity);
        UnitTestUtils.customAssertEquals(m.multiply(identity), m);
        UnitTestUtils.customAssertEquals(identity.multiply(mInv), mInv);
        UnitTestUtils.customAssertEquals(m2.multiply(identity), m2);
        try {
            m.multiply(new BlockFieldMatrix<Fraction>(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testSeveralBlocks() {
        FieldMatrix<Fraction> m =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), 37, 41);
        for (int i = 0; i < m.getRowDimension(); ++i) {
            for (int j = 0; j < m.getColumnDimension(); ++j) {
                m.setEntry(i, j, new Fraction(i * 11 + j, 11));
            }
        }

        FieldMatrix<Fraction> mT = m.transpose();
        assertEquals(m.getRowDimension(), mT.getColumnDimension());
        assertEquals(m.getColumnDimension(), mT.getRowDimension());
        for (int i = 0; i < mT.getRowDimension(); ++i) {
            for (int j = 0; j < mT.getColumnDimension(); ++j) {
                assertEquals(m.getEntry(j, i), mT.getEntry(i, j));
            }
        }

        FieldMatrix<Fraction> mPm = m.add(m);
        for (int i = 0; i < mPm.getRowDimension(); ++i) {
            for (int j = 0; j < mPm.getColumnDimension(); ++j) {
                assertEquals(m.getEntry(i, j).multiply(new Fraction(2)), mPm.getEntry(i, j));
            }
        }

        FieldMatrix<Fraction> mPmMm = mPm.subtract(m);
        for (int i = 0; i < mPmMm.getRowDimension(); ++i) {
            for (int j = 0; j < mPmMm.getColumnDimension(); ++j) {
                assertEquals(m.getEntry(i, j), mPmMm.getEntry(i, j));
            }
        }

        FieldMatrix<Fraction> mTm = mT.multiply(m);
        for (int i = 0; i < mTm.getRowDimension(); ++i) {
            for (int j = 0; j < mTm.getColumnDimension(); ++j) {
                Fraction sum = Fraction.ZERO;
                for (int k = 0; k < mT.getColumnDimension(); ++k) {
                    sum = sum.add(new Fraction(k * 11 + i, 11).multiply(new Fraction(k * 11 + j, 11)));
                }
                assertEquals(sum, mTm.getEntry(i, j));
            }
        }

        FieldMatrix<Fraction> mmT = m.multiply(mT);
        for (int i = 0; i < mmT.getRowDimension(); ++i) {
            for (int j = 0; j < mmT.getColumnDimension(); ++j) {
                Fraction sum = Fraction.ZERO;
                for (int k = 0; k < m.getColumnDimension(); ++k) {
                    sum = sum.add(new Fraction(i * 11 + k, 11).multiply(new Fraction(j * 11 + k, 11)));
                }
                assertEquals(sum, mmT.getEntry(i, j));
            }
        }

        FieldMatrix<Fraction> sub1 = m.getSubMatrix(2, 9, 5, 20);
        for (int i = 0; i < sub1.getRowDimension(); ++i) {
            for (int j = 0; j < sub1.getColumnDimension(); ++j) {
                assertEquals(new Fraction((i + 2) * 11 + (j + 5), 11), sub1.getEntry(i, j));
            }
        }

        FieldMatrix<Fraction> sub2 = m.getSubMatrix(10, 12, 3, 40);
        for (int i = 0; i < sub2.getRowDimension(); ++i) {
            for (int j = 0; j < sub2.getColumnDimension(); ++j) {
                assertEquals(new Fraction((i + 10) * 11 + (j + 3), 11), sub2.getEntry(i, j));
            }
        }

        FieldMatrix<Fraction> sub3 = m.getSubMatrix(30, 34, 0, 5);
        for (int i = 0; i < sub3.getRowDimension(); ++i) {
            for (int j = 0; j < sub3.getColumnDimension(); ++j) {
                assertEquals(new Fraction((i + 30) * 11 + (j + 0), 11), sub3.getEntry(i, j));
            }
        }

        FieldMatrix<Fraction> sub4 = m.getSubMatrix(30, 32, 32, 35);
        for (int i = 0; i < sub4.getRowDimension(); ++i) {
            for (int j = 0; j < sub4.getColumnDimension(); ++j) {
                assertEquals(new Fraction((i + 30) * 11 + (j + 32), 11), sub4.getEntry(i, j));
            }
        }

    }

    //Additional Test for BlockFieldMatrix<Fraction>Test.testMultiply

    private Fraction[][] d3 = new Fraction[][] {
            {new Fraction(1),new Fraction(2),new Fraction(3),new Fraction(4)},
            {new Fraction(5),new Fraction(6),new Fraction(7),new Fraction(8)}
    };
    private Fraction[][] d4 = new Fraction[][] {
            {new Fraction(1)},
            {new Fraction(2)},
            {new Fraction(3)},
            {new Fraction(4)}
    };
    private Fraction[][] d5 = new Fraction[][] {{new Fraction(30)},{new Fraction(70)}};

    @Test
    void testMultiply2() {
       FieldMatrix<Fraction> m3 = new BlockFieldMatrix<Fraction>(d3);
       FieldMatrix<Fraction> m4 = new BlockFieldMatrix<Fraction>(d4);
       FieldMatrix<Fraction> m5 = new BlockFieldMatrix<Fraction>(d5);
       UnitTestUtils.customAssertEquals(m3.multiply(m4), m5);
    }

    @Test
    void testMultiplyTransposedBlockRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xfaa1594a49a1359el);
        final FieldMatrixChangingVisitor<Binary64> randomSetter = new DefaultFieldMatrixChangingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public Binary64 visit(final int row, final int column, final Binary64 value) {
                return new Binary64(randomGenerator.nextDouble());
            }
        };
        final FieldMatrixPreservingVisitor<Binary64> zeroChecker = new DefaultFieldMatrixPreservingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public void visit(final int row, final int column, final Binary64 value) {
                assertEquals(0.0, value.doubleValue(), 1.0e-15);
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final BlockFieldMatrix<Binary64> a = new BlockFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final BlockFieldMatrix<Binary64> b = new BlockFieldMatrix<>(Binary64Field.getInstance(), interm, cols);
                    b.walkInOptimizedOrder(randomSetter);
                    a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).walkInOptimizedOrder(zeroChecker);
                }
            }
        }
    }

    @Test
    void testMultiplyTransposedArray2DRowRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xac2d0185fc69670bl);
        final FieldMatrixChangingVisitor<Binary64> randomSetter = new DefaultFieldMatrixChangingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public Binary64 visit(final int row, final int column, final Binary64 value) {
                return new Binary64(randomGenerator.nextDouble());
            }
        };
        final FieldMatrixPreservingVisitor<Binary64> zeroChecker = new DefaultFieldMatrixPreservingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public void visit(final int row, final int column, final Binary64 value) {
                assertEquals(0.0, value.doubleValue(), 1.0e-15);
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final BlockFieldMatrix<Binary64> a = new BlockFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final Array2DRowFieldMatrix<Binary64> b = new Array2DRowFieldMatrix<>(Binary64Field.getInstance(), interm, cols);
                    b.walkInOptimizedOrder(randomSetter);
                    a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).walkInOptimizedOrder(zeroChecker);
                }
            }
        }
    }

    @Test
    void testMultiplyTransposedWrongDimensions() {
        try {
            new BlockFieldMatrix<Binary64>(Binary64Field.getInstance(), 2, 3).
            multiplyTransposed(new BlockFieldMatrix<Binary64>(Binary64Field.getInstance(), 3, 2));
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
        final FieldMatrixChangingVisitor<Binary64> randomSetter = new DefaultFieldMatrixChangingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public Binary64 visit(final int row, final int column, final Binary64 value) {
                return new Binary64(randomGenerator.nextDouble());
            }
        };
        final FieldMatrixPreservingVisitor<Binary64> zeroChecker = new DefaultFieldMatrixPreservingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public void visit(final int row, final int column, final Binary64 value) {
                assertEquals(0.0, value.doubleValue(), 2.0e-14);
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final BlockFieldMatrix<Binary64> a = new BlockFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final BlockFieldMatrix<Binary64> b = new BlockFieldMatrix<>(Binary64Field.getInstance(), rows, interm);
                    b.walkInOptimizedOrder(randomSetter);
                    a.transposeMultiply(b).subtract(a.transpose().multiply(b)).walkInOptimizedOrder(zeroChecker);
                }
            }
        }
    }

    @Test
    void testTransposeMultiplyArray2DRowRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xac2d0185fc69670bl);
        final FieldMatrixChangingVisitor<Binary64> randomSetter = new DefaultFieldMatrixChangingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public Binary64 visit(final int row, final int column, final Binary64 value) {
                return new Binary64(randomGenerator.nextDouble());
            }
        };
        final FieldMatrixPreservingVisitor<Binary64> zeroChecker = new DefaultFieldMatrixPreservingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public void visit(final int row, final int column, final Binary64 value) {
                assertEquals(0.0, value.doubleValue(), 1.0e-15);
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final BlockFieldMatrix<Binary64> a = new BlockFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final Array2DRowFieldMatrix<Binary64> b = new Array2DRowFieldMatrix<>(Binary64Field.getInstance(), rows, interm);
                    b.walkInOptimizedOrder(randomSetter);
                    a.transposeMultiply(b).subtract(a.transpose().multiply(b)).walkInOptimizedOrder(zeroChecker);
                }
            }
        }
    }

    @Test
    void testTransposeMultiplyWrongDimensions() {
        try {
            new BlockFieldMatrix<Binary64>(Binary64Field.getInstance(), 2, 3).
            transposeMultiply(new BlockFieldMatrix<Binary64>(Binary64Field.getInstance(), 3, 2));
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(id);
        assertEquals(new Fraction(3),m.getTrace());
        m = new BlockFieldMatrix<Fraction>(testData2);
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        UnitTestUtils.customAssertEquals(new BlockFieldMatrix<Fraction>(testDataPlus2),
                                         m.scalarAdd(new Fraction(2)));
    }

    /** test operate */
    @Test
    void testOperate() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(id);
        UnitTestUtils.customAssertEquals(testVector, m.operate(testVector));
        UnitTestUtils.customAssertEquals(testVector, m.operate(new ArrayFieldVector<Fraction>(testVector)).toArray());
        m = new BlockFieldMatrix<Fraction>(bigSingular);
        try {
            m.operate(testVector);
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testOperateLarge() {
        int p = (11 * BlockFieldMatrix.BLOCK_SIZE) / 10;
        int q = (11 * BlockFieldMatrix.BLOCK_SIZE) / 10;
        int r =  BlockFieldMatrix.BLOCK_SIZE / 2;
        Random random = new Random(111007463902334l);
        FieldMatrix<Fraction> m1 = createRandomMatrix(random, p, q);
        FieldMatrix<Fraction> m2 = createRandomMatrix(random, q, r);
        FieldMatrix<Fraction> m1m2 = m1.multiply(m2);
        for (int i = 0; i < r; ++i) {
            UnitTestUtils.customAssertEquals(m1m2.getColumn(i), m1.operate(m2.getColumn(i)));
        }
    }

    @Test
    void testOperatePremultiplyLarge() {
        int p = (11 * BlockFieldMatrix.BLOCK_SIZE) / 10;
        int q = (11 * BlockFieldMatrix.BLOCK_SIZE) / 10;
        int r =  BlockFieldMatrix.BLOCK_SIZE / 2;
        Random random = new Random(111007463902334l);
        FieldMatrix<Fraction> m1 = createRandomMatrix(random, p, q);
        FieldMatrix<Fraction> m2 = createRandomMatrix(random, q, r);
        FieldMatrix<Fraction> m1m2 = m1.multiply(m2);
        for (int i = 0; i < p; ++i) {
            UnitTestUtils.customAssertEquals(m1m2.getRow(i), m2.preMultiply(m1.getRow(i)));
        }
    }

    /** test issue MATH-209 */
    @Test
    void testMath209() {
        FieldMatrix<Fraction> a = new BlockFieldMatrix<Fraction>(new Fraction[][] {
                { new Fraction(1), new Fraction(2) },
                { new Fraction(3), new Fraction(4) },
                { new Fraction(5), new Fraction(6) }
        });
        Fraction[] b = a.operate(new Fraction[] { new Fraction(1), new Fraction(1) });
        assertEquals(a.getRowDimension(), b.length);
        assertEquals( new Fraction(3), b[0]);
        assertEquals( new Fraction(7), b[1]);
        assertEquals(new Fraction(11), b[2]);
    }

    /** test transpose */
    @Test
    void testTranspose() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        FieldMatrix<Fraction> mIT = new FieldLUDecomposition<Fraction>(m).getSolver().getInverse().transpose();
        FieldMatrix<Fraction> mTI = new FieldLUDecomposition<Fraction>(m.transpose()).getSolver().getInverse();
        UnitTestUtils.customAssertEquals(mIT, mTI);
        m = new BlockFieldMatrix<Fraction>(testData2);
        FieldMatrix<Fraction> mt = new BlockFieldMatrix<Fraction>(testData2T);
        UnitTestUtils.customAssertEquals(mt, m.transpose());
    }

    /** test preMultiply by vector */
    @Test
    void testPremultiplyVector() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        UnitTestUtils.customAssertEquals(m.preMultiply(testVector), preMultTest);
        UnitTestUtils.customAssertEquals(m.preMultiply(new ArrayFieldVector<Fraction>(testVector).toArray()),
                                         preMultTest);
        m = new BlockFieldMatrix<Fraction>(bigSingular);
        try {
            m.preMultiply(testVector);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testPremultiply() {
        FieldMatrix<Fraction> m3 = new BlockFieldMatrix<Fraction>(d3);
        FieldMatrix<Fraction> m4 = new BlockFieldMatrix<Fraction>(d4);
        FieldMatrix<Fraction> m5 = new BlockFieldMatrix<Fraction>(d5);
        UnitTestUtils.customAssertEquals(m4.preMultiply(m3), m5);

        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> mInv = new BlockFieldMatrix<Fraction>(testDataInv);
        BlockFieldMatrix<Fraction> identity = new BlockFieldMatrix<Fraction>(id);
        UnitTestUtils.customAssertEquals(m.preMultiply(mInv), identity);
        UnitTestUtils.customAssertEquals(mInv.preMultiply(m), identity);
        UnitTestUtils.customAssertEquals(m.preMultiply(identity), m);
        UnitTestUtils.customAssertEquals(identity.preMultiply(mInv), mInv);
        try {
            m.preMultiply(new BlockFieldMatrix<Fraction>(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testGetVectors() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        UnitTestUtils.customAssertEquals(m.getRow(0), testDataRow1);
        UnitTestUtils.customAssertEquals(m.getColumn(2), testDataCol3);
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        assertEquals(m.getEntry(0,1),new Fraction(2));
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
        Fraction[][] matrixData = {
                {new Fraction(1),new Fraction(2),new Fraction(3)},
                {new Fraction(2),new Fraction(5),new Fraction(3)}
        };
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(matrixData);
        // One more with three rows, two columns
        Fraction[][] matrixData2 = {
                {new Fraction(1),new Fraction(2)},
                {new Fraction(2),new Fraction(5)},
                {new Fraction(1), new Fraction(7)}
        };
        FieldMatrix<Fraction> n = new BlockFieldMatrix<Fraction>(matrixData2);
        // Now multiply m by n
        FieldMatrix<Fraction> p = m.multiply(n);
        assertEquals(2, p.getRowDimension());
        assertEquals(2, p.getColumnDimension());
        // Invert p
        FieldMatrix<Fraction> pInverse = new FieldLUDecomposition<Fraction>(p).getSolver().getInverse();
        assertEquals(2, pInverse.getRowDimension());
        assertEquals(2, pInverse.getColumnDimension());

        // Solve example
        Fraction[][] coefficientsData = {
                {new Fraction(2), new Fraction(3), new Fraction(-2)},
                {new Fraction(-1), new Fraction(7), new Fraction(6)},
                {new Fraction(4), new Fraction(-3), new Fraction(-5)}
        };
        FieldMatrix<Fraction> coefficients = new BlockFieldMatrix<Fraction>(coefficientsData);
        Fraction[] constants = {
            new Fraction(1), new Fraction(-2), new Fraction(1)
        };
        Fraction[] solution;
        solution = new FieldLUDecomposition<Fraction>(coefficients)
            .getSolver()
            .solve(new ArrayFieldVector<Fraction>(constants, false)).toArray();
        assertEquals(new Fraction(2).multiply(solution[0]).
                     add(new Fraction(3).multiply(solution[1])).
                     subtract(new Fraction(2).multiply(solution[2])),
                     constants[0]);
        assertEquals(new Fraction(-1).multiply(solution[0]).
                     add(new Fraction(7).multiply(solution[1])).
                     add(new Fraction(6).multiply(solution[2])),
                     constants[1]);
        assertEquals(new Fraction(4).multiply(solution[0]).
                     subtract(new Fraction(3).multiply(solution[1])).
                     subtract(new Fraction(5).multiply(solution[2])),
                     constants[2]);

    }

    // test submatrix accessors
    @Test
    void testGetSubMatrix() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
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

    private void checkGetSubMatrix(FieldMatrix<Fraction> m, Fraction[][] reference,
                                   int startRow, int endRow, int startColumn, int endColumn) {
        try {
            FieldMatrix<Fraction> sub = m.getSubMatrix(startRow, endRow, startColumn, endColumn);
            if (reference != null) {
                assertEquals(new BlockFieldMatrix<Fraction>(reference), sub);
            } else {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentException"
                     + " or MathIllegalArgumentException or MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (reference != null) {
                throw e;
            }
        }
    }

    private void checkGetSubMatrix(FieldMatrix<Fraction> m, Fraction[][] reference,
                                   int[] selectedRows, int[] selectedColumns) {
        try {
            FieldMatrix<Fraction> sub = m.getSubMatrix(selectedRows, selectedColumns);
            if (reference != null) {
                assertEquals(new BlockFieldMatrix<Fraction>(reference), sub);
            } else {
                fail("Expecting MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (reference != null) {
                throw e;
            }
        }
    }

    @Test
    void testGetSetMatrixLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        FieldMatrix<Fraction> sub =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n - 4, n - 4).scalarAdd(new Fraction(1));

        m.setSubMatrix(sub.getData(), 2, 2);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if ((i < 2) || (i > n - 3) || (j < 2) || (j > n - 3)) {
                    assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        assertEquals(sub, m.getSubMatrix(2, n - 3, 2, n - 3));
    }

    @Test
    void testCopySubMatrix() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
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
        checkCopy(m, null, new int[] {}, new int[] { 0 });
        checkCopy(m, null, new int[] { 0 }, new int[] { 4 });
    }

    private void checkCopy(FieldMatrix<Fraction> m, Fraction[][] reference,
                           int startRow, int endRow, int startColumn, int endColumn) {
        try {
            Fraction[][] sub = (reference == null) ?
                             new Fraction[1][1] :
                             new Fraction[reference.length][reference[0].length];
            m.copySubMatrix(startRow, endRow, startColumn, endColumn, sub);
            if (reference != null) {
                assertEquals(new BlockFieldMatrix<Fraction>(reference), new BlockFieldMatrix<Fraction>(sub));
            } else {
                fail("Expecting MathIllegalArgumentException or MathIllegalArgumentException or MathIllegalArgumentException");
            }
        } catch (MathIllegalArgumentException e) {
            if (reference != null) {
                throw e;
            }
        }
    }

    private void checkCopy(FieldMatrix<Fraction> m, Fraction[][] reference,
                           int[] selectedRows, int[] selectedColumns) {
        try {
            Fraction[][] sub = (reference == null) ?
                    new Fraction[1][1] :
                    new Fraction[reference.length][reference[0].length];
            m.copySubMatrix(selectedRows, selectedColumns, sub);
            if (reference != null) {
                assertEquals(new BlockFieldMatrix<Fraction>(reference), new BlockFieldMatrix<Fraction>(sub));
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
        FieldMatrix<Fraction> m     = new BlockFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mRow0 = new BlockFieldMatrix<Fraction>(subRow0);
        FieldMatrix<Fraction> mRow3 = new BlockFieldMatrix<Fraction>(subRow3);
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mRow3 = new BlockFieldMatrix<Fraction>(subRow3);
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
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        FieldMatrix<Fraction> sub =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), 1, n).scalarAdd(new Fraction(1));

        m.setRowMatrix(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        assertEquals(sub, m.getRowMatrix(2));

    }

    @Test
    void testGetColumnMatrix() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mColumn1 = new BlockFieldMatrix<Fraction>(subColumn1);
        FieldMatrix<Fraction> mColumn3 = new BlockFieldMatrix<Fraction>(subColumn3);
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mColumn3 = new BlockFieldMatrix<Fraction>(subColumn3);
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
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        FieldMatrix<Fraction> sub =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, 1).scalarAdd(new Fraction(1));

        m.setColumnMatrix(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        assertEquals(sub, m.getColumnMatrix(2));

    }

    @Test
    void testGetRowVector() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mRow0 = new ArrayFieldVector<Fraction>(subRow0[0]);
        FieldVector<Fraction> mRow3 = new ArrayFieldVector<Fraction>(subRow3[0]);
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mRow3 = new ArrayFieldVector<Fraction>(subRow3[0]);
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
            m.setRowVector(0, new ArrayFieldVector<Fraction>(FractionField.getInstance(), 5));
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetSetRowVectorLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        FieldVector<Fraction> sub = new ArrayFieldVector<Fraction>(n, new Fraction(1));

        m.setRowVector(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        assertEquals(sub, m.getRowVector(2));

    }

    @Test
    void testGetColumnVector() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mColumn1 = columnToVector(subColumn1);
        FieldVector<Fraction> mColumn3 = columnToVector(subColumn3);
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mColumn3 = columnToVector(subColumn3);
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
            m.setColumnVector(0, new ArrayFieldVector<Fraction>(FractionField.getInstance(), 5));
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetSetColumnVectorLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        FieldVector<Fraction> sub = new ArrayFieldVector<Fraction>(n, new Fraction(1));

        m.setColumnVector(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        assertEquals(sub, m.getColumnVector(2));

    }

    private FieldVector<Fraction> columnToVector(Fraction[][] column) {
        Fraction[] data = new Fraction[column.length];
        for (int i = 0; i < data.length; ++i) {
            data[i] = column[i][0];
        }
        return new ArrayFieldVector<Fraction>(data, false);
    }

    @Test
    void testGetRow() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
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
            m.setRow(0, new Fraction[5]);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetSetRowLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        Fraction[] sub = new Fraction[n];
        Arrays.fill(sub, new Fraction(1));

        m.setRow(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        checkArrays(sub, m.getRow(2));

    }

    @Test
    void testGetColumn() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        Fraction[] mColumn1 = columnToArray(subColumn1);
        Fraction[] mColumn3 = columnToArray(subColumn3);
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        Fraction[] mColumn3 = columnToArray(subColumn3);
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
            m.setColumn(0, new Fraction[5]);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testGetSetColumnLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        Fraction[] sub = new Fraction[n];
        Arrays.fill(sub, new Fraction(1));

        m.setColumn(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        checkArrays(sub, m.getColumn(2));

    }

    private Fraction[] columnToArray(Fraction[][] column) {
        Fraction[] data = new Fraction[column.length];
        for (int i = 0; i < data.length; ++i) {
            data[i] = column[i][0];
        }
        return data;
    }

    private void checkArrays(Fraction[] expected, Fraction[] actual) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(expected[i], actual[i]);
        }
    }

    @Test
    void testEqualsAndHashCode() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> m1 = (BlockFieldMatrix<Fraction>) m.copy();
        BlockFieldMatrix<Fraction> mt = (BlockFieldMatrix<Fraction>) m.transpose();
        assertTrue(m.hashCode() != mt.hashCode());
        assertEquals(m.hashCode(), m1.hashCode());
        assertEquals(m, m);
        assertEquals(m, m1);
        assertNotEquals(null, m);
        assertNotEquals(m, mt);
        assertNotEquals(m, new BlockFieldMatrix<Fraction>(bigSingular));
    }

    @Test
    void testToString() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        assertEquals("BlockFieldMatrix{{1,2,3},{2,5,3},{1,0,8}}", m.toString());
    }

    @Test
    void testSetSubMatrix() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        m.setSubMatrix(detData2,1,1);
        FieldMatrix<Fraction> expected = new BlockFieldMatrix<Fraction>
            (new Fraction[][] {{new Fraction(1),new Fraction(2),new Fraction(3)},{new Fraction(2),new Fraction(1),new Fraction(3)},{new Fraction(1),new Fraction(2),new Fraction(4)}});
        assertEquals(expected, m);

        m.setSubMatrix(detData2,0,0);
        expected = new BlockFieldMatrix<Fraction>
            (new Fraction[][] {{new Fraction(1),new Fraction(3),new Fraction(3)},{new Fraction(2),new Fraction(4),new Fraction(3)},{new Fraction(1),new Fraction(2),new Fraction(4)}});
        assertEquals(expected, m);

        m.setSubMatrix(testDataPlus2,0,0);
        expected = new BlockFieldMatrix<Fraction>
            (new Fraction[][] {{new Fraction(3),new Fraction(4),new Fraction(5)},{new Fraction(4),new Fraction(7),new Fraction(5)},{new Fraction(3),new Fraction(2),new Fraction(10)}});
        assertEquals(expected, m);

        // javadoc example
        BlockFieldMatrix<Fraction> matrix =
            new BlockFieldMatrix<Fraction>(new Fraction[][] {
                    {new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4)},
                    {new Fraction(5), new Fraction(6), new Fraction(7), new Fraction(8)},
                    {new Fraction(9), new Fraction(0), new Fraction(1) , new Fraction(2)}
            });
        matrix.setSubMatrix(new Fraction[][] {
                {new Fraction(3), new Fraction(4)},
                {new Fraction(5), new Fraction(6)}
        }, 1, 1);
        expected =
            new BlockFieldMatrix<Fraction>(new Fraction[][] {
                    {new Fraction(1), new Fraction(2), new Fraction(3),new Fraction(4)},
                    {new Fraction(5), new Fraction(3), new Fraction(4), new Fraction(8)},
                    {new Fraction(9), new Fraction(5) ,new Fraction(6), new Fraction(2)}
            });
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
            m.setSubMatrix(new Fraction[][] {{new Fraction(1)}, {new Fraction(2), new Fraction(3)}}, 0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        // empty
        try {
            m.setSubMatrix(new Fraction[][] {{}}, 0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
    }

    @Test
    void testWalk() {
        int rows    = 150;
        int columns = 75;

        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInRowOrder(new SetVisitor());
        GetVisitor getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInRowOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            assertEquals(new Fraction(0), m.getEntry(i, 0));
            assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            assertEquals(new Fraction(0), m.getEntry(0, j));
            assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInColumnOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInColumnOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            assertEquals(new Fraction(0), m.getEntry(i, 0));
            assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            assertEquals(new Fraction(0), m.getEntry(0, j));
            assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            assertEquals(new Fraction(0), m.getEntry(i, 0));
            assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            assertEquals(new Fraction(0), m.getEntry(0, j));
            assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            assertEquals(new Fraction(0), m.getEntry(i, 0));
            assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            assertEquals(new Fraction(0), m.getEntry(0, j));
            assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

    }

    @Test
    void testSerial()  {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        assertEquals(m,UnitTestUtils.serializeAndRecover(m));
    }

    private static class SetVisitor extends DefaultFieldMatrixChangingVisitor<Fraction> {
        public SetVisitor() {
            super(Fraction.ZERO);
        }
        @Override
        public Fraction visit(int i, int j, Fraction value) {
            return new Fraction(i * 11 + j, 11);
        }
    }

    private static class GetVisitor extends DefaultFieldMatrixPreservingVisitor<Fraction> {
        private int count;
        public GetVisitor() {
            super(Fraction.ZERO);
            count = 0;
        }
        @Override
        public void visit(int i, int j, Fraction value) {
            ++count;
            assertEquals(new Fraction(i * 11 + j, 11), value);
        }
        public int getCount() {
            return count;
        }
    }

    private BlockFieldMatrix<Fraction> createRandomMatrix(Random r, int rows, int columns) {
        BlockFieldMatrix<Fraction> m =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                int p = r.nextInt(20) - 10;
                int q = r.nextInt(20) - 10;
                if (q == 0) {
                    q = 1;
                }
                m.setEntry(i, j, new Fraction(p, q));
            }
        }
        return m;
    }
}

