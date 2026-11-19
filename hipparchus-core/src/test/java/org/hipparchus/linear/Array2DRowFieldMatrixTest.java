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
import org.hipparchus.fraction.Fraction;
import org.hipparchus.fraction.FractionField;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link Array2DRowFieldMatrix} class.
 *
 */

final class Array2DRowFieldMatrixTest {

    // 3 x 3 identity matrix
    protected Fraction[][] id = { {new Fraction(1),new Fraction(0),new Fraction(0)}, {new Fraction(0),new Fraction(1),new Fraction(0)}, {new Fraction(0),new Fraction(0),new Fraction(1)} };

    // Test data for group operations
    protected Fraction[][] testData = { {new Fraction(1),new Fraction(2),new Fraction(3)}, {new Fraction(2),new Fraction(5),new Fraction(3)}, {new Fraction(1),new Fraction(0),new Fraction(8)} };
    protected Fraction[][] testDataLU = {{new Fraction(2), new Fraction(5), new Fraction(3)}, {new Fraction(1, 2), new Fraction(-5, 2), new Fraction(13, 2)}, {new Fraction(1, 2), new Fraction(1, 5), new Fraction(1, 5)}};
    protected Fraction[][] testDataPlus2 = { {new Fraction(3),new Fraction(4),new Fraction(5)}, {new Fraction(4),new Fraction(7),new Fraction(5)}, {new Fraction(3),new Fraction(2),new Fraction(10)} };
    protected Fraction[][] testDataMinus = { {new Fraction(-1),new Fraction(-2),new Fraction(-3)}, {new Fraction(-2),new Fraction(-5),new Fraction(-3)},
       {new Fraction(-1),new Fraction(0),new Fraction(-8)} };
    protected Fraction[] testDataRow1 = {new Fraction(1),new Fraction(2),new Fraction(3)};
    protected Fraction[] testDataCol3 = {new Fraction(3),new Fraction(3),new Fraction(8)};
    protected Fraction[][] testDataInv =
        { {new Fraction(-40),new Fraction(16),new Fraction(9)}, {new Fraction(13),new Fraction(-5),new Fraction(-3)}, {new Fraction(5),new Fraction(-2),new Fraction(-1)} };
    protected Fraction[] preMultTest = {new Fraction(8),new Fraction(12),new Fraction(33)};
    protected Fraction[][] testData2 ={ {new Fraction(1),new Fraction(2),new Fraction(3)}, {new Fraction(2),new Fraction(5),new Fraction(3)}};
    protected Fraction[][] testData2T = { {new Fraction(1),new Fraction(2)}, {new Fraction(2),new Fraction(5)}, {new Fraction(3),new Fraction(3)}};
    protected Fraction[][] testDataPlusInv =
        { {new Fraction(-39),new Fraction(18),new Fraction(12)}, {new Fraction(15),new Fraction(0),new Fraction(0)}, {new Fraction(6),new Fraction(-2),new Fraction(7)} };

    // lu decomposition tests
    protected Fraction[][] luData = { {new Fraction(2),new Fraction(3),new Fraction(3)}, {new Fraction(0),new Fraction(5),new Fraction(7)}, {new Fraction(6),new Fraction(9),new Fraction(8)} };
    protected Fraction[][] luDataLUDecomposition = { {new Fraction(6),new Fraction(9),new Fraction(8)}, {new Fraction(0),new Fraction(5),new Fraction(7)},
            {new Fraction(1, 3),new Fraction(0),new Fraction(1, 3)} };

    // singular matrices
    protected Fraction[][] singular = { {new Fraction(2),new Fraction(3)}, {new Fraction(2),new Fraction(3)} };
    protected Fraction[][] bigSingular = {{new Fraction(1),new Fraction(2),new Fraction(3),new Fraction(4)}, {new Fraction(2),new Fraction(5),new Fraction(3),new Fraction(4)},
        {new Fraction(7),new Fraction(3),new Fraction(256),new Fraction(1930)}, {new Fraction(3),new Fraction(7),new Fraction(6),new Fraction(8)}}; // 4th row = 1st + 2nd
    protected Fraction[][] detData = { {new Fraction(1),new Fraction(2),new Fraction(3)}, {new Fraction(4),new Fraction(5),new Fraction(6)}, {new Fraction(7),new Fraction(8),new Fraction(10)} };
    protected Fraction[][] detData2 = { {new Fraction(1), new Fraction(3)}, {new Fraction(2), new Fraction(4)}};

    // vectors
    protected Fraction[] testVector = {new Fraction(1),new Fraction(2),new Fraction(3)};
    protected Fraction[] testVector2 = {new Fraction(1),new Fraction(2),new Fraction(3),new Fraction(4)};

    // submatrix accessor tests
    protected Fraction[][] subTestData = {{new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4)}, {new Fraction(3, 2), new Fraction(5, 2), new Fraction(7, 2), new Fraction(9, 2)},
            {new Fraction(2), new Fraction(4), new Fraction(6), new Fraction(8)}, {new Fraction(4), new Fraction(5), new Fraction(6), new Fraction(7)}};
    // array selections
    protected Fraction[][] subRows02Cols13 = { {new Fraction(2), new Fraction(4)}, {new Fraction(4), new Fraction(8)}};
    protected Fraction[][] subRows03Cols12 = { {new Fraction(2), new Fraction(3)}, {new Fraction(5), new Fraction(6)}};
    protected Fraction[][] subRows03Cols123 = { {new Fraction(2), new Fraction(3), new Fraction(4)} , {new Fraction(5), new Fraction(6), new Fraction(7)}};
    // effective permutations
    protected Fraction[][] subRows20Cols123 = { {new Fraction(4), new Fraction(6), new Fraction(8)} , {new Fraction(2), new Fraction(3), new Fraction(4)}};
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
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(testData2);
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
        Array2DRowFieldMatrix<Fraction> m1 = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(m1.getData());
        assertEquals(m2,m1);
        Array2DRowFieldMatrix<Fraction> m3 = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m4 = new Array2DRowFieldMatrix<Fraction>(m3.getData(), false);
        assertEquals(m4,m3);
    }

    /** test add */
    @Test
    void testAdd() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> mInv = new Array2DRowFieldMatrix<Fraction>(testDataInv);
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
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(testData2);
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
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(testDataInv);
        UnitTestUtils.customAssertEquals(m.subtract(m2), m2.scalarMultiply(new Fraction(-1)).add(m));
        try {
            m.subtract(new Array2DRowFieldMatrix<Fraction>(testData2));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    /** test multiply */
    @Test
    void testMultiply() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> mInv = new Array2DRowFieldMatrix<Fraction>(testDataInv);
        Array2DRowFieldMatrix<Fraction> identity = new Array2DRowFieldMatrix<Fraction>(id);
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(testData2);
        UnitTestUtils.customAssertEquals(m.multiply(mInv), identity);
        UnitTestUtils.customAssertEquals(mInv.multiply(m), identity);
        UnitTestUtils.customAssertEquals(m.multiply(identity), m);
        UnitTestUtils.customAssertEquals(identity.multiply(mInv), mInv);
        UnitTestUtils.customAssertEquals(m2.multiply(identity), m2);
        try {
            m.multiply(new Array2DRowFieldMatrix<Fraction>(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    //Additional Test for Array2DRowFieldMatrix<Fraction>Test.testMultiply

    private final Fraction[][] d3 = new Fraction[][] {{new Fraction(1),new Fraction(2),new Fraction(3),new Fraction(4)},{new Fraction(5),new Fraction(6),new Fraction(7),new Fraction(8)}};
    private final Fraction[][] d4 = new Fraction[][] {{new Fraction(1)},{new Fraction(2)},{new Fraction(3)},{new Fraction(4)}};
    private final Fraction[][] d5 = new Fraction[][] {{new Fraction(30)},{new Fraction(70)}};

    @Test
    void testMultiply2() {
       FieldMatrix<Fraction> m3 = new Array2DRowFieldMatrix<Fraction>(d3);
       FieldMatrix<Fraction> m4 = new Array2DRowFieldMatrix<Fraction>(d4);
       FieldMatrix<Fraction> m5 = new Array2DRowFieldMatrix<Fraction>(d5);
       UnitTestUtils.customAssertEquals(m3.multiply(m4), m5);
    }

    @Test
    void testMultiplyTransposedArray2DRowRealMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0xdeff3d383a112763l);
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
                final Array2DRowFieldMatrix<Binary64> a = new Array2DRowFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
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
    void testMultiplyTransposedBlockFieldMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x463e54fb50b900fel);
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
                final Array2DRowFieldMatrix<Binary64> a = new Array2DRowFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
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
    void testMultiplyTransposedWrongDimensions() {
        try {
            new Array2DRowFieldMatrix<Binary64>(Binary64Field.getInstance(), 2, 3).
            multiplyTransposed(new Array2DRowFieldMatrix<Binary64>(Binary64Field.getInstance(), 3, 2));
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
                final Array2DRowFieldMatrix<Binary64> a = new Array2DRowFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
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
    void testTransposeMultiplyBlockFieldMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x463e54fb50b900fel);
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
                final Array2DRowFieldMatrix<Binary64> a = new Array2DRowFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
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
    void testTransposeMultiplyWrongDimensions() {
        try {
            new Array2DRowFieldMatrix<Binary64>(Binary64Field.getInstance(), 2, 3).
            transposeMultiply(new Array2DRowFieldMatrix<Binary64>(Binary64Field.getInstance(), 3, 2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(2, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(3, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testPower() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        FieldMatrix<Fraction> mInv = new Array2DRowFieldMatrix<Fraction>(testDataInv);
        FieldMatrix<Fraction> mPlusInv = new Array2DRowFieldMatrix<Fraction>(testDataPlusInv);
        FieldMatrix<Fraction> identity = new Array2DRowFieldMatrix<Fraction>(id);

        UnitTestUtils.customAssertEquals(m.power(0), identity);
        UnitTestUtils.customAssertEquals(mInv.power(0), identity);
        UnitTestUtils.customAssertEquals(mPlusInv.power(0), identity);

        UnitTestUtils.customAssertEquals(m.power(1), m);
        UnitTestUtils.customAssertEquals(mInv.power(1), mInv);
        UnitTestUtils.customAssertEquals(mPlusInv.power(1), mPlusInv);

        FieldMatrix<Fraction> C1 = m.copy();
        FieldMatrix<Fraction> C2 = mInv.copy();
        FieldMatrix<Fraction> C3 = mPlusInv.copy();

        // stop at 5 to avoid overflow
        for (int i = 2; i <= 5; ++i) {
            C1 = C1.multiply(m);
            C2 = C2.multiply(mInv);
            C3 = C3.multiply(mPlusInv);

            UnitTestUtils.customAssertEquals(m.power(i), C1);
            UnitTestUtils.customAssertEquals(mInv.power(i), C2);
            UnitTestUtils.customAssertEquals(mPlusInv.power(i), C3);
        }

        try {
            FieldMatrix<Fraction> mNotSquare = new Array2DRowFieldMatrix<Fraction>(testData2T);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(id);
        assertEquals(new Fraction(3),m.getTrace(),"identity trace");
        m = new Array2DRowFieldMatrix<Fraction>(testData2);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        UnitTestUtils.customAssertEquals(new Array2DRowFieldMatrix<Fraction>(testDataPlus2), m.scalarAdd(new Fraction(2)));
    }

    /** test operate */
    @Test
    void testOperate() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(id);
        UnitTestUtils.customAssertEquals(testVector, m.operate(testVector));
        UnitTestUtils.customAssertEquals(testVector, m.operate(new ArrayFieldVector<Fraction>(testVector)).toArray());
        m = new Array2DRowFieldMatrix<Fraction>(bigSingular);
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
        FieldMatrix<Fraction> a = new Array2DRowFieldMatrix<Fraction>(new Fraction[][] {
                { new Fraction(1), new Fraction(2) }, { new Fraction(3), new Fraction(4) }, { new Fraction(5), new Fraction(6) }
        }, false);
        Fraction[] b = a.operate(new Fraction[] { new Fraction(1), new Fraction(1) });
        assertEquals(a.getRowDimension(), b.length);
        assertEquals( new Fraction(3), b[0]);
        assertEquals( new Fraction(7), b[1]);
        assertEquals(new Fraction(11), b[2]);
    }

    /** test transpose */
    @Test
    void testTranspose() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        FieldMatrix<Fraction> mIT = new FieldLUDecomposition<Fraction>(m).getSolver().getInverse().transpose();
        FieldMatrix<Fraction> mTI = new FieldLUDecomposition<Fraction>(m.transpose()).getSolver().getInverse();
        UnitTestUtils.customAssertEquals(mIT, mTI);
        m = new Array2DRowFieldMatrix<Fraction>(testData2);
        FieldMatrix<Fraction> mt = new Array2DRowFieldMatrix<Fraction>(testData2T);
        UnitTestUtils.customAssertEquals(mt, m.transpose());
    }

    /** test preMultiply by vector */
    @Test
    void testPremultiplyVector() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        UnitTestUtils.customAssertEquals(m.preMultiply(testVector), preMultTest);
        UnitTestUtils.customAssertEquals(m.preMultiply(new ArrayFieldVector<Fraction>(testVector).toArray()),
                                         preMultTest);
        m = new Array2DRowFieldMatrix<Fraction>(bigSingular);
        try {
            m.preMultiply(testVector);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testPremultiply() {
        FieldMatrix<Fraction> m3 = new Array2DRowFieldMatrix<Fraction>(d3);
        FieldMatrix<Fraction> m4 = new Array2DRowFieldMatrix<Fraction>(d4);
        FieldMatrix<Fraction> m5 = new Array2DRowFieldMatrix<Fraction>(d5);
        UnitTestUtils.customAssertEquals(m4.preMultiply(m3), m5);

        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> mInv = new Array2DRowFieldMatrix<Fraction>(testDataInv);
        Array2DRowFieldMatrix<Fraction> identity = new Array2DRowFieldMatrix<Fraction>(id);
        UnitTestUtils.customAssertEquals(m.preMultiply(mInv), identity);
        UnitTestUtils.customAssertEquals(mInv.preMultiply(m), identity);
        UnitTestUtils.customAssertEquals(m.preMultiply(identity), m);
        UnitTestUtils.customAssertEquals(identity.preMultiply(mInv), mInv);
        try {
            m.preMultiply(new Array2DRowFieldMatrix<Fraction>(bigSingular));
            fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    void testGetVectors() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        assertEquals(m.getEntry(0,1), new Fraction(2), "get entry");
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(matrixData);
        // One more with three rows, two columns
        Fraction[][] matrixData2 = {
                {new Fraction(1),new Fraction(2)},
                {new Fraction(2),new Fraction(5)},
                {new Fraction(1), new Fraction(7)}
        };
        FieldMatrix<Fraction> n = new Array2DRowFieldMatrix<Fraction>(matrixData2);
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
        FieldMatrix<Fraction> coefficients = new Array2DRowFieldMatrix<Fraction>(coefficientsData);
        Fraction[] constants = {
            new Fraction(1), new Fraction(-2), new Fraction(1)
        };
        Fraction[] solution;
        solution = new FieldLUDecomposition<Fraction>(coefficients)
            .getSolver()
            .solve(new ArrayFieldVector<Fraction>(constants, false)).toArray();
        assertEquals(new Fraction(2).multiply(solution[0]).
                     add(new Fraction(3).multiply(solution[1])).
                     subtract(new Fraction(2).multiply(solution[2])), constants[0]);
        assertEquals(new Fraction(-1).multiply(solution[0]).
                     add(new Fraction(7).multiply(solution[1])).
                     add(new Fraction(6).multiply(solution[2])), constants[1]);
        assertEquals(new Fraction(4).multiply(solution[0]).
                     subtract(new Fraction(3).multiply(solution[1])).
                     subtract(new Fraction(5).multiply(solution[2])), constants[2]);

    }

    // test submatrix accessors
    @Test
    void testGetSubMatrix() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
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
                assertEquals(new Array2DRowFieldMatrix<Fraction>(reference), sub);
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
                assertEquals(new Array2DRowFieldMatrix<Fraction>(reference), sub);
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

    @Test
    void testCopySubMatrix() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
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

    private void checkCopy(FieldMatrix<Fraction> m, Fraction[][] reference,
                           int startRow, int endRow, int startColumn, int endColumn) {
        try {
            Fraction[][] sub = (reference == null) ?
                             new Fraction[1][1] :
                             new Fraction[reference.length][reference[0].length];
            m.copySubMatrix(startRow, endRow, startColumn, endColumn, sub);
            if (reference != null) {
                assertEquals(new Array2DRowFieldMatrix<Fraction>(reference), new Array2DRowFieldMatrix<Fraction>(sub));
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
                assertEquals(new Array2DRowFieldMatrix<Fraction>(reference), new Array2DRowFieldMatrix<Fraction>(sub));
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mRow0 = new Array2DRowFieldMatrix<Fraction>(subRow0);
        FieldMatrix<Fraction> mRow3 = new Array2DRowFieldMatrix<Fraction>(subRow3);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mRow3 = new Array2DRowFieldMatrix<Fraction>(subRow3);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mColumn1 = new Array2DRowFieldMatrix<Fraction>(subColumn1);
        FieldMatrix<Fraction> mColumn3 = new Array2DRowFieldMatrix<Fraction>(subColumn3);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mColumn3 = new Array2DRowFieldMatrix<Fraction>(subColumn3);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mRow0 = new ArrayFieldVector<Fraction>(subRow0[0]);
        FieldVector<Fraction> mRow3 = new ArrayFieldVector<Fraction>(subRow3[0]);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
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
    void testGetColumnVector() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mColumn1 = columnToVector(subColumn1);
        FieldVector<Fraction> mColumn3 = columnToVector(subColumn3);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
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

    private FieldVector<Fraction> columnToVector(Fraction[][] column) {
        Fraction[] data = new Fraction[column.length];
        for (int i = 0; i < data.length; ++i) {
            data[i] = column[i][0];
        }
        return new ArrayFieldVector<Fraction>(data, false);
    }

    @Test
    void testGetRow() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        assertNotSame(subRow3[0][0], m.getRow(0)[0]);
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
    void testGetColumn() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        Fraction[] mColumn3 = columnToArray(subColumn3);
        assertNotSame(mColumn3[0], m.getColumn(1)[0]);
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
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m1 = (Array2DRowFieldMatrix<Fraction>) m.copy();
        Array2DRowFieldMatrix<Fraction> mt = (Array2DRowFieldMatrix<Fraction>) m.transpose();
        assertTrue(m.hashCode() != mt.hashCode());
        assertEquals(m.hashCode(), m1.hashCode());
        assertEquals(m, m);
        assertEquals(m, m1);
        assertNotEquals(null, m);
        assertNotEquals(m, mt);
        assertNotEquals(m, new Array2DRowFieldMatrix<Fraction>(bigSingular));
    }

    @Test
    void testToString() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        assertEquals("Array2DRowFieldMatrix{{1,2,3},{2,5,3},{1,0,8}}", m.toString());
        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance());
        assertEquals("Array2DRowFieldMatrix{}", m.toString());
    }

    @Test
    void testSetSubMatrix() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        m.setSubMatrix(detData2,1,1);
        FieldMatrix<Fraction> expected = new Array2DRowFieldMatrix<Fraction>
            (new Fraction[][] {
                    {new Fraction(1),new Fraction(2),new Fraction(3)},
                    {new Fraction(2),new Fraction(1),new Fraction(3)},
                    {new Fraction(1),new Fraction(2),new Fraction(4)}
             });
        assertEquals(expected, m);

        m.setSubMatrix(detData2,0,0);
        expected = new Array2DRowFieldMatrix<Fraction>
            (new Fraction[][] {
                    {new Fraction(1),new Fraction(3),new Fraction(3)},
                    {new Fraction(2),new Fraction(4),new Fraction(3)},
                    {new Fraction(1),new Fraction(2),new Fraction(4)}
             });
        assertEquals(expected, m);

        m.setSubMatrix(testDataPlus2,0,0);
        expected = new Array2DRowFieldMatrix<Fraction>
            (new Fraction[][] {
                    {new Fraction(3),new Fraction(4),new Fraction(5)},
                    {new Fraction(4),new Fraction(7),new Fraction(5)},
                    {new Fraction(3),new Fraction(2),new Fraction(10)}
             });
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
            m.setSubMatrix(null, 1, 1);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException e) {
            // expected
        }
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance());
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

        FieldMatrix<Fraction> m =
            new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInRowOrder(new SetVisitor());
        GetVisitor getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
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

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInColumnOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
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

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
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

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor);
        assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
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
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        assertEquals(m,UnitTestUtils.serializeAndRecover(m));
    }

    private static class SetVisitor extends DefaultFieldMatrixChangingVisitor<Fraction> {
        public SetVisitor() {
            super(Fraction.ZERO);
        }
        @Override
        public Fraction visit(int i, int j, Fraction value) {
            return new Fraction(i * 1024 + j, 1024);
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
            assertEquals(new Fraction(i * 1024 + j, 1024), value);
        }
        public int getCount() {
            return count;
        }
    }

    //--------------- -----------------Protected methods

    /** extracts the l  and u matrices from compact lu representation */
    protected void splitLU(FieldMatrix<Fraction> lu,
                           Fraction[][] lowerData,
                           Fraction[][] upperData) {
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
                    upperData[i][j] = Fraction.ZERO;
                } else if (i == j) {
                    lowerData[i][j] = Fraction.ONE;
                    upperData[i][j] = lu.getEntry(i, j);
                } else {
                    lowerData[i][j] = Fraction.ZERO;
                    upperData[i][j] = lu.getEntry(i, j);
                }
            }
        }
    }

    /** Returns the result of applying the given row permutation to the matrix */
    protected FieldMatrix<Fraction> permuteRows(FieldMatrix<Fraction> matrix, int[] permutation) {
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
        Fraction[][] out = new Fraction[m][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out[i][j] = matrix.getEntry(permutation[i], j);
            }
        }
        return new Array2DRowFieldMatrix<Fraction>(out);
    }
}
