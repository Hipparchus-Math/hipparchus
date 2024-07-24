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

import org.hipparchus.Field;
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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link SparseFieldMatrix} class.
 *
 */
public class SparseFieldMatrixTest {
    // 3 x 3 identity matrix
    protected Fraction[][] id = { {new Fraction(1), new Fraction(0), new Fraction(0) }, { new Fraction(0), new Fraction(1), new Fraction(0) }, { new Fraction(0), new Fraction(0), new Fraction(1) } };
    // Test data for group operations
    protected Fraction[][] testData = { { new Fraction(1), new Fraction(2), new Fraction(3) }, { new Fraction(2), new Fraction(5), new Fraction(3) },
            { new Fraction(1), new Fraction(0), new Fraction(8) } };
    protected Fraction[][] testDataLU = null;
    protected Fraction[][] testDataPlus2 = { { new Fraction(3), new Fraction(4), new Fraction(5) }, { new Fraction(4), new Fraction(7), new Fraction(5) },
            { new Fraction(3), new Fraction(2), new Fraction(10) } };
    protected Fraction[][] testDataMinus = { { new Fraction(-1), new Fraction(-2), new Fraction(-3) },
            { new Fraction(-2), new Fraction(-5), new Fraction(-3) }, { new Fraction(-1), new Fraction(0), new Fraction(-8) } };
    protected Fraction[] testDataRow1 = { new Fraction(1), new Fraction(2), new Fraction(3) };
    protected Fraction[] testDataCol3 = { new Fraction(3), new Fraction(3), new Fraction(8) };
    protected Fraction[][] testDataInv = { { new Fraction(-40), new Fraction(16), new Fraction(9) }, { new Fraction(13), new Fraction(-5), new Fraction(-3) },
            { new Fraction(5), new Fraction(-2), new Fraction(-1) } };
    protected Fraction[] preMultTest = { new Fraction(8), new Fraction(12), new Fraction(33) };
    protected Fraction[][] testData2 = { { new Fraction(1), new Fraction(2), new Fraction(3) }, { new Fraction(2), new Fraction(5), new Fraction(3) } };
    protected Fraction[][] testData2T = { { new Fraction(1), new Fraction(2) }, { new Fraction(2), new Fraction(5) }, { new Fraction(3), new Fraction(3) } };
    protected Fraction[][] testDataPlusInv = { { new Fraction(-39), new Fraction(18), new Fraction(12) },
            { new Fraction(15), new Fraction(0), new Fraction(0) }, { new Fraction(6), new Fraction(-2), new Fraction(7) } };

    // lu decomposition tests
    protected Fraction[][] luData = { { new Fraction(2), new Fraction(3), new Fraction(3) }, { new Fraction(0), new Fraction(5), new Fraction(7) }, { new Fraction(6), new Fraction(9), new Fraction(8) } };
    protected Fraction[][] luDataLUDecomposition = null;

    // singular matrices
    protected Fraction[][] singular = { { new Fraction(2), new Fraction(3) }, { new Fraction(2), new Fraction(3) } };
    protected Fraction[][] bigSingular = { { new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4) },
            { new Fraction(2), new Fraction(5), new Fraction(3), new Fraction(4) }, { new Fraction(7), new Fraction(3), new Fraction(256), new Fraction(1930) }, { new Fraction(3), new Fraction(7), new Fraction(6), new Fraction(8) } }; // 4th

    // row
    // =
    // 1st
    // +
    // 2nd
    protected Fraction[][] detData = { { new Fraction(1), new Fraction(2), new Fraction(3) }, { new Fraction(4), new Fraction(5), new Fraction(6) },
            { new Fraction(7), new Fraction(8), new Fraction(10) } };
    protected Fraction[][] detData2 = { { new Fraction(1), new Fraction(3) }, { new Fraction(2), new Fraction(4) } };

    // vectors
    protected Fraction[] testVector = { new Fraction(1), new Fraction(2), new Fraction(3) };
    protected Fraction[] testVector2 = { new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4) };

    // submatrix accessor tests
    protected Fraction[][] subTestData = null;

    // array selections
    protected Fraction[][] subRows02Cols13 = { {new Fraction(2), new Fraction(4) }, { new Fraction(4), new Fraction(8) } };
    protected Fraction[][] subRows03Cols12 = { { new Fraction(2), new Fraction(3) }, { new Fraction(5), new Fraction(6) } };
    protected Fraction[][] subRows03Cols123 = { { new Fraction(2), new Fraction(3), new Fraction(4) }, { new Fraction(5), new Fraction(6), new Fraction(7) } };

    // effective permutations
    protected Fraction[][] subRows20Cols123 = { { new Fraction(4), new Fraction(6), new Fraction(8) }, { new Fraction(2), new Fraction(3), new Fraction(4) } };
    protected Fraction[][] subRows31Cols31 = null;

    // contiguous ranges
    protected Fraction[][] subRows01Cols23 = null;
    protected Fraction[][] subRows23Cols00 = { { new Fraction(2) }, { new Fraction(4) } };
    protected Fraction[][] subRows00Cols33 = { { new Fraction(4) } };

    // row matrices
    protected Fraction[][] subRow0 = { { new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4) } };
    protected Fraction[][] subRow3 = { { new Fraction(4), new Fraction(5), new Fraction(6), new Fraction(7) } };

    // column matrices
    protected Fraction[][] subColumn1 = null;
    protected Fraction[][] subColumn3 = null;

    // tolerances
    protected double entryTolerance = 10E-16;
    protected double normTolerance = 10E-14;
    protected Field<Fraction> field = FractionField.getInstance();

    public SparseFieldMatrixTest() {
            testDataLU = new Fraction[][]{ { new Fraction(2), new Fraction(5), new Fraction(3) }, { new Fraction(.5d), new Fraction(-2.5d), new Fraction(6.5d) },
                    { new Fraction(0.5d), new Fraction(0.2d), new Fraction(.2d) } };
            luDataLUDecomposition = new Fraction[][]{ { new Fraction(6), new Fraction(9), new Fraction(8) },
                { new Fraction(0), new Fraction(5), new Fraction(7) }, { new Fraction(0.33333333333333), new Fraction(0), new Fraction(0.33333333333333) } };
            subTestData = new Fraction [][]{ { new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4) },
                    { new Fraction(1.5), new Fraction(2.5), new Fraction(3.5), new Fraction(4.5) }, { new Fraction(2), new Fraction(4), new Fraction(6), new Fraction(8) }, { new Fraction(4), new Fraction(5), new Fraction(6), new Fraction(7) } };
            subRows31Cols31 = new Fraction[][]{ { new Fraction(7), new Fraction(5) }, { new Fraction(4.5), new Fraction(2.5) } };
            subRows01Cols23 = new Fraction[][]{ { new Fraction(3), new Fraction(4) }, { new Fraction(3.5), new Fraction(4.5) } };
            subColumn1 = new Fraction [][]{ { new Fraction(2) }, { new Fraction(2.5) }, { new Fraction(4) }, { new Fraction(5) } };
            subColumn3 = new Fraction[][]{ { new Fraction(4) }, { new Fraction(4.5) }, { new Fraction(8) }, { new Fraction(7) } };
    }

    /** test dimensions */
    @Test
    void testDimensions() {
        SparseFieldMatrix<Fraction> m = createSparseMatrix(testData);
        SparseFieldMatrix<Fraction> m2 = createSparseMatrix(testData2);
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
        SparseFieldMatrix<Fraction> m1 = createSparseMatrix(testData);
        FieldMatrix<Fraction> m2 = m1.copy();
        assertEquals(m1.getClass(), m2.getClass());
        assertEquals((m2), m1);
        SparseFieldMatrix<Fraction> m3 = createSparseMatrix(testData);
        FieldMatrix<Fraction> m4 = m3.copy();
        assertEquals(m3.getClass(), m4.getClass());
        assertEquals((m4), m3);
    }

    /** test add */
    @Test
    void testAdd() {
        SparseFieldMatrix<Fraction> m = createSparseMatrix(testData);
        SparseFieldMatrix<Fraction> mInv = createSparseMatrix(testDataInv);
        SparseFieldMatrix<Fraction> mDataPlusInv = createSparseMatrix(testDataPlusInv);
        FieldMatrix<Fraction> mPlusMInv = m.add(mInv);
        for (int row = 0; row < m.getRowDimension(); row++) {
            for (int col = 0; col < m.getColumnDimension(); col++) {
                assertEquals(mDataPlusInv.getEntry(row, col).doubleValue(), mPlusMInv.getEntry(row, col).doubleValue(),
                    entryTolerance,
                    "sum entry entry");
            }
        }
    }

    /** test add failure */
    @Test
    void testAddFail() {
        SparseFieldMatrix<Fraction> m = createSparseMatrix(testData);
        SparseFieldMatrix<Fraction> m2 = createSparseMatrix(testData2);
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
        SparseFieldMatrix<Fraction> m = createSparseMatrix(testData);
        SparseFieldMatrix<Fraction> n = createSparseMatrix(testDataInv);
        customAssertClose("m-n = m + -n", m.subtract(n),
                          n.scalarMultiply(new Fraction(-1)).add(m), entryTolerance);
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
        SparseFieldMatrix<Fraction> m = createSparseMatrix(testData);
        SparseFieldMatrix<Fraction> mInv = createSparseMatrix(testDataInv);
        SparseFieldMatrix<Fraction> identity = createSparseMatrix(id);
        SparseFieldMatrix<Fraction> m2 = createSparseMatrix(testData2);
        customAssertClose("inverse multiply", m.multiply(mInv), identity,
                          entryTolerance);
        customAssertClose("inverse multiply", m.multiply(new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), testDataInv)), identity,
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

    private Fraction[][] d3 = new Fraction[][] { { new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4) }, { new Fraction(5), new Fraction(6), new Fraction(7), new Fraction(8) } };
    private Fraction[][] d4 = new Fraction[][] { { new Fraction(1) }, { new Fraction(2) }, { new Fraction(3) }, { new Fraction(4) } };
    private Fraction[][] d5 = new Fraction[][] { { new Fraction(30) }, { new Fraction(70) } };

    @Test
    void testMultiply2() {
        FieldMatrix<Fraction> m3 = createSparseMatrix(d3);
        FieldMatrix<Fraction> m4 = createSparseMatrix(d4);
        FieldMatrix<Fraction> m5 = createSparseMatrix(d5);
        customAssertClose("m3*m4=m5", m3.multiply(m4), m5, entryTolerance);
    }

    @Test
    void testMultiplyTransposedSparseFieldMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x5f31d5645cf821efl);
        final FieldMatrixChangingVisitor<Binary64> randomSetter = new DefaultFieldMatrixChangingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public Binary64 visit(final int row, final int column, final Binary64 value) {
                return new Binary64(randomGenerator.nextDouble());
            }
        };
        final FieldMatrixPreservingVisitor<Binary64> zeroChecker = new DefaultFieldMatrixPreservingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public void visit(final int row, final int column, final Binary64 value) {
                assertEquals(0.0, value.doubleValue(), 3.0e-14);
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final SparseFieldMatrix<Binary64> a = new SparseFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final SparseFieldMatrix<Binary64> b = new SparseFieldMatrix<>(Binary64Field.getInstance(), interm, cols);
                    b.walkInOptimizedOrder(randomSetter);
                   a.multiplyTransposed(b).subtract(a.multiply(b.transpose())).walkInOptimizedOrder(zeroChecker);
                }
            }
        }
    }

    @Test
    void testMultiplyTransposedWrongDimensions() {
        try {
            new SparseFieldMatrix<>(Binary64Field.getInstance(), 2, 3).
            multiplyTransposed(new SparseFieldMatrix<>(Binary64Field.getInstance(), 3, 2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testTransposeMultiplySparseFieldMatrix() {
        RandomGenerator randomGenerator = new Well1024a(0x5f31d5645cf821efl);
        final FieldMatrixChangingVisitor<Binary64> randomSetter = new DefaultFieldMatrixChangingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public Binary64 visit(final int row, final int column, final Binary64 value) {
                return new Binary64(randomGenerator.nextDouble());
            }
        };
        final FieldMatrixPreservingVisitor<Binary64> zeroChecker = new DefaultFieldMatrixPreservingVisitor<Binary64>(Binary64Field.getInstance().getZero()) {
            public void visit(final int row, final int column, final Binary64 value) {
                assertEquals(0.0, value.doubleValue(), 3.0e-14);
            }
        };
        for (int rows = 1; rows <= 64; rows += 7) {
            for (int cols = 1; cols <= 64; cols += 7) {
                final SparseFieldMatrix<Binary64> a = new SparseFieldMatrix<>(Binary64Field.getInstance(), rows, cols);
                a.walkInOptimizedOrder(randomSetter);
                for (int interm = 1; interm <= 64; interm += 7) {
                    final SparseFieldMatrix<Binary64> b = new SparseFieldMatrix<>(Binary64Field.getInstance(), rows, interm);
                    b.walkInOptimizedOrder(randomSetter);
                   a.transposeMultiply(b).subtract(a.transpose().multiply(b)).walkInOptimizedOrder(zeroChecker);
                }
            }
        }
    }

    @Test
    void testTransposeMultiplyWrongDimensions() {
        try {
            new SparseFieldMatrix<>(Binary64Field.getInstance(), 2, 3).
            transposeMultiply(new SparseFieldMatrix<>(Binary64Field.getInstance(), 3, 2));
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
        FieldMatrix<Fraction> m = createSparseMatrix(id);
        assertEquals(3d, m.getTrace().doubleValue(), entryTolerance, "identity trace");
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
        FieldMatrix<Fraction> m = createSparseMatrix(testData);
        customAssertClose("scalar add", createSparseMatrix(testDataPlus2),
                          m.scalarAdd(new Fraction(2)), entryTolerance);
    }

    /** test operate */
    @Test
    void testOperate() {
        FieldMatrix<Fraction> m = createSparseMatrix(id);
        customAssertClose("identity operate", testVector, m.operate(testVector),
                          entryTolerance);
        customAssertClose("identity operate", testVector, m.operate(
                new ArrayFieldVector<Fraction>(testVector)).toArray(), entryTolerance);
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
        FieldMatrix<Fraction> a = createSparseMatrix(new Fraction[][] {
                { new Fraction(1), new Fraction(2) }, { new Fraction(3), new Fraction(4) }, { new Fraction(5), new Fraction(6) } });
        Fraction[] b = a.operate(new Fraction[] { new Fraction(1), new Fraction(1) });
        assertEquals(a.getRowDimension(), b.length);
        assertEquals(3.0, b[0].doubleValue(), 1.0e-12);
        assertEquals(7.0, b[1].doubleValue(), 1.0e-12);
        assertEquals(11.0, b[2].doubleValue(), 1.0e-12);
    }

    /** test transpose */
    @Test
    void testTranspose() {
        FieldMatrix<Fraction> m = createSparseMatrix(testData);
        FieldMatrix<Fraction> mIT = new FieldLUDecomposition<Fraction>(m).getSolver().getInverse().transpose();
        FieldMatrix<Fraction> mTI = new FieldLUDecomposition<Fraction>(m.transpose()).getSolver().getInverse();
        customAssertClose("inverse-transpose", mIT, mTI, normTolerance);
        m = createSparseMatrix(testData2);
        FieldMatrix<Fraction> mt = createSparseMatrix(testData2T);
        customAssertClose("transpose", mt, m.transpose(), normTolerance);
    }

    /** test preMultiply by vector */
    @Test
    void testPremultiplyVector() {
        FieldMatrix<Fraction> m = createSparseMatrix(testData);
        customAssertClose("premultiply", m.preMultiply(testVector), preMultTest,
                          normTolerance);
        customAssertClose("premultiply", m.preMultiply(
            new ArrayFieldVector<Fraction>(testVector).toArray()), preMultTest, normTolerance);
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
        FieldMatrix<Fraction> m3 = createSparseMatrix(d3);
        FieldMatrix<Fraction> m4 = createSparseMatrix(d4);
        FieldMatrix<Fraction> m5 = createSparseMatrix(d5);
        customAssertClose("m3*m4=m5", m4.preMultiply(m3), m5, entryTolerance);

        SparseFieldMatrix<Fraction> m = createSparseMatrix(testData);
        SparseFieldMatrix<Fraction> mInv = createSparseMatrix(testDataInv);
        SparseFieldMatrix<Fraction> identity = createSparseMatrix(id);
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
        FieldMatrix<Fraction> m = createSparseMatrix(testData);
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
        FieldMatrix<Fraction> m = createSparseMatrix(testData);
        assertEquals(2d, m.getEntry(0, 1).doubleValue(), entryTolerance, "get entry");
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
        Fraction[][] matrixData = { { new Fraction(1), new Fraction(2), new Fraction(3) }, { new Fraction(2), new Fraction(5), new Fraction(3) } };
        FieldMatrix<Fraction> m = createSparseMatrix(matrixData);
        // One more with three rows, two columns
        Fraction[][] matrixData2 = { { new Fraction(1), new Fraction(2) }, { new Fraction(2), new Fraction(5) }, { new Fraction(1), new Fraction(7) } };
        FieldMatrix<Fraction> n = createSparseMatrix(matrixData2);
        // Now multiply m by n
        FieldMatrix<Fraction> p = m.multiply(n);
        assertEquals(2, p.getRowDimension());
        assertEquals(2, p.getColumnDimension());
        // Invert p
        FieldMatrix<Fraction> pInverse = new FieldLUDecomposition<Fraction>(p).getSolver().getInverse();
        assertEquals(2, pInverse.getRowDimension());
        assertEquals(2, pInverse.getColumnDimension());

        // Solve example
        Fraction[][] coefficientsData = { { new Fraction(2), new Fraction(3), new Fraction(-2) }, { new Fraction(-1), new Fraction(7), new Fraction(6) },
                { new Fraction(4), new Fraction(-3), new Fraction(-5) } };
        FieldMatrix<Fraction> coefficients = createSparseMatrix(coefficientsData);
        Fraction[] constants = { new Fraction(1), new Fraction(-2), new Fraction(1) };
        Fraction[] solution;
        solution = new FieldLUDecomposition<Fraction>(coefficients)
            .getSolver()
            .solve(new ArrayFieldVector<Fraction>(constants, false)).toArray();
        assertEquals((new Fraction(2).multiply((solution[0])).add(new Fraction(3).multiply(solution[1])).subtract(new Fraction(2).multiply(solution[2]))).doubleValue(),
                constants[0].doubleValue(), 1E-12);
        assertEquals(((new Fraction(-1).multiply(solution[0])).add(new Fraction(7).multiply(solution[1])).add(new Fraction(6).multiply(solution[2]))).doubleValue(),
                constants[1].doubleValue(), 1E-12);
        assertEquals(((new Fraction(4).multiply(solution[0])).subtract(new Fraction(3).multiply( solution[1])).subtract(new Fraction(5).multiply(solution[2]))).doubleValue(),
                constants[2].doubleValue(), 1E-12);

    }

    // test submatrix accessors
    @Test
    void testSubMatrix() {
        FieldMatrix<Fraction> m = createSparseMatrix(subTestData);
        FieldMatrix<Fraction> mRows23Cols00 = createSparseMatrix(subRows23Cols00);
        FieldMatrix<Fraction> mRows00Cols33 = createSparseMatrix(subRows00Cols33);
        FieldMatrix<Fraction> mRows01Cols23 = createSparseMatrix(subRows01Cols23);
        FieldMatrix<Fraction> mRows02Cols13 = createSparseMatrix(subRows02Cols13);
        FieldMatrix<Fraction> mRows03Cols12 = createSparseMatrix(subRows03Cols12);
        FieldMatrix<Fraction> mRows03Cols123 = createSparseMatrix(subRows03Cols123);
        FieldMatrix<Fraction> mRows20Cols123 = createSparseMatrix(subRows20Cols123);
        FieldMatrix<Fraction> mRows31Cols31 = createSparseMatrix(subRows31Cols31);
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
        FieldMatrix<Fraction> m = createSparseMatrix(subTestData);
        FieldMatrix<Fraction> mRow0 = createSparseMatrix(subRow0);
        FieldMatrix<Fraction> mRow3 = createSparseMatrix(subRow3);
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
        FieldMatrix<Fraction> m = createSparseMatrix(subTestData);
        FieldMatrix<Fraction> mColumn1 = createSparseMatrix(subColumn1);
        FieldMatrix<Fraction> mColumn3 = createSparseMatrix(subColumn3);
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
        FieldMatrix<Fraction> m = createSparseMatrix(subTestData);
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
    void testGetColumnVector() {
        FieldMatrix<Fraction> m = createSparseMatrix(subTestData);
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

    private FieldVector<Fraction> columnToVector(Fraction[][] column) {
        Fraction[] data = new Fraction[column.length];
        for (int i = 0; i < data.length; ++i) {
            data[i] = column[i][0];
        }
        return new ArrayFieldVector<Fraction>(data, false);
    }

    @Test
    void testEqualsAndHashCode() {
        SparseFieldMatrix<Fraction> m = createSparseMatrix(testData);
        SparseFieldMatrix<Fraction> m1 = (SparseFieldMatrix<Fraction>) m.copy();
        SparseFieldMatrix<Fraction> mt = (SparseFieldMatrix<Fraction>) m.transpose();
        assertTrue(m.hashCode() != mt.hashCode());
        assertEquals(m.hashCode(), m1.hashCode());
        assertEquals(m, m);
        assertEquals(m, m1);
        assertNotEquals(null, m);
        assertNotEquals(m, mt);
        assertNotEquals(m, createSparseMatrix(bigSingular));
    }

    /* Disable for now
    @Test
    public void testToString() {
        SparseFieldMatrix<Fraction> m = createSparseMatrix(testData);
        Assertions.assertEquals("SparseFieldMatrix<Fraction>{{1.0,2.0,3.0},{2.0,5.0,3.0},{1.0,0.0,8.0}}",
            m.toString());
        m = new SparseFieldMatrix<Fraction>(field, 1, 1);
        Assertions.assertEquals("SparseFieldMatrix<Fraction>{{0.0}}", m.toString());
    }
    */

    @Test
    void testSetSubMatrix() {
        SparseFieldMatrix<Fraction> m = createSparseMatrix(testData);
        m.setSubMatrix(detData2, 1, 1);
        FieldMatrix<Fraction> expected = createSparseMatrix(new Fraction[][] {
                { new Fraction(1), new Fraction(2), new Fraction(3) }, { new Fraction(2), new Fraction(1), new Fraction(3) }, { new Fraction(1), new Fraction(2), new Fraction(4) } });
        assertEquals(expected, m);

        m.setSubMatrix(detData2, 0, 0);
        expected = createSparseMatrix(new Fraction[][] {
                { new Fraction(1), new Fraction(3), new Fraction(3) }, { new Fraction(2), new Fraction(4), new Fraction(3) }, { new Fraction(1), new Fraction(2), new Fraction(4) } });
        assertEquals(expected, m);

        m.setSubMatrix(testDataPlus2, 0, 0);
        expected = createSparseMatrix(new Fraction[][] {
                { new Fraction(3), new Fraction(4), new Fraction(5) }, { new Fraction(4), new Fraction(7), new Fraction(5) }, { new Fraction(3), new Fraction(2), new Fraction(10) } });
        assertEquals(expected, m);

        // javadoc example
        SparseFieldMatrix<Fraction> matrix =
            createSparseMatrix(new Fraction[][] {
        { new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4) }, { new Fraction(5), new Fraction(6), new Fraction(7), new Fraction(8) }, { new Fraction(9), new Fraction(0), new Fraction(1), new Fraction(2) } });
        matrix.setSubMatrix(new Fraction[][] { { new Fraction(3), new Fraction(4) }, { new Fraction(5), new Fraction(6) } }, 1, 1);
        expected = createSparseMatrix(new Fraction[][] {
                { new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4) }, { new Fraction(5), new Fraction(3), new Fraction(4), new Fraction(8) }, { new Fraction(9), new Fraction(5), new Fraction(6), new Fraction(2) } });
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
            new SparseFieldMatrix<Fraction>(field, 0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        // ragged
        try {
            m.setSubMatrix(new Fraction[][] { { new Fraction(1) }, { new Fraction(2), new Fraction(3) } }, 0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        // empty
        try {
            m.setSubMatrix(new Fraction[][] { {} }, 0, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
    }

    // --------------- -----------------Protected methods

    /** verifies that two matrices are close (1-norm) */
    protected void customAssertClose(String msg, FieldMatrix<Fraction> m, FieldMatrix<Fraction> n,
                                     double tolerance) {
        for(int i=0; i < m.getRowDimension(); i++){
            for(int j=0; j < m.getColumnDimension(); j++){
                assertEquals(m.getEntry(i,j).doubleValue(), n.getEntry(i,j).doubleValue(), tolerance, msg);
            }

        }
    }

    /** verifies that two vectors are close (sup norm) */
    protected void customAssertClose(String msg, Fraction[] m, Fraction[] n,
                                     double tolerance) {
        if (m.length != n.length) {
            fail("vectors not same length");
        }
        for (int i = 0; i < m.length; i++) {
            assertEquals(m[i].doubleValue(), n[i].doubleValue(),
                    tolerance,
                    msg + " " + i + " elements differ");
        }
    }

    private SparseFieldMatrix<Fraction> createSparseMatrix(Fraction[][] data) {
        SparseFieldMatrix<Fraction> matrix = new SparseFieldMatrix<Fraction>(field, data.length, data[0].length);
        for (int row = 0; row < data.length; row++) {
            for (int col = 0; col < data[row].length; col++) {
                matrix.setEntry(row, col, data[row][col]);
            }
        }
        return matrix;
    }
}
