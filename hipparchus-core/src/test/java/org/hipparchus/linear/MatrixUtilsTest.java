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

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.List;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.fraction.BigFraction;
import org.hipparchus.fraction.Fraction;
import org.hipparchus.fraction.FractionField;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test cases for the {@link MatrixUtils} class.
 *
 */

public final class MatrixUtilsTest {

    protected double[][] testData = { {1d,2d,3d}, {2d,5d,3d}, {1d,0d,8d} };
    protected double[][] testData3x3Singular = { { 1, 4, 7, }, { 2, 5, 8, }, { 3, 6, 9, } };
    protected double[][] testData3x4 = { { 12, -51, 4, 1 }, { 6, 167, -68, 2 }, { -4, 24, -41, 3 } };
    protected double[][] nullMatrix = null;
    protected double[] row = {1,2,3};
    protected BigDecimal[] bigRow =
        {new BigDecimal(1),new BigDecimal(2),new BigDecimal(3)};
    protected String[] stringRow = {"1", "2", "3"};
    protected Fraction[] fractionRow =
        {new Fraction(1),new Fraction(2),new Fraction(3)};
    protected double[][] rowMatrix = {{1,2,3}};
    protected BigDecimal[][] bigRowMatrix =
        {{new BigDecimal(1), new BigDecimal(2), new BigDecimal(3)}};
    protected String[][] stringRowMatrix = {{"1", "2", "3"}};
    protected Fraction[][] fractionRowMatrix =
        {{new Fraction(1), new Fraction(2), new Fraction(3)}};
    protected double[] col = {0,4,6};
    protected BigDecimal[] bigCol =
        {new BigDecimal(0),new BigDecimal(4),new BigDecimal(6)};
    protected String[] stringCol = {"0","4","6"};
    protected Fraction[] fractionCol =
        {new Fraction(0),new Fraction(4),new Fraction(6)};
    protected double[] nullDoubleArray = null;
    protected double[][] colMatrix = {{0},{4},{6}};
    protected BigDecimal[][] bigColMatrix =
        {{new BigDecimal(0)},{new BigDecimal(4)},{new BigDecimal(6)}};
    protected String[][] stringColMatrix = {{"0"}, {"4"}, {"6"}};
    protected Fraction[][] fractionColMatrix =
        {{new Fraction(0)},{new Fraction(4)},{new Fraction(6)}};

    @Test
    public void testCreateRealMatrix() {
        Assert.assertEquals(new BlockRealMatrix(testData),
                MatrixUtils.createRealMatrix(testData));
        try {
            MatrixUtils.createRealMatrix(new double[][] {{1}, {1,2}});  // ragged
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createRealMatrix(new double[][] {{}, {}});  // no columns
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createRealMatrix(null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testcreateFieldMatrix() {
        Assert.assertEquals(new Array2DRowFieldMatrix<Fraction>(asFraction(testData)),
                     MatrixUtils.createFieldMatrix(asFraction(testData)));
        Assert.assertEquals(new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), fractionColMatrix),
                     MatrixUtils.createFieldMatrix(fractionColMatrix));
        try {
            MatrixUtils.createFieldMatrix(asFraction(new double[][] {{1}, {1,2}}));  // ragged
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createFieldMatrix(asFraction(new double[][] {{}, {}}));  // no columns
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createFieldMatrix((Fraction[][])null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testCreateRealVector() {
        RealVector v1 = MatrixUtils.createRealVector(new double[] { 0.0, 1.0, 2.0, 3.0 });
        Assert.assertEquals(4, v1.getDimension());
        for (int i = 0; i < v1.getDimension(); ++i) {
            Assert.assertEquals(i, v1.getEntry(i), 1.0e-15);
        }
        RealVector v2 = MatrixUtils.createRealVector(7);
        Assert.assertEquals(7, v2.getDimension());
        for (int i = 0; i < v2.getDimension(); ++i) {
            Assert.assertEquals(0.0, v2.getEntry(i), 1.0e-15);
        }
    }

    @Test
    public void testCreateFieldVector() {
        FieldVector<Binary64> v1 = MatrixUtils.createFieldVector(new Binary64[] {
            new Binary64(0.0), new Binary64(1.0), new Binary64(2.0), new Binary64(3.0)
        });
        Assert.assertEquals(4, v1.getDimension());
        for (int i = 0; i < v1.getDimension(); ++i) {
            Assert.assertEquals(i, v1.getEntry(i).getReal(), 1.0e-15);
        }
        FieldVector<Binary64> v2 = MatrixUtils.createFieldVector(Binary64Field.getInstance(), 7);
        Assert.assertEquals(7, v2.getDimension());
        for (int i = 0; i < v2.getDimension(); ++i) {
            Assert.assertEquals(0.0, v2.getEntry(i).getReal(), 1.0e-15);
        }
    }

    @Test
    public void testCreateRowRealMatrix() {
        Assert.assertEquals(MatrixUtils.createRowRealMatrix(row),
                     new BlockRealMatrix(rowMatrix));
        try {
            MatrixUtils.createRowRealMatrix(new double[] {});  // empty
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createRowRealMatrix(null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testCreateRowFieldMatrix() {
        Assert.assertEquals(MatrixUtils.createRowFieldMatrix(asFraction(row)),
                     new Array2DRowFieldMatrix<Fraction>(asFraction(rowMatrix)));
        Assert.assertEquals(MatrixUtils.createRowFieldMatrix(fractionRow),
                     new Array2DRowFieldMatrix<Fraction>(fractionRowMatrix));
        try {
            MatrixUtils.createRowFieldMatrix(new Fraction[] {});  // empty
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createRowFieldMatrix((Fraction[]) null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testCreateColumnRealMatrix() {
        Assert.assertEquals(MatrixUtils.createColumnRealMatrix(col),
                     new BlockRealMatrix(colMatrix));
        try {
            MatrixUtils.createColumnRealMatrix(new double[] {});  // empty
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createColumnRealMatrix(null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testCreateColumnFieldMatrix() {
        Assert.assertEquals(MatrixUtils.createColumnFieldMatrix(asFraction(col)),
                     new Array2DRowFieldMatrix<Fraction>(asFraction(colMatrix)));
        Assert.assertEquals(MatrixUtils.createColumnFieldMatrix(fractionCol),
                     new Array2DRowFieldMatrix<Fraction>(fractionColMatrix));

        try {
            MatrixUtils.createColumnFieldMatrix(new Fraction[] {});  // empty
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createColumnFieldMatrix((Fraction[]) null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    /**
     * Verifies that the matrix is an identity matrix
     */
    protected void checkIdentityMatrix(RealMatrix m) {
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j =0; j < m.getColumnDimension(); j++) {
                if (i == j) {
                    Assert.assertEquals(m.getEntry(i, j), 1d, 0);
                } else {
                    Assert.assertEquals(m.getEntry(i, j), 0d, 0);
                }
            }
        }
    }

    @Test
    public void testCreateIdentityMatrix() {
        checkIdentityMatrix(MatrixUtils.createRealIdentityMatrix(3));
        checkIdentityMatrix(MatrixUtils.createRealIdentityMatrix(2));
        checkIdentityMatrix(MatrixUtils.createRealIdentityMatrix(1));
        try {
            MatrixUtils.createRealIdentityMatrix(0);
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    /**
     * Verifies that the matrix is an identity matrix
     */
    protected void checkIdentityFieldMatrix(FieldMatrix<Fraction> m) {
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j =0; j < m.getColumnDimension(); j++) {
                if (i == j) {
                    Assert.assertEquals(m.getEntry(i, j), Fraction.ONE);
                } else {
                    Assert.assertEquals(m.getEntry(i, j), Fraction.ZERO);
                }
            }
        }
    }

    @Test
    public void testcreateFieldIdentityMatrix() {
        checkIdentityFieldMatrix(MatrixUtils.createFieldIdentityMatrix(FractionField.getInstance(), 3));
        checkIdentityFieldMatrix(MatrixUtils.createFieldIdentityMatrix(FractionField.getInstance(), 2));
        checkIdentityFieldMatrix(MatrixUtils.createFieldIdentityMatrix(FractionField.getInstance(), 1));
        try {
            MatrixUtils.createRealIdentityMatrix(0);
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testBigFractionConverter() {
        BigFraction[][] bfData = {
                { new BigFraction(1), new BigFraction(2), new BigFraction(3) },
                { new BigFraction(2), new BigFraction(5), new BigFraction(3) },
                { new BigFraction(1), new BigFraction(0), new BigFraction(8) }
        };
        FieldMatrix<BigFraction> m = new Array2DRowFieldMatrix<BigFraction>(bfData, false);
        RealMatrix converted = MatrixUtils.bigFractionMatrixToRealMatrix(m);
        RealMatrix reference = new Array2DRowRealMatrix(testData, false);
        Assert.assertEquals(0.0, converted.subtract(reference).getNorm1(), 0.0);
    }

    @Test
    public void testFractionConverter() {
        Fraction[][] fData = {
                { new Fraction(1), new Fraction(2), new Fraction(3) },
                { new Fraction(2), new Fraction(5), new Fraction(3) },
                { new Fraction(1), new Fraction(0), new Fraction(8) }
        };
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(fData, false);
        RealMatrix converted = MatrixUtils.fractionMatrixToRealMatrix(m);
        RealMatrix reference = new Array2DRowRealMatrix(testData, false);
        Assert.assertEquals(0.0, converted.subtract(reference).getNorm1(), 0.0);
    }

    public static final Fraction[][] asFraction(double[][] data) {
        Fraction d[][] = new Fraction[data.length][];
        try {
            for (int i = 0; i < data.length; ++i) {
                double[] dataI = data[i];
                Fraction[] dI  = new Fraction[dataI.length];
                for (int j = 0; j < dataI.length; ++j) {
                    dI[j] = new Fraction(dataI[j]);
                }
                d[i] = dI;
            }
        } catch (MathIllegalStateException fce) {
            Assert.fail(fce.getMessage());
        }
        return d;
    }

    public static final Fraction[] asFraction(double[] data) {
        Fraction d[] = new Fraction[data.length];
        try {
            for (int i = 0; i < data.length; ++i) {
                d[i] = new Fraction(data[i]);
            }
        } catch (MathIllegalStateException fce) {
            Assert.fail(fce.getMessage());
        }
        return d;
    }

    @Test
    public void testSolveLowerTriangularSystem(){
        RealMatrix rm = new Array2DRowRealMatrix(
                new double[][] { {2,0,0,0 }, { 1,1,0,0 }, { 3,3,3,0 }, { 3,3,3,4 } },
                       false);
        RealVector b = new ArrayRealVector(new double[] { 2,3,4,8 }, false);
        MatrixUtils.solveLowerTriangularSystem(rm, b);
        UnitTestUtils.assertEquals( new double[]{1,2,-1.66666666666667, 1.0}  , b.toArray() , 1.0e-12);
    }


    /*
     * Taken from R manual http://stat.ethz.ch/R-manual/R-patched/library/base/html/backsolve.html
     */
    @Test
    public void testSolveUpperTriangularSystem(){
        RealMatrix rm = new Array2DRowRealMatrix(
                new double[][] { {1,2,3 }, { 0,1,1 }, { 0,0,2 } },
                       false);
        RealVector b = new ArrayRealVector(new double[] { 8,4,2 }, false);
        MatrixUtils.solveUpperTriangularSystem(rm, b);
        UnitTestUtils.assertEquals( new double[]{-1,3,1}  , b.toArray() , 1.0e-12);
    }

    /**
     * This test should probably be replaced by one that could show
     * whether this algorithm can sometimes perform better (precision- or
     * performance-wise) than the direct inversion of the whole matrix.
     */
    @Test
    public void testBlockInverse() {
        final double[][] data = {
            { -1, 0, 123, 4 },
            { -56, 78.9, -0.1, -23.4 },
            { 5.67, 8, -9, 1011 },
            { 12, 345, -67.8, 9 },
        };

        final RealMatrix m = new Array2DRowRealMatrix(data);
        final int len = data.length;
        final double tol = 1e-14;

        for (int splitIndex = 0; splitIndex < 3; splitIndex++) {
            final RealMatrix mInv = MatrixUtils.blockInverse(m, splitIndex);
            final RealMatrix id = m.multiply(mInv);

            // Check that we recovered the identity matrix.
            for (int i = 0; i < len; i++) {
                for (int j = 0; j < len; j++) {
                    final double entry = id.getEntry(i, j);
                    if (i == j) {
                        Assert.assertEquals("[" + i + "][" + j + "]",
                                            1, entry, tol);
                    } else {
                        Assert.assertEquals("[" + i + "][" + j + "]",
                                            0, entry, tol);
                    }
                }
            }
        }
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testBlockInverseNonInvertible() {
        final double[][] data = {
            { -1, 0, 123, 4 },
            { -56, 78.9, -0.1, -23.4 },
            { 5.67, 8, -9, 1011 },
            { 5.67, 8, -9, 1011 },
        };

        MatrixUtils.blockInverse(new Array2DRowRealMatrix(data), 2);
    }

    @Test
    public void testIsSymmetric() {
        final double eps = Math.ulp(1d);

        final double[][] dataSym = {
            { 1, 2, 3 },
            { 2, 2, 5 },
            { 3, 5, 6 },
        };
        Assert.assertTrue(MatrixUtils.isSymmetric(MatrixUtils.createRealMatrix(dataSym), eps));

        final double[][] dataNonSym = {
            { 1, 2, -3 },
            { 2, 2, 5 },
            { 3, 5, 6 },
        };
        Assert.assertFalse(MatrixUtils.isSymmetric(MatrixUtils.createRealMatrix(dataNonSym), eps));
    }

    @Test
    public void testIsSymmetricTolerance() {
        final double eps = 1e-4;

        final double[][] dataSym1 = {
            { 1,   1, 1.00009 },
            { 1,   1, 1       },
            { 1.0, 1, 1       },
        };
        Assert.assertTrue(MatrixUtils.isSymmetric(MatrixUtils.createRealMatrix(dataSym1), eps));
        final double[][] dataSym2 = {
            { 1,   1, 0.99990 },
            { 1,   1, 1       },
            { 1.0, 1, 1       },
        };
        Assert.assertTrue(MatrixUtils.isSymmetric(MatrixUtils.createRealMatrix(dataSym2), eps));

        final double[][] dataNonSym1 = {
            { 1,   1, 1.00011 },
            { 1,   1, 1       },
            { 1.0, 1, 1       },
        };
        Assert.assertFalse(MatrixUtils.isSymmetric(MatrixUtils.createRealMatrix(dataNonSym1), eps));
        final double[][] dataNonSym2 = {
            { 1,   1, 0.99989 },
            { 1,   1, 1       },
            { 1.0, 1, 1       },
        };
        Assert.assertFalse(MatrixUtils.isSymmetric(MatrixUtils.createRealMatrix(dataNonSym2), eps));
    }

    @Test
    public void testCheckSymmetric1() {
        final double[][] dataSym = {
            { 1, 2, 3 },
            { 2, 2, 5 },
            { 3, 5, 6 },
        };
        MatrixUtils.checkSymmetric(MatrixUtils.createRealMatrix(dataSym), Math.ulp(1d));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testCheckSymmetric2() {
        final double[][] dataNonSym = {
            { 1, 2, -3 },
            { 2, 2, 5 },
            { 3, 5, 6 },
        };
        MatrixUtils.checkSymmetric(MatrixUtils.createRealMatrix(dataNonSym), Math.ulp(1d));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testInverseSingular() {
        RealMatrix m = MatrixUtils.createRealMatrix(testData3x3Singular);
        MatrixUtils.inverse(m);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testInverseNonSquare() {
        RealMatrix m = MatrixUtils.createRealMatrix(testData3x4);
        MatrixUtils.inverse(m);
    }

    @Test
    public void testInverseDiagonalMatrix() {
        final double[] data = { 1, 2, 3 };
        final RealMatrix m = new DiagonalMatrix(data);
        final RealMatrix inverse = MatrixUtils.inverse(m);

        final RealMatrix result = m.multiply(inverse);
        UnitTestUtils.assertEquals("MatrixUtils.inverse() returns wrong result",
                MatrixUtils.createRealIdentityMatrix(data.length), result, Math.ulp(1d));
    }

    @Test
    public void testInverseRealMatrix() {
        RealMatrix m = MatrixUtils.createRealMatrix(testData);
        final RealMatrix inverse = MatrixUtils.inverse(m);

        final RealMatrix result = m.multiply(inverse);
        UnitTestUtils.assertEquals("MatrixUtils.inverse() returns wrong result",
                MatrixUtils.createRealIdentityMatrix(testData.length), result, 1e-12);
    }

    @Test
    public void testMatrixExponentialNonSquare() {
        double[][] exponentArr = {
                {0.0001, 0.001},
                {0.001, -0.0001},
                {0.001, -0.0001}
        };
        RealMatrix exponent = MatrixUtils.createRealMatrix(exponentArr);

        try {
            MatrixUtils.matrixExponential(exponent);  // ragged
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testMatrixExponential3() {
        double[][] exponentArr = {
                {0.0001, 0.001},
                {0.001, -0.0001}
        };
        RealMatrix exponent = MatrixUtils.createRealMatrix(exponentArr);

        double[][] expectedResultArr = {
                {1.00010050501688, 0.00100000016833332},
                {0.00100000016833332, 0.999900504983209}
        };
        RealMatrix expectedResult = MatrixUtils.createRealMatrix(expectedResultArr);

        UnitTestUtils.assertEquals("matrixExponential pade3 incorrect result",
                expectedResult, MatrixUtils.matrixExponential(exponent), 32.0 * Math.ulp(1.0));
    }


    @Test
    public void testMatrixExponential5() {
        double[][] exponentArr = {
                {0.1, 0.1},
                {0.001, -0.1}
        };
        RealMatrix exponent = MatrixUtils.createRealMatrix(exponentArr);

        double[][] expectedResultArr = {
                {1.10522267021001, 0.100168418362112},
                {0.00100168418362112, 0.904885833485786}
        };
        RealMatrix expectedResult = MatrixUtils.createRealMatrix(expectedResultArr);

        UnitTestUtils.assertEquals("matrixExponential pade5 incorrect result",
                expectedResult, MatrixUtils.matrixExponential(exponent), 2.0 * Math.ulp(1.0));
    }

    @Test
    public void testMatrixExponential7() {
        double[][] exponentArr = {
                {0.5, 0.1},
                {0.001, -0.5}
        };
        RealMatrix exponent = MatrixUtils.createRealMatrix(exponentArr);

        double[][] expectedResultArr = {
                {1.64878192423569, 0.104220769814317},
                {0.00104220769814317, 0.606574226092523}
        };
        RealMatrix expectedResult = MatrixUtils.createRealMatrix(expectedResultArr);

        UnitTestUtils.assertEquals("matrixExponential pade7 incorrect result",
                expectedResult, MatrixUtils.matrixExponential(exponent), 32.0 * Math.ulp(1.0));
    }

    @Test
    public void testMatrixExponential9() {
        double[][] exponentArr = {
                {1.8, 0.3},
                {0.001, -0.9}
        };
        RealMatrix exponent = MatrixUtils.createRealMatrix(exponentArr);

        double[][] expectedResultArr = {
                {6.05008743087114, 0.627036746099251},
                {0.00209012248699751, 0.406756715977872}
        };
        RealMatrix expectedResult = MatrixUtils.createRealMatrix(expectedResultArr);

        UnitTestUtils.assertEquals("matrixExponential pade9 incorrect result",
                expectedResult, MatrixUtils.matrixExponential(exponent), 16.0 * Math.ulp(1.0));
    }

    @Test
    public void testMatrixExponential13() {
        double[][] exponentArr1 = {
                {3.4, 1.2},
                {0.001, -0.9}
        };
        RealMatrix exponent1 = MatrixUtils.createRealMatrix(exponentArr1);

        double[][] expectedResultArr1 = {
                {29.9705442872504, 8.2499077972773},
                {0.00687492316439775, 0.408374680340048}
        };
        RealMatrix expectedResult1 = MatrixUtils.createRealMatrix(expectedResultArr1);

        UnitTestUtils.assertEquals("matrixExponential pade13-1 incorrect result",
                expectedResult1, MatrixUtils.matrixExponential(exponent1), 16.0 * Math.ulp(30.0));


        double[][] exponentArr2 = {
                {1.0, 1e5},
                {0.001, -1.0}
        };
        RealMatrix exponent2 = MatrixUtils.createRealMatrix(exponentArr2);

        double[][] expectedResultArr2 = {
                {12728.3536593144, 115190017.08756},
                {1.1519001708756, 10424.5533175632}
        };
        RealMatrix expectedResult2 = MatrixUtils.createRealMatrix(expectedResultArr2);

        UnitTestUtils.assertEquals("matrixExponential pade13-2 incorrect result",
                expectedResult2, MatrixUtils.matrixExponential(exponent2), 65536.0 * Math.ulp(1e8));


        double[][] exponentArr3 = {
                {-1e4, 1e4},
                {1.0, -1.0}
        };
        RealMatrix exponent3 = MatrixUtils.createRealMatrix(exponentArr3);

        double[][] expectedResultArr3 = {
                {9.99900009999e-05, 0.999900009999},
                {9.99900009999e-05, 0.999900009999}
        };
        RealMatrix expectedResult3 = MatrixUtils.createRealMatrix(expectedResultArr3);

        UnitTestUtils.assertEquals("matrixExponential pade13-3 incorrect result",
                expectedResult3, MatrixUtils.matrixExponential(exponent3), 4096.0 * Math.ulp(1.0));
    }

    @Test
    public void testOrthonormalize1() {

        final List<RealVector> basis =
                        MatrixUtils.orthonormalize(Arrays.asList(new ArrayRealVector(new double[] {  1, 2, 2 }),
                                                                 new ArrayRealVector(new double[] { -1, 0, 2 }),
                                                                 new ArrayRealVector(new double[] {  0, 0, 1 })),
                                                   Precision.EPSILON, DependentVectorsHandler.GENERATE_EXCEPTION);
        Assert.assertEquals(3, basis.size());
        checkBasis(basis);
        checkVector(basis.get(0),  1.0 / 3.0,  2.0 / 3.0, 2.0 / 3.0);
        checkVector(basis.get(1), -2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0);
        checkVector(basis.get(2),  2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0);

    }

    @Test
    public void testOrthonormalize2() {

        final List<RealVector> basis =
                        MatrixUtils.orthonormalize(Arrays.asList(new ArrayRealVector(new double[] { 3, 1 }),
                                                                 new ArrayRealVector(new double[] { 2, 2 })),
                                                   Precision.EPSILON, DependentVectorsHandler.GENERATE_EXCEPTION);
        final double s10 = FastMath.sqrt(10);
        Assert.assertEquals(2, basis.size());
        checkBasis(basis);
        checkVector(basis.get(0),  3 / s10,  1 / s10);
        checkVector(basis.get(1), -1 / s10,  3 / s10);

    }

    @Test
    public void testOrthonormalize3() {

        final double small = 1.0e-12;
        final List<RealVector> basis =
                        MatrixUtils.orthonormalize(Arrays.asList(new ArrayRealVector(new double[] { 1, small, small }),
                                                                 new ArrayRealVector(new double[] { 1, small, 0     }),
                                                                 new ArrayRealVector(new double[] { 1, 0,     small })),
                                                   Precision.EPSILON, DependentVectorsHandler.GENERATE_EXCEPTION);
        Assert.assertEquals(3, basis.size());
        checkBasis(basis);
        checkVector(basis.get(0), 1,  small, small);
        checkVector(basis.get(1), 0,  0,     -1   );
        checkVector(basis.get(2), 0, -1,      0   );

    }

    @Test
    public void testOrthonormalizeIncompleteBasis() {

        final double small = 1.0e-12;
        final List<RealVector> basis =
                        MatrixUtils.orthonormalize(Arrays.asList(new ArrayRealVector(new double[] { 1, small, small }),
                                                                 new ArrayRealVector(new double[] { 1, small, 0     })),
                                                   Precision.EPSILON, DependentVectorsHandler.GENERATE_EXCEPTION);
        Assert.assertEquals(2, basis.size());
        checkBasis(basis);
        checkVector(basis.get(0), 1,  small, small);
        checkVector(basis.get(1), 0,  0,     -1   );

    }

    @Test
    public void testOrthonormalizeDependent() {
        final double small = 1.0e-12;
        try {
            MatrixUtils.orthonormalize(Arrays.asList(new ArrayRealVector(new double[] { 1, small, small }),
                                                     new ArrayRealVector(new double[] { 1, small, small })),
                                       Precision.EPSILON, DependentVectorsHandler.GENERATE_EXCEPTION);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.ZERO_NORM, miae.getSpecifier());
        }
    }

    @Test
    public void testOrthonormalizeDependentGenerateException() {
        try {
            MatrixUtils.orthonormalize(Arrays.asList(new ArrayRealVector(new double[] { 2, 3, 0 }),
                                                     new ArrayRealVector(new double[] { 2, 7, 0 }),
                                                     new ArrayRealVector(new double[] { 4, 5, 0 }),
                                                     new ArrayRealVector(new double[] { 0, 0, 1 })),
                                       7 * Precision.EPSILON, DependentVectorsHandler.GENERATE_EXCEPTION);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.ZERO_NORM, miae.getSpecifier());
        }

    }

    @Test
    public void testOrthonormalizeDependentAddZeroNorm() {
        List<RealVector> basis = MatrixUtils.orthonormalize(Arrays.asList(new ArrayRealVector(new double[] { 2, 3, 0 }),
                                                                          new ArrayRealVector(new double[] { 2, 7, 0 }),
                                                                          new ArrayRealVector(new double[] { 4, 5, 0 }),
                                                                          new ArrayRealVector(new double[] { 0, 0, 1 })),
                                                            7 * Precision.EPSILON, DependentVectorsHandler.ADD_ZERO_VECTOR);
        Assert.assertEquals(4, basis.size());
        Assert.assertEquals(0, basis.get(2).getEntry(0), 1.0e-15);
        Assert.assertEquals(0, basis.get(2).getEntry(1), 1.0e-15);
        Assert.assertEquals(0, basis.get(2).getEntry(2), 1.0e-15);
    }

    @Test
    public void testOrthonormalizeDependentReduceToSpan() {
        List<RealVector> basis = MatrixUtils.orthonormalize(Arrays.asList(new ArrayRealVector(new double[] { 2, 3, 0 }),
                                                                          new ArrayRealVector(new double[] { 2, 7, 0 }),
                                                                          new ArrayRealVector(new double[] { 4, 5, 0 }),
                                                                          new ArrayRealVector(new double[] { 0, 0, 1 })),
                                                            7 * Precision.EPSILON, DependentVectorsHandler.REDUCE_BASE_TO_SPAN);
        Assert.assertEquals(3, basis.size());
        Assert.assertEquals(0, basis.get(2).getEntry(0), 1.0e-15);
        Assert.assertEquals(0, basis.get(2).getEntry(1), 1.0e-15);
        Assert.assertEquals(1, basis.get(2).getEntry(2), 1.0e-15);
    }

    @Test
    public void testFieldOrthonormalize1() {
        doTestOrthonormalize1(Binary64Field.getInstance());
    }

    @Test
    public void testFieldOrthonormalize2() {
        doTestOrthonormalize2(Binary64Field.getInstance());
    }

    @Test
    public void testFieldOrthonormalize3() {
        doTestOrthonormalize3(Binary64Field.getInstance());
    }

    @Test
    public void testFieldOrthonormalizeIncompleteBasis() {
        doTestOrthonormalizeIncompleteBasis(Binary64Field.getInstance());
    }

    @Test
    public void testFieldOrthonormalizeDependent() {
        doTestOrthonormalizeDependent(Binary64Field.getInstance());
    }

    @Test
    public void testFieldOrthonormalizeDependentGenerateException() {
        doTestOrthonormalizeDependentGenerateException(Binary64Field.getInstance());
    }

    @Test
    public void testFieldOrthonormalizeDependentAddZeroNorm() {
        doTestOrthonormalizeDependentAddZeroNorm(Binary64Field.getInstance());
    }

    @Test
    public void testFieldOrthonormalizeDependentReduceToSpan() {
        doTestOrthonormalizeDependentReduceToSpan(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestOrthonormalize1(final Field<T> field) {

        final List<FieldVector<T>> basis =
                        MatrixUtils.orthonormalize(field,
                                                   Arrays.asList(convert(field, new ArrayRealVector(new double[] {  1, 2, 2 })),
                                                                 convert(field, new ArrayRealVector(new double[] { -1, 0, 2 })),
                                                                 convert(field, new ArrayRealVector(new double[] {  0, 0, 1 }))),
                                                   field.getZero().newInstance(Precision.EPSILON),
                                                   DependentVectorsHandler.GENERATE_EXCEPTION);
        Assert.assertEquals(3, basis.size());
        checkBasis(field, basis);
        checkVector(basis.get(0),  1.0 / 3.0,  2.0 / 3.0, 2.0 / 3.0);
        checkVector(basis.get(1), -2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0);
        checkVector(basis.get(2),  2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0);

    }

    private <T extends CalculusFieldElement<T>> void doTestOrthonormalize2(final Field<T> field) {

        final List<FieldVector<T>> basis =
                        MatrixUtils.orthonormalize(field,
                                                   Arrays.asList(convert(field, new ArrayRealVector(new double[] { 3, 1 })),
                                                                 convert(field, new ArrayRealVector(new double[] { 2, 2 }))),
                                                   field.getZero().newInstance(Precision.EPSILON),
                                                   DependentVectorsHandler.GENERATE_EXCEPTION);
        final double s10 = FastMath.sqrt(10);
        Assert.assertEquals(2, basis.size());
        checkBasis(field, basis);
        checkVector(basis.get(0),  3 / s10,  1 / s10);
        checkVector(basis.get(1), -1 / s10,  3 / s10);

    }

    private <T extends CalculusFieldElement<T>> void doTestOrthonormalize3(final Field<T> field) {

        final double small = 1.0e-12;
        final List<FieldVector<T>> basis =
                        MatrixUtils.orthonormalize(field,
                                                   Arrays.asList(convert(field, new ArrayRealVector(new double[] { 1, small, small })),
                                                                 convert(field, new ArrayRealVector(new double[] { 1, small, 0     })),
                                                                 convert(field, new ArrayRealVector(new double[] { 1, 0,     small }))),
                                                   field.getZero().newInstance(Precision.EPSILON),
                                                   DependentVectorsHandler.GENERATE_EXCEPTION);
        Assert.assertEquals(3, basis.size());
        checkBasis(field, basis);
        checkVector(basis.get(0), 1,  small, small);
        checkVector(basis.get(1), 0,  0,     -1   );
        checkVector(basis.get(2), 0, -1,      0   );

    }

    private <T extends CalculusFieldElement<T>> void doTestOrthonormalizeIncompleteBasis(final Field<T> field) {

        final double small = 1.0e-12;
        final List<FieldVector<T>> basis =
                        MatrixUtils.orthonormalize(field,
                                                   Arrays.asList(convert(field, new ArrayRealVector(new double[] { 1, small, small })),
                                                                 convert(field, new ArrayRealVector(new double[] { 1, small, 0     }))),
                                                   field.getZero().newInstance(Precision.EPSILON),
                                                   DependentVectorsHandler.GENERATE_EXCEPTION);
        Assert.assertEquals(2, basis.size());
        checkBasis(field, basis);
        checkVector(basis.get(0), 1,  small, small);
        checkVector(basis.get(1), 0,  0,     -1   );

    }

    private <T extends CalculusFieldElement<T>> void doTestOrthonormalizeDependent(final Field<T> field) {
        final double small = 1.0e-12;
        try {
            MatrixUtils.orthonormalize(field,
                                       Arrays.asList(convert(field, new ArrayRealVector(new double[] { 1, small, small })),
                                                     convert(field, new ArrayRealVector(new double[] { 1, small, small }))),
                                       field.getZero().newInstance(Precision.EPSILON),
                                       DependentVectorsHandler.GENERATE_EXCEPTION);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.ZERO_NORM, miae.getSpecifier());
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestOrthonormalizeDependentGenerateException(final Field<T> field) {
        try {
            MatrixUtils.orthonormalize(field,
                                       Arrays.asList(convert(field, new ArrayRealVector(new double[] { 2, 3, 0 })),
                                                     convert(field, new ArrayRealVector(new double[] { 2, 7, 0 })),
                                                     convert(field, new ArrayRealVector(new double[] { 4, 5, 0 })),
                                                     convert(field, new ArrayRealVector(new double[] { 0, 0, 1 }))),
                                       field.getZero().newInstance(7 * Precision.EPSILON),
                                       DependentVectorsHandler.GENERATE_EXCEPTION);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.ZERO_NORM, miae.getSpecifier());
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestOrthonormalizeDependentAddZeroNorm(final Field<T> field) {
        List<FieldVector<T>> basis = MatrixUtils.orthonormalize(field,
                                                                Arrays.asList(convert(field, new ArrayRealVector(new double[] { 2, 3, 0 })),
                                                                              convert(field, new ArrayRealVector(new double[] { 2, 7, 0 })),
                                                                              convert(field, new ArrayRealVector(new double[] { 4, 5, 0 })),
                                                                              convert(field, new ArrayRealVector(new double[] { 0, 0, 1 }))),
                                                                field.getZero().newInstance(7 * Precision.EPSILON),
                                                                DependentVectorsHandler.ADD_ZERO_VECTOR);
       Assert.assertEquals(4, basis.size());
       Assert.assertEquals(0, basis.get(2).getEntry(0).getReal(), 1.0e-15);
       Assert.assertEquals(0, basis.get(2).getEntry(1).getReal(), 1.0e-15);
       Assert.assertEquals(0, basis.get(2).getEntry(2).getReal(), 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestOrthonormalizeDependentReduceToSpan(final Field<T> field) {
        List<FieldVector<T>> basis = MatrixUtils.orthonormalize(field,
                                                                Arrays.asList(convert(field, new ArrayRealVector(new double[] { 2, 3, 0 })),
                                                                              convert(field, new ArrayRealVector(new double[] { 2, 7, 0 })),
                                                                              convert(field, new ArrayRealVector(new double[] { 4, 5, 0 })),
                                                                              convert(field, new ArrayRealVector(new double[] { 0, 0, 1 }))),
                                                                field.getZero().newInstance(7 * Precision.EPSILON),
                                                                DependentVectorsHandler.REDUCE_BASE_TO_SPAN);
        Assert.assertEquals(3, basis.size());
        Assert.assertEquals(0, basis.get(2).getEntry(0).getReal(), 1.0e-15);
        Assert.assertEquals(0, basis.get(2).getEntry(1).getReal(), 1.0e-15);
        Assert.assertEquals(1, basis.get(2).getEntry(2).getReal(), 1.0e-15);
    }

    private void checkVector(final RealVector v, double... p) {
        Assert.assertEquals(p.length, v.getDimension());
        for (int i = 0; i < p.length; ++i) {
            Assert.assertEquals(p[i], v.getEntry(i), 1.0e-15);
        }
    }

    private void checkBasis(final List<RealVector> basis) {
        for (int i = 0; i < basis.size(); ++i) {
            for (int j = i; j < basis.size(); ++j) {
                Assert.assertEquals(i == j ? 1.0 : 0.0, basis.get(i).dotProduct(basis.get(j)), 1.0e-12);
            }
        }        
    }

    private <T extends CalculusFieldElement<T>> void checkVector(final FieldVector<T> v, double... p) {
        Assert.assertEquals(p.length, v.getDimension());
        for (int i = 0; i < p.length; ++i) {
            Assert.assertEquals(p[i], v.getEntry(i).getReal(), 1.0e-15);
        }
    }

    private <T extends CalculusFieldElement<T>> void checkBasis(final Field<T> field, final List<FieldVector<T>> basis) {
        for (int i = 0; i < basis.size(); ++i) {
            for (int j = i; j < basis.size(); ++j) {
                Assert.assertEquals(i == j ? 1.0 : 0.0, basis.get(i).dotProduct(basis.get(j)).getReal(), 1.0e-12);
            }
        }        
    }

    private <T extends CalculusFieldElement<T>> FieldVector<T> convert(final Field<T> field, final RealVector v) {
        ArrayFieldVector<T> c = new ArrayFieldVector<T>(v.getDimension(), field.getZero());
        for (int k = 0; k < v.getDimension(); ++k) {
            c.setEntry(k, field.getZero().newInstance(v.getEntry(k)));
        }
        return c;
    }

}
