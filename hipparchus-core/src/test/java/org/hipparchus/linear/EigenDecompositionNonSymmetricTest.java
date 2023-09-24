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

import static org.junit.Assert.assertNotNull;

import java.util.Arrays;
import java.util.Random;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexComparator;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class EigenDecompositionNonSymmetricTest {

    private Complex[] refValues;
    private RealMatrix matrix;

    /** test dimensions */
    @Test
    public void testDimensions() {
        final int m = matrix.getRowDimension();
        EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(matrix);
        Assert.assertEquals(m, ed.getV().getRowDimension());
        Assert.assertEquals(m, ed.getV().getColumnDimension());
        Assert.assertEquals(m, ed.getD().getColumnDimension());
        Assert.assertEquals(m, ed.getD().getColumnDimension());
        Assert.assertEquals(m, ed.getVInv().getRowDimension());
        Assert.assertEquals(m, ed.getVInv().getColumnDimension());
    }

    /** test eigenvalues */
    @Test
    public void testEigenvalues() {
        EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(matrix);
        Complex[] eigenValues = ed.getEigenvalues();
        Assert.assertEquals(refValues.length, eigenValues.length);
        for (int i = 0; i < refValues.length; ++i) {
            Assert.assertEquals(refValues[i].getRealPart(),      eigenValues[i].getRealPart(), 3.0e-15);
            Assert.assertEquals(refValues[i].getImaginaryPart(), eigenValues[i].getImaginaryPart(), 3.0e-15);
        }
    }

    @Test
    public void testNonSymmetric() {
        // Vandermonde matrix V(x;i,j) = x_i^{n - j} with x = (-1,-2,3,4)
        double[][] vData = { { -1.0, 1.0, -1.0, 1.0 },
                             { -8.0, 4.0, -2.0, 1.0 },
                             { 27.0, 9.0,  3.0, 1.0 },
                             { 64.0, 16.0, 4.0, 1.0 } };
        checkNonSymmetricMatrix(MatrixUtils.createRealMatrix(vData));

        RealMatrix randMatrix = MatrixUtils.createRealMatrix(new double[][] {
                {0,  1,     0,     0},
                {1,  0,     2.e-7, 0},
                {0, -2.e-7, 0,     1},
                {0,  0,     1,     0}
        });
        checkNonSymmetricMatrix(randMatrix);

        // from http://eigen.tuxfamily.org/dox/classEigen_1_1RealSchur.html
        double[][] randData2 = {
                {  0.680, -0.3300, -0.2700, -0.717, -0.687,  0.0259 },
                { -0.211,  0.5360,  0.0268,  0.214, -0.198,  0.6780 },
                {  0.566, -0.4440,  0.9040, -0.967, -0.740,  0.2250 },
                {  0.597,  0.1080,  0.8320, -0.514, -0.782, -0.4080 },
                {  0.823, -0.0452,  0.2710, -0.726,  0.998,  0.2750 },
                { -0.605,  0.2580,  0.4350,  0.608, -0.563,  0.0486 }
        };
        checkNonSymmetricMatrix(MatrixUtils.createRealMatrix(randData2));
    }

    @Test
    public void testRandomNonSymmetricMatrix() {
        for (int run = 0; run < 100; run++) {
            RandomGenerator r = new Well1024a(0x171956baefeac83el);

            // matrix size
            int size = r.nextInt(20) + 4;

            double[][] data = new double[size][size];
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    data[i][j] = r.nextInt(100);
                }
            }

            RealMatrix m = MatrixUtils.createRealMatrix(data);
            checkNonSymmetricMatrix(m);
        }
    }

    /**
     * Tests the porting of a bugfix in Jama-1.0.3 (from changelog):
     *
     *  Patched hqr2 method in Jama.EigenvalueDecomposition to avoid infinite loop;
     *  Thanks Frederic Devernay &lt;frederic.devernay@m4x.org&gt;
     */
    @Test
    public void testMath1051() {
        double[][] data = {
                {0,0,0,0,0},
                {0,0,0,0,1},
                {0,0,0,1,0},
                {1,1,0,0,1},
                {1,0,1,0,1}
        };

        RealMatrix m = MatrixUtils.createRealMatrix(data);
        checkNonSymmetricMatrix(m);
    }

    @Test
    public void testNormalDistributionNonSymmetricMatrix() {
        for (int run = 0; run < 100; run++) {
            final RandomGenerator r = new Well1024a(0x511d8551a5641ea2l);
            final RandomDataGenerator gen = new RandomDataGenerator(100);

            // matrix size
            int size = r.nextInt(20) + 4;

            double[][] data = new double[size][size];
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    data[i][j] = gen.nextNormal(0.0, r.nextDouble() * 5);
                }
            }

            RealMatrix m = MatrixUtils.createRealMatrix(data);
            checkNonSymmetricMatrix(m);
        }
    }

    @Test
    public void testMath848() {
        double[][] data = {
                { 0.1849449280, -0.0646971046,  0.0774755812, -0.0969651755, -0.0692648806,  0.3282344352, -0.0177423074,  0.2063136340},
                {-0.0742700134, -0.0289063030, -0.0017269460, -0.0375550146, -0.0487737922, -0.2616837868, -0.0821201295, -0.2530000167},
                { 0.2549910127,  0.0995733692, -0.0009718388,  0.0149282808,  0.1791878897, -0.0823182816,  0.0582629256,  0.3219545182},
                {-0.0694747557, -0.1880649148, -0.2740630911,  0.0720096468, -0.1800836914, -0.3518996425,  0.2486747833,  0.6257938167},
                { 0.0536360918, -0.1339297778,  0.2241579764, -0.0195327484, -0.0054103808,  0.0347564518,  0.5120802482, -0.0329902864},
                {-0.5933332356, -0.2488721082,  0.2357173629,  0.0177285473,  0.0856630593, -0.3567126300, -0.1600668126, -0.1010899621},
                {-0.0514349819, -0.0854319435,  0.1125050061,  0.0063453560, -0.2250000688, -0.2209343090,  0.1964623477, -0.1512329924},
                { 0.0197395947, -0.1997170581, -0.1425959019, -0.2749477910, -0.0969467073,  0.0603688520, -0.2826905192,  0.1794315473}};
        RealMatrix m = MatrixUtils.createRealMatrix(data);
        checkNonSymmetricMatrix(m);
    }

    /**
     * Checks that the eigen decomposition of a general (non-symmetric) matrix is valid by
     * checking: A*V = V*D
     */
    private void checkNonSymmetricMatrix(final RealMatrix m) {
        try {
            EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(m, 1.0e-15);

            RealMatrix d = ed.getD();
            RealMatrix v = ed.getV();

            RealMatrix x = m.multiply(v);
            RealMatrix y = v.multiply(d);

            double diffNorm = x.subtract(y).getNorm1();
            Assert.assertTrue("The norm of (X-Y) is too large: " + diffNorm + ", matrix=" + m.toString(),
                    x.subtract(y).getNorm1() < 1000 * Precision.EPSILON * FastMath.max(x.getNorm1(), y.getNorm1()));

            RealMatrix invV = new LUDecomposition(v).getSolver().getInverse();
            double norm = v.multiply(d).multiply(invV).subtract(m).getNorm1();
            Assert.assertEquals(0.0, norm, 1.0e-10);
        } catch (Exception e) {
            Assert.fail("Failed to create EigenDecomposition for matrix " + m.toString() + ", ex=" + e.toString());
        }
    }

    /** test eigenvectors */
    @Test
    public void testEigenvectors() {
        EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(matrix);
        FieldMatrix<Complex> cMatrix = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),
                                                                     matrix.getRowDimension(),
                                                                     matrix.getColumnDimension());
        for (int i = 0; i < matrix.getRowDimension(); ++i) {
            for (int j = 0; j < matrix.getColumnDimension(); ++j) {
                cMatrix.setEntry(i, j, new Complex(matrix.getEntry(j, i)));
            }
        }

        for (int i = 0; i < matrix.getRowDimension(); ++i) {
            final Complex              lambda = ed.getEigenvalue(i);
            final FieldVector<Complex> v      = ed.getEigenvector(i);
            final FieldVector<Complex> mV     = cMatrix.operate(v);
            for (int k = 0; k < v.getDimension(); ++k) {
                Assert.assertEquals(0, mV.getEntry(k).subtract(v.getEntry(k).multiply(lambda)).norm(), 1.0e-13);
            }
        }
    }

    /** test A = VDV⁻¹ */
    @Test
    public void testAEqualVDVInv() {
        EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(matrix);
        RealMatrix v  = ed.getV();
        RealMatrix d  = ed.getD();
        RealMatrix vI = MatrixUtils.inverse(v);
        double norm = v.multiply(d).multiply(vI).subtract(matrix).getNorm1();
        Assert.assertEquals(0, norm, 6.0e-13);
    }

    /**
     * Matrix with eigenvalues {2, 0, 12}
     */
    @Test
    public void testDistinctEigenvalues() {
        RealMatrix distinct = MatrixUtils.createRealMatrix(new double[][] {
                {3, 1, -4},
                {1, 3, -4},
                {-4, -4, 8}
        });
        EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(distinct);
        checkEigenValues((new Complex[] { new Complex( 2), new Complex( 0), new Complex(12) }), ed, 1E-12);
        checkEigenVector((new Complex[] { new Complex( 1), new Complex(-1), new Complex( 0) }), ed, 1E-12);
        checkEigenVector((new Complex[] { new Complex( 1), new Complex( 1), new Complex( 1) }), ed, 1E-12);
        checkEigenVector((new Complex[] { new Complex(-1), new Complex(-1), new Complex( 2) }), ed, 1E-12);
    }

    /**
     * Verifies operation on indefinite matrix
     */
    @Test
    public void testZeroDivide() {
        RealMatrix indefinite = MatrixUtils.createRealMatrix(new double [][] {
                { 0.0, 1.0, -1.0 },
                { 1.0, 1.0, 0.0 },
                { -1.0,0.0, 1.0 }
        });
        EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(indefinite);
        checkEigenValues((new Complex[] { new Complex(2), new Complex(1), new Complex(-1) }), ed, 1E-12);
        Complex isqrt3 = new Complex(1 / FastMath.sqrt(3.0));
        checkEigenVector((new Complex[] {isqrt3, isqrt3, isqrt3.negate() }), ed, 1E-12);
        Complex isqrt2 = new Complex(1 / FastMath.sqrt(2.0));
        checkEigenVector((new Complex[] { new Complex(0.0), isqrt2.negate(), isqrt2.negate() }), ed, 1E-12);
        Complex isqrt6 = new Complex(1 / FastMath.sqrt(6.0));
        checkEigenVector((new Complex[] { isqrt6.multiply(2), isqrt6.negate(), isqrt6 }), ed, 1E-12);
    }

    /** test eigenvalues for a big matrix. */
    @Test
    public void testBigMatrix() {
        Random r = new Random(17748333525117l);
        double[] bigValues = new double[200];
        for (int i = 0; i < bigValues.length; ++i) {
            bigValues[i] = 2 * r.nextDouble() - 1;
        }
        Arrays.sort(bigValues);
        EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(createTestMatrix(r, bigValues));
        Complex[] eigenValues = ed.getEigenvalues();
        Arrays.sort(eigenValues, new ComplexComparator());
        Assert.assertEquals(bigValues.length, eigenValues.length);
        for (int i = 0; i < bigValues.length; ++i) {
            Assert.assertEquals(bigValues[i], eigenValues[i].getRealPart(), 2.0e-14);
        }
    }

    /**
     * Verifies operation on very small values.
     * Matrix with eigenvalues {2e-100, 0, 12e-100}
     */
    @Test
    public void testTinyValues() {
        final double tiny = 1.0e-100;
        RealMatrix distinct = MatrixUtils.createRealMatrix(new double[][] {
                {3, 1, -4},
                {1, 3, -4},
                {-4, -4, 8}
        });
        distinct = distinct.scalarMultiply(tiny);

        final EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(distinct);
        checkEigenValues(new Complex[] { new Complex(2).multiply(tiny), new Complex(0).multiply(tiny), new Complex(12).multiply(tiny) },
                         ed, 1e-12 * tiny);
        checkEigenVector(new Complex[] { new Complex( 1), new Complex(-1), new Complex(0) }, ed, 1e-12);
        checkEigenVector(new Complex[] { new Complex( 1), new Complex( 1), new Complex(1) }, ed, 1e-12);
        checkEigenVector(new Complex[] { new Complex(-1), new Complex(-1), new Complex(2) }, ed, 1e-12);
    }

    @Test
    public void testDeterminantWithCompleEigenValues() {
        final RealMatrix m = MatrixUtils.createRealMatrix(new double[][] { { -3, -1.5, -3 }, { 0, -1, 0 }, { 1, 0, 0 } });
        EigenDecompositionNonSymmetric decomposition = new EigenDecompositionNonSymmetric(m);
        Assert.assertEquals(-3.0, decomposition.getDeterminant().getRealPart(),      1.0e-15);
        Assert.assertEquals( 0.0, decomposition.getDeterminant().getImaginaryPart(), 1.0e-15);
    }

    @Test
    public void testReal() {
        // AA = [1 2;1 -3];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 1, 2 }, { 1, -3 }
        });

        EigenDecompositionNonSymmetric ed = new EigenDecompositionNonSymmetric(A);

        assertNotNull(ed.getD());
        assertNotNull(ed.getV());

        final double s2 = FastMath.sqrt(2);
        final double s3 = FastMath.sqrt(3);
        RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
            { s2 * s3 - 1.0, 0 }, { 0,  -1 - s2 * s3 }
        });

        RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
            { (s2 + 2 * s3) / 5, (s3 - 3 * s2)     /  6 },
            { (2 * s2 - s3) / 5, (4 * s3 + 3 * s2) / 12 }
        });

        UnitTestUtils.assertEquals("D", D_expected, ed.getD(), 1.0e-15);
        UnitTestUtils.assertEquals("V", V_expected, ed.getV(), 1.0e-15);

        // checking definition of the decomposition A = V*D*inv(V)
        UnitTestUtils.assertEquals("A", A,
                                   ed.getV().multiply(ed.getD()).multiply(MatrixUtils.inverse(ed.getV())),
                                   2.0e-15);

    }

    @Test
    public void testImaginary() {
        // AA = [3 -2;4 -1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 3, -2 }, { 4, -1 }
        });

        EigenDecompositionNonSymmetric ordEig = new EigenDecompositionNonSymmetric(A);

        assertNotNull(ordEig.getD());
        assertNotNull(ordEig.getV());

        RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 1, 2 }, { -2, 1 }
        });

        RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
            { -0.5, 0.5 }, { 0, 1 }
        });

        UnitTestUtils.assertEquals("D", D_expected, ordEig.getD(), 1.0e-15);
        UnitTestUtils.assertEquals("V", V_expected, ordEig.getV(), 1.0e-15);

        // checking definition of the decomposition A = V*D*inv(V)
        UnitTestUtils.assertEquals("A", A, ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV())), 1.0e-15);

        // checking A*V = D*V
        FieldMatrix<Complex> ac = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), 2, 2);
        for (int i = 0; i < ac.getRowDimension(); ++i) {
            for (int j = 0; j < ac.getColumnDimension(); ++j) {
                ac.setEntry(i, j, new Complex(A.getEntry(i, j)));
            }
        }
        for (int i = 0; i < ordEig.getEigenvalues().length; ++i) {
            final Complex              li  = ordEig.getEigenvalue(i);
            final FieldVector<Complex> ei  = ordEig.getEigenvector(i);
            final FieldVector<Complex> aei = ac.operate(ei);
            for (int j = 0; j < ei.getDimension(); ++j) {
                Assert.assertEquals(aei.getEntry(j).getRealPart(),      li.multiply(ei.getEntry(j)).getRealPart(),      1.0e-10);
                Assert.assertEquals(aei.getEntry(j).getImaginaryPart(), li.multiply(ei.getEntry(j)).getImaginaryPart(), 1.0e-10);
            }
        }

    }

    @Test
    public void testImaginary33() {
        // AA = [3 -2 0;4 -1 0;1 1 1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 3, -2, 0 }, { 4, -1, 0 }, { 1, 1, 1 }
        });

        EigenDecompositionNonSymmetric ordEig = new EigenDecompositionNonSymmetric(A);

        assertNotNull(ordEig.getD());
        assertNotNull(ordEig.getV());

        RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 1, 2, 0 }, { -2, 1, 0}, {0, 0, 1 }
        });

        final double a = FastMath.sqrt(8.0 / 17.0);
        final double b = FastMath.sqrt(17.0) / 2.0;
        RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 0,        a, 0 },
            { a,        a, 0 },
            { a, -0.5 * a, b }
        });
        UnitTestUtils.assertEquals("D", D_expected, ordEig.getD(), 1.0e-15);
        UnitTestUtils.assertEquals("V", V_expected, ordEig.getV(), 1.0e-15);

        // checking definition of the decomposition A = V*D*inv(V)
        UnitTestUtils.assertEquals("A", A,
                                   ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV())),
                                   8.0e-15);

        // checking A*V = D*V
        FieldMatrix<Complex> ac = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), 3, 3);
        for (int i = 0; i < ac.getRowDimension(); ++i) {
            for (int j = 0; j < ac.getColumnDimension(); ++j) {
                ac.setEntry(i, j, new Complex(A.getEntry(i, j)));
            }
        }
        for (int i = 0; i < ordEig.getEigenvalues().length; ++i) {
            final Complex              li  = ordEig.getEigenvalue(i);
            final FieldVector<Complex> ei  = ordEig.getEigenvector(i);
            final FieldVector<Complex> aei = ac.operate(ei);
            for (int j = 0; j < ei.getDimension(); ++j) {
                Assert.assertEquals(aei.getEntry(j).getRealPart(),      li.multiply(ei.getEntry(j)).getRealPart(),      1.0e-10);
                Assert.assertEquals(aei.getEntry(j).getImaginaryPart(), li.multiply(ei.getEntry(j)).getImaginaryPart(), 1.0e-10);
            }
        }

    }

    @Test
    public void testImaginaryNullEigenvalue() {
        // AA = [3 -2 0;4 -1 0;3 -2 0];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 3, -2, 0 }, { 4, -1, 0 }, { 3, -2, 0 }
        });

        EigenDecompositionNonSymmetric ordEig = new EigenDecompositionNonSymmetric(A);

        assertNotNull(ordEig.getD());
        assertNotNull(ordEig.getV());

        RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 1, 2, 0 }, { -2, 1, 0 }, { 0, 0, 0 }
        });

        final double a  = FastMath.sqrt(11.0 / 50.0);
        final double s2 = FastMath.sqrt(2.0);
        RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
            {      -a, -2 * a,  0.0 },
            {  -3 * a,     -a,  0.0 },
            {      -a, -2 * a,   s2 }
        });

        UnitTestUtils.assertEquals("D", D_expected, ordEig.getD(), 1.0e-15);
        UnitTestUtils.assertEquals("V", V_expected, ordEig.getV(), 2.0e-15);

        // checking definition of the decomposition A = V*D*inv(V)
        UnitTestUtils.assertEquals("A", A,
                                   ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV())),
                                   1.0e-14);

        // checking A*V = D*V
        FieldMatrix<Complex> ac = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), 3, 3);
        for (int i = 0; i < ac.getRowDimension(); ++i) {
            for (int j = 0; j < ac.getColumnDimension(); ++j) {
                ac.setEntry(i, j, new Complex(A.getEntry(i, j)));
            }
        }

        for (int i = 0; i < ordEig.getEigenvalues().length; ++i) {
            final Complex              li  = ordEig.getEigenvalue(i);
            final FieldVector<Complex> ei  = ordEig.getEigenvector(i);
            final FieldVector<Complex> aei = ac.operate(ei);
            for (int j = 0; j < ei.getDimension(); ++j) {
                Assert.assertEquals(aei.getEntry(j).getRealPart(),      li.multiply(ei.getEntry(j)).getRealPart(),      1.0e-10);
                Assert.assertEquals(aei.getEntry(j).getImaginaryPart(), li.multiply(ei.getEntry(j)).getImaginaryPart(), 1.0e-10);
            }
        }

    }

    /**
     * Verifies that the given EigenDecomposition has eigenvalues equivalent to
     * the targetValues, ignoring the order of the values and allowing
     * values to differ by tolerance.
     */
    protected void checkEigenValues(Complex[] targetValues,
                                    EigenDecompositionNonSymmetric ed, double tolerance) {
        Complex[] observed = ed.getEigenvalues();
        for (int i = 0; i < observed.length; i++) {
            Assert.assertTrue(isIncludedValue(observed[i], targetValues, tolerance));
            Assert.assertTrue(isIncludedValue(targetValues[i], observed, tolerance));
        }
    }


    /**
     * Returns true iff there is an entry within tolerance of value in
     * searchArray.
     */
    private boolean isIncludedValue(Complex value, Complex[] searchArray, double tolerance) {
       boolean found = false;
       int i = 0;
       while (!found && i < searchArray.length) {
           if (value.subtract(searchArray[i]).norm() < tolerance) {
               found = true;
           }
           i++;
       }
       return found;
    }

    /**
     * Returns true iff eigenVector is a scalar multiple of one of the columns
     * of ed.getV().  Does not try linear combinations - i.e., should only be
     * used to find vectors in one-dimensional eigenspaces.
     */
    protected void checkEigenVector(Complex[] eigenVector, EigenDecompositionNonSymmetric ed, double tolerance) {
        Assert.assertTrue(isIncludedColumn(eigenVector, ed.getV(), tolerance));
    }

    /**
     * Returns true iff there is a column that is a scalar multiple of column
     * in searchMatrix (modulo tolerance)
     */
    private boolean isIncludedColumn(Complex[] column, RealMatrix searchMatrix, double tolerance) {
        boolean found = false;
        int i = 0;
        while (!found && i < searchMatrix.getColumnDimension()) {
            Complex multiplier = Complex.ONE;
            boolean matching = true;
            int j = 0;
            while (matching && j < searchMatrix.getRowDimension()) {
                double colEntry = searchMatrix.getEntry(j, i);
                // Use the first entry where both are non-zero as scalar
                if (multiplier.subtract(1).norm() <= FastMath.ulp(1.0) && FastMath.abs(colEntry) > 1E-14
                        && column[j].norm() > 1e-14) {
                    multiplier = column[j].reciprocal().multiply(colEntry);
                }
                if (column[j].multiply(multiplier).subtract(colEntry).norm() > tolerance) {
                    matching = false;
                }
                j++;
            }
            found = matching;
            i++;
        }
        return found;
    }

    @Before
    public void setUp() {
        double[] real = {
        		0.001, 1.000, 1.001, 2.001, 2.002, 2.003
        };
        refValues = new Complex[real.length];
        for (int i = 0; i < refValues.length; ++i) {
            refValues[i] = new Complex(real[i]);
        }
        matrix = createTestMatrix(new Random(35992629946426l), real);
    }

    @After
    public void tearDown() {
        refValues = null;
        matrix    = null;
    }

    static RealMatrix createTestMatrix(final Random r, final double[] eigenValues) {
        final int n = eigenValues.length;
        final RealMatrix v = createOrthogonalMatrix(r, n);
        final RealMatrix d = MatrixUtils.createRealDiagonalMatrix(eigenValues);
        return v.multiply(d).multiplyTransposed(v);
    }

    public static RealMatrix createOrthogonalMatrix(final Random r, final int size) {

        final double[][] data = new double[size][size];

        for (int i = 0; i < size; ++i) {
            final double[] dataI = data[i];
            double norm2 = 0;
            do {

                // generate randomly row I
                for (int j = 0; j < size; ++j) {
                    dataI[j] = 2 * r.nextDouble() - 1;
                }

                // project the row in the subspace orthogonal to previous rows
                for (int k = 0; k < i; ++k) {
                    final double[] dataK = data[k];
                    double dotProduct = 0;
                    for (int j = 0; j < size; ++j) {
                        dotProduct += dataI[j] * dataK[j];
                    }
                    for (int j = 0; j < size; ++j) {
                        dataI[j] -= dotProduct * dataK[j];
                    }
                }

                // normalize the row
                norm2 = 0;
                for (final double dataIJ : dataI) {
                    norm2 += dataIJ * dataIJ;
                }
                final double inv = 1.0 / FastMath.sqrt(norm2);
                for (int j = 0; j < size; ++j) {
                    dataI[j] *= inv;
                }

            } while (norm2 * size < 0.01);
        }

        return MatrixUtils.createRealMatrix(data);

    }
}
