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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

public class EigenDecompositionSymmetricTest {

    private double[] refValues;
    private RealMatrix matrix;

    @Test
    void testDimension1() {
        RealMatrix matrix =
            MatrixUtils.createRealMatrix(new double[][] { { 1.5 } });
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        assertEquals(1.5, ed.getEigenvalue(0), 1.0e-15);
    }

    @Test
    void testDimension2() {
        RealMatrix matrix =
            MatrixUtils.createRealMatrix(new double[][] {
                    { 59.0, 12.0 },
                    { 12.0, 66.0 }
            });
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        assertEquals(75.0, ed.getEigenvalue(0), 1.0e-15);
        assertEquals(50.0, ed.getEigenvalue(1), 1.0e-15);
    }

    @Test
    void testDimension3() {
        RealMatrix matrix =
            MatrixUtils.createRealMatrix(new double[][] {
                                   {  39632.0, -4824.0, -16560.0 },
                                   {  -4824.0,  8693.0,   7920.0 },
                                   { -16560.0,  7920.0,  17300.0 }
                               });
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        assertEquals(50000.0, ed.getEigenvalue(0), 3.0e-11);
        assertEquals(12500.0, ed.getEigenvalue(1), 3.0e-11);
        assertEquals( 3125.0, ed.getEigenvalue(2), 3.0e-11);
    }

    @Test
    void testDimension3MultipleRoot() {
        RealMatrix matrix =
            MatrixUtils.createRealMatrix(new double[][] {
                    {  5,   10,   15 },
                    { 10,   20,   30 },
                    { 15,   30,   45 }
            });
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        assertEquals(70.0, ed.getEigenvalue(0), 3.0e-11);
        assertEquals(0.0,  ed.getEigenvalue(1), 3.0e-11);
        assertEquals(0.0,  ed.getEigenvalue(2), 3.0e-11);
        assertEquals(matrix.getRowDimension(),    ed.getSolver().getRowDimension());
        assertEquals(matrix.getColumnDimension(), ed.getSolver().getColumnDimension());
    }

    @Test
    void testDimension4WithSplit() {
        RealMatrix matrix =
            MatrixUtils.createRealMatrix(new double[][] {
                                   {  0.784, -0.288,  0.000,  0.000 },
                                   { -0.288,  0.616,  0.000,  0.000 },
                                   {  0.000,  0.000,  0.164, -0.048 },
                                   {  0.000,  0.000, -0.048,  0.136 }
                               });
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        assertEquals(1.0, ed.getEigenvalue(0), 1.0e-15);
        assertEquals(0.4, ed.getEigenvalue(1), 1.0e-15);
        assertEquals(0.2, ed.getEigenvalue(2), 1.0e-15);
        assertEquals(0.1, ed.getEigenvalue(3), 1.0e-15);
        assertEquals(matrix.getRowDimension(),    ed.getSolver().getRowDimension());
        assertEquals(matrix.getColumnDimension(), ed.getSolver().getColumnDimension());
    }

    @Test
    void testDimension4WithoutSplit() {
        RealMatrix matrix =
            MatrixUtils.createRealMatrix(new double[][] {
                                   {  0.5608, -0.2016,  0.1152, -0.2976 },
                                   { -0.2016,  0.4432, -0.2304,  0.1152 },
                                   {  0.1152, -0.2304,  0.3088, -0.1344 },
                                   { -0.2976,  0.1152, -0.1344,  0.3872 }
                               });
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        assertEquals(1.0, ed.getEigenvalue(0), 1.0e-15);
        assertEquals(0.4, ed.getEigenvalue(1), 1.0e-15);
        assertEquals(0.2, ed.getEigenvalue(2), 1.0e-15);
        assertEquals(0.1, ed.getEigenvalue(3), 1.0e-15);
        assertEquals(matrix.getRowDimension(),    ed.getSolver().getRowDimension());
        assertEquals(matrix.getColumnDimension(), ed.getSolver().getColumnDimension());
    }

    @Test
    void testMath308() {

        double[] mainTridiagonal = {
            22.330154644539597, 46.65485522478641, 17.393672330044705, 54.46687435351116, 80.17800767709437
        };
        double[] secondaryTridiagonal = {
            13.04450406501361, -5.977590941539671, 2.9040909856707517, 7.1570352792841225
        };

        // the reference values have been computed using routine DSTEMR
        // from the fortran library LAPACK version 3.2.1
        double[] refEigenValues = {
            82.044413207204002, 53.456697699894512, 52.536278520113882, 18.847969733754262, 14.138204224043099
        };
        RealVector[] refEigenVectors = {
            new ArrayRealVector(new double[] { -0.000462690386766, -0.002118073109055,  0.011530080757413,  0.252322434584915,  0.967572088232592 }),
            new ArrayRealVector(new double[] {  0.314647769490148,  0.750806415553905, -0.167700312025760, -0.537092972407375,  0.143854968127780 }),
            new ArrayRealVector(new double[] {  0.222368839324646,  0.514921891363332, -0.021377019336614,  0.801196801016305, -0.207446991247740 }),
            new ArrayRealVector(new double[] { -0.713933751051495,  0.190582113553930, -0.671410443368332,  0.056056055955050, -0.006541576993581 }),
            new ArrayRealVector(new double[] { -0.584677060845929,  0.367177264979103,  0.721453187784497, -0.052971054621812,  0.005740715188257 })
        };

        EigenDecompositionSymmetric decomposition;
        decomposition = new EigenDecompositionSymmetric(mainTridiagonal, secondaryTridiagonal);

        double[] eigenValues = decomposition.getEigenvalues();
        for (int i = 0; i < refEigenValues.length; ++i) {
            assertEquals(refEigenValues[i], eigenValues[i], 1.0e-5);
            assertEquals(0, refEigenVectors[i].subtract(decomposition.getEigenvector(i)).getNorm(), 2.0e-7);
        }

        assertEquals(mainTridiagonal.length, decomposition.getSolver().getRowDimension());
        assertEquals(mainTridiagonal.length, decomposition.getSolver().getColumnDimension());
    }

    @Test
    void testMathpbx02() {

        double[] mainTridiagonal = {
              7484.860960227216, 18405.28129035345, 13855.225609560746,
             10016.708722343366, 559.8117399576674, 6750.190788301587,
                71.21428769782159
        };
        double[] secondaryTridiagonal = {
             -4175.088570476366,1975.7955858241994,5193.178422374075,
              1995.286659169179,75.34535882933804,-234.0808002076056
        };

        // the reference values have been computed using routine DSTEMR
        // from the fortran library LAPACK version 3.2.1
        double[] refEigenValues = {
                20654.744890306974412,16828.208208485466457,
                6893.155912634994820,6757.083016675340332,
                5887.799885688558788,64.309089923240379,
                57.992628792736340
        };
        RealVector[] refEigenVectors = {
                new ArrayRealVector(new double[] {-0.270356342026904, 0.852811091326997, 0.399639490702077, 0.198794657813990, 0.019739323307666, 0.000106983022327, -0.000001216636321}),
                new ArrayRealVector(new double[] {0.179995273578326,-0.402807848153042,0.701870993525734,0.555058211014888,0.068079148898236,0.000509139115227,-0.000007112235617}),
                new ArrayRealVector(new double[] {-0.399582721284727,-0.056629954519333,-0.514406488522827,0.711168164518580,0.225548081276367,0.125943999652923,-0.004321507456014}),
                new ArrayRealVector(new double[] {0.058515721572821,0.010200130057739,0.063516274916536,-0.090696087449378,-0.017148420432597,0.991318870265707,-0.034707338554096}),
                new ArrayRealVector(new double[] {0.855205995537564,0.327134656629775,-0.265382397060548,0.282690729026706,0.105736068025572,-0.009138126622039,0.000367751821196}),
                new ArrayRealVector(new double[] {-0.002913069901144,-0.005177515777101,0.041906334478672,-0.109315918416258,0.436192305456741,0.026307315639535,0.891797507436344}),
                new ArrayRealVector(new double[] {-0.005738311176435,-0.010207611670378,0.082662420517928,-0.215733886094368,0.861606487840411,-0.025478530652759,-0.451080697503958})
        };

        // the following line triggers the exception
        EigenDecompositionSymmetric decomposition;
        decomposition = new EigenDecompositionSymmetric(mainTridiagonal, secondaryTridiagonal);

        double[] eigenValues = decomposition.getEigenvalues();
        for (int i = 0; i < refEigenValues.length; ++i) {
            assertEquals(refEigenValues[i], eigenValues[i], 1.0e-3);
            if (refEigenVectors[i].dotProduct(decomposition.getEigenvector(i)) < 0) {
                assertEquals(0, refEigenVectors[i].add(decomposition.getEigenvector(i)).getNorm(), 1.0e-5);
            } else {
                assertEquals(0, refEigenVectors[i].subtract(decomposition.getEigenvector(i)).getNorm(), 1.0e-5);
            }
        }

    }

    @Test
    void testMathpbx03() {

        double[] mainTridiagonal = {
            1809.0978259647177,3395.4763425956166,1832.1894584712693,3804.364873592377,
            806.0482458637571,2403.656427234185,28.48691431556015
        };
        double[] secondaryTridiagonal = {
            -656.8932064545833,-469.30804108920734,-1021.7714889369421,
            -1152.540497328983,-939.9765163817368,-12.885877015422391
        };

        // the reference values have been computed using routine DSTEMR
        // from the fortran library LAPACK version 3.2.1
        double[] refEigenValues = {
            4603.121913685183245,3691.195818048970978,2743.442955402465032,1657.596442107321764,
            1336.797819095331306,30.129865209677519,17.035352085224986
        };

        RealVector[] refEigenVectors = {
            new ArrayRealVector(new double[] {-0.036249830202337,0.154184732411519,-0.346016328392363,0.867540105133093,-0.294483395433451,0.125854235969548,-0.000354507444044}),
            new ArrayRealVector(new double[] {-0.318654191697157,0.912992309960507,-0.129270874079777,-0.184150038178035,0.096521712579439,-0.070468788536461,0.000247918177736}),
            new ArrayRealVector(new double[] {-0.051394668681147,0.073102235876933,0.173502042943743,-0.188311980310942,-0.327158794289386,0.905206581432676,-0.004296342252659}),
            new ArrayRealVector(new double[] {0.838150199198361,0.193305209055716,-0.457341242126146,-0.166933875895419,0.094512811358535,0.119062381338757,-0.000941755685226}),
            new ArrayRealVector(new double[] {0.438071395458547,0.314969169786246,0.768480630802146,0.227919171600705,-0.193317045298647,-0.170305467485594,0.001677380536009}),
            new ArrayRealVector(new double[] {-0.003726503878741,-0.010091946369146,-0.067152015137611,-0.113798146542187,-0.313123000097908,-0.118940107954918,0.932862311396062}),
            new ArrayRealVector(new double[] {0.009373003194332,0.025570377559400,0.170955836081348,0.291954519805750,0.807824267665706,0.320108347088646,0.360202112392266}),
        };

        // the following line triggers the exception
        EigenDecompositionSymmetric decomposition;
        decomposition = new EigenDecompositionSymmetric(mainTridiagonal, secondaryTridiagonal);

        double[] eigenValues = decomposition.getEigenvalues();
        for (int i = 0; i < refEigenValues.length; ++i) {
            assertEquals(refEigenValues[i], eigenValues[i], 1.0e-4);
            if (refEigenVectors[i].dotProduct(decomposition.getEigenvector(i)) < 0) {
                assertEquals(0, refEigenVectors[i].add(decomposition.getEigenvector(i)).getNorm(), 1.0e-5);
            } else {
                assertEquals(0, refEigenVectors[i].subtract(decomposition.getEigenvector(i)).getNorm(), 1.0e-5);
            }
        }

    }

    /** test a matrix already in tridiagonal form. */
    @Test
    void testTridiagonal() {
        Random r = new Random(4366663527842l);
        double[] ref = new double[30];
        for (int i = 0; i < ref.length; ++i) {
            if (i < 5) {
                ref[i] = 2 * r.nextDouble() - 1;
            } else {
                ref[i] = 0.0001 * r.nextDouble() + 6;
            }
        }
        Arrays.sort(ref);
        TriDiagonalTransformer t =
            new TriDiagonalTransformer(createTestMatrix(r, ref));
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(t.getMainDiagonalRef(), t.getSecondaryDiagonalRef());
        double[] eigenValues = ed.getEigenvalues();
        assertEquals(ref.length, eigenValues.length);
        for (int i = 0; i < ref.length; ++i) {
            assertEquals(ref[ref.length - i - 1], eigenValues[i], 2.0e-14);
        }

    }

    @Test
    void testGitHubIssue30() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {
            { 23473.684554963584, 4273.093076392109 },
            { 4273.093076392048,  4462.13956661408  }
        });
        EigenDecompositionSymmetric ed = new EigenDecompositionSymmetric(m, 1.0e-13, true);
        assertEquals(0.0,
                            ed.getV().multiply(ed.getD()).multiply(ed.getVT()).subtract(m).getNorm1(),
                            1.0e-10);
    }

    /** test dimensions */
    @Test
    void testDimensions() {
        final int m = matrix.getRowDimension();
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        assertEquals(m, ed.getV().getRowDimension());
        assertEquals(m, ed.getV().getColumnDimension());
        assertEquals(m, ed.getD().getColumnDimension());
        assertEquals(m, ed.getD().getColumnDimension());
        assertEquals(m, ed.getVT().getRowDimension());
        assertEquals(m, ed.getVT().getColumnDimension());
    }

    /** test eigenvalues */
    @Test
    void testEigenvalues() {
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        double[] eigenValues = ed.getEigenvalues();
        assertEquals(refValues.length, eigenValues.length);
        for (int i = 0; i < refValues.length; ++i) {
            assertEquals(refValues[i], eigenValues[i], 3.0e-15);
        }
    }

    @Test
    void testSymmetric() {
        RealMatrix symmetric = MatrixUtils.createRealMatrix(new double[][] {
                {4, 1, 1},
                {1, 2, 3},
                {1, 3, 6}
        });

        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(symmetric);

        DiagonalMatrix d  = ed.getD();
        RealMatrix     v  = ed.getV();
        RealMatrix     vT = ed.getVT();

        double norm = v.multiply(d).multiply(vT).subtract(symmetric).getNorm1();
        assertEquals(0, norm, 6.0e-13);
    }

    @Test
    void testSquareRoot() {
        final double[][] data = {
            { 33, 24,  7 },
            { 24, 57, 11 },
            {  7, 11,  9 }
        };

        final EigenDecompositionSymmetric dec = new EigenDecompositionSymmetric(MatrixUtils.createRealMatrix(data));
        final RealMatrix sqrtM = dec.getSquareRoot();

        // Reconstruct initial matrix.
        final RealMatrix m = sqrtM.multiply(sqrtM);

        final int dim = data.length;
        for (int r = 0; r < dim; r++) {
            for (int c = 0; c < dim; c++) {
                assertEquals(data[r][c], m.getEntry(r, c), 1e-13, "m[" + r + "][" + c + "]");
            }
        }
    }

    @Test
    void testSquareRootNonSymmetric() {
        assertThrows(MathRuntimeException.class, () -> {
            final double[][] data = {
                {1, 2, 4},
                {2, 3, 5},
                {11, 5, 9}
            };

            final EigenDecompositionSymmetric dec = new EigenDecompositionSymmetric(MatrixUtils.createRealMatrix(data));
            @SuppressWarnings("unused")
            final RealMatrix sqrtM = dec.getSquareRoot();
        });
    }

    @Test
    void testSquareRootNonPositiveDefinite() {
        assertThrows(MathRuntimeException.class, () -> {
            final double[][] data = {
                {1, 2, 4},
                {2, 3, 5},
                {4, 5, -9}
            };

            final EigenDecompositionSymmetric dec = new EigenDecompositionSymmetric(MatrixUtils.createRealMatrix(data));
            @SuppressWarnings("unused")
            final RealMatrix sqrtM = dec.getSquareRoot();
        });
    }

    @Test
    void testIsNonSingularTinyOutOfOrderEigenvalue() {
        try {
            new EigenDecompositionSymmetric(MatrixUtils.createRealMatrix(new double[][] {
                { 1e-13, 0 },
                { 1,     1 },
            }));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NON_SYMMETRIC_MATRIX, miae.getSpecifier());
        }
    }

    /** test eigenvectors */
    @Test
    void testEigenvectors() {
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        for (int i = 0; i < matrix.getRowDimension(); ++i) {
            double lambda = ed.getEigenvalue(i);
            RealVector v  = ed.getEigenvector(i);
            RealVector mV = matrix.operate(v);
            assertEquals(0, mV.subtract(v.mapMultiplyToSelf(lambda)).getNorm(), 1.0e-13);
        }
    }

    /** test A = VDVt */
    @Test
    void testAEqualVDVt() {
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(matrix);
        RealMatrix     v  = ed.getV();
        DiagonalMatrix d  = ed.getD();
        RealMatrix     vT = ed.getVT();
        double norm = v.multiply(d).multiply(vT).subtract(matrix).getNorm1();
        assertEquals(0, norm, 6.0e-13);
    }

    /** test that V is orthogonal */
    @Test
    void testVOrthogonal() {
        RealMatrix v = new EigenDecompositionSymmetric(matrix).getV();
        RealMatrix vTv = v.transposeMultiply(v);
        RealMatrix id  = MatrixUtils.createRealIdentityMatrix(vTv.getRowDimension());
        assertEquals(0, vTv.subtract(id).getNorm1(), 2.0e-13);
    }

    /** test diagonal matrix */
    @Test
    void testDiagonal() {
        double[] diagonal = new double[] { -3.0, -2.0, 2.0, 5.0 };
        RealMatrix m = MatrixUtils.createRealDiagonalMatrix(diagonal);
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(m);
        assertEquals(diagonal[0], ed.getEigenvalue(3), 2.0e-15);
        assertEquals(diagonal[1], ed.getEigenvalue(2), 2.0e-15);
        assertEquals(diagonal[2], ed.getEigenvalue(1), 2.0e-15);
        assertEquals(diagonal[3], ed.getEigenvalue(0), 2.0e-15);
    }

    /**
     * Matrix with eigenvalues {8, -1, -1}
     */
    @Test
    void testRepeatedEigenvalue() {
        RealMatrix repeated = MatrixUtils.createRealMatrix(new double[][] {
                {3,  2,  4},
                {2,  0,  2},
                {4,  2,  3}
        });
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(repeated);
        checkEigenValues((new double[] {8, -1, -1}), ed, 1E-12);
        checkEigenVector((new double[] {2, 1, 2}), ed, 1E-12);
    }

    /**
     * Matrix with eigenvalues {2, 0, 12}
     */
    @Test
    void testDistinctEigenvalues() {
        RealMatrix distinct = MatrixUtils.createRealMatrix(new double[][] {
                {3, 1, -4},
                {1, 3, -4},
                {-4, -4, 8}
        });
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(distinct);
        checkEigenValues((new double[] {2, 0, 12}), ed, 1E-12);
        checkEigenVector((new double[] {1, -1, 0}), ed, 1E-12);
        checkEigenVector((new double[] {1, 1, 1}), ed, 1E-12);
        checkEigenVector((new double[] {-1, -1, 2}), ed, 1E-12);
    }

    /**
     * Verifies operation on indefinite matrix
     */
    @Test
    void testZeroDivide() {
        RealMatrix indefinite = MatrixUtils.createRealMatrix(new double [][] {
                { 0.0, 1.0, -1.0 },
                { 1.0, 1.0, 0.0 },
                { -1.0,0.0, 1.0 }
        });
        EigenDecompositionSymmetric ed;
        ed = new EigenDecompositionSymmetric(indefinite);
        checkEigenValues((new double[] {2, 1, -1}), ed, 1E-12);
        double isqrt3 = 1/FastMath.sqrt(3.0);
        checkEigenVector((new double[] {isqrt3,isqrt3,-isqrt3}), ed, 1E-12);
        double isqrt2 = 1/FastMath.sqrt(2.0);
        checkEigenVector((new double[] {0.0,-isqrt2,-isqrt2}), ed, 1E-12);
        double isqrt6 = 1/FastMath.sqrt(6.0);
        checkEigenVector((new double[] {2*isqrt6,-isqrt6,isqrt6}), ed, 1E-12);
    }

    /**
     * Verifies operation on very small values.
     * Matrix with eigenvalues {2e-100, 0, 12e-100}
     */
    @Test
    void testTinyValues() {
        final double tiny = 1e-100;
        RealMatrix distinct = MatrixUtils.createRealMatrix(new double[][] {
                {3, 1, -4},
                {1, 3, -4},
                {-4, -4, 8}
        });
        distinct = distinct.scalarMultiply(tiny);

        final EigenDecompositionSymmetric ed = new EigenDecompositionSymmetric(distinct);
        checkEigenValues(MathArrays.scale(tiny, new double[] {2, 0, 12}), ed, 1e-12 * tiny);
        checkEigenVector(new double[] {1, -1, 0}, ed, 1e-12);
        checkEigenVector(new double[] {1, 1, 1}, ed, 1e-12);
        checkEigenVector(new double[] {-1, -1, 2}, ed, 1e-12);
    }

    /**
     * Verifies that a custom epsilon value is used when testing for singular
     */
    @Test
    void testCustomEpsilon() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(new double[][] {
            {  1.76208738E-13, -9.37625373E-13, -1.94760551E-12, -2.56572222E-11, -7.30093964E-11, -1.98340808E-09 },
            { -9.37625373E-13,  5.00812620E-12,  1.06017205E-11,  1.40431472E-10,  3.62452521E-10,  1.05830167E-08 },
            { -1.94760551E-12,  1.06017205E-11,  3.15658331E-11,  2.32155752E-09, -1.53067748E-09,  2.23110293E-08 },
            { -2.56572222E-11,  1.40431472E-10,  2.32155752E-09,  8.81161492E-07, -8.70304198E-07,  2.93564832E-07 },
            { -7.30093964E-11,  3.62452521E-10, -1.53067748E-09, -8.70304198E-07,  9.42413982E-07,  7.81029359E-07 },
            { -1.98340808E-09,  1.05830167E-08,  2.23110293E-08,  2.93564832E-07,  7.81029359E-07,  2.23721205E-05 }
        });

        final EigenDecompositionSymmetric defaultEd = new EigenDecompositionSymmetric(matrix);
        assertFalse(defaultEd.getSolver().isNonSingular());

        final double customEpsilon = 1e-20;
        final EigenDecompositionSymmetric customEd = new EigenDecompositionSymmetric(matrix, customEpsilon, false);
        assertTrue(customEd.getSolver().isNonSingular());
        assertEquals(customEpsilon, customEd.getEpsilon(), 1.0e-25);
    }

    @Test
    void testSingular() {
        final RealMatrix m = MatrixUtils.createRealMatrix(new double[][] { { 1, 0 }, { 0, 0 } });
        DecompositionSolver solver = new EigenDecompositionSymmetric(m).getSolver();
        try {
            solver.solve(new ArrayRealVector(new double[] { 1.0, 1.0 }));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.SINGULAR_MATRIX, miae.getSpecifier());
        }
        try {
            solver.solve(new Array2DRowRealMatrix(new double[][] { { 1.0, 1.0 }, { 1.0, 2.0 } }));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.SINGULAR_MATRIX, miae.getSpecifier());
        }
    }

    @Test
    void testIncreasingOrder() {

        final RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {
            { 81.0, 63.0, 55.0, 49.0, 0.0, 0.0},
            { 63.0, 82.0, 80.0, 69.0, 0.0, 0.0},
            { 55.0, 80.0, 92.0, 75.0, 0.0, 0.0},
            { 49.0, 69.0, 75.0, 73.0, 0.0, 0.0},
            {  0.0,  0.0,  0.0,  0.0, 0.0, 0.0},
            {  0.0,  0.0,  0.0,  0.0, 0.0, 0.0}
        });

        EigenDecompositionSymmetric ed = new EigenDecompositionSymmetric(m,
                                                                         EigenDecompositionSymmetric.DEFAULT_EPSILON,
                                                                         false);

        assertEquals(  0.0,               ed.getD().getEntry(0, 0), 1.0e-14);
        assertEquals(  0.0,               ed.getD().getEntry(1, 1), 1.0e-14);
        assertEquals(  4.793253233672134, ed.getD().getEntry(2, 2), 1.0e-14);
        assertEquals(  7.429392849462756, ed.getD().getEntry(3, 3), 1.0e-14);
        assertEquals( 36.43053556404571,  ed.getD().getEntry(4, 4), 1.0e-14);
        assertEquals(279.34681835281935,  ed.getD().getEntry(5, 5), 1.0e-14);

    }

    /**
     * Verifies that the given EigenDecomposition has eigenvalues equivalent to
     * the targetValues, ignoring the order of the values and allowing
     * values to differ by tolerance.
     */
    protected void checkEigenValues(double[] targetValues,
                                    EigenDecompositionSymmetric ed, double tolerance) {
        double[] observed = ed.getEigenvalues();
        for (int i = 0; i < observed.length; i++) {
            assertTrue(isIncludedValue(observed[i], targetValues, tolerance));
            assertTrue(isIncludedValue(targetValues[i], observed, tolerance));
        }
    }


    /**
     * Returns true iff there is an entry within tolerance of value in
     * searchArray.
     */
    private boolean isIncludedValue(double value, double[] searchArray,
                                    double tolerance) {
       boolean found = false;
       int i = 0;
       while (!found && i < searchArray.length) {
           if (FastMath.abs(value - searchArray[i]) < tolerance) {
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
    protected void checkEigenVector(double[] eigenVector,
                                    EigenDecompositionSymmetric ed, double tolerance) {
        assertTrue(isIncludedColumn(eigenVector, ed.getV(), tolerance));
    }

    /**
     * Returns true iff there is a column that is a scalar multiple of column
     * in searchMatrix (modulo tolerance)
     */
    private boolean isIncludedColumn(double[] column, RealMatrix searchMatrix,
                                     double tolerance) {
        boolean found = false;
        int i = 0;
        while (!found && i < searchMatrix.getColumnDimension()) {
            double multiplier = 1.0;
            boolean matching = true;
            int j = 0;
            while (matching && j < searchMatrix.getRowDimension()) {
                double colEntry = searchMatrix.getEntry(j, i);
                // Use the first entry where both are non-zero as scalar
                if (FastMath.abs(multiplier - 1.0) <= FastMath.ulp(1.0) && FastMath.abs(colEntry) > 1E-14
                        && FastMath.abs(column[j]) > 1e-14) {
                    multiplier = colEntry / column[j];
                }
                if (FastMath.abs(column[j] * multiplier - colEntry) > tolerance) {
                    matching = false;
                }
                j++;
            }
            found = matching;
            i++;
        }
        return found;
    }

    @BeforeEach
    void setUp() {
        refValues = new double[] {
                2.003, 2.002, 2.001, 1.001, 1.000, 0.001
        };
        matrix = createTestMatrix(new Random(35992629946426l), refValues);
    }

    @AfterEach
    void tearDown() {
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
