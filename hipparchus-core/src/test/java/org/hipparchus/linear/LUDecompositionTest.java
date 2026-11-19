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
import org.hipparchus.fraction.Fraction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class LUDecompositionTest {
    private double[][] testData = {
            { 1.0, 2.0, 3.0},
            { 2.0, 5.0, 3.0},
            { 1.0, 0.0, 8.0}
    };
    private double[][] testDataMinus = {
            { -1.0, -2.0, -3.0},
            { -2.0, -5.0, -3.0},
            { -1.0,  0.0, -8.0}
    };
    private double[][] luData = {
            { 2.0, 3.0, 3.0 },
            { 0.0, 5.0, 7.0 },
            { 6.0, 9.0, 8.0 }
    };

    // singular matrices
    private double[][] singular = {
            { 2.0, 3.0 },
            { 2.0, 3.0 }
    };
    private double[][] bigSingular = {
            { 1.0, 2.0,   3.0,    4.0 },
            { 2.0, 5.0,   3.0,    4.0 },
            { 7.0, 3.0, 256.0, 1930.0 },
            { 3.0, 7.0,   6.0,    8.0 }
    }; // 4th row = 1st + 2nd

    private static final double entryTolerance = 10e-16;

    private static final double normTolerance = 10e-14;

    /** test dimensions */
    @Test
    void testDimensions() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData);
        LUDecomposition LU = new LUDecomposition(matrix);
        assertEquals(testData.length, LU.getL().getRowDimension());
        assertEquals(testData.length, LU.getL().getColumnDimension());
        assertEquals(testData.length, LU.getU().getRowDimension());
        assertEquals(testData.length, LU.getU().getColumnDimension());
        assertEquals(testData.length, LU.getP().getRowDimension());
        assertEquals(testData.length, LU.getP().getColumnDimension());

    }

    /** test non-square matrix */
    @Test
    void testNonSquare() {
        try {
            new LUDecomposition(MatrixUtils.createRealMatrix(new double[3][2]));
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ime) {
            assertEquals(LocalizedCoreFormats.NON_SQUARE_MATRIX, ime.getSpecifier());
        }
    }

    /** test PA = LU */
    @Test
    void testPAEqualLU() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData);
        LUDecomposition lu = new LUDecomposition(matrix);
        RealMatrix l = lu.getL();
        RealMatrix u = lu.getU();
        RealMatrix p = lu.getP();
        double norm = l.multiply(u).subtract(p.multiply(matrix)).getNorm1();
        assertEquals(0, norm, normTolerance);

        matrix = MatrixUtils.createRealMatrix(testDataMinus);
        lu = new LUDecomposition(matrix);
        l = lu.getL();
        u = lu.getU();
        p = lu.getP();
        norm = l.multiply(u).subtract(p.multiply(matrix)).getNorm1();
        assertEquals(0, norm, normTolerance);

        matrix = MatrixUtils.createRealIdentityMatrix(17);
        lu = new LUDecomposition(matrix);
        l = lu.getL();
        u = lu.getU();
        p = lu.getP();
        norm = l.multiply(u).subtract(p.multiply(matrix)).getNorm1();
        assertEquals(0, norm, normTolerance);

        matrix = MatrixUtils.createRealMatrix(singular);
        lu = new LUDecomposition(matrix);
        assertFalse(lu.getSolver().isNonSingular());
        assertNull(lu.getL());
        assertNull(lu.getU());
        assertNull(lu.getP());

        matrix = MatrixUtils.createRealMatrix(bigSingular);
        lu = new LUDecomposition(matrix);
        assertFalse(lu.getSolver().isNonSingular());
        assertNull(lu.getL());
        assertNull(lu.getU());
        assertNull(lu.getP());

    }

    /** test that L is lower triangular with unit diagonal */
    @Test
    void testLLowerTriangular() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData);
        RealMatrix l = new LUDecomposition(matrix).getL();
        for (int i = 0; i < l.getRowDimension(); i++) {
            assertEquals(1, l.getEntry(i, i), entryTolerance);
            for (int j = i + 1; j < l.getColumnDimension(); j++) {
                assertEquals(0, l.getEntry(i, j), entryTolerance);
            }
        }
    }

    /** test that U is upper triangular */
    @Test
    void testUUpperTriangular() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData);
        RealMatrix u = new LUDecomposition(matrix).getU();
        for (int i = 0; i < u.getRowDimension(); i++) {
            for (int j = 0; j < i; j++) {
                assertEquals(0, u.getEntry(i, j), entryTolerance);
            }
        }
    }

    /** test that P is a permutation matrix */
    @Test
    void testPPermutation() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData);
        RealMatrix p   = new LUDecomposition(matrix).getP();

        RealMatrix ppT = p.multiplyTransposed(p);
        RealMatrix id  = MatrixUtils.createRealIdentityMatrix(p.getRowDimension());
        assertEquals(0, ppT.subtract(id).getNorm1(), normTolerance);

        for (int i = 0; i < p.getRowDimension(); i++) {
            int zeroCount  = 0;
            int oneCount   = 0;
            int otherCount = 0;
            for (int j = 0; j < p.getColumnDimension(); j++) {
                final double e = p.getEntry(i, j);
                if (e == 0) {
                    ++zeroCount;
                } else if (e == 1) {
                    ++oneCount;
                } else {
                    ++otherCount;
                }
            }
            assertEquals(p.getColumnDimension() - 1, zeroCount);
            assertEquals(1, oneCount);
            assertEquals(0, otherCount);
        }

        for (int j = 0; j < p.getColumnDimension(); j++) {
            int zeroCount  = 0;
            int oneCount   = 0;
            int otherCount = 0;
            for (int i = 0; i < p.getRowDimension(); i++) {
                final double e = p.getEntry(i, j);
                if (e == 0) {
                    ++zeroCount;
                } else if (e == 1) {
                    ++oneCount;
                } else {
                    ++otherCount;
                }
            }
            assertEquals(p.getRowDimension() - 1, zeroCount);
            assertEquals(1, oneCount);
            assertEquals(0, otherCount);
        }

    }

    /** test singular */
    @Test
    void testSingular() {
        final RealMatrix m = MatrixUtils.createRealMatrix(testData);
        LUDecomposition lu = new LUDecomposition(m);
        assertTrue(lu.getSolver().isNonSingular());
        assertEquals(-1.0, lu.getDeterminant(), 1.0e-15);
        lu = new LUDecomposition(m.getSubMatrix(0, 1, 0, 1));
        assertTrue(lu.getSolver().isNonSingular());
        assertEquals(+1.0, lu.getDeterminant(), 1.0e-15);
        lu = new LUDecomposition(MatrixUtils.createRealMatrix(singular));
        assertFalse(lu.getSolver().isNonSingular());
        assertEquals(0.0, lu.getDeterminant(), 1.0e-15);
        lu = new LUDecomposition(MatrixUtils.createRealMatrix(bigSingular));
        assertFalse(lu.getSolver().isNonSingular());
        assertEquals(0.0, lu.getDeterminant(), 1.0e-15);
        try {
            lu.getSolver().solve(new ArrayRealVector(new double[] { 1, 1, 1, 1 }));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.SINGULAR_MATRIX, miae.getSpecifier());
        }
    }

    /** test matrices values */
    @Test
    void testMatricesValues1() {
       LUDecomposition lu =
            new LUDecomposition(MatrixUtils.createRealMatrix(testData));
        RealMatrix lRef = MatrixUtils.createRealMatrix(new double[][] {
                { 1.0, 0.0, 0.0 },
                { 0.5, 1.0, 0.0 },
                { 0.5, 0.2, 1.0 }
        });
        RealMatrix uRef = MatrixUtils.createRealMatrix(new double[][] {
                { 2.0,  5.0, 3.0 },
                { 0.0, -2.5, 6.5 },
                { 0.0,  0.0, 0.2 }
        });
        RealMatrix pRef = MatrixUtils.createRealMatrix(new double[][] {
                { 0.0, 1.0, 0.0 },
                { 0.0, 0.0, 1.0 },
                { 1.0, 0.0, 0.0 }
        });
        int[] pivotRef = { 1, 2, 0 };

        // check values against known references
        RealMatrix l = lu.getL();
        assertEquals(0, l.subtract(lRef).getNorm1(), 1.0e-13);
        RealMatrix u = lu.getU();
        assertEquals(0, u.subtract(uRef).getNorm1(), 1.0e-13);
        RealMatrix p = lu.getP();
        assertEquals(0, p.subtract(pRef).getNorm1(), 1.0e-13);
        int[] pivot = lu.getPivot();
        for (int i = 0; i < pivotRef.length; ++i) {
            assertEquals(pivotRef[i], pivot[i]);
        }

        // check the same cached instance is returned the second time
        assertSame(l, lu.getL());
        assertSame(u, lu.getU());
        assertSame(p, lu.getP());

    }

    /** test matrices values */
    @Test
    void testMatricesValues2() {
       LUDecomposition lu =
            new LUDecomposition(MatrixUtils.createRealMatrix(luData));
        RealMatrix lRef = MatrixUtils.createRealMatrix(new double[][] {
                {    1.0,    0.0, 0.0 },
                {    0.0,    1.0, 0.0 },
                { 1.0 / 3.0, 0.0, 1.0 }
        });
        RealMatrix uRef = MatrixUtils.createRealMatrix(new double[][] {
                { 6.0, 9.0,    8.0    },
                { 0.0, 5.0,    7.0    },
                { 0.0, 0.0, 1.0 / 3.0 }
        });
        RealMatrix pRef = MatrixUtils.createRealMatrix(new double[][] {
                { 0.0, 0.0, 1.0 },
                { 0.0, 1.0, 0.0 },
                { 1.0, 0.0, 0.0 }
        });
        int[] pivotRef = { 2, 1, 0 };

        // check values against known references
        RealMatrix l = lu.getL();
        assertEquals(0, l.subtract(lRef).getNorm1(), 1.0e-13);
        RealMatrix u = lu.getU();
        assertEquals(0, u.subtract(uRef).getNorm1(), 1.0e-13);
        RealMatrix p = lu.getP();
        assertEquals(0, p.subtract(pRef).getNorm1(), 1.0e-13);
        int[] pivot = lu.getPivot();
        for (int i = 0; i < pivotRef.length; ++i) {
            assertEquals(pivotRef[i], pivot[i]);
        }

        // check the same cached instance is returned the second time
        assertSame(l, lu.getL());
        assertSame(u, lu.getU());
        assertSame(p, lu.getP());
    }

    @Test
    void testSolve() {
        LUDecomposition lu =
                        new LUDecomposition(new Array2DRowRealMatrix(testData));
        DecompositionSolver solver = lu.getSolver();
        RealVector solution = solver.solve(new ArrayRealVector(new double[] {
            new Fraction(1, 2).doubleValue(), new Fraction(2, 3).doubleValue(), new Fraction(3,4).doubleValue()
        }));
        assertEquals(testData.length, solution.getDimension());
        assertEquals(new Fraction(-31, 12).doubleValue(), solution.getEntry(0), 1.0e-14);
        assertEquals(new Fraction( 11, 12).doubleValue(), solution.getEntry(1), 1.0e-14);
        assertEquals(new Fraction(  5, 12).doubleValue(), solution.getEntry(2), 1.0e-14);
        assertEquals(testData.length,    solver.getRowDimension());
        assertEquals(testData[0].length, solver.getColumnDimension());
        try {
            solver.solve(new ArrayRealVector(new double[] { 1, 1 }));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
    }

}
