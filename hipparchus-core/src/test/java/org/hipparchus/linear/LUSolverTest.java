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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class LUSolverTest {
    private double[][] testData = {
            { 1.0, 2.0, 3.0},
            { 2.0, 5.0, 3.0},
            { 1.0, 0.0, 8.0}
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

    /** test threshold impact */
    @Test
    void testThreshold() {
        final RealMatrix matrix = MatrixUtils.createRealMatrix(new double[][] {
                                                       { 1.0, 2.0, 3.0},
                                                       { 2.0, 5.0, 3.0},
                                                       { 4.000001, 9.0, 9.0}
                                                     });
        assertFalse(new LUDecomposer(1.0e-5).decompose(matrix).isNonSingular());
        assertTrue(new LUDecomposer(1.0e-10).decompose(matrix).isNonSingular());
    }

    /** test singular */
    @Test
    void testSingular() {
        DecompositionSolver solver =
            new LUDecomposition(MatrixUtils.createRealMatrix(testData)).getSolver();
        assertTrue(solver.isNonSingular());
        solver = new LUDecomposition(MatrixUtils.createRealMatrix(singular)).getSolver();
        assertFalse(solver.isNonSingular());
        solver = new LUDecomposition(MatrixUtils.createRealMatrix(bigSingular)).getSolver();
        assertFalse(solver.isNonSingular());
    }

    /** test solve dimension errors */
    @Test
    void testSolveDimensionErrors() {
        DecompositionSolver solver =
            new LUDecomposition(MatrixUtils.createRealMatrix(testData)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[2][2]);
        try {
            solver.solve(b);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException iae) {
            // expected behavior
        }
        try {
            solver.solve(b.getColumnVector(0));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException iae) {
            // expected behavior
        }
        try {
            solver.solve(new ArrayRealVectorTest.RealVectorTestImpl(b.getColumn(0)));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException iae) {
            // expected behavior
        }
    }

    /** test solve singularity errors */
    @Test
    void testSolveSingularityErrors() {
        DecompositionSolver solver =
            new LUDecomposition(MatrixUtils.createRealMatrix(singular)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[2][2]);
        try {
            solver.solve(b);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException ime) {
            // expected behavior
        }
        try {
            solver.solve(b.getColumnVector(0));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException ime) {
            // expected behavior
        }
        try {
            solver.solve(new ArrayRealVectorTest.RealVectorTestImpl(b.getColumn(0)));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException ime) {
            // expected behavior
        }
    }

    /** test solve */
    @Test
    void testSolve() {
        DecompositionSolver solver =
            new LUDecomposition(MatrixUtils.createRealMatrix(testData)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[][] {
                { 1, 0 }, { 2, -5 }, { 3, 1 }
        });
        RealMatrix xRef = MatrixUtils.createRealMatrix(new double[][] {
                { 19, -71 }, { -6, 22 }, { -2, 9 }
        });

        // using RealMatrix
        assertEquals(0, solver.solve(b).subtract(xRef).getNorm1(), 1.0e-13);

        // using ArrayRealVector
        for (int i = 0; i < b.getColumnDimension(); ++i) {
            assertEquals(0,
                         solver.solve(b.getColumnVector(i)).subtract(xRef.getColumnVector(i)).getNorm(),
                         1.0e-13);
        }

        // using RealVector with an alternate implementation
        for (int i = 0; i < b.getColumnDimension(); ++i) {
            ArrayRealVectorTest.RealVectorTestImpl v =
                new ArrayRealVectorTest.RealVectorTestImpl(b.getColumn(i));
            assertEquals(0,
                         solver.solve(v).subtract(xRef.getColumnVector(i)).getNorm(),
                         1.0e-13);
        }
    }

    /** test determinant */
    @Test
    void testDeterminant() {
        assertEquals( -1, getDeterminant(MatrixUtils.createRealMatrix(testData)), 1.0e-15);
        assertEquals(-10, getDeterminant(MatrixUtils.createRealMatrix(luData)), 1.0e-14);
        assertEquals(  0, getDeterminant(MatrixUtils.createRealMatrix(singular)), 1.0e-17);
        assertEquals(  0, getDeterminant(MatrixUtils.createRealMatrix(bigSingular)), 1.0e-10);
    }

    private double getDeterminant(RealMatrix m) {
        return new LUDecomposition(m).getDeterminant();
    }
}
