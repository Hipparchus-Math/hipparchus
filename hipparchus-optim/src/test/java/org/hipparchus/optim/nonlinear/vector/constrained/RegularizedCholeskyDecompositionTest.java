/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.DecompositionSolver;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class RegularizedCholeskyDecompositionTest {

    @Test
    public void testNonSquare() {
        try {
            new RegularizedCholeskyDecomposition(MatrixUtils.createRealMatrix(6, 5));
            Assertions.fail("an exception should have been thrown");
        }  catch (final MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.NON_SQUARE_MATRIX, e.getSpecifier());
        }
    }

    @Test
    public void testNotSymmetric() {
        try {
            final RealMatrix matrix = MatrixUtils.createRealIdentityMatrix(4);
            matrix.setEntry(1, 2, 2.5);
            new RegularizedCholeskyDecomposition(matrix);
            Assertions.fail("an exception should have been thrown");
        }  catch (final MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.NON_SYMMETRIC_MATRIX, e.getSpecifier());
        }
    }

    @Test
    public void testGetters() {
        final RealMatrix matrix = MatrixUtils.createRealMatrix(new double[][] {
                {  1.0,  1.0,  1.0,  1.0 },
                {  1.0,  5.0,  5.0,  5.0 },
                {  1.0,  5.0, 14.0, 14.0 },
                {  1.0,  5.0, 14.0, 15.0 }
        });
        RegularizedCholeskyDecomposition decomposition = new RegularizedCholeskyDecomposition(matrix);


        Assertions.assertEquals(36.0, decomposition.getDeterminant(),     1.0e-12);

        final RealMatrix l = decomposition.getL();
        Assertions.assertEquals(1.0, l.getEntry(0, 0), 1.0e-12);
        Assertions.assertEquals(1.0, l.getEntry(1, 0), 1.0e-12);
        Assertions.assertEquals(2.0, l.getEntry(1, 1), 1.0e-12);
        Assertions.assertEquals(1.0, l.getEntry(2, 0), 1.0e-12);
        Assertions.assertEquals(2.0, l.getEntry(2, 1), 1.0e-12);
        Assertions.assertEquals(3.0, l.getEntry(2, 2), 1.0e-12);
        Assertions.assertEquals(1.0, l.getEntry(3, 0), 1.0e-12);
        Assertions.assertEquals(2.0, l.getEntry(3, 1), 1.0e-12);
        Assertions.assertEquals(3.0, l.getEntry(3, 2), 1.0e-12);
        Assertions.assertEquals(1.0, l.getEntry(3, 3), 1.0e-12);
        Assertions.assertSame(l, decomposition.getL());

        final RealMatrix lt = decomposition.getLT();
        Assertions.assertEquals(1.0, lt.getEntry(0, 0), 1.0e-12);
        Assertions.assertEquals(1.0, lt.getEntry(0, 1), 1.0e-12);
        Assertions.assertEquals(2.0, lt.getEntry(1, 1), 1.0e-12);
        Assertions.assertEquals(1.0, lt.getEntry(0, 2), 1.0e-12);
        Assertions.assertEquals(2.0, lt.getEntry(1, 2), 1.0e-12);
        Assertions.assertEquals(3.0, lt.getEntry(2, 2), 1.0e-12);
        Assertions.assertEquals(1.0, lt.getEntry(0, 3), 1.0e-12);
        Assertions.assertEquals(2.0, lt.getEntry(1, 3), 1.0e-12);
        Assertions.assertEquals(3.0, lt.getEntry(2, 3), 1.0e-12);
        Assertions.assertEquals(1.0, lt.getEntry(3, 3), 1.0e-12);
        Assertions.assertSame(lt, decomposition.getLT());

    }

    @Test
    public void testDimensionsMismatchV() {
        try {
            final RealMatrix matrix = MatrixUtils.createRealIdentityMatrix(4);
            new RegularizedCholeskyDecomposition(matrix).getSolver().solve(MatrixUtils.createRealVector(5));
            Assertions.fail("an exception should have been thrown");
        }  catch (final MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, e.getSpecifier());
        }
    }

    @Test
    public void testSolveV() {
        final RealMatrix matrix = MatrixUtils.createRealMatrix(new double[][] {
                {  1.0,  1.0,  1.0,  1.0 },
                {  1.0,  5.0,  5.0,  5.0 },
                {  1.0,  5.0, 14.0, 14.0 },
                {  1.0,  5.0, 14.0, 15.0 }
        });
        DecompositionSolver solver = new RegularizedCholeskyDecomposition(matrix).getSolver();

        RealVector resultV = solver.solve(MatrixUtils.createRealVector(new double[] { 1.0, 1.0, -1.0, 2.0 }));
        Assertions.assertEquals(1.0,         resultV.getEntry(0), 1.0e-12);
        Assertions.assertEquals(2.0   / 9.0, resultV.getEntry(1), 1.0e-12);
        Assertions.assertEquals(-29.0 / 9.0, resultV.getEntry(2), 1.0e-12);
        Assertions.assertEquals(3.0,         resultV.getEntry(3), 1.0e-12);

    }

    @Test
    public void testDimensionsMismatchM() {
        try {
            final RealMatrix matrix = MatrixUtils.createRealIdentityMatrix(4);
            new RegularizedCholeskyDecomposition(matrix).getSolver().solve(MatrixUtils.createRealIdentityMatrix(5));
            Assertions.fail("an exception should have been thrown");
        }  catch (final MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, e.getSpecifier());
        }
    }

    @Test
    public void testSolveM() {
        final RealMatrix matrix = MatrixUtils.createRealMatrix(new double[][] {
                {  1.0,  1.0,  1.0,  1.0 },
                {  1.0,  5.0,  5.0,  5.0 },
                {  1.0,  5.0, 14.0, 14.0 },
                {  1.0,  5.0, 14.0, 15.0 }
        });
        DecompositionSolver solver = new RegularizedCholeskyDecomposition(matrix).getSolver();
        Assertions.assertEquals(4, solver.getRowDimension());
        Assertions.assertEquals(4, solver.getColumnDimension());

        RealMatrix resultM = solver.getInverse();
        System.out.println(resultM);
        Assertions.assertEquals( 1.25,       resultM.getEntry(0, 0), 1.0e-12);
        Assertions.assertEquals(-0.25,       resultM.getEntry(0, 1), 1.0e-12);
        Assertions.assertEquals( 0.00,       resultM.getEntry(0, 2), 1.0e-12);
        Assertions.assertEquals( 0.00,       resultM.getEntry(0, 3), 1.0e-12);
        Assertions.assertEquals(-0.25,       resultM.getEntry(1, 0), 1.0e-12);
        Assertions.assertEquals(13.0 / 36.0, resultM.getEntry(1, 1), 1.0e-12);
        Assertions.assertEquals(-1.0 /  9.0, resultM.getEntry(1, 2), 1.0e-12);
        Assertions.assertEquals( 0.00,       resultM.getEntry(1, 3), 1.0e-12);
        Assertions.assertEquals( 0.0,        resultM.getEntry(2, 0), 1.0e-12);
        Assertions.assertEquals(-1.0 /  9.0, resultM.getEntry(2, 1), 1.0e-12);
        Assertions.assertEquals(10.0 /  9.0, resultM.getEntry(2, 2), 1.0e-12);
        Assertions.assertEquals(-1.00,       resultM.getEntry(2, 3), 1.0e-12);
        Assertions.assertEquals( 0.00,       resultM.getEntry(3, 0), 1.0e-12);
        Assertions.assertEquals( 0.00,       resultM.getEntry(3, 1), 1.0e-12);
        Assertions.assertEquals(-1.00,       resultM.getEntry(3, 2), 1.0e-12);
        Assertions.assertEquals( 1.00,       resultM.getEntry(3, 3), 1.0e-12);

    }

}
