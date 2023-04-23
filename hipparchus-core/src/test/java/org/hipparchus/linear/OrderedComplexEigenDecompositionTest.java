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
package org.hipparchus.linear;

import org.junit.Assert;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.complex.Complex;
import org.junit.Test;

public class OrderedComplexEigenDecompositionTest {

    private final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
        { 3, -2 }, { 4, -1 } });

    @Test
    public void testDefinition() {
        OrderedComplexEigenDecomposition eigenDecomp = new OrderedComplexEigenDecomposition(A);

        FieldMatrix<Complex> A = MatrixUtils.createFieldMatrix(new Complex[][] {
            { new Complex(3, 0), new Complex(-2, 0) },
            { new Complex(4, 0), new Complex(-1, 0) } });

        double threshold = 1.0e-15;

        // testing AV = lamba V - [0]
        compareVectors(A.operate(eigenDecomp.getEigenvector(0)),
                       eigenDecomp.getEigenvector(0).mapMultiply(eigenDecomp.getEigenvalues()[0]),
                       threshold);

        // testing AV = lamba V - [1]
        compareVectors(A.operate(eigenDecomp.getEigenvector(1)),
                       eigenDecomp.getEigenvector(1).mapMultiply(eigenDecomp.getEigenvalues()[1]),
                       threshold);

        // checking definition of the decomposition A*V = V*D
        compareMatrices(A.multiply(eigenDecomp.getV()),
                        eigenDecomp.getV().multiply(eigenDecomp.getD()),
                        threshold);

    }

    @Test
    public void testIssue173() {

        OrderedComplexEigenDecomposition eigenDecomp = new OrderedComplexEigenDecomposition(A);

        UnitTestUtils.assertEquals(eigenDecomp.getEigenvalues()[0], eigenDecomp.getD().getEntry(0, 0), 1.0e-10);
        UnitTestUtils.assertEquals(eigenDecomp.getEigenvalues()[1], eigenDecomp.getD().getEntry(1, 1), 1.0e-10);

        UnitTestUtils.assertEquals(eigenDecomp.getEigenvector(0).getEntry(0), eigenDecomp.getV().getEntry(0, 0), 1.0e-10);
        UnitTestUtils.assertEquals(eigenDecomp.getEigenvector(0).getEntry(1), eigenDecomp.getV().getEntry(1, 0), 1.0e-10);
        UnitTestUtils.assertEquals(eigenDecomp.getEigenvector(1).getEntry(0), eigenDecomp.getV().getEntry(0, 1), 1.0e-10);
        UnitTestUtils.assertEquals(eigenDecomp.getEigenvector(1).getEntry(1), eigenDecomp.getV().getEntry(1, 0), 1.0e-10);

    }

    @Test
    public void testEqualEigenValues() {
        double[][] d = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(d);
        ComplexEigenDecomposition ed = new OrderedComplexEigenDecomposition(matrix);
        for (Complex z : ed.getEigenvalues()) {
            Assert.assertEquals(1.0, z.getReal(),      1.0e-15);
            Assert.assertEquals(0.0, z.getImaginary(), 1.0e-15);
        }
    }

    @Test
    public void testIssue111A() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {
            {-0.0,0.0000002463,-0.000000462,0.0,-0.2843196976,-0.0,-0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {-0.0000002463,-0.0,0.0000004577,0.0,-0.0,-0.2843196976,-0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {0.000000462,-0.0000004577,-0.0,0.0,-0.0,-0.0,-0.2843196976,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {-0.0000004577,-0.000000462,-0.0000002463,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {0.0,0.0,0.0,0.0,-0.0000000015,0.0000005751,-0.0000015871,0.0,0.0,0.0,0.0,-0.0000104075,0.0,0.0},
            {0.0,0.0,0.0,0.0,-0.000000421,-0.0000000048,0.0000013528,0.0,0.0,0.0,0.0,0.0,-0.000007717,0.0},
            {0.0,0.0,0.0,0.0,0.0000005375,-0.00000062,0.0000000063,0.0,0.0,0.0,0.0,0.0,0.0,-0.0000035509},
            {-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0000002463,-0.000000462,0.0000004577,0.0,0.0,0.0},
            {0.0,-1.0,0.0,0.0,0.0,0.0,0.0,-0.0000002463,0.0,0.0000004577,0.000000462,0.0,0.0,0.0},
            {0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.000000462,-0.0000004577,0.0,0.0000002463,0.0,0.0,0.0},
            {0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.2843196976,0.0,0.0,0.0,0.0000000015,0.000000421,-0.0000005375},
            {0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.2843196976,0.0,0.0,-0.0000005751,0.0000000048,0.00000062},
            {0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.2843196976,0.0,0.0000015871,-0.0000013528,-0.0000000063}
        });
        
        Assert.assertNotNull(new OrderedComplexEigenDecomposition(m,
                                                                  1.0e-6,
                                                                  ComplexEigenDecomposition.DEFAULT_EPSILON,
                                                                  ComplexEigenDecomposition.DEFAULT_EPSILON_AV_VD_CHECK));
        
    }


    @Test
    public void testIssue111B() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {
            {-0.0,0.0000002463,-0.000000462,0.0,-0.1095029805,-0.0,-0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {-0.0000002463,-0.0,0.0000004577,0.0,-0.0,-0.1095029805,-0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {0.000000462,-0.0000004577,-0.0,0.0,-0.0,-0.0,-0.1095029805,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {-0.0000004577,-0.000000462,-0.0000002463,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {0.0,0.0,0.0,0.0,-0.0000000015,0.0000005751,-0.0000015871,0.0,0.0,0.0,0.0,-0.0000104075,0.0,0.0},
            {0.0,0.0,0.0,0.0,-0.000000421,-0.0000000048,0.0000013528,0.0,0.0,0.0,0.0,0.0,-0.000007717,0.0},
            {0.0,0.0,0.0,0.0,0.0000005375,-0.00000062,0.0000000063,0.0,0.0,0.0,0.0,0.0,0.0,-0.0000035509},
            {-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0000002463,-0.000000462,0.0000004577,0.0,0.0,0.0},
            {0.0,-1.0,0.0,0.0,0.0,0.0,0.0,-0.0000002463,0.0,0.0000004577,0.000000462,0.0,0.0,0.0},
            {0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.000000462,-0.0000004577,0.0,0.0000002463,0.0,0.0,0.0},
            {0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.1095029805,0.0,0.0,0.0,0.0000000015,0.000000421,-0.0000005375},
            {0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.1095029805,0.0,0.0,-0.0000005751,0.0000000048,0.00000062},
            {0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.1095029805,0.0,0.0000015871,-0.0000013528,-0.0000000063}
        });
        
        Assert.assertNotNull(new OrderedComplexEigenDecomposition(m,
                                                                  1.0e-6,
                                                                  ComplexEigenDecomposition.DEFAULT_EPSILON,
                                                                  ComplexEigenDecomposition.DEFAULT_EPSILON_AV_VD_CHECK));
        
    }

    @Test
    public void testIssue111C() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {
            {-0.0,0.0000002463,-0.000000462,0.0,-0.109502989,-0.0,-0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {-0.0000002463,-0.0,0.0000004577,0.0,-0.0,-0.109502989,-0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {0.000000462,-0.0000004577,-0.0,0.0,-0.0,-0.0,-0.109502989,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {-0.0000004577,-0.000000462,-0.0000002463,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {0.0,0.0,0.0,0.0,-0.0000000015,0.0000005751,-0.0000015871,0.0,0.0,0.0,0.0,-0.0000104075,0.0,0.0},
            {0.0,0.0,0.0,0.0,-0.000000421,-0.0000000048,0.0000013528,0.0,0.0,0.0,0.0,0.0,-0.000007717,0.0},
            {0.0,0.0,0.0,0.0,0.0000005375,-0.00000062,0.0000000063,0.0,0.0,0.0,0.0,0.0,0.0,-0.0000035509},
            {-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0000002463,-0.000000462,0.0000004577,0.0,0.0,0.0},
            {0.0,-1.0,0.0,0.0,0.0,0.0,0.0,-0.0000002463,0.0,0.0000004577,0.000000462,0.0,0.0,0.0},
            {0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.000000462,-0.0000004577,0.0,0.0000002463,0.0,0.0,0.0},
            {0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
            {0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.109502989,0.0,0.0,0.0,0.0000000015,0.000000421,-0.0000005375},
            {0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.109502989,0.0,0.0,-0.0000005751,0.0000000048,0.00000062},
            {0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.109502989,0.0,0.0000015871,-0.0000013528,-0.0000000063}
        });
        
        Assert.assertNotNull(new OrderedComplexEigenDecomposition(m,
                                                                  1.0e-6,
                                                                  ComplexEigenDecomposition.DEFAULT_EPSILON,
                                                                  ComplexEigenDecomposition.DEFAULT_EPSILON_AV_VD_CHECK));
        
    }

    /**
     * Verify two complex vectors.
     * @param vector1 first vector
     * @param vector2 second vector
     * @param threshold threshold
     */
    private static void compareVectors(final FieldVector<Complex> vector1,
                                       final FieldVector<Complex> vector2,
                                       final double threshold) {

        // Verify size
        Assert.assertEquals(vector1.getDimension(), vector2.getDimension());

        // Loop on vector entries
        for (int index = 0; index < vector1.getDimension(); index++) {
            final Complex complex1 = vector1.getEntry(index);
            final Complex complex2 = vector2.getEntry(index);
            Assert.assertEquals(complex1.getReal(),      complex2.getReal(),      threshold);
            Assert.assertEquals(complex1.getImaginary(), complex2.getImaginary(), threshold);
        }

    }

    /**
     * Verify two complex matrices.
     * @param vector1 first vector
     * @param vector2 second vector
     * @param threshold threshold
     */
    private static void compareMatrices(final FieldMatrix<Complex> matrix1,
                                        final FieldMatrix<Complex> matrix2,
                                        final double threshold) {

        // Verify size
        Assert.assertEquals(matrix1.getRowDimension(),    matrix2.getRowDimension());
        Assert.assertEquals(matrix1.getColumnDimension(), matrix2.getColumnDimension());

        // Loop on matrices entries
        for (int row = 0; row < matrix1.getRowDimension(); row++) {
            for (int column = 0; column < matrix1.getColumnDimension(); column++) {
                final Complex complex1 = matrix1.getEntry(row, column);
                final Complex complex2 = matrix2.getEntry(row, column);
                Assert.assertEquals(complex1.getReal(),      complex2.getReal(),      threshold);
                Assert.assertEquals(complex1.getImaginary(), complex2.getImaginary(), threshold);
            }

        }

    }

}
