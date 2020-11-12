/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.linear;

import org.junit.Assert;

import org.hipparchus.complex.Complex;
import org.junit.Test;

public class ComplexEigenDecompositionTest {

    @Test
    public void testGetEigenValues() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        Complex ev1 = eigenDecomp.getEigenvalues()[0];
        Complex ev2 = eigenDecomp.getEigenvalues()[1];
        Assert.assertEquals(new Complex(1, 2), ev1);
        Assert.assertEquals(new Complex(1, -2), ev2);
    }

    @Test
    public void testHasComplexEigenValues() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        Assert.assertTrue(eigenDecomp.hasComplexEigenvalues());
    }

    @Test
    public void testGetDeterminant() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        Assert.assertEquals(5, eigenDecomp.getDeterminant(), 1E-10d);
    }

    @Test
    public void testGetEigenVectors() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        FieldVector<Complex> ev1 = eigenDecomp.getEigenvector(0);
        FieldVector<Complex> ev2 = eigenDecomp.getEigenvector(1);
        Assert.assertEquals(new Complex(.5, .5), ev1.getEntry(0));
        Assert.assertEquals(new Complex(1), ev1.getEntry(1));
        Assert.assertEquals(new Complex(.5, -.5), ev2.getEntry(0));
        Assert.assertEquals(new Complex(1), ev2.getEntry(1));
    }

    @Test
    public void testGetV() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        FieldMatrix<Complex> V = eigenDecomp.getV();
        Assert.assertEquals(new Complex(.5, .5), V.getEntry(0, 0));
        Assert.assertEquals(new Complex(.5, -.5), V.getEntry(0, 1));
        Assert.assertEquals(new Complex(1), V.getEntry(1, 0));
        Assert.assertEquals(new Complex(1), V.getEntry(1, 1));
    }

    @Test
    public void testGetVT() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        FieldMatrix<Complex> V = eigenDecomp.getVT();
        Assert.assertEquals(new Complex(.5, .5), V.getEntry(0, 0));
        Assert.assertEquals(new Complex(.5, -.5), V.getEntry(1, 0));
        Assert.assertEquals(new Complex(1), V.getEntry(0, 1));
        Assert.assertEquals(new Complex(1), V.getEntry(1, 1));
    }

    @Test
    public void testGetD() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        FieldMatrix<Complex> D = eigenDecomp.getD();
        Assert.assertEquals(new Complex(1, 2), D.getEntry(0, 0));
        Assert.assertEquals(new Complex(0), D.getEntry(0, 1));
        Assert.assertEquals(new Complex(0), D.getEntry(0, 1));
        Assert.assertEquals(new Complex(1, -2), D.getEntry(1, 1));
    }

    @Test
    public void testIssue99() {

        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 1.0, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 }
        });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);

        Assert.assertEquals(3, eigenDecomp.getEigenvalues().length);
        for (Complex z : eigenDecomp.getEigenvalues()) {
            Assert.assertEquals(Complex.ONE, z);
        }

        FieldVector<Complex> v0 = eigenDecomp.getEigenvector(0);
        Assert.assertEquals(Complex.ONE,  v0.getEntry(0));
        Assert.assertEquals(Complex.ZERO, v0.getEntry(1));
        Assert.assertEquals(Complex.ZERO, v0.getEntry(2));

        FieldVector<Complex> v1 = eigenDecomp.getEigenvector(1);
        Assert.assertEquals(Complex.ZERO, v1.getEntry(0));
        Assert.assertEquals(Complex.ONE,  v1.getEntry(1));
        Assert.assertEquals(Complex.ZERO, v1.getEntry(2));

        FieldVector<Complex> v2 = eigenDecomp.getEigenvector(2);
        Assert.assertEquals(Complex.ZERO, v2.getEntry(0));
        Assert.assertEquals(Complex.ZERO, v2.getEntry(1));
        Assert.assertEquals(Complex.ONE,  v2.getEntry(2));

    }

    @Test
    public void testDefinition() {
        final RealMatrix aR = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(aR);

        FieldMatrix<Complex> aC = MatrixUtils.createFieldMatrix(new Complex[][] {
            { new Complex(3, 0), new Complex(-2, 0) },
            { new Complex(4, 0), new Complex(-1, 0) } });

        // testing AV = lamba V - [0]
        Assert.assertEquals(aC.operate(eigenDecomp.getEigenvector(0)),
                            eigenDecomp.getEigenvector(0).mapMultiply(eigenDecomp.getEigenvalues()[0]));

        // testing AV = lamba V - [1]
        Assert.assertEquals(aC.operate(eigenDecomp.getEigenvector(1)),
                            eigenDecomp.getEigenvector(1).mapMultiply(eigenDecomp.getEigenvalues()[1]));

        // checking definition of the decomposition A*V = V*D
        Assert.assertEquals(aC.multiply(eigenDecomp.getV()),
                            eigenDecomp.getV().multiply(eigenDecomp.getD()));
    }

}
