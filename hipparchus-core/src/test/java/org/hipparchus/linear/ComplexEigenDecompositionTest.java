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

import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.Assert;
import org.junit.Test;

public class ComplexEigenDecompositionTest {

    @Test
    public void testNonSquare() {
        try {
            new ComplexEigenDecomposition(MatrixUtils.createRealMatrix(2, 3));
            Assert.fail("an axception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.NON_SQUARE_MATRIX, miae.getSpecifier());
        }
    }

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
        checkVector(eigenDecomp.getEigenvector(0), buildVector(new Complex(.5,  .5), new Complex(1)), 1.0e-15);
        checkVector(eigenDecomp.getEigenvector(1), buildVector(new Complex(.5, -.5), new Complex(1)), 1.0e-15);
    }

    @Test
    public void testEigenValuesAndVectors() {
        final RealMatrix aR = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        final ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(aR);
        final FieldMatrix<Complex> aC = toComplex(aR);
        for (int i = 0; i < aR.getRowDimension(); ++i) {
            final Complex              lambda = eigenDecomp.getEigenvalues()[i];
            final FieldVector<Complex> u      = eigenDecomp.getEigenvector(i);
            final FieldVector<Complex> v      = aC.operate(u);
            checkVector(v, u.mapMultiplyToSelf(lambda), 1.0e-12);
        }
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

        checkVector(buildVector(Complex.ONE,  Complex.ZERO, Complex.ZERO), eigenDecomp.getEigenvector(0), 1.0e-12);
        checkVector(buildVector(Complex.ZERO, Complex.ONE,  Complex.ZERO), eigenDecomp.getEigenvector(1), 1.0e-12);
        checkVector(buildVector(Complex.ZERO, Complex.ZERO, Complex.ONE), eigenDecomp.getEigenvector(2), 1.0e-12);

    }

    @Test
    public void testDefinition() {

        final RealMatrix aR = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(aR);
        FieldMatrix<Complex> aC = toComplex(aR);

        // testing AV = lamba V - [0]
        checkVector(aC.operate(eigenDecomp.getEigenvector(0)),
                    eigenDecomp.getEigenvector(0).mapMultiply(eigenDecomp.getEigenvalues()[0]),
                    1.0e-12);

        // testing AV = lamba V - [1]
        checkVector(aC.operate(eigenDecomp.getEigenvector(1)),
                    eigenDecomp.getEigenvector(1).mapMultiply(eigenDecomp.getEigenvalues()[1]),
                    1.0e-12);

        // checking definition of the decomposition A*V = V*D
        checkMatrix(aC.multiply(eigenDecomp.getV()),
                    eigenDecomp.getV().multiply(eigenDecomp.getD()),
                    1.0e-12);
    }

    private FieldMatrix<Complex> toComplex(final RealMatrix m) {
        FieldMatrix<Complex> c = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),
                                                               m.getRowDimension(),
                                                               m.getColumnDimension());
        for (int i = 0; i < m.getRowDimension(); ++i) {
            for (int j = 0; j < m.getColumnDimension(); ++j) {
                c.setEntry(i, j, new Complex(m.getEntry(i, j)));
            }
        }

        return c;

    }

    private FieldVector<Complex> buildVector(final Complex... vi) {
        return new ArrayFieldVector<>(vi);
    }

    private void checkVector(final FieldVector<Complex> v, final FieldVector<Complex> reference, final double tol) {
        Assert.assertEquals(reference.getDimension(), v.getDimension());
        for (int i = 0; i < reference.getDimension(); ++i) {
            Assert.assertEquals("" + (reference.getEntry(i).getReal() - v.getEntry(i).getReal()),
                                reference.getEntry(i).getReal(), v.getEntry(i).getReal(), tol);
            Assert.assertEquals("" + (reference.getEntry(i).getImaginary() - v.getEntry(i).getImaginary()),
                                reference.getEntry(i).getImaginary(), v.getEntry(i).getImaginary(), tol);
        }
    }

    private void checkMatrix(final FieldMatrix<Complex> m, final FieldMatrix<Complex> reference, final double tol) {
        Assert.assertEquals(reference.getRowDimension(), m.getRowDimension());
        Assert.assertEquals(reference.getColumnDimension(), m.getColumnDimension());
        for (int i = 0; i < reference.getRowDimension(); ++i) {
            for (int j = 0; j < reference.getColumnDimension(); ++j) {
                Assert.assertEquals("" + (reference.getEntry(i, j).getReal() - m.getEntry(i, j).getReal()),
                                    reference.getEntry(i, j).getReal(), m.getEntry(i, j).getReal(), tol);
                Assert.assertEquals("" + (reference.getEntry(i, j).getImaginary() - m.getEntry(i, j).getImaginary()),
                                    reference.getEntry(i, j).getImaginary(), m.getEntry(i, j).getImaginary(), tol);
            }
        }
    }

}
