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

import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class ComplexEigenDecompositionTest {

    @Test
    void testNonSquare() {
        try {
            new ComplexEigenDecomposition(MatrixUtils.createRealMatrix(2, 3), 1.0e-5, 1.0e-12, 1.0e-6);
            fail("an axception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NON_SQUARE_MATRIX, miae.getSpecifier());
        }
    }

    @Test
    void testRealEigenValues() {
        final RealMatrix m = MatrixUtils.createRealMatrix(new double[][] { { 2, 0 }, { 0, 3 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(m);
        assertFalse(eigenDecomp.hasComplexEigenvalues());
    }

    @Test
    void testGetEigenValues() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        Complex ev1 = eigenDecomp.getEigenvalues()[0];
        Complex ev2 = eigenDecomp.getEigenvalues()[1];
        assertEquals(new Complex(1, +2), ev1);
        assertEquals(new Complex(1, -2), ev2);
    }

    @Test
    void testHasComplexEigenValues() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        assertTrue(eigenDecomp.hasComplexEigenvalues());
    }

    @Test
    void testGetDeterminant() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        assertEquals(5, eigenDecomp.getDeterminant(), 1.0e-12);
    }

    @Test
    void testGetEigenVectors() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        checkScaledVector(eigenDecomp.getEigenvector(0), buildVector(new Complex(1), new Complex(1, -1)), 1.0e-15);
        checkScaledVector(eigenDecomp.getEigenvector(1), buildVector(new Complex(1), new Complex(1, +1)), 1.0e-15);
    }

    @Test
    void testEigenValuesAndVectors() {
        final RealMatrix aR = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        final ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(aR);
        final FieldMatrix<Complex> aC = toComplex(aR);
        for (int i = 0; i < aR.getRowDimension(); ++i) {
            final Complex              lambda = eigenDecomp.getEigenvalues()[i];
            final FieldVector<Complex> u      = eigenDecomp.getEigenvector(i);
            final FieldVector<Complex> v      = aC.operate(u);
            checkScaledVector(v, u.mapMultiplyToSelf(lambda), 1.0e-12);
        }
    }

    @Test
    void testGetV() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        FieldMatrix<Complex> V = eigenDecomp.getV();        
        assertEquals(0.0, new Complex(.5, .5).subtract(V.getEntry(0, 0)).norm(), 1.0e-15);
        assertEquals(0.0, new Complex(.5, -.5).subtract(V.getEntry(0, 1)).norm(), 1.0e-15);
        assertEquals(0.0, new Complex(1).subtract(V.getEntry(1, 0)).norm(), 1.0e-15);
        assertEquals(0.0, new Complex(1).subtract(V.getEntry(1, 1)).norm(), 1.0e-15);
    }

    @Test
    void testGetVT() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        FieldMatrix<Complex> V = eigenDecomp.getVT();
        assertEquals(0.0, new Complex(.5, .5).subtract(V.getEntry(0, 0)).norm(), 1.0e-15);
        assertEquals(0.0, new Complex(.5, -.5).subtract(V.getEntry(1, 0)).norm(), 1.0e-15);
        assertEquals(0.0, new Complex(1).subtract(V.getEntry(0, 1)).norm(), 1.0e-15);
        assertEquals(0.0, new Complex(1).subtract(V.getEntry(1, 1)).norm(), 1.0e-15);
    }

    @Test
    void testGetD() {
        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);
        FieldMatrix<Complex> D = eigenDecomp.getD();
        assertEquals(0.0, new Complex(1, +2).subtract(D.getEntry(0, 0)).norm(), 1.0e-15);
        assertEquals(0.0, new Complex(0).subtract(D.getEntry(0, 1)).norm(), 1.0e-15);
        assertEquals(0.0, new Complex(0).subtract(D.getEntry(0, 1)).norm(), 1.0e-15);
        assertEquals(0.0, new Complex(1, -2).subtract(D.getEntry(1, 1)).norm(), 1.0e-15);
    }

    @Test
    void testEqualEigenvalues() {

        final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 1.0, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 }
        });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(A);

        assertEquals(3, eigenDecomp.getEigenvalues().length);
        for (Complex z : eigenDecomp.getEigenvalues()) {
            assertEquals(Complex.ONE, z);
        }

        checkScaledVector(buildVector(Complex.ONE,  Complex.ZERO, Complex.ZERO), eigenDecomp.getEigenvector(0), 1.0e-12);
        checkScaledVector(buildVector(Complex.ZERO, Complex.ONE,  Complex.ZERO), eigenDecomp.getEigenvector(1), 1.0e-12);
        checkScaledVector(buildVector(Complex.ZERO, Complex.ZERO, Complex.ONE),  eigenDecomp.getEigenvector(2), 1.0e-12);

    }

    @Test
    void testDefinition() {

        final RealMatrix aR = MatrixUtils.createRealMatrix(new double[][] { { 3, -2 }, { 4, -1 } });
        ComplexEigenDecomposition eigenDecomp = new ComplexEigenDecomposition(aR);
        FieldMatrix<Complex> aC = toComplex(aR);

        // testing AV = lamba V - [0]
        checkScaledVector(aC.operate(eigenDecomp.getEigenvector(0)),
                    eigenDecomp.getEigenvector(0).mapMultiply(eigenDecomp.getEigenvalues()[0]),
                    1.0e-12);

        // testing AV = lamba V - [1]
        checkScaledVector(aC.operate(eigenDecomp.getEigenvector(1)),
                    eigenDecomp.getEigenvector(1).mapMultiply(eigenDecomp.getEigenvalues()[1]),
                    1.0e-12);

        // checking definition of the decomposition A*V = V*D
        checkMatrix(aC.multiply(eigenDecomp.getV()),
                    eigenDecomp.getV().multiply(eigenDecomp.getD()),
                    1.0e-12);
    }

    @Test
    void testIssue249() {

        // the characteristic polynomial of this matrix is (1-λ)³,
        // so the matrix has a single eigen value λ=1 and algebraic multiplicity 3.
        // for this eigen value, there are only two eigen vectors: v₁ = (0, 1, 0) and v₂ = (0, 0, 1)
        // as the geometric multiplicity is only 2
        final RealMatrix matrix = new Array2DRowRealMatrix(new double[][] {
            {  1.0, 0.0, 0.0 },
            { -2.0, 1.0, 0.0 },
            {  0.0, 0.0, 1.0 }
        });

        final ComplexEigenDecomposition ced = new OrderedComplexEigenDecomposition(matrix,
                  ComplexEigenDecomposition.DEFAULT_EIGENVECTORS_EQUALITY,
                  ComplexEigenDecomposition.DEFAULT_EPSILON,
                  ComplexEigenDecomposition.DEFAULT_EPSILON_AV_VD_CHECK,
                  (c1, c2) -> Double.compare(c2.norm(), c1.norm()));

        assertEquals(3, ced.getEigenvalues().length);
        assertEquals(1.0, ced.getEigenvector(0).dotProduct(ced.getEigenvector(0)).norm(), 1.0e-15);
        assertEquals(1.0, ced.getEigenvector(1).dotProduct(ced.getEigenvector(1)).norm(), 1.0e-15);
        assertEquals(0.0, ced.getEigenvector(2).dotProduct(ced.getEigenvector(2)).norm(), 1.0e-15);

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

    private void checkScaledVector(final FieldVector<Complex> v, final FieldVector<Complex> reference, final double tol) {

        assertEquals(reference.getDimension(), v.getDimension());

        // find the global scaling factor, using the maximum reference component
        Complex scale = Complex.NaN;
        double  norm  = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < reference.getDimension(); i++) {
            final Complex ri = reference.getEntry(i);
            final double  ni = FastMath.hypot(ri.getReal(), ri.getImaginary());
            if (ni > norm) {
                scale = ri.divide(v.getEntry(i));
                norm  = ni;
            }
        }

        // check vector, applying the scaling factor
        for (int i = 0; i < reference.getDimension(); ++i) {
            final Complex ri = reference.getEntry(i);
            final Complex vi = v.getEntry(i);
            final Complex si = vi.multiply(scale);
            assertEquals(ri.getReal(), si.getReal(), tol, "" + (ri.getReal() - si.getReal()));
            assertEquals(ri.getImaginary(), si.getImaginary(), tol, "" + (ri.getImaginary() - si.getImaginary()));
        }

    }

    private void checkMatrix(final FieldMatrix<Complex> m, final FieldMatrix<Complex> reference, final double tol) {
        assertEquals(reference.getRowDimension(), m.getRowDimension());
        assertEquals(reference.getColumnDimension(), m.getColumnDimension());
        for (int i = 0; i < reference.getRowDimension(); ++i) {
            for (int j = 0; j < reference.getColumnDimension(); ++j) {
                assertEquals(reference.getEntry(i, j).getReal(), m.getEntry(i, j).getReal(), tol, "" + (reference.getEntry(i, j).getReal() - m.getEntry(i, j).getReal()));
                assertEquals(reference.getEntry(i, j).getImaginary(), m.getEntry(i, j).getImaginary(), tol, "" + (reference.getEntry(i, j).getImaginary() - m.getEntry(i, j).getImaginary()));
            }
        }
    }

}
