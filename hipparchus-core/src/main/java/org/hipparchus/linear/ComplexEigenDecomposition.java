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
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Given a matrix A, it computes a complex eigen decomposition A = VDV^{T}. It
 * checks the definition in runtime using AV = VD.
 *
 * Complex Eigen Decomposition, it differs from the EigenDecomposition since it
 * computes the eigen vectors as complex eigen vectors (if applicable).
 *
 * Compute complex eigen values from the schur transform. Compute complex eigen
 * vectors based on eigen values and the inverse iteration method.
 *
 * see: https://en.wikipedia.org/wiki/Inverse_iteration
 * https://en.wikiversity.org/wiki/Shifted_inverse_iteration
 * http://www.robots.ox.ac.uk/~sjrob/Teaching/EngComp/ecl4.pdf
 * http://www.math.ohiou.edu/courses/math3600/lecture16.pdf
 *
 */
public class ComplexEigenDecomposition {

    /** Internally used epsilon criteria. */
    private static final double EPSILON = 1e-12;
    /** Internally used epsilon criteria for equals. */
    private static final double EPSILON_EQUALS = 1e-6;
    /** complex eigenvalues. */
    private Complex[] eigenvalues;
    /** Eigenvectors. */
    private ArrayFieldVector<Complex>[] eigenvectors;
    /** Cached value of V. */
    private FieldMatrix<Complex> V;
    /** Cached value of D. */
    private FieldMatrix<Complex> D;

    /**
     * Constructor for decomposition.
     *
     * @param matrix
     *            real matrix.
     */
    public ComplexEigenDecomposition(final RealMatrix matrix) {

        if (!matrix.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   matrix.getRowDimension(), matrix.getColumnDimension());
        }

        // computing the eigen values
        findEigenValues(matrix);
        // computing the eigen vectors
        findEigenVectors(convertToFieldComplex(matrix));

        // V
        final int m = eigenvectors.length;
        V = MatrixUtils.createFieldMatrix(ComplexField.getInstance(), m, m);
        for (int k = 0; k < m; ++k) {
            V.setColumnVector(k, eigenvectors[k]);
        }

        // D
        D = MatrixUtils.createFieldDiagonalMatrix(eigenvalues);

        checkDefinition(matrix);
    }

    /**
     * Getter of the eigen values.
     *
     * @return igen values.
     */
    public Complex[] getEigenvalues() {
        return eigenvalues.clone();
    }

    /**
     * Getter of the eigen vectors.
     *
     * @param i
     *            which eigen vector.
     * @return eigen vector.
     */
    public FieldVector<Complex> getEigenvector(final int i) {
        return eigenvectors[i].copy();
    }

    /**
     * Confirm if there are complex eigen values.
     *
     * @return true if there are complex eigen values.
     */
    public boolean hasComplexEigenvalues() {
        for (int i = 0; i < eigenvalues.length; i++) {
            if (!Precision.equals(eigenvalues[i].getImaginary(), 0.0, EPSILON)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Computes the determinant.
     *
     * @return the determinant.
     */
    public double getDeterminant() {
        Complex determinant = new Complex(1, 0);
        for (Complex lambda : eigenvalues) {
            determinant = determinant.multiply(lambda);
        }
        return determinant.getReal();
    }

    /**
     * Getter V.
     *
     * @return V.
     */
    public FieldMatrix<Complex> getV() {
        return V;
    }

    /**
     * Getter D.
     *
     * @return D.
     */
    public FieldMatrix<Complex> getD() {
        return D;
    }

    /**
     * Getter VT.
     *
     * @return VT.
     */
    public FieldMatrix<Complex> getVT() {
        return V.transpose();
    }

    /**
     * Compute eigen values using the Schur transform.
     *
     * @param matrix
     *            real matrix to compute eigen values.
     */
    protected void findEigenValues(final RealMatrix matrix) {
        final SchurTransformer schurTransform = new SchurTransformer(matrix);
        final double[][] matT = schurTransform.getT().getData();

        eigenvalues = new Complex[matT.length];

        for (int i = 0; i < eigenvalues.length; i++) {
            if (i == (eigenvalues.length - 1) || Precision.equals(matT[i + 1][i], 0.0, EPSILON)) {
                eigenvalues[i] = new Complex(matT[i][i]);
            } else {
                final double x = matT[i + 1][i + 1];
                final double p = 0.5 * (matT[i][i] - x);
                final double z = FastMath.sqrt(FastMath.abs(p * p + matT[i + 1][i] * matT[i][i + 1]));
                eigenvalues[i] = new Complex(x + p, z);
                eigenvalues[i + 1] = new Complex(x + p, -z);
                i++;
            }
        }

    }

    /**
     * Compute the eigen vectors using the inverse power method.
     *
     * @param matrix
     *            real matrix to compute eigen vectors.
     */
    @SuppressWarnings("unchecked")
    protected void findEigenVectors(final FieldMatrix<Complex> matrix) {
        // number of eigen values/vectors
        int n = eigenvalues.length;

        // identity
        final FieldMatrix<Complex> ident =
                        MatrixUtils.createFieldIdentityMatrix(ComplexField.getInstance(), n);

        // eigen vectors
        eigenvectors = (ArrayFieldVector<Complex>[]) new ArrayFieldVector[n];

        // computing eigen vector based on eigen values and inverse iteration
        for (int i = 0; i < eigenvalues.length; i++) {
            Complex eigenValue = eigenvalues[i];

            // mu multiplied by Identity
            // muI
            FieldMatrix<Complex> eigv = ident.scalarMultiply(eigenValue.add(EPSILON));

            // A-muI
            FieldMatrix<Complex> Aeigv = (FieldMatrix<Complex>) matrix.subtract(eigv);

            // finding inverse of (A - muI)
            // (A - muI) (A - muI)^{-1} = I
            FieldLUDecomposition<Complex> luDecomp = new FieldLUDecomposition<>(Aeigv);
            FieldMatrix<Complex> inv_Aeigv = luDecomp.getSolver().getInverse();

            // starting with a unitary vector
            Complex[] unityVector = new Complex[n];
            for (int k = 0; k < n; k++) {
                unityVector[k] = new Complex(1, 0);
            }
            FieldVector<Complex> eigenVector = MatrixUtils.createFieldVector(unityVector);

            for (int k = 0; k < 2; k++) {
                FieldVector<Complex> eigenVector_new = inv_Aeigv.operate(eigenVector);
                // normalizing
                Complex norm = getNormInf(eigenVector_new);
                eigenVector = eigenVector_new.mapDivide(norm);
            }

            eigenVector = eigenVector.mapAdd(Complex.ZERO);

            eigenvectors[i] = new ArrayFieldVector<>(eigenVector.toArray());
        }
    }

    /**
     * Compute the infinity norm of the a given vector.
     *
     * @param vector
     *            vector.
     * @return infinity norm.
     */
    private Complex getNormInf(FieldVector<Complex> vector) {
        Complex norm = null;
        for (int i = 0; i < vector.getDimension(); i++) {
            if (norm == null) {
                norm = vector.getEntry(i);
            } else if (norm.abs() < vector.getEntry(i).abs()) {
                norm = vector.getEntry(i);
            }
        }
        return norm;
    }

    /**
     * Check definition of the decomposition in runtime.
     *
     * @param matrix
     *            matrix to be decomposed.
     */
    protected void checkDefinition(final RealMatrix matrix) {
        FieldMatrix<Complex> matrixC = convertToFieldComplex(matrix);

        // checking definition of the decomposition
        // testing A*V = V*D
        FieldMatrix<Complex> AV = matrixC.multiply(getV());
        FieldMatrix<Complex> VD = getV().multiply(getD());
        if (!equalsWithPrecision(AV, VD, EPSILON_EQUALS)) {
            throw new MathRuntimeException(LocalizedCoreFormats.FAILED_DECOMPOSITION,
                                           matrix.getRowDimension(), matrix.getColumnDimension());

        }

    }

    /**
     * Helper method that checks with two matrix is equals taking into account a
     * given precision.
     *
     * @param matrix1 first matrix to compare
     * @param matrix2 second matrix to compare
     * @param tolerance tolerance on matrices entries
     * @return true is matrices entries are equal within tolerance,
     * false otherwise
     */
    private boolean equalsWithPrecision(final FieldMatrix<Complex> matrix1,
                                        final FieldMatrix<Complex> matrix2, final double tolerance) {
        boolean toRet = true;
        for (int i = 0; i < matrix1.getRowDimension(); i++) {
            for (int j = 0; j < matrix1.getColumnDimension(); j++) {
                Complex c1 = matrix1.getEntry(i, j);
                Complex c2 = matrix2.getEntry(i, j);
                if (c1.add(c2.negate()).abs() > tolerance) {
                    toRet = false;
                    break;
                }
            }
        }
        return toRet;
    }

    /**
     * It converts a real matrix into a complex field matrix.
     *
     * @param matrix
     *            real matrix.
     * @return complex matrix.
     */
    private FieldMatrix<Complex> convertToFieldComplex(RealMatrix matrix) {
        final FieldMatrix<Complex> toRet =
                        MatrixUtils.createFieldIdentityMatrix(ComplexField.getInstance(),
                                                              matrix.getRowDimension());
        for (int i = 0; i < toRet.getRowDimension(); i++) {
            for (int j = 0; j < toRet.getColumnDimension(); j++) {
                toRet.setEntry(i, j, new Complex(matrix.getEntry(i, j)));
            }
        }
        return toRet;
    }
}
