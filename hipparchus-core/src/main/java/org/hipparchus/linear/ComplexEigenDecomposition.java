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

import java.lang.reflect.Array;

import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Given a matrix A, it computes a complex eigen decomposition AV = VD.
 *
 * <p>
 * Complex Eigen Decomposition differs from the {@link EigenDecompositionSymmetric} since it
 * computes the eigen vectors as complex eigen vectors (if applicable).
 * </p>
 *
 * <p>
 * Beware that in the complex case, you do not always have \(V \times V^{T} = I\) or even a
 * diagonal matrix, even if the eigenvectors that form the columns of the V
 * matrix are independent. On example is the square matrix
 * \[
 * A = \left(\begin{matrix}
 * 3 &amp; -2\\
 * 4 &amp; -1
 * \end{matrix}\right)
 * \]
 * which has two conjugate eigenvalues \(\lambda_1=1+2i\) and \(\lambda_2=1-2i\)
 * with associated eigenvectors \(v_1^T = (1, 1-i)\) and \(v_2^T = (1, 1+i)\).
 * \[
 * V\timesV^T = \left(\begin{matrix}
 * 2 &amp; 2\\
 * 2 &amp; 0
 * \end{matrix}\right)
 * \]
 * which is not the identity matrix. Therefore, despite \(A \times V = V \times D\),
 * \(A \ne V \times D \time V^T\), which would hold for real eigendecomposition.
 * </p>
 * <p>
 * Also note that for consistency with Wolfram langage
 * <a href="https://reference.wolfram.com/language/ref/Eigenvectors.html">eigenvectors</a>,
 * we add zero vectors when the geometric multiplicity of the eigenvalue is smaller than
 * its algebraic multiplicity (hence the regular eigenvector matrix should be non-square).
 * With these additional null vectors, the eigenvectors matrix becomes square. This happens
 * for example with the square matrix
 * \[
 * A = \left(\begin{matrix}
 *  1 &amp; 0 &amp; 0\\
 * -2 &amp; 1 &amp; 0\\
 *  0 &amp; 0 &amp; 1
 * \end{matrix}\right)
 * \]
 * Its characteristic polynomial is \((1-\lambda)^3\), hence is has one eigen value \(\lambda=1\)
 * with algebraic multiplicity 3. However, this eigenvalue leads to only two eigenvectors
 * \(v_1=(0, 1, 0)\) and \(v_2=(0, 0, 1)\), hence its geometric multiplicity is only 2, not 3.
 * So we add a third zero vector \(v_3=(0, 0, 0)\), in the same way Wolfram language does.
 * </p>
 *
 * Compute complex eigen values from the Schur transform. Compute complex eigen
 * vectors based on eigen values and the inverse iteration method.
 *
 * see: https://en.wikipedia.org/wiki/Inverse_iteration
 * https://en.wikiversity.org/wiki/Shifted_inverse_iteration
 * http://www.robots.ox.ac.uk/~sjrob/Teaching/EngComp/ecl4.pdf
 * http://www.math.ohiou.edu/courses/math3600/lecture16.pdf
 *
 */
public class ComplexEigenDecomposition {

    /** Default threshold below which eigenvectors are considered equal. */
    public static final double DEFAULT_EIGENVECTORS_EQUALITY = 1.0e-5;
    /** Default value to use for internal epsilon. */
    public static final double DEFAULT_EPSILON = 1e-12;
    /** Internally used epsilon criteria for final AV=VD check. */
    public static final double DEFAULT_EPSILON_AV_VD_CHECK = 1e-6;
    /** Maximum number of inverse iterations. */
    private static final int MAX_ITER = 10;
    /** complex eigenvalues. */
    private Complex[] eigenvalues;
    /** Eigenvectors. */
    private FieldVector<Complex>[] eigenvectors;
    /** Cached value of V. */
    private FieldMatrix<Complex> V;
    /** Cached value of D. */
    private FieldMatrix<Complex> D;
    /** Internally used threshold below which eigenvectors are considered equal. */
    private final double eigenVectorsEquality;
    /** Internally used epsilon criteria. */
    private final double epsilon;
    /** Internally used epsilon criteria for final AV=VD check. */
    private final double epsilonAVVDCheck;

    /**
     * Constructor for decomposition.
     * <p>
     * This constructor uses the default values {@link #DEFAULT_EIGENVECTORS_EQUALITY},
     * {@link #DEFAULT_EPSILON} and {@link #DEFAULT_EPSILON_AV_VD_CHECK}
     * </p>
     * @param matrix
     *            real matrix.
     */
    public ComplexEigenDecomposition(final RealMatrix matrix) {
        this(matrix, DEFAULT_EIGENVECTORS_EQUALITY,
             DEFAULT_EPSILON, DEFAULT_EPSILON_AV_VD_CHECK);
    }

    /**
     * Constructor for decomposition.
     * <p>
     * The {@code eigenVectorsEquality} threshold is used to ensure the L∞-normalized
     * eigenvectors found using inverse iteration are different from each other.
     * if \(min(|e_i-e_j|,|e_i+e_j|)\) is smaller than this threshold, the algorithm
     * considers it has found again an already known vector, so it drops it and attempts
     * a new inverse iteration with a different start vector. This value should be
     * much larger than {@code epsilon} which is used for convergence
     * </p>
     * @param matrix real matrix.
     * @param eigenVectorsEquality threshold below which eigenvectors are considered equal
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @param epsilonAVVDCheck Epsilon criteria for final AV=VD check
     * @since 1.8
     */
    public ComplexEigenDecomposition(final RealMatrix matrix, final double eigenVectorsEquality,
                                     final double epsilon, final double epsilonAVVDCheck) {

        if (!matrix.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   matrix.getRowDimension(), matrix.getColumnDimension());
        }
        this.eigenVectorsEquality = eigenVectorsEquality;
        this.epsilon              = epsilon;
        this.epsilonAVVDCheck     = epsilonAVVDCheck;

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
     * @return eigen values.
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

    /** Reset eigenvalues and eigen vectors from matrices.
     * <p>
     * This method is intended to be called by sub-classes (mainly {@link OrderedComplexEigenDecomposition})
     * that reorder the matrices elements. It rebuild the eigenvalues and eigen vectors arrays
     * from the D and V matrices.
     * </p>
     * @since 2.1
     */
    protected void matricesToEigenArrays() {
        for (int i = 0; i < eigenvalues.length; ++i) {
            eigenvalues[i] = D.getEntry(i, i);
        }
        for (int i = 0; i < eigenvectors.length; ++i) {
            for (int j = 0; j < eigenvectors[i].getDimension(); ++j) {
                eigenvectors[i].setEntry(j, V.getEntry(j, i));
            }
        }
    }

    /**
     * Confirm if there are complex eigen values.
     *
     * @return true if there are complex eigen values.
     */
    public boolean hasComplexEigenvalues() {
        for (int i = 0; i < eigenvalues.length; i++) {
            if (!Precision.equals(eigenvalues[i].getImaginary(), 0.0, epsilon)) {
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
            if (i == (eigenvalues.length - 1) || Precision.equals(matT[i + 1][i], 0.0, epsilon)) {
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

        // eigen vectors
        eigenvectors = (FieldVector<Complex>[]) Array.newInstance(FieldVector.class, n);

        // computing eigen vector based on eigen values and inverse iteration
        for (int i = 0; i < eigenvalues.length; i++) {

            // shifted non-singular matrix matrix A-(λ+ε)I that is close to the singular matrix A-λI
            Complex mu = eigenvalues[i].add(epsilon);
            final FieldMatrix<Complex> shifted = matrix.copy();
            for (int k = 0; k < matrix.getColumnDimension(); ++k) {
                shifted.setEntry(k, k, shifted.getEntry(k, k).subtract(mu));
            }

            // solver for linear system (A - (λ+ε)I) Bₖ₊₁ = Bₖ
            FieldDecompositionSolver<Complex> solver = new FieldQRDecomposition<>(shifted).getSolver();

            // loop over possible start vectors
            for (int p = 0; eigenvectors[i] == null && p < matrix.getColumnDimension(); ++p) {

                // find a vector to start iterations
                FieldVector<Complex> b = findStart(p);

                if (getNorm(b).norm() > Precision.SAFE_MIN) {
                    // start vector is a good candidate for inverse iteration

                    // perform inverse iteration
                    double delta = Double.POSITIVE_INFINITY;
                    for (int k = 0; delta > epsilon && k < MAX_ITER; k++) {

                        // solve (A - (λ+ε)) Bₖ₊₁ = Bₖ
                        final FieldVector<Complex> bNext = solver.solve(b);

                        // normalize according to L∞ norm
                        normalize(bNext);

                        // compute convergence criterion, comparing Bₖ and both ±Bₖ₊₁
                        // as iterations sometimes flip between two opposite vectors
                        delta = separation(b, bNext);

                        // prepare next iteration
                        b = bNext;

                    }

                    // check we have not found again an already known vector
                    for (int j = 0; b != null && j < i; ++j) {
                        if (separation(eigenvectors[j], b) <= eigenVectorsEquality) {
                            // the selected start vector leads us to found a known vector again,
                            // we must try another start
                            b = null;
                        }
                    }
                    eigenvectors[i] = b;

                }
            }

            if (eigenvectors[i] == null) {
                // for consistency with Wolfram langage
                // https://reference.wolfram.com/language/ref/Eigenvectors.html
                // we add zero vectors when the geometric multiplicity of the eigenvalue
                // is smaller than its algebraic multiplicity (hence the regular eigenvector
                // matrix should be non-square). With these additional null vectors, the
                // eigenvectors matrix becomes square
                eigenvectors[i] = MatrixUtils.createFieldVector(ComplexField.getInstance(), n);
            }

        }
    }

    /** Find a start vector orthogonal to all already found normalized eigenvectors.
     * @param index index of the vector
     * @return start vector
     */
    private FieldVector<Complex> findStart(final int index) {

        // create vector
        final FieldVector<Complex> start =
                        MatrixUtils.createFieldVector(ComplexField.getInstance(),
                                                      eigenvalues.length);

        // initialize with a canonical vector
        start.setEntry(index, Complex.ONE);

        return start;

    }

    /**
     * Compute the L∞ norm of the a given vector.
     *
     * @param vector
     *            vector.
     * @return L∞ norm.
     */
    private Complex getNorm(FieldVector<Complex> vector) {
        double  normR = 0;
        Complex normC = Complex.ZERO;
        for (int i = 0; i < vector.getDimension(); i++) {
            final Complex ci = vector.getEntry(i);
            final double  ni = FastMath.hypot(ci.getReal(), ci.getImaginary());
            if (ni > normR) {
                normR = ni;
                normC = ci;
            }
        }
        return normC;
    }

    /** Normalize a vector with respect to L∞ norm.
     * @param v vector to normalized
     */
    private void normalize(final FieldVector<Complex> v) {
        final Complex invNorm = getNorm(v).reciprocal();
        for (int j = 0; j < v.getDimension(); ++j) {
            v.setEntry(j, v.getEntry(j).multiply(invNorm));
        }
    }

    /** Compute the separation between two normalized vectors (which may be in opposite directions).
     * @param v1 first normalized vector
     * @param v2 second normalized vector
     * @return min (|v1 - v2|, |v1+v2|)
     */
    private double separation(final FieldVector<Complex> v1, final FieldVector<Complex> v2) {
        double deltaPlus  = 0;
        double deltaMinus = 0;
        for (int j = 0; j < v1.getDimension(); ++j) {
            final Complex bCurrj = v1.getEntry(j);
            final Complex bNextj = v2.getEntry(j);
            deltaPlus  = FastMath.max(deltaPlus,
                                      FastMath.hypot(bNextj.getReal()      + bCurrj.getReal(),
                                                     bNextj.getImaginary() + bCurrj.getImaginary()));
            deltaMinus = FastMath.max(deltaMinus,
                                      FastMath.hypot(bNextj.getReal()      - bCurrj.getReal(),
                                                     bNextj.getImaginary() - bCurrj.getImaginary()));
        }
        return FastMath.min(deltaPlus, deltaMinus);
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
        if (!equalsWithPrecision(AV, VD, epsilonAVVDCheck)) {
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
                if (c1.add(c2.negate()).norm() > tolerance) {
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
