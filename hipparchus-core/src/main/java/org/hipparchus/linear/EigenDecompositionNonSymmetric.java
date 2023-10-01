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

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Calculates the eigen decomposition of a non-symmetric real matrix.
 * <p>
 * The eigen decomposition of matrix A is a set of two matrices:
 * \(V\) and \(D\) such that \(A V = V D\) where $\(A\),
 * \(V\) and \(D\) are all \(m \times m\) matrices.
 * <p>
 * This class is similar in spirit to the {@code EigenvalueDecomposition}
 * class from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a>
 * library, with the following changes:
 * <ul>
 *   <li>a {@link #getVInv() getVInv} method has been added,</li>
 *   <li>z {@link #getEigenvalue(int) getEigenvalue} method to pick up a
 *       single eigenvalue has been added,</li>
 *   <li>a {@link #getEigenvector(int) getEigenvector} method to pick up a
 *       single eigenvector has been added,</li>
 *   <li>a {@link #getDeterminant() getDeterminant} method has been added.</li>
 * </ul>
 * <p>
 * This class supports non-symmetric matrices, which have complex eigenvalues.
 * Support for symmetric matrices is provided by {@link EigenDecompositionSymmetric}.
 * </p>
 * <p>
 * As \(A\) is not symmetric, then the eigenvalue matrix \(D\) is block diagonal with
 * the real eigenvalues in 1-by-1 blocks and any complex eigenvalues, \(\lambda \pm i \mu\),
 * in 2-by-2 blocks:
 * </p>
 * <p>
 * \[
 *   \begin{bmatrix}
 *    \lambda &amp; \mu\\
 *    -\mu    &amp; \lambda
 *   \end{bmatrix}
 * \]
 * </p>
 * <p>
 * The columns of \(V\) represent the eigenvectors in the sense that \(A V = V D\),
 * i.e. {@code A.multiply(V)} equals {@code V.multiply(D)}.
 * The matrix \(V\) may be badly conditioned, or even singular, so the validity of the
 * equation \(A = V D V^{-1}\) depends upon the condition of \(V\).
 * </p>
 * <p>
 * This implementation is based on the paper by A. Drubrulle, R.S. Martin and
 * J.H. Wilkinson "The Implicit QL Algorithm" in Wilksinson and Reinsch (1971)
 * Handbook for automatic computation, vol. 2, Linear algebra, Springer-Verlag,
 * New-York.
 *
 * @see <a href="http://mathworld.wolfram.com/EigenDecomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix">Wikipedia</a>
 * @since 3.0
 */
public class EigenDecompositionNonSymmetric {
    /** Default epsilon value to use for internal epsilon **/
    public static final double DEFAULT_EPSILON = 1e-12;
    /** Internally used epsilon criteria. */
    private final double epsilon;
    /** eigenvalues. */
    private Complex[] eigenvalues;
    /** Eigenvectors. */
    private List<FieldVector<Complex>> eigenvectors;
    /** Cached value of \(V\). */
    private RealMatrix cachedV;
    /** Cached value of \(D\). */
    private RealMatrix cachedD;
    /** Cached value of \(V^{-1}\). */
    private RealMatrix cachedVInv;

    /** Calculates the eigen decomposition of the given real matrix.
     * @param matrix Matrix to decompose.
     * @throws MathIllegalStateException if the algorithm fails to converge.
     * @throws MathRuntimeException if the decomposition of a general matrix
     * results in a matrix with zero norm
     */
    public EigenDecompositionNonSymmetric(final RealMatrix matrix) {
        this(matrix, DEFAULT_EPSILON);
    }

    /** Calculates the eigen decomposition of the given real matrix.
     * @param matrix Matrix to decompose.
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @throws MathIllegalStateException if the algorithm fails to converge.
     * @throws MathRuntimeException if the decomposition of a general matrix
     * results in a matrix with zero norm
     */
    public EigenDecompositionNonSymmetric(final RealMatrix matrix, double epsilon)
        throws MathRuntimeException {
        this.epsilon = epsilon;
        findEigenVectorsFromSchur(transformToSchur(matrix));
    }

    /**
     * Gets the matrix V of the decomposition.
     * V is a matrix whose columns hold either the real or the
     * imaginary part of eigenvectors.
     *
     * @return the V matrix.
     */
    public RealMatrix getV() {

        if (cachedV == null) {
            final int m = eigenvectors.size();
            cachedV = MatrixUtils.createRealMatrix(m, m);
            for (int k = 0; k < m; ++k) {
                final FieldVector<Complex> ek = eigenvectors.get(k);
                if (eigenvalues[k].getImaginaryPart() >= 0) {
                    // either it is a real eigenvalue, or it is the first of two conjugate eigenvalues,
                    // we pick up the real part of the eigenvector
                    for (int l = 0; l < m; ++l) {
                        cachedV.setEntry(l, k, ek.getEntry(l).getRealPart());
                    }
                } else {
                    // second of two conjugate eigenvalues,
                    // we pick up the imaginary part of the eigenvector
                    for (int l = 0; l < m; ++l) {
                        cachedV.setEntry(l, k, -ek.getEntry(l).getImaginaryPart());
                    }
                }
            }
        }

        // return the cached matrix
        return cachedV;

    }

    /**
     * Gets the block diagonal matrix D of the decomposition.
     * D is a block diagonal matrix.
     * Real eigenvalues are on the diagonal while complex values are on
     * 2x2 blocks { {real +imaginary}, {-imaginary, real} }.
     *
     * @return the D matrix.
     */
    public RealMatrix getD() {

        if (cachedD == null) {
            // cache the matrix for subsequent calls
            cachedD = MatrixUtils.createRealMatrix(eigenvalues.length, eigenvalues.length);
            for (int i = 0; i < eigenvalues.length; ++i) {
                cachedD.setEntry(i, i, eigenvalues[i].getRealPart());
                if (Precision.compareTo(eigenvalues[i].getImaginaryPart(), 0.0, epsilon) > 0) {
                    cachedD.setEntry(i, i + 1, eigenvalues[i].getImaginaryPart());
                } else if (Precision.compareTo(eigenvalues[i].getImaginaryPart(), 0.0, epsilon) < 0) {
                    cachedD.setEntry(i, i - 1, eigenvalues[i].getImaginaryPart());
                }
            }
        }
        return cachedD;
    }

    /**
     * Get's the value for epsilon which is used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     *
     * @return the epsilon value.
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * Gets the inverse of the matrix V of the decomposition.
     *
     * @return the inverse of the V matrix.
     */
    public RealMatrix getVInv() {

        if (cachedVInv == null) {
            cachedVInv = MatrixUtils.inverse(getV(), epsilon);
        }

        // return the cached matrix
        return cachedVInv;
    }

    /**
     * Gets a copy of the eigenvalues of the original matrix.
     *
     * @return a copy of the eigenvalues of the original matrix.
     *
     * @see #getD()
     * @see #getEigenvalue(int)
     */
    public Complex[] getEigenvalues() {
        return eigenvalues.clone();
    }

    /**
     * Returns the i<sup>th</sup> eigenvalue of the original matrix.
     *
     * @param i index of the eigenvalue (counting from 0)
     * @return i<sup>th</sup> eigenvalue of the original matrix.
     *
     * @see #getD()
     * @see #getEigenvalues()
     */
    public Complex getEigenvalue(final int i) {
        return eigenvalues[i];
    }

    /**
     * Gets a copy of the i<sup>th</sup> eigenvector of the original matrix.
     * <p>
     * Note that if the the i<sup>th</sup> is complex this method will throw
     * an exception.
     * </p>
     * @param i Index of the eigenvector (counting from 0).
     * @return a copy of the i<sup>th</sup> eigenvector of the original matrix.
     * @see #getD()
     */
    public FieldVector<Complex> getEigenvector(final int i) {
        return eigenvectors.get(i).copy();
    }

    /**
     * Computes the determinant of the matrix.
     *
     * @return the determinant of the matrix.
     */
    public Complex getDeterminant() {
        Complex determinant = Complex.ONE;
        for (int i = 0; i < eigenvalues.length; ++i) {
            determinant = determinant.multiply(eigenvalues[i]);
        }
        return determinant;
    }

    /**
     * Transforms the matrix to Schur form and calculates the eigenvalues.
     *
     * @param matrix Matrix to transform.
     * @return the {@link SchurTransformer Schur transform} for this matrix
     */
    private SchurTransformer transformToSchur(final RealMatrix matrix) {
        final SchurTransformer schurTransform = new SchurTransformer(matrix, epsilon);
        final double[][] matT = schurTransform.getT().getData();
        final double norm = matrix.getNorm1();

        eigenvalues = new Complex[matT.length];

        int i = 0;
        while (i < eigenvalues.length) {
            if (i == (eigenvalues.length - 1) ||
                Precision.equals(matT[i + 1][i], 0.0, norm * epsilon)) {
                eigenvalues[i] = new Complex(matT[i][i]);
                i++;
            } else {
                final double x   = matT[i + 1][i + 1];
                final double p   = 0.5 * (matT[i][i] - x);
                final double z   = FastMath.sqrt(FastMath.abs(p * p + matT[i + 1][i] * matT[i][i + 1]));
                eigenvalues[i++] = new Complex(x + p, +z);
                eigenvalues[i++] = new Complex(x + p, -z);
            }
        }
        return schurTransform;
    }

    /**
     * Performs a division of two complex numbers.
     *
     * @param xr real part of the first number
     * @param xi imaginary part of the first number
     * @param yr real part of the second number
     * @param yi imaginary part of the second number
     * @return result of the complex division
     */
    private Complex cdiv(final double xr, final double xi,
                         final double yr, final double yi) {
        return new Complex(xr, xi).divide(new Complex(yr, yi));
    }

    /**
     * Find eigenvectors from a matrix transformed to Schur form.
     *
     * @param schur the schur transformation of the matrix
     * @throws MathRuntimeException if the Schur form has a norm of zero
     */
    private void findEigenVectorsFromSchur(final SchurTransformer schur)
        throws MathRuntimeException {
        final double[][] matrixT = schur.getT().getData();
        final double[][] matrixP = schur.getP().getData();

        final int n = matrixT.length;

        // compute matrix norm
        double norm = 0.0;
        for (int i = 0; i < n; i++) {
           for (int j = FastMath.max(i - 1, 0); j < n; j++) {
               norm += FastMath.abs(matrixT[i][j]);
           }
        }

        // we can not handle a matrix with zero norm
        if (norm == 0.0) {
           throw new MathRuntimeException(LocalizedCoreFormats.ZERO_NORM);
        }

        // Backsubstitute to find vectors of upper triangular form

        double r = 0.0;
        double s = 0.0;
        double z = 0.0;

        for (int idx = n - 1; idx >= 0; idx--) {
            double p = eigenvalues[idx].getRealPart();
            double q = eigenvalues[idx].getImaginaryPart();

            if (Precision.equals(q, 0.0)) {
                // Real vector
                int l = idx;
                matrixT[idx][idx] = 1.0;
                for (int i = idx - 1; i >= 0; i--) {
                    double w = matrixT[i][i] - p;
                    r = 0.0;
                    for (int j = l; j <= idx; j++) {
                        r += matrixT[i][j] * matrixT[j][idx];
                    }
                    if (Precision.compareTo(eigenvalues[i].getImaginaryPart(), 0.0, epsilon) < 0) {
                        z = w;
                        s = r;
                    } else {
                        l = i;
                        if (Precision.equals(eigenvalues[i].getImaginaryPart(), 0.0)) {
                            if (w != 0.0) {
                                matrixT[i][idx] = -r / w;
                            } else {
                                matrixT[i][idx] = -r / (Precision.EPSILON * norm);
                            }
                        } else {
                            // Solve real equations
                            double x = matrixT[i][i + 1];
                            double y = matrixT[i + 1][i];
                            q = (eigenvalues[i].getRealPart() - p) * (eigenvalues[i].getRealPart() - p) +
                                eigenvalues[i].getImaginaryPart() * eigenvalues[i].getImaginaryPart();
                            double t = (x * s - z * r) / q;
                            matrixT[i][idx] = t;
                            if (FastMath.abs(x) > FastMath.abs(z)) {
                                matrixT[i + 1][idx] = (-r - w * t) / x;
                            } else {
                                matrixT[i + 1][idx] = (-s - y * t) / z;
                            }
                        }

                        // Overflow control
                        double t = FastMath.abs(matrixT[i][idx]);
                        if ((Precision.EPSILON * t) * t > 1) {
                            for (int j = i; j <= idx; j++) {
                                matrixT[j][idx] /= t;
                            }
                        }
                    }
                }
            } else if (q < 0.0) {
                // Complex vector
                int l = idx - 1;

                // Last vector component imaginary so matrix is triangular
                if (FastMath.abs(matrixT[idx][idx - 1]) > FastMath.abs(matrixT[idx - 1][idx])) {
                    matrixT[idx - 1][idx - 1] = q / matrixT[idx][idx - 1];
                    matrixT[idx - 1][idx]     = -(matrixT[idx][idx] - p) / matrixT[idx][idx - 1];
                } else {
                    final Complex result = cdiv(0.0, -matrixT[idx - 1][idx],
                                                matrixT[idx - 1][idx - 1] - p, q);
                    matrixT[idx - 1][idx - 1] = result.getReal();
                    matrixT[idx - 1][idx]     = result.getImaginary();
                }

                matrixT[idx][idx - 1] = 0.0;
                matrixT[idx][idx]     = 1.0;

                for (int i = idx - 2; i >= 0; i--) {
                    double ra = 0.0;
                    double sa = 0.0;
                    for (int j = l; j <= idx; j++) {
                        ra += matrixT[i][j] * matrixT[j][idx - 1];
                        sa += matrixT[i][j] * matrixT[j][idx];
                    }
                    double w = matrixT[i][i] - p;

                    if (Precision.compareTo(eigenvalues[i].getImaginaryPart(), 0.0, epsilon) < 0) {
                        z = w;
                        r = ra;
                        s = sa;
                    } else {
                        l = i;
                        if (Precision.equals(eigenvalues[i].getImaginaryPart(), 0.0)) {
                            final Complex c = cdiv(-ra, -sa, w, q);
                            matrixT[i][idx - 1] = c.getReal();
                            matrixT[i][idx] = c.getImaginary();
                        } else {
                            // Solve complex equations
                            double x = matrixT[i][i + 1];
                            double y = matrixT[i + 1][i];
                            double vr = (eigenvalues[i].getRealPart() - p) * (eigenvalues[i].getRealPart() - p) +
                                        eigenvalues[i].getImaginaryPart() * eigenvalues[i].getImaginaryPart() - q * q;
                            final double vi = (eigenvalues[i].getRealPart() - p) * 2.0 * q;
                            if (Precision.equals(vr, 0.0) && Precision.equals(vi, 0.0)) {
                                vr = Precision.EPSILON * norm *
                                     (FastMath.abs(w) + FastMath.abs(q) + FastMath.abs(x) +
                                      FastMath.abs(y) + FastMath.abs(z));
                            }
                            final Complex c     = cdiv(x * r - z * ra + q * sa,
                                                       x * s - z * sa - q * ra, vr, vi);
                            matrixT[i][idx - 1] = c.getReal();
                            matrixT[i][idx]     = c.getImaginary();

                            if (FastMath.abs(x) > (FastMath.abs(z) + FastMath.abs(q))) {
                                matrixT[i + 1][idx - 1] = (-ra - w * matrixT[i][idx - 1] +
                                                           q * matrixT[i][idx]) / x;
                                matrixT[i + 1][idx]     = (-sa - w * matrixT[i][idx] -
                                                           q * matrixT[i][idx - 1]) / x;
                            } else {
                                final Complex c2        = cdiv(-r - y * matrixT[i][idx - 1],
                                                               -s - y * matrixT[i][idx], z, q);
                                matrixT[i + 1][idx - 1] = c2.getReal();
                                matrixT[i + 1][idx]     = c2.getImaginary();
                            }
                        }

                        // Overflow control
                        double t = FastMath.max(FastMath.abs(matrixT[i][idx - 1]),
                                                FastMath.abs(matrixT[i][idx]));
                        if ((Precision.EPSILON * t) * t > 1) {
                            for (int j = i; j <= idx; j++) {
                                matrixT[j][idx - 1] /= t;
                                matrixT[j][idx] /= t;
                            }
                        }
                    }
                }
            }
        }

        // Back transformation to get eigenvectors of original matrix
        for (int j = n - 1; j >= 0; j--) {
            for (int i = 0; i <= n - 1; i++) {
                z = 0.0;
                for (int k = 0; k <= FastMath.min(j, n - 1); k++) {
                    z += matrixP[i][k] * matrixT[k][j];
                }
                matrixP[i][j] = z;
            }
        }

        eigenvectors = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            FieldVector<Complex> ei = new ArrayFieldVector<>(ComplexField.getInstance(), n);
            for (int j = 0; j < n; j++) {
                if (Precision.compareTo(eigenvalues[i].getImaginaryPart(), 0.0, epsilon) > 0) {
                    ei.setEntry(j, new Complex(matrixP[j][i], +matrixP[j][i + 1]));
                } else if (Precision.compareTo(eigenvalues[i].getImaginaryPart(), 0.0, epsilon) < 0) {
                    ei.setEntry(j, new Complex(matrixP[j][i - 1], -matrixP[j][i]));
                } else {
                    ei.setEntry(j, new Complex(matrixP[j][i]));
                }
            }
            eigenvectors.add(ei);
        }
    }
}
