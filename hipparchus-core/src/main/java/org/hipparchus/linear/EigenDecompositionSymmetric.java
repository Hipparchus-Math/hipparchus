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
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Calculates the eigen decomposition of a symmetric real matrix.
 * <p>
 * The eigen decomposition of matrix A is a set of two matrices:
 * \(V\) and \(D\) such that \(A V = V D\) where $\(A\),
 * \(V\) and \(D\) are all \(m \times m\) matrices.
 * <p>
 * This class is similar in spirit to the {@code EigenvalueDecomposition}
 * class from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a>
 * library, with the following changes:
 * </p>
 * <ul>
 *   <li>a {@link #getVT() getVt} method has been added,</li>
 *   <li>a {@link #getEigenvalue(int) getEigenvalue} method to pick up a
 *       single eigenvalue has been added,</li>
 *   <li>a {@link #getEigenvector(int) getEigenvector} method to pick up a
 *       single eigenvector has been added,</li>
 *   <li>a {@link #getDeterminant() getDeterminant} method has been added.</li>
 *   <li>a {@link #getSolver() getSolver} method has been added.</li>
 * </ul>
 * <p>
 * As \(A\) is symmetric, then \(A = V D V^T\) where the eigenvalue matrix \(D\)
 * is diagonal and the eigenvector matrix \(V\) is orthogonal, i.e.
 * {@code A = V.multiply(D.multiply(V.transpose()))} and
 * {@code V.multiply(V.transpose())} equals the identity matrix.
 * </p>
 * <p>
 * The columns of \(V\) represent the eigenvectors in the sense that \(A V = V D\),
 * i.e. {@code A.multiply(V)} equals {@code V.multiply(D)}.
 * The matrix \(V\) may be badly conditioned, or even singular, so the validity of the
 * equation \(A = V D V^{-1}\) depends upon the condition of \(V\).
 * </p>
 * This implementation is based on the paper by A. Drubrulle, R.S. Martin and
 * J.H. Wilkinson "The Implicit QL Algorithm" in Wilksinson and Reinsch (1971)
 * Handbook for automatic computation, vol. 2, Linear algebra, Springer-Verlag,
 * New-York.
 *
 * @see <a href="http://mathworld.wolfram.com/EigenDecomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix">Wikipedia</a>
 */
public class EigenDecompositionSymmetric {

    /** Default epsilon value to use for internal epsilon **/
    public static final double DEFAULT_EPSILON = 1e-12;

    /** Maximum number of iterations accepted in the implicit QL transformation */
    private static final byte MAX_ITER = 30;

    /** Internally used epsilon criteria. */
    private final double epsilon;

    /** Eigenvalues. */
    private double[] eigenvalues;

    /** Eigenvectors. */
    private ArrayRealVector[] eigenvectors;

    /** Cached value of V. */
    private RealMatrix cachedV;

    /** Cached value of D. */
    private DiagonalMatrix cachedD;

    /** Cached value of Vt. */
    private RealMatrix cachedVt;

    /**
     * Calculates the eigen decomposition of the given symmetric real matrix.
     * <p>
     * This constructor uses the {@link #DEFAULT_EPSILON default epsilon} and
     * decreasing order for eigenvalues.
     * </p>
     * @param matrix Matrix to decompose.
     * @throws MathIllegalStateException if the algorithm fails to converge.
     * @throws MathRuntimeException if the decomposition of a general matrix
     * results in a matrix with zero norm
     */
    public EigenDecompositionSymmetric(final RealMatrix matrix) {
        this(matrix, DEFAULT_EPSILON, true);
    }

    /**
     * Calculates the eigen decomposition of the given real matrix.
     * <p>
     * Supports decomposition of a general matrix since 3.1.
     *
     * @param matrix Matrix to decompose.
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @param decreasing if true, eigenvalues will be sorted in decreasing order
     * @throws MathIllegalStateException if the algorithm fails to converge.
     * @throws MathRuntimeException if the decomposition of a general matrix
     * results in a matrix with zero norm
     * @since 3.0
     */
    public EigenDecompositionSymmetric(final RealMatrix matrix,
                                       final double epsilon, final boolean decreasing)
        throws MathRuntimeException {

        this.epsilon = epsilon;
        MatrixUtils.checkSymmetric(matrix, epsilon);

        // transform the matrix to tridiagonal
        final TriDiagonalTransformer transformer = new TriDiagonalTransformer(matrix);

        findEigenVectors(transformer.getMainDiagonalRef(),
                         transformer.getSecondaryDiagonalRef(),
                         transformer.getQ().getData(),
                         decreasing);

    }

    /**
     * Calculates the eigen decomposition of the symmetric tridiagonal matrix.
     * <p>
     * The Householder matrix is assumed to be the identity matrix.
     * </p>
     * <p>
     * This constructor uses the {@link #DEFAULT_EPSILON default epsilon} and
     * decreasing order for eigenvalues.
     * </p>
     * @param main Main diagonal of the symmetric tridiagonal form.
     * @param secondary Secondary of the tridiagonal form.
     * @throws MathIllegalStateException if the algorithm fails to converge.
     */
    public EigenDecompositionSymmetric(final double[] main, final double[] secondary) {
        this(main, secondary, DEFAULT_EPSILON, true);
    }

    /**
     * Calculates the eigen decomposition of the symmetric tridiagonal
     * matrix.  The Householder matrix is assumed to be the identity matrix.
     *
     * @param main Main diagonal of the symmetric tridiagonal form.
     * @param secondary Secondary of the tridiagonal form.
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @param decreasing if true, eigenvalues will be sorted in decreasing order
     * @throws MathIllegalStateException if the algorithm fails to converge.
     * @since 3.0
     */
    public EigenDecompositionSymmetric(final double[] main, final double[] secondary,
                                       final double epsilon, final boolean decreasing) {
        this.epsilon = epsilon;
        final int size = main.length;
        final double[][] z = new double[size][size];
        for (int i = 0; i < size; i++) {
            z[i][i] = 1.0;
        }
        findEigenVectors(main.clone(), secondary.clone(), z, decreasing);
    }

    /**
     * Gets the matrix V of the decomposition.
     * V is an orthogonal matrix, i.e. its transpose is also its inverse.
     * The columns of V are the eigenvectors of the original matrix.
     * No assumption is made about the orientation of the system axes formed
     * by the columns of V (e.g. in a 3-dimension space, V can form a left-
     * or right-handed system).
     *
     * @return the V matrix.
     */
    public RealMatrix getV() {

        if (cachedV == null) {
            final int m = eigenvectors.length;
            cachedV = MatrixUtils.createRealMatrix(m, m);
            for (int k = 0; k < m; ++k) {
                cachedV.setColumnVector(k, eigenvectors[k]);
            }
        }
        // return the cached matrix
        return cachedV;
    }

    /**
     * Gets the diagonal matrix D of the decomposition.
     * D is a diagonal matrix.
     * @return the D matrix.
     *
     * @see #getEigenvalues()
      */
    public DiagonalMatrix getD() {

        if (cachedD == null) {
            // cache the matrix for subsequent calls
            cachedD = new DiagonalMatrix(eigenvalues);
        }

        return cachedD;

    }

    /**
     * Get's the value for epsilon which is used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     *
     * @return the epsilon value.
     */
    public double getEpsilon() { return epsilon; }

    /**
     * Gets the transpose of the matrix V of the decomposition.
     * V is an orthogonal matrix, i.e. its transpose is also its inverse.
     * The columns of V are the eigenvectors of the original matrix.
     * No assumption is made about the orientation of the system axes formed
     * by the columns of V (e.g. in a 3-dimension space, V can form a left-
     * or right-handed system).
     *
     * @return the transpose of the V matrix.
     */
    public RealMatrix getVT() {

        if (cachedVt == null) {
            final int m = eigenvectors.length;
            cachedVt = MatrixUtils.createRealMatrix(m, m);
            for (int k = 0; k < m; ++k) {
                cachedVt.setRowVector(k, eigenvectors[k]);
            }
        }

        // return the cached matrix
        return cachedVt;
    }

    /**
     * Gets a copy of the eigenvalues of the original matrix.
     *
     * @return a copy of the eigenvalues of the original matrix.
     *
     * @see #getD()
     * @see #getEigenvalue(int)
     */
    public double[] getEigenvalues() {
        return eigenvalues.clone();
    }

    /**
     * Returns the i<sup>th</sup> eigenvalue of the original matrix.
     *
     * @param i index of the eigenvalue (counting from 0)
     * @return real part of the i<sup>th</sup> eigenvalue of the original
     * matrix.
     *
     * @see #getD()
     * @see #getEigenvalues()
     */
    public double getEigenvalue(final int i) {
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
    public RealVector getEigenvector(final int i) {
        return eigenvectors[i].copy();
    }

    /**
     * Computes the determinant of the matrix.
     *
     * @return the determinant of the matrix.
     */
    public double getDeterminant() {
        double determinant = 1;
        for (int i = 0; i < eigenvalues.length; ++i) {
            determinant *= eigenvalues[i];
        }
        return determinant;
    }

    /**
     * Computes the square-root of the matrix.
     * This implementation assumes that the matrix is positive definite.
     *
     * @return the square-root of the matrix.
     * @throws MathRuntimeException if the matrix is not
     * symmetric or not positive definite.
     */
    public RealMatrix getSquareRoot() {

        final double[] sqrtEigenValues = new double[eigenvalues.length];
        for (int i = 0; i < eigenvalues.length; i++) {
            final double eigen = eigenvalues[i];
            if (eigen <= 0) {
                throw new MathRuntimeException(LocalizedCoreFormats.UNSUPPORTED_OPERATION);
            }
            sqrtEigenValues[i] = FastMath.sqrt(eigen);
        }
        final RealMatrix sqrtEigen = MatrixUtils.createRealDiagonalMatrix(sqrtEigenValues);
        final RealMatrix v = getV();
        final RealMatrix vT = getVT();

        return v.multiply(sqrtEigen).multiply(vT);

    }

    /** Gets a solver for finding the \(A \times X = B\) solution in exact linear sense.
     * @return a solver
     */
    public DecompositionSolver getSolver() {
        return new Solver();
    }

    /** Specialized solver. */
    private class Solver implements DecompositionSolver {

        /**
         * Solves the linear equation \(A \times X = B\)for symmetric matrices A.
         * <p>
         * This method only finds exact linear solutions, i.e. solutions for
         * which ||A &times; X - B|| is exactly 0.
         * </p>
         *
         * @param b Right-hand side of the equation \(A \times X = B\).
         * @return a Vector X that minimizes the 2-norm of \(A \times X - B\).
         *
         * @throws MathIllegalArgumentException if the matrices dimensions do not match.
         * @throws MathIllegalArgumentException if the decomposed matrix is singular.
         */
        @Override
        public RealVector solve(final RealVector b) {
            if (!isNonSingular()) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.SINGULAR_MATRIX);
            }

            final int m = eigenvalues.length;
            if (b.getDimension() != m) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                       b.getDimension(), m);
            }

            final double[] bp = new double[m];
            for (int i = 0; i < m; ++i) {
                final ArrayRealVector v = eigenvectors[i];
                final double[] vData = v.getDataRef();
                final double s = v.dotProduct(b) / eigenvalues[i];
                for (int j = 0; j < m; ++j) {
                    bp[j] += s * vData[j];
                }
            }

            return new ArrayRealVector(bp, false);
        }

        /** {@inheritDoc} */
        @Override
        public RealMatrix solve(RealMatrix b) {

            if (!isNonSingular()) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.SINGULAR_MATRIX);
            }

            final int m = eigenvalues.length;
            if (b.getRowDimension() != m) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                       b.getRowDimension(), m);
            }

            final int nColB = b.getColumnDimension();
            final double[][] bp = new double[m][nColB];
            final double[] tmpCol = new double[m];
            for (int k = 0; k < nColB; ++k) {
                for (int i = 0; i < m; ++i) {
                    tmpCol[i] = b.getEntry(i, k);
                    bp[i][k]  = 0;
                }
                for (int i = 0; i < m; ++i) {
                    final ArrayRealVector v = eigenvectors[i];
                    final double[] vData = v.getDataRef();
                    double s = 0;
                    for (int j = 0; j < m; ++j) {
                        s += v.getEntry(j) * tmpCol[j];
                    }
                    s /= eigenvalues[i];
                    for (int j = 0; j < m; ++j) {
                        bp[j][k] += s * vData[j];
                    }
                }
            }

            return new Array2DRowRealMatrix(bp, false);

        }

        /**
         * Checks whether the decomposed matrix is non-singular.
         *
         * @return true if the decomposed matrix is non-singular.
         */
        @Override
        public boolean isNonSingular() {
            double largestEigenvalueNorm = 0.0;
            // Looping over all values (in case they are not sorted in decreasing
            // order of their norm).
            for (int i = 0; i < eigenvalues.length; ++i) {
                largestEigenvalueNorm = FastMath.max(largestEigenvalueNorm, FastMath.abs(eigenvalues[i]));
            }
            // Corner case: zero matrix, all exactly 0 eigenvalues
            if (largestEigenvalueNorm == 0.0) {
                return false;
            }
            for (int i = 0; i < eigenvalues.length; ++i) {
                // Looking for eigenvalues that are 0, where we consider anything much much smaller
                // than the largest eigenvalue to be effectively 0.
                if (Precision.equals(FastMath.abs(eigenvalues[i]) / largestEigenvalueNorm, 0, epsilon)) {
                    return false;
                }
            }
            return true;
        }

        /**
         * Get the inverse of the decomposed matrix.
         *
         * @return the inverse matrix.
         * @throws MathIllegalArgumentException if the decomposed matrix is singular.
         */
        @Override
        public RealMatrix getInverse() {
            if (!isNonSingular()) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.SINGULAR_MATRIX);
            }

            final int m = eigenvalues.length;
            final double[][] invData = new double[m][m];

            for (int i = 0; i < m; ++i) {
                final double[] invI = invData[i];
                for (int j = 0; j < m; ++j) {
                    double invIJ = 0;
                    for (int k = 0; k < m; ++k) {
                        final double[] vK = eigenvectors[k].getDataRef();
                        invIJ += vK[i] * vK[j] / eigenvalues[k];
                    }
                    invI[j] = invIJ;
                }
            }
            return MatrixUtils.createRealMatrix(invData);
        }

        /** {@inheritDoc} */
        @Override
        public int getRowDimension() {
            return eigenvalues.length;
        }

        /** {@inheritDoc} */
        @Override
        public int getColumnDimension() {
            return eigenvalues.length;
        }

    }

    /**
     * Find eigenvalues and eigenvectors (Dubrulle et al., 1971)
     * @param main main diagonal of the tridiagonal matrix
     * @param secondary secondary diagonal of the tridiagonal matrix
     * @param householderMatrix Householder matrix of the transformation
     * @param decreasing if true, eigenvalues will be sorted in decreasing order
     * to tridiagonal form.
     */
    private void findEigenVectors(final double[] main, final double[] secondary,
                                  final double[][] householderMatrix, final boolean decreasing) {
        final double[][]z = householderMatrix.clone();
        final int n = main.length;
        eigenvalues = new double[n];
        final double[] e = new double[n];
        for (int i = 0; i < n - 1; i++) {
            eigenvalues[i] = main[i];
            e[i] = secondary[i];
        }
        eigenvalues[n - 1] = main[n - 1];
        e[n - 1] = 0;

        // Determine the largest main and secondary value in absolute term.
        double maxAbsoluteValue = 0;
        for (int i = 0; i < n; i++) {
            if (FastMath.abs(eigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = FastMath.abs(eigenvalues[i]);
            }
            if (FastMath.abs(e[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = FastMath.abs(e[i]);
            }
        }
        // Make null any main and secondary value too small to be significant
        if (maxAbsoluteValue != 0) {
            for (int i=0; i < n; i++) {
                if (FastMath.abs(eigenvalues[i]) <= Precision.EPSILON * maxAbsoluteValue) {
                    eigenvalues[i] = 0;
                }
                if (FastMath.abs(e[i]) <= Precision.EPSILON * maxAbsoluteValue) {
                    e[i]=0;
                }
            }
        }

        for (int j = 0; j < n; j++) {
            int its = 0;
            int m;
            do {
                for (m = j; m < n - 1; m++) {
                    double delta = FastMath.abs(eigenvalues[m]) +
                        FastMath.abs(eigenvalues[m + 1]);
                    if (FastMath.abs(e[m]) + delta == delta) {
                        break;
                    }
                }
                if (m != j) {
                    if (its == MAX_ITER) {
                        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED,
                                                            MAX_ITER);
                    }
                    its++;
                    double q = (eigenvalues[j + 1] - eigenvalues[j]) / (2 * e[j]);
                    double t = FastMath.sqrt(1 + q * q);
                    if (q < 0.0) {
                        q = eigenvalues[m] - eigenvalues[j] + e[j] / (q - t);
                    } else {
                        q = eigenvalues[m] - eigenvalues[j] + e[j] / (q + t);
                    }
                    double u = 0.0;
                    double s = 1.0;
                    double c = 1.0;
                    int i;
                    for (i = m - 1; i >= j; i--) {
                        double p = s * e[i];
                        double h = c * e[i];
                        if (FastMath.abs(p) >= FastMath.abs(q)) {
                            c = q / p;
                            t = FastMath.sqrt(c * c + 1.0);
                            e[i + 1] = p * t;
                            s = 1.0 / t;
                            c *= s;
                        } else {
                            s = p / q;
                            t = FastMath.sqrt(s * s + 1.0);
                            e[i + 1] = q * t;
                            c = 1.0 / t;
                            s *= c;
                        }
                        if (e[i + 1] == 0.0) {
                            eigenvalues[i + 1] -= u;
                            e[m] = 0.0;
                            break;
                        }
                        q = eigenvalues[i + 1] - u;
                        t = (eigenvalues[i] - q) * s + 2.0 * c * h;
                        u = s * t;
                        eigenvalues[i + 1] = q + u;
                        q = c * t - h;
                        for (int ia = 0; ia < n; ia++) {
                            p = z[ia][i + 1];
                            z[ia][i + 1] = s * z[ia][i] + c * p;
                            z[ia][i] = c * z[ia][i] - s * p;
                        }
                    }
                    if (t == 0.0 && i >= j) {
                        continue;
                    }
                    eigenvalues[j] -= u;
                    e[j] = q;
                    e[m] = 0.0;
                }
            } while (m != j);
        }

        // Sort the eigen values (and vectors) in desired order
        for (int i = 0; i < n; i++) {
            int k = i;
            double p = eigenvalues[i];
            for (int j = i + 1; j < n; j++) {
                if (eigenvalues[j] > p == decreasing) {
                    k = j;
                    p = eigenvalues[j];
                }
            }
            if (k != i) {
                eigenvalues[k] = eigenvalues[i];
                eigenvalues[i] = p;
                for (int j = 0; j < n; j++) {
                    p = z[j][i];
                    z[j][i] = z[j][k];
                    z[j][k] = p;
                }
            }
        }

        // Determine the largest eigen value in absolute term.
        maxAbsoluteValue = 0;
        for (int i = 0; i < n; i++) {
            if (FastMath.abs(eigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = FastMath.abs(eigenvalues[i]);
            }
        }
        // Make null any eigen value too small to be significant
        if (maxAbsoluteValue != 0.0) {
            for (int i=0; i < n; i++) {
                if (FastMath.abs(eigenvalues[i]) < Precision.EPSILON * maxAbsoluteValue) {
                    eigenvalues[i] = 0;
                }
            }
        }
        eigenvectors = new ArrayRealVector[n];
        for (int i = 0; i < n; i++) {
            eigenvectors[i] = new ArrayRealVector(n);
            for (int j = 0; j < n; j++) {
                eigenvectors[i].setEntry(j, z[j][i]);
            }
        }
    }

}
