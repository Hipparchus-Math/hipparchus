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
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Class transforming a general complex matrix to Schur form.
 * <p>A m &times; m complex matrix A can be written as the product of three matrices: A = P
 * &times; T &times; P<sup>H</sup> with P a unitary matrix and T an upper triangular
 * matrix. Both P and T are m &times; m matrices.</p>
 * <p>Transformation to Schur form is often not a goal by itself, but it is an
 * intermediate step in more general decomposition algorithms like
 * {@link ComplexEigenDecomposition eigen decomposition}. This class is therefore
 * intended for expert use. As a consequence of this explicitly limited scope,
 * many methods directly returns references to internal arrays, not copies.</p>
 *
 * @see <a href="http://mathworld.wolfram.com/SchurDecomposition.html">Schur Decomposition - MathWorld</a>
 */
public class ComplexSchurTransformer {
    /** Maximum allowed iterations for convergence of the transformation. */
    private static final int MAX_ITERATIONS = 100;

    /** P matrix. */
    private final Complex[][] matrixP;
    /** T matrix. */
    private final Complex[][] matrixT;
    /** Cached value of P. */
    private FieldMatrix<Complex> cachedP;
    /** Cached value of T. */
    private FieldMatrix<Complex> cachedT;
    /** Cached value of Pt. */
    private FieldMatrix<Complex> cachedPt;

    /** Epsilon criteria. */
    private final double epsilon;

    /**
     * Build the transformation to Schur form of a general complex matrix.
     *
     * @param matrix matrix to transform
     * @throws MathIllegalArgumentException if the matrix is not square
     */
    public ComplexSchurTransformer(final FieldMatrix<Complex> matrix) {
        this(matrix, Precision.EPSILON);
    }

    /**
     * Build the transformation to Schur form of a general complex matrix.
     *
     * @param matrix matrix to transform
     * @param epsilon convergence criteria
     * @throws MathIllegalArgumentException if the matrix is not square
     */
    public ComplexSchurTransformer(final FieldMatrix<Complex> matrix, final double epsilon) {
        if (!matrix.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   matrix.getRowDimension(), matrix.getColumnDimension());
        }

        this.epsilon = epsilon;
        final int n = matrix.getRowDimension();
        matrixP = new Complex[n][n];
        matrixT = new Complex[n][n];

        // Initialize T with the input matrix and P with Identity
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrixT[i][j] = matrix.getEntry(i, j);
                matrixP[i][j] = (i == j) ? Complex.ONE : Complex.ZERO;
            }
        }

        // 1. Reduce to Hessenberg form
        transformToHessenberg();

        // 2. Reduce Hessenberg to Schur form (Triangular)
        transformToSchur();
    }

    /**
     * Returns the unitary matrix P of the transform.
     * <p>P is a unitary matrix, i.e. its inverse is also its conjugate transpose.</p>
     *
     * @return the P matrix
     */
    public FieldMatrix<Complex> getP() {
        if (cachedP == null) {
            cachedP = MatrixUtils.createFieldMatrix(matrixP);
        }
        return cachedP;
    }

    /**
     * Returns the conjugate transpose of the matrix P of the transform.
     * <p>P is a unitary matrix, i.e. its inverse is also its conjugate transpose.</p>
     *
     * @return the conjugate transpose of the P matrix
     */
    public FieldMatrix<Complex> getPT() {
        if (cachedPt == null) {
            cachedPt = getP().transpose().map(Complex::conjugate);
        }
        return cachedPt;
    }

    /**
     * Returns the triangular Schur matrix T of the transform.
     *
     * @return the T matrix
     */
    public FieldMatrix<Complex> getT() {
        if (cachedT == null) {
            cachedT = MatrixUtils.createFieldMatrix(matrixT);
        }
        return cachedT;
    }

    /**
     * Transform the matrix to Hessenberg form using Householder reflections.
     */
    private void transformToHessenberg() {
        final int n = matrixT.length;

        for (int k = 0; k < n - 2; k++) {
            // Find Householder vector to eliminate T[k+2..n][k]
            double normSq = 0;
            for (int i = k + 1; i < n; i++) {
                normSq += matrixT[i][k].norm() * matrixT[i][k].norm();
            }

            if (normSq == 0) {
                continue;
            }

            double normX = FastMath.sqrt(normSq);
            Complex pivot = matrixT[k + 1][k];

            // Compute alpha = - (pivot / |pivot|) * normX
            Complex alpha;
            if (pivot.norm() != 0) {
                alpha = pivot.divide(pivot.norm()).multiply(normX).negate();
            } else {
                alpha = new Complex(-normX);
            }

            // u = x - alpha * e1
            // Copy u to a separate array to avoid corrupting it during the Left Multiply
            int uLen = n - (k + 1);
            Complex[] u = new Complex[uLen];
            u[0] = pivot.subtract(alpha);
            for (int i = 1; i < uLen; i++) {
                u[i] = matrixT[k + 1 + i][k];
            }

            double uNormSq = normSq - pivot.norm() * pivot.norm() + u[0].norm() * u[0].norm();
            if (uNormSq == 0) {
                continue;
            }
            double beta = 2.0 / uNormSq;

            // Apply P = I - beta * u * u^H

            // 1. Left Multiply: T = P * T
            // Affects rows k+1 .. n
            for (int j = k + 1; j < n; j++) {
                Complex dot = Complex.ZERO;
                for (int i = 0; i < uLen; i++) {
                    dot = dot.add(u[i].conjugate().multiply(matrixT[k + 1 + i][j]));
                }
                dot = dot.multiply(beta);

                for (int i = 0; i < uLen; i++) {
                    matrixT[k + 1 + i][j] = matrixT[k + 1 + i][j].subtract(dot.multiply(u[i]));
                }
            }

            // 2. Right Multiply: T = T * P^H
            // Affects columns k+1 .. n
            for (int i = 0; i < n; i++) {
                Complex dot = Complex.ZERO;
                for (int j = 0; j < uLen; j++) {
                    dot = dot.add(matrixT[i][k + 1 + j].multiply(u[j]));
                }
                dot = dot.multiply(beta);

                for (int j = 0; j < uLen; j++) {
                    matrixT[i][k + 1 + j] = matrixT[i][k + 1 + j].subtract(dot.multiply(u[j].conjugate()));
                }
            }

            // 3. Update P: P = P * P^H
            for (int i = 0; i < n; i++) {
                Complex dot = Complex.ZERO;
                for (int j = 0; j < uLen; j++) {
                    dot = dot.add(matrixP[i][k + 1 + j].multiply(u[j]));
                }
                dot = dot.multiply(beta);

                for (int j = 0; j < uLen; j++) {
                    matrixP[i][k + 1 + j] = matrixP[i][k + 1 + j].subtract(dot.multiply(u[j].conjugate()));
                }
            }

            // Restore Hessenberg form manually in column k
            matrixT[k + 1][k] = alpha;
            for (int i = k + 2; i < n; i++) {
                matrixT[i][k] = Complex.ZERO;
            }
        }
    }

    /**
     * Transform Hessenberg matrix to Schur form using Single-Shift QR Algorithm.
     */
    private void transformToSchur() {
        final int n = matrixT.length;
        final double norm = getNorm();
        final ComplexShiftInfo shift = new ComplexShiftInfo();

        int iteration = 0;
        int iu = n - 1;
        while (iu > 0) {
            final int il = findSmallSubDiagonalElement(iu, norm);

            if (il == iu) {
                // Convergence: Add accumulated shift back to eigenvalue
                matrixT[iu][iu] = matrixT[iu][iu].add(shift.exShift);
                iu--;
                iteration = 0;
            } else {
                computeShift(iu, iteration, shift);
                if (++iteration > MAX_ITERATIONS) {
                    throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED,
                                                        MAX_ITERATIONS);
                }

                performQRStep(il, iu, shift);
            }
        }
        if (n > 0) {
            matrixT[0][0] = matrixT[0][0].add(shift.exShift);
        }
    }

    private double getNorm() {
        double norm = 0.0;
        for (final Complex[] complexes : matrixT) {
            for (int j = 0; j < matrixT.length; j++) {
                norm += complexes[j].norm();
            }
        }
        return norm;
    }

    private int findSmallSubDiagonalElement(final int startIdx, final double norm) {
        int l = startIdx;
        while (l > 0) {
            double s = matrixT[l - 1][l - 1].norm() + matrixT[l][l].norm();
            if (s == 0.0) {
                s = norm;
            }
            if (matrixT[l][l - 1].norm() < epsilon * s) {
                break;
            }
            l--;
        }
        return l;
    }

    private void computeShift(final int idx, final int iteration,
        final ComplexShiftInfo shift) {
      // If iterations are high, use ad-hoc (exceptional) shifts to break cycles.
      // We check iterations 10, 20, and 30 to inject shifts if convergence is stalling.
      if (iteration == 10 || iteration == 20 || iteration == 30) {
        // Accumulate the previous shift into the total shift
        shift.exShift = shift.exShift.add(shift.x);

        // Calculate the magnitude of the sub-diagonal elements to determine the scale
        // Safe check for bounds is included.
        double norm1 = matrixT[idx][idx - 1].norm();
        double norm2 = (idx > 1) ? matrixT[idx - 1][idx - 2].norm() : 0.0;
        double s = norm1 + norm2;

        // Vary the shift strategy to ensure we don't get stuck in a specific cycle.
        // - Iteration 10 & 20: Use the classic empirical shift factor 0.75
        // - Iteration 30: Use 0.964, a heuristic value (from JAMA/EISPACK) known to break
        //   cycles that are stable under 0.75.
        double alpha = (iteration == 30) ? 0.964 : 0.75;

        shift.x = new Complex(alpha * s);

        // Subtract this new exceptional shift from the diagonal to stabilize the matrix
        for (int i = 0; i <= idx; i++) {
          matrixT[i][i] = matrixT[i][i].subtract(shift.x);
        }
        return;
      }

      // Standard Wilkinson Shift
      // Computes the eigenvalues of the bottom 2x2 block and selects the one closest to T[idx][idx]
      Complex a = matrixT[idx - 1][idx - 1];
      Complex b = matrixT[idx - 1][idx];
      Complex c = matrixT[idx][idx - 1];
      Complex d = matrixT[idx][idx];

      Complex tr = a.add(d);
      Complex det = a.multiply(d).subtract(b.multiply(c));

      // Roots of characteristic equation: x^2 - tr*x + det = 0
      Complex discriminant = tr.multiply(tr).subtract(det.multiply(4));
      Complex sqrtDisc = discriminant.sqrt();

      Complex root1 = tr.add(sqrtDisc).multiply(0.5);
      Complex root2 = tr.subtract(sqrtDisc).multiply(0.5);

      // Choose the root closer to d (the bottom-right element) to maximize convergence
      if (root1.subtract(d).norm() < root2.subtract(d).norm()) {
        shift.x = root1;
      } else {
        shift.x = root2;
      }
    }

    private void performQRStep(final int il, final int iu, final ComplexShiftInfo shift) {
        Complex x = matrixT[il][il].subtract(shift.x);
        Complex y = matrixT[il + 1][il];

        for (int k = il; k < iu; k++) {
            if (k > il) {
                x = matrixT[k][k - 1];
                y = matrixT[k + 1][k - 1];
            }

            double absX = x.norm();
            double absY = y.norm();

            if (absY == 0) {
                continue;
            }

            double norm = FastMath.hypot(absX, absY);
            if (norm == 0) {
                continue;
            }

            double c = absX / norm;
            Complex alpha = (absX == 0) ? Complex.ONE : x.divide(absX);
            Complex s = alpha.multiply(y.conjugate()).divide(norm);

            // Apply G to T (from left) -> Row ops
            // IMPORTANT: If k > il, we must update the column k-1 (the bulge) as well.
            int startCol = (k == il) ? k : k - 1;
            for (int j = startCol; j < matrixT.length; j++) {
                Complex t1 = matrixT[k][j];
                Complex t2 = matrixT[k + 1][j];
                matrixT[k][j] = t1.multiply(c).add(s.multiply(t2));
                matrixT[k + 1][j] = t1.multiply(s.conjugate()).negate().add(t2.multiply(c));
            }

            // Apply G^H to T (from right) -> Col ops
            for (int i = 0; i <= FastMath.min(iu, k + 2); i++) {
                Complex t1 = matrixT[i][k];
                Complex t2 = matrixT[i][k + 1];
                matrixT[i][k] = t1.multiply(c).add(t2.multiply(s.conjugate()));
                matrixT[i][k + 1] = t1.multiply(s).negate().add(t2.multiply(c));
            }

            // Accumulate P
            for (int i = 0; i < matrixP.length; i++) {
                Complex p1 = matrixP[i][k];
                Complex p2 = matrixP[i][k + 1];
                matrixP[i][k] = p1.multiply(c).add(p2.multiply(s.conjugate()));
                matrixP[i][k + 1] = p1.multiply(s).negate().add(p2.multiply(c));
            }

            if (k > il) {
                matrixT[k + 1][k - 1] = Complex.ZERO;
            }
        }
    }

    private static class ComplexShiftInfo {

        /** Current shift. */
        private Complex x = Complex.ZERO;

        /** Accumulated shift. */
        private Complex exShift = Complex.ZERO;

    }

    public static Complex[] getEigenvalues(FieldMatrix<Complex> matrix) {
      // 1. Perform the Schur Decomposition
      ComplexSchurTransformer transformer = new ComplexSchurTransformer(matrix);

      // 2. Get the upper triangular matrix T
      FieldMatrix<Complex> T = transformer.getT();

      // 3. Extract the diagonal entries
      int n = T.getRowDimension();
      Complex[] eigenvalues = new Complex[n];

      for (int i = 0; i < n; i++) {
          eigenvalues[i] = T.getEntry(i, i);
      }

      return eigenvalues;
  }
}
