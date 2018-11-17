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
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;

/**
 *
 * This solver computes the solution using the following approach:
 *
 * 1. Compute the Hamiltonian matrix 2. Extract its complex eigen vectors (not
 * the best solution, a better solution would be ordered Schur transformation)
 * 3. Approximate the initial solution given by 2 using the Kleinman algorithm
 * (an iterative method)
 */
public class RiccatiEquationSolverImpl implements RiccatiEquationSolver {

    /** Internally used maximum iterations. */
    private static final int MAX_ITERATIONS = 100;

    /** Internally used epsilon criteria. */
    private static final double EPSILON = 1e-8;

    /** The solution of the algebraic Riccati equation. */
    private final RealMatrix P;

    /** The computed K. */
    private final RealMatrix K;

    /**
     * Constructor of the solver. A and B should be compatible. B and R must be
     * multiplicative compatible. A and Q must be multiplicative compatible. R
     * must be invertible.
     *
     * @param A state transition matrix
     * @param B control multipliers matrix
     * @param Q state cost matrix
     * @param R control cost matrix
     */
    public RiccatiEquationSolverImpl(final RealMatrix A, final RealMatrix B,
                                     final RealMatrix Q, final RealMatrix R) {

        // checking A
        if (!A.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   A.getRowDimension(), A.getColumnDimension());
        }
        if (A.getColumnDimension() != B.getRowDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   A.getRowDimension(), B.getRowDimension());
        }
        MatrixUtils.checkMultiplicationCompatible(B, R);
        MatrixUtils.checkMultiplicationCompatible(A, Q);

        // checking R
        final SingularValueDecomposition svd = new SingularValueDecomposition(R);
        if (!svd.getSolver().isNonSingular()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SINGULAR_MATRIX);
        }

        final RealMatrix R_inv = svd.getSolver().getInverse();

        P = computeP(A, B, Q, R, R_inv, MAX_ITERATIONS, EPSILON);

        K = R_inv.multiplyTransposed(B).multiply(P);
    }

    /**
     * Compute an initial stable solution and then applies the Kleinman
     * algorithm to approximate it using an EPSILON.
     *
     * @param A state transition matrix
     * @param B control multipliers matrix
     * @param Q state cost matrix
     * @param R control cost matrix
     * @param R_inv inverse of matrix R
     * @param maxIterations maximum number of iterations
     * @param epsilon epsilon to be used
     * @return matrix P, solution of the algebraic Riccati equation
     */
    private RealMatrix computeP(final RealMatrix A, final RealMatrix B,
                                final RealMatrix Q, final RealMatrix R, final RealMatrix R_inv,
                                final int maxIterations, final double epsilon) {
        final RealMatrix P_ = computeInitialP(A, B, Q, R_inv);
        return approximateP(A, B, Q, R, R_inv, P_, maxIterations, epsilon);
    }

    /**
     * Compute initial P using the Hamiltonian and the ordered eigen values
     * decomposition.
     *
     * @param A state transition matrix
     * @param B control multipliers matrix
     * @param Q state cost matrix
     * @param R_inv inverse of matrix R
     * @return initial solution
     */
    private RealMatrix computeInitialP(final RealMatrix A, final RealMatrix B,
                                       final RealMatrix Q, final RealMatrix R_inv) {
        final RealMatrix B_tran = B.transpose();

        // computing the Hamiltonian Matrix
        final RealMatrix m11 = A;
        final RealMatrix m12 = B.multiply(R_inv).multiply(B_tran).scalarMultiply(-1).scalarAdd(0);
        final RealMatrix m21 = Q.scalarMultiply(-1).scalarAdd(0);
        final RealMatrix m22 = A.transpose().scalarMultiply(-1).scalarAdd(0);
        if (m11.getRowDimension() != m12.getRowDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   m11.getRowDimension(), m12.getRowDimension());
        }
        if (m21.getRowDimension() != m22.getRowDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   m21.getRowDimension(), m22.getRowDimension());
        }
        if (m11.getColumnDimension() != m21.getColumnDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   m11.getColumnDimension(), m21.getColumnDimension());
        }
        if (m21.getColumnDimension() != m22.getColumnDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   m21.getColumnDimension(), m22.getColumnDimension());
        }

        // defining M
        final RealMatrix m = MatrixUtils.createRealMatrix(m11.getRowDimension() + m21.getRowDimension(),
                                                          m11.getColumnDimension() + m12.getColumnDimension());
        // defining submatrixes
        m.setSubMatrix(m11.getData(), 0, 0);
        m.setSubMatrix(m12.getData(), 0, m11.getColumnDimension());
        m.setSubMatrix(m21.getData(), m11.getRowDimension(), 0);
        m.setSubMatrix(m22.getData(), m11.getRowDimension(), m11.getColumnDimension());

        // eigen decomposition
        // numerically bad, but it is used for the initial stable solution for
        // the
        // Kleinman Algorithm
        // it must be ordered in order to work with submatrices
        final OrderedComplexEigenDecomposition eigenDecomposition = new OrderedComplexEigenDecomposition(m);
        final FieldMatrix<Complex> u = eigenDecomposition.getV();

        // solving linear system
        // P = U_21*U_11^-1
        final FieldMatrix<Complex> u11 = u.getSubMatrix(0,
                                                        m11.getRowDimension() - 1, 0, m11.getColumnDimension() - 1);
        final FieldMatrix<Complex> u21 = u.getSubMatrix(m11.getRowDimension(),
                                                        2 * m11.getRowDimension() - 1, 0, m11.getColumnDimension() - 1);

        final FieldDecompositionSolver<Complex> solver = new FieldLUDecomposition<Complex>(u11).getSolver();

        if (!solver.isNonSingular()) {
            throw new MathRuntimeException(LocalizedCoreFormats.SINGULAR_MATRIX);
        }

        // solving U_11^{-1}
        FieldMatrix<Complex> u11_inv = solver.getInverse();

        // P = U_21*U_11^-1
        FieldMatrix<Complex> p = u21.multiply(u11_inv);

        // converting to realmatrix - ignoring precision errors in imaginary
        // components
        return convertToRealMatrix(p, Double.MAX_VALUE);

    }

    /**
     * Applies the Kleinman's algorithm.
     *
     * @param A state transition matrix
     * @param B control multipliers matrix
     * @param Q state cost matrix
     * @param R control cost matrix
     * @param R_inv inverse of matrix R
     * @param initialP initial solution
     * @param maxIterations maximum number of iterations allowed
     * @param epsilon convergence threshold
     * @return improved solution
     */
    private RealMatrix approximateP(final RealMatrix A, final RealMatrix B,
                                    final RealMatrix Q, final RealMatrix R, final RealMatrix R_inv,
                                    final RealMatrix initialP, final int maxIterations, final double epsilon) {
        RealMatrix K_ = null;
        RealMatrix P_ = initialP;

        double error = 1;
        int i = 1;
        while (error > epsilon) {
            K_ = P_.multiply(B).multiply(R_inv).scalarMultiply(-1);

            // X = AA+BB*K1';
            final RealMatrix X = A.add(B.multiplyTransposed(K_));
            // Y = -K1*RR*K1' - QQ;
            final RealMatrix Y = K_.multiply(R).multiplyTransposed(K_).scalarMultiply(-1).subtract(Q);

            final Array2DRowRealMatrix X_ = (Array2DRowRealMatrix) X.transpose();
            final Array2DRowRealMatrix Y_ = (Array2DRowRealMatrix) Y;
            final Array2DRowRealMatrix eyeX =
                            (Array2DRowRealMatrix) MatrixUtils.createRealIdentityMatrix(X_.getRowDimension());

            // X1=kron(X',eye(size(X))) + kron(eye(size(X)),X');
            final RealMatrix X__ = X_.kroneckerProduct(eyeX).add(eyeX.kroneckerProduct(X_));
            // Y1=reshape(Y,prod(size(Y)),1); %%stack
            final RealMatrix Y__ = Y_.stack();

            // PX = inv(X1)*Y1;
            // sensitive to numerical erros
            // final RealMatrix PX = MatrixUtils.inverse(X__).multiply(Y__);
            DecompositionSolver solver = new LUDecomposition(X__).getSolver();
            if (!solver.isNonSingular()) {
                throw new MathRuntimeException(LocalizedCoreFormats.SINGULAR_MATRIX);
            }
            final RealMatrix PX = solver.solve(Y__);

            // P = reshape(PX,sqrt(length(PX)),sqrt(length(PX))); %%unstack
            final RealMatrix P__ = ((Array2DRowRealMatrix) PX).unstackSquare();

            // aerror = norm(P - P1);
            final RealMatrix diff = P__.subtract(P_);
            // calculationg l2 norm
            final SingularValueDecomposition svd = new SingularValueDecomposition(diff);
            error = svd.getNorm();

            P_ = P__;
            i++;
            if (i > maxIterations) {
                throw new MathRuntimeException(LocalizedCoreFormats.CONVERGENCE_FAILED);
            }
        }

        return P_;
    }

    /** {inheritDoc} */
    @Override
    public RealMatrix getP() {
        return P;
    }

    /** {inheritDoc} */
    @Override
    public RealMatrix getK() {
        return K;
    }

    /**
     * Converts a given complex matrix into a real matrix taking into account a
     * precision for the imaginary components.
     *
     * @param matrix complex field matrix
     * @param tolerance tolerance on the imaginary part
     * @return real matrix.
     */
    private RealMatrix convertToRealMatrix(FieldMatrix<Complex> matrix, double tolerance) {
        final RealMatrix toRet = MatrixUtils.createRealMatrix(matrix.getRowDimension(), matrix.getRowDimension());
        for (int i = 0; i < toRet.getRowDimension(); i++) {
            for (int j = 0; j < toRet.getColumnDimension(); j++) {
                Complex c = matrix.getEntry(i, j);
                if (c.getImaginary() != 0 && FastMath.abs(c.getImaginary()) > tolerance) {
                    throw new MathRuntimeException(LocalizedCoreFormats.COMPLEX_CANNOT_BE_CONSIDERED_A_REAL_NUMBER,
                                                   c.getReal(), c.getImaginary());
                }
                toRet.setEntry(i, j, c.getReal());
            }
        }
        return toRet;
    }

}
