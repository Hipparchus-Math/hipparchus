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
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS361 (TP361) – Nonlinear problem with 5 variables and 6 nonlinear inequality constraints.
 * Based on Hock & Schittkowski test problem TP361.
 */
public class HS361Test {

    private static final int DIM = 5;
    private static final int NUM_INEQUALITIES = 6;

    // Coefficient arrays A, B, C, D (Fortran DATA blocks, 1-based → 0-based).
    private static final double[] A = {
        -0.8720288849e7, // A(1)
         0.1505125253e6, // A(2)
        -0.1566950325e3, // A(3)
         0.4764703222e6, // A(4)
         0.7294828271e6  // A(5)
    };

    private static final double[] B = {
        -0.145421402e6,
         0.29311506e4,
        -0.40427932e2,
         0.5106192e4,
         0.1571136e5
    };

    private static final double[] C_COEFF = {
        -0.1550111084e6,
         0.436053352e4,
         0.129492344e2,
         0.10236884e5,
         0.13176786e5
    };

    private static final double[] D = {
        -0.3266695104e6,
         0.739068412e4,
        -0.278986976e2,
         0.16643076e5,
         0.30988146e5
    };

    /**
     * Objective function for HS361.
     *
     * Fortran (MODE 2):
     *
     *   FX = A(1)
     *   DO I = 2,5
     *       FX = FX + A(I) * X(I)
     *   END DO
     *   FX = X(1) * FX - 0.24345D+5
     *   FX = -FX
     *
     * Definendo S = A(1) + Σ_{i=2..5} A(i)*X(i), abbiamo:
     *
     *   FX = -( X(1) * S - 24345 ) = 24345 - X(1) * S
     *
     * Quindi la funzione da minimizzare è esattamente:
     *   F(x) = 24345.0 - x1 * ( A1 + Σ_{i=2..5} Ai * Xi )
     */
    private static class HS361Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector X) {
            final double x1 = X.getEntry(0);

            // S = A(1) + Σ_{i=2..5} A(i) * X(i)
            double sumAX = A[0];
            for (int i = 1; i < DIM; i++) {
                sumAX += A[i] * X.getEntry(i);
            }

            // FX = 24345 - x1 * S
            return 24345.0 - x1 * sumAX;
        }

        @Override
        public RealVector gradient(RealVector X) {
            final double x1 = X.getEntry(0);
            final double[] g = new double[DIM];

            // S = A(1) + Σ_{i=2..5} A(i) * X(i)
            double sumAX = A[0];
            for (int i = 1; i < DIM; i++) {
                sumAX += A[i] * X.getEntry(i);
            }

            // From F(x) = 24345 - x1 * S:
            // ∂F/∂x1 = -S
            g[0] = -sumAX;

            // For i = 2..5 (Java 1..4): ∂F/∂x_i = -x1 * A(i)
            for (int i = 1; i < DIM; i++) {
                g[i] = -x1 * A[i];
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector X) {
            // Hessian: H_1i = H_i1 = -A(i) per i=2..5, tutto il resto 0.
            final RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);

            for (int i = 1; i < DIM; i++) {
                H.setEntry(0, i, -A[i]); // ∂²F / (∂x1 ∂xi)
                H.setEntry(i, 0, -A[i]); // ∂²F / (∂xi ∂x1)
            }

            // Tutti gli altri termini sono 0.
            return H;
        }
    }

    /**
     * Nonlinear inequality constraints:
     *
     * Fortran (MODE 4):
     *
     *   H(1) = X(1) * ( B(1) + Σ_{i=2..5} B(i)*X(i) )
     *   H(2) = X(1) * ( C(1) + Σ_{i=2..5} C(i)*X(i) )
     *   H(3) = X(1) * ( D(1) + Σ_{i=2..5} D(i)*X(i) )
     *
     *   G(1) =  H(1)
     *   G(2) =  H(2)
     *   G(3) =  H(3)
     *   G(4) = 0.294D+5 - H(1) = 29400 - H(1)
     *   G(5) = 0.294D+5 - H(2) = 29400 - H(2)
     *   G(6) = 0.2772D+6 - H(3) = 277200 - H(3)
     *
     * Tutte le G(j) sono del tipo G(j) ≥ 0.
     */
    private static class HS361Ineq extends InequalityConstraint {

        HS361Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQUALITIES])); // RHS = 0 per tutte
        }

        private double[] calculateH(RealVector X) {
            final double x1 = X.getEntry(0);
            double[] h = new double[3];

            double sumB = B[0];
            double sumC = C_COEFF[0];
            double sumD = D[0];

            for (int i = 1; i < DIM; i++) {
                final double xi = X.getEntry(i);
                sumB += B[i] * xi;
                sumC += C_COEFF[i] * xi;
                sumD += D[i] * xi;
            }

            h[0] = x1 * sumB;
            h[1] = x1 * sumC;
            h[2] = x1 * sumD;

            return h;
        }

        private RealMatrix calculateJacobianH(RealVector X) {
            final double x1 = X.getEntry(0);
            final RealMatrix HH = new Array2DRowRealMatrix(3, DIM);

            // Colonna 1: dH_j/dx1 = coeff(1) + Σ_{i=2..5} coeff(i)*X(i)
            double sumB = B[0];
            double sumC = C_COEFF[0];
            double sumD = D[0];

            for (int i = 1; i < DIM; i++) {
                final double xi = X.getEntry(i);
                sumB += B[i] * xi;
                sumC += C_COEFF[i] * xi;
                sumD += D[i] * xi;
            }

            HH.setEntry(0, 0, sumB); // dH1/dx1
            HH.setEntry(1, 0, sumC); // dH2/dx1
            HH.setEntry(2, 0, sumD); // dH3/dx1

            // Colonne 2..5: dH_j/dx_i = coeff_j(i) * X(1)
            final double[][] coeffs = { B, C_COEFF, D };
            for (int j = 0; j < 3; j++) {
                final double[] Cj = coeffs[j];
                for (int i = 1; i < DIM; i++) {
                    HH.setEntry(j, i, Cj[i] * x1);
                }
            }

            return HH;
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector X) {
            final double[] h = calculateH(X);
            final double[] g = new double[NUM_INEQUALITIES];

            // G(1..3) = H(j) >= 0
            g[0] = h[0];
            g[1] = h[1];
            g[2] = h[2];

            // G(4..6) = C_j - H(j) >= 0
            g[3] = 29400.0  - h[0]; // 0.294D+5
            g[4] = 29400.0  - h[1]; // 0.294D+5
            g[5] = 277200.0 - h[2]; // 0.2772D+6

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(RealVector X) {
            final RealMatrix HH = calculateJacobianH(X); // 3x5
            final RealMatrix GG = new Array2DRowRealMatrix(NUM_INEQUALITIES, DIM);

            // G(1..3): grad(Gj) = grad(Hj)
            for (int j = 0; j < 3; j++) {
                GG.setRowVector(j, HH.getRowVector(j));
            }

            // G(4..6): grad(Gj) = -grad(H_{j-3})
            for (int j = 0; j < 3; j++) {
                GG.setRowVector(j + 3, HH.getRowVector(j).mapMultiply(-1.0));
            }

            return GG;
        }
    }

    // Fortran initial guess (MODE 1)
    private static final double[] X_START = {
        2.52,  // 0.252D+1
        2.0,   // 0.2D+1
        37.5,  // 0.375D+2
        9.25,  // 0.925D+1
        6.8    // 0.68D+1
    };

    // Fortran FEX
    private static final double F_EXPECTED = -0.77641212e6;

    // Fortran XEX
    private static final double[] X_EXPECTED = {
        0.68128605,
        2.4,
        20.0,
        9.3,
        7.0
    };

    @Test
    public void testHS361Optimization() {
        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        final HS361Obj objective = new HS361Obj();
        final HS361Ineq inequality = new HS361Ineq();

      
        // Bounds: traduzione diretta da LXL/LXU, XL, XU.
        final double[] lowerBounds = {
            0.0,                     // XL(1) = 0.0, LXL(1)=TRUE
            1.2,                     // XL(2) = 1.2, LXL(2)=TRUE
            20.0,                    // XL(3) = 20.0, LXL(3)=TRUE
            9.0,                     // XL(4) = 9.0, LXL(4)=TRUE
            Double.NEGATIVE_INFINITY // LXL(5)=FALSE → nessun lower bound
        };

        final double[] upperBounds = {
            Double.POSITIVE_INFINITY, // LXU(1)=FALSE → nessun upper bound
            2.4,                      // XU(2) = 2.4
            60.0,                     // XU(3) = 60.0
            9.3,                      // XU(4) = 9.3
            7.0                       // XU(5) = 7.0
        };
         
        final SimpleBounds bounds = new SimpleBounds(lowerBounds, upperBounds);
        final LagrangeSolution sol = opt.optimize(
                HSProblemTestUtils.newCentralDifferenceOption(),
                new InitialGuess(X_START),
                new ObjectiveFunction(objective),
                inequality,
                bounds
        );

        HSProblemTestUtils.assertExpectedObjective(F_EXPECTED, sol);
    }
}