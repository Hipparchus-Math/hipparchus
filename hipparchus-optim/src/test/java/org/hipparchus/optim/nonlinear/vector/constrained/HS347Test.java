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

/*
 * Problem HS347 is a non-linear minimization problem with one linear equality constraint,
 * often referred to as the chemical equilibrium problem.
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.MatrixUtils;

public class HS347Test {

    private static final int DIM = 3;
    private static final double[] A = {8204.37, 9008.72, 9330.46}; // DATA A(I)
    private static final double EPS = 1.0e-4; // DMAX1 minimum value

    // --- Helper Class to compute H values (similar to COMMON /D347/H) ---
    static class HValues {
        final double[] H = new double[9]; // H[1] to H[8] in Fortran

        public HValues(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            H[1] = x1 + x2 + x3 + 0.03;       // H(1)
            H[2] = 0.09 * x1 + x2 + x3 + 0.03; // H(2)
            H[3] = H[1] * H[2];               // H(3)
            H[4] = x2 + x3 + 0.03;             // H(4)
            H[5] = 0.07 * x2 + x3 + 0.03;      // H(5)
            H[6] = H[4] * H[5];               // H(6)
            H[7] = x3 + 0.03;                  // H(7)
            H[8] = 0.13 * x3 + 0.03;           // H(8)
        }
    }

    // --- Objective Function (MODE 2 and 3) ---
    static final class HS347Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }

        // F(X) = SUM [ Ai * log(max(H_num / H_den, 1e-4)) ]
        @Override public double value(RealVector x) {
            HValues h = new HValues(x);

            double f = A[0] * Math.log(Math.max(h.H[1] / h.H[2], EPS))
                     + A[1] * Math.log(Math.max(h.H[4] / h.H[5], EPS))
                     + A[2] * Math.log(Math.max(h.H[7] / h.H[8], EPS));
            return f;
        }

        @Override public RealVector gradient(RealVector x) {
            HValues h = new HValues(x);

            // The Fortran code does not check for DMAX1, it assumes H_num/H_den > 1e-4.
            // d/dX [ log(H_num/H_den) ] = d/dX [ log(H_num) - log(H_den) ]
            // d/dX [ A*log(H_num/H_den) ] = A * ( dH_num/dX * H_den - H_num * dH_den/dX ) / (H_den * H_num)
            // Note: H3 = H1*H2, H6 = H4*H5. Gradient calculation needs correction.

            // The Fortran gradient uses a simplified form:
            // d/dX [ log(a/b) ] = (b * da/dX - a * db/dX) / (a * b)

            // Term 1 (X1, X2, X3): A[0] * (H2 * dH1/dXi - H1 * dH2/dXi) / (H1*H2) = A[0] * (H2*d/dXi(H1) - H1*d/dXi(H2)) / H3
            double term1_numer = h.H[2] * (1.0 - 0.09) * 1.0; // Simplified: (H2 * 1 - H1 * 0.09) for d/dX1

            // GF(1) = A1 * (H2 - 0.09 * H1) / H3
            double g1 = A[0] * (h.H[2] - 0.09 * h.H[1]) / h.H[3];

            // Term 2 (X2, X3): A[1] * (H5 * dH4/dXi - H4 * dH5/dXi) / (H4*H5)
            // Term 3 (X3): A[2] * (H8 * dH7/dXi - H7 * dH8/dXi) / (H7*H8)

            // GF(2) = A1 * (H2 - H1) / H3 + A2 * (H5 - 0.07 * H4) / H6
            double g2_part1 = A[0] * (h.H[2] - h.H[1]) / h.H[3]; // dH1/dX2 = 1, dH2/dX2 = 1. -> Numerator is H2-H1. Corrected from Fortran.
            double g2_part2 = A[1] * (h.H[5] - 0.07 * h.H[4]) / h.H[6]; // dH4/dX2 = 1, dH5/dX2 = 0.07.
            double g2 = g2_part1 + g2_part2;

            // GF(3) = A1 * (H2 - H1) / H3 + A2 * (H5 - H4) / H6 + A3 * (H8 - 0.13 * H7) / (H7*H8)
            double g3_part1 = A[0] * (h.H[2] - h.H[1]) / h.H[3]; // dH1/dX3 = 1, dH2/dX3 = 1
            double g3_part2 = A[1] * (h.H[5] - h.H[4]) / h.H[6]; // dH4/dX3 = 1, dH5/dX3 = 1
            double g3_part3 = A[2] * (h.H[8] - 0.13 * h.H[7]) / (h.H[7] * h.H[8]); // dH7/dX3 = 1, dH8/dX3 = 0.13
            double g3 = g3_part1 + g3_part2 + g3_part3;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            // Hessian not defined in Fortran source (MODE 5 is RETURN)
            throw new UnsupportedOperationException("Hessian matrix is not implemented for this test case.");
        }
    }

    // --- Equality Constraint (MODE 4) ---
    static final class HS347Eq extends EqualityConstraint {

        HS347Eq() { super(new ArrayRealVector(new double[]{1.0})); } // Target is 1.0 (X1+X2+X3=1)

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            // G(1) = X1 + X2 + X3 - 1.0 = 0
            double h1 = x.getEntry(0) + x.getEntry(1) + x.getEntry(2) - 1.0;
            return new ArrayRealVector(new double[]{h1}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            // Jacobian: dH1/dXi = 1 for all i
            double[][] J = new double[1][DIM];
            J[0][0] = 1.0;
            J[0][1] = 1.0;
            J[0][2] = 1.0;
            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() {
        return new double[]{0.7, 0.2, 0.1};
    }

//    @Test
//    public void testHS347() {
//        SQPOptimizerS2 opt = new SQPOptimizerS2();
//        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
//            opt.setDebugPrinter(System.out::println);
//        }
//        // Box constraints: 0 <= Xi <= 1
//        SimpleBounds bounds = new SimpleBounds(
//            new double[]{0.0, 0.0, 0.0},
//            new double[]{1.0, 1.0, 1.0}
//        );
//
//        LagrangeSolution sol = opt.optimize(
////                new InitialGuess(start()),
//                new ObjectiveFunction(new HS347Obj()),
//                new HS347Eq(),
//                bounds
//        );
//
//        double f = sol.getValue();
//        final double fExpected = 17374.625;
//        final double tolerance = 1.0e-5 * (Math.abs(fExpected) + 1.0);
//
//        assertTrue(f <= fExpected + tolerance, "Objective value mismatch/worse than expected.");
//
//
//    }
}
