/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
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

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * HS367 (TP367) – 7-variable nonlinear constrained problem.
 *
 * N = 7
 * - 3 inequality constraints (G1, G2, G3 >= 0 in the original formulation)
 * - 2 equality constraints (G4 = 0, G5 = 0)
 *
 * Decision variables: x1..x7
 *
 * Objective:
 *   f(x) =
 *     -5*x1
 *     -5*x2
 *     -4*x3
 *     - x1*x3
 *     -6*x4
 *     -5*x5 / (1 + x5)
 *     -8*x6 / (1 + x6)
 *     -10 * ( 1 - 2*exp(-x7) + exp(-2*x7) )
 *
 * Bounds:
 *   0 <= xi, no upper bounds.
 *
 * Expected reference:
 *   FEX = -0.37412960D+2  (used as an upper bound: FEX >= f, LEX = .FALSE.)
 */
public class HS367Test {

    private static final int DIM = 7;
    private static final int NUM_INEQ = 3; // G1, G2, G3
    private static final int NUM_EQ   = 2; // G4, G5

    // -------------------------------------------------------------------------
    // Objective function
    // -------------------------------------------------------------------------
    private static class HS367Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);

            double term1 = -5.0 * x1;
            double term2 = -5.0 * x2;
            double term3 = -4.0 * x3;
            double term4 = -x1 * x3;
            double term5 = -6.0 * x4;
            double term6 = -5.0 * x5 / (1.0 + x5);
            double term7 = -8.0 * x6 / (1.0 + x6);
            double term8 = -10.0 * (1.0 - 2.0 * FastMath.exp(-x7) + FastMath.exp(-2.0 * x7));

            return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8;
        }

        @Override
        public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);

            double[] g = new double[DIM];

            // df/dx1 = -5 - x3
            g[0] = -5.0 - x3;

            // df/dx2 = -5
            g[1] = -5.0;

            // df/dx3 = -4 - x1
            g[2] = -4.0 - x1;

            // df/dx4 = -6
            g[3] = -6.0;

            // df/dx5 = d/dx5 [-5*x5/(1+x5)] = -5 / (1+x5)^2
            g[4] = -5.0 / FastMath.pow(1.0 + x5, 2.0);

            // df/dx6 = d/dx6 [-8*x6/(1+x6)] = -8 / (1+x6)^2
            g[5] = -8.0 / FastMath.pow(1.0 + x6, 2.0);

            // df/dx7 = -20*(exp(-x7) - exp(-2*x7))
            double e1 = FastMath.exp(-x7);
            double e2 = FastMath.exp(-2.0 * x7);
            g[6] = -20.0 * (e1 - e2);

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);

            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);

            // Cross derivative between x1 and x3: d²/dx1dx3 of (-x1*x3) = -1
            H.setEntry(0, 2, -1.0);
            H.setEntry(2, 0, -1.0);

            // Second derivative wrt x5: d²/dx5²[-5*x5/(1+x5)] = 10/(1+x5)^3
            double denom5 = FastMath.pow(1.0 + x5, 3.0);
            H.setEntry(4, 4, 10.0 / denom5);

            // Second derivative wrt x6: d²/dx6²[-8*x6/(1+x6)] = 16/(1+x6)^3
            double denom6 = FastMath.pow(1.0 + x6, 3.0);
            H.setEntry(5, 5, 16.0 / denom6);

            // Second derivative wrt x7:
            // df/dx7 = -20*(e^{-x7} - e^{-2x7})
            // d²f/dx7² = 20*(e^{-x7} - 2*e^{-2x7})
            double e1 = FastMath.exp(-x7);
            double e2 = FastMath.exp(-2.0 * x7);
            H.setEntry(6, 6, 20.0 * (e1 - 2.0 * e2));

            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Equality constraints (G4, G5)
    // -------------------------------------------------------------------------
    private static class HS367Eq extends EqualityConstraint {

        HS367Eq() {
            // Both constraints have right-hand side 0
            super(new ArrayRealVector(new double[NUM_EQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        /**
         * Equality constraints:
         *
         * G4 = 2*x4 + x5 + 0.8*x6 + x7 - 5 = 0
         * G5 = x2^2 + x3^2 + x5^2 + x6^2 - 5 = 0
         */
        @Override
        public RealVector value(RealVector x) {
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);

            double g4 = 2.0 * x4 + x5 + 0.8 * x6 + x7 - 5.0;
            double g5 = x2 * x2 + x3 * x3 + x5 * x5 + x6 * x6 - 5.0;

            return new ArrayRealVector(new double[] { g4, g5 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);

            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);

            // G4 gradient = [0, 0, 0, 2, 1, 0.8, 1]
            J.setEntry(0, 0, 0.0);
            J.setEntry(0, 1, 0.0);
            J.setEntry(0, 2, 0.0);
            J.setEntry(0, 3, 2.0);
            J.setEntry(0, 4, 1.0);
            J.setEntry(0, 5, 0.8);
            J.setEntry(0, 6, 1.0);

            // G5 gradient = [0, 2*x2, 2*x3, 0, 2*x5, 2*x6, 0]
            J.setEntry(1, 0, 0.0);
            J.setEntry(1, 1, 2.0 * x2);
            J.setEntry(1, 2, 2.0 * x3);
            J.setEntry(1, 3, 0.0);
            J.setEntry(1, 4, 2.0 * x5);
            J.setEntry(1, 5, 2.0 * x6);
            J.setEntry(1, 6, 0.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints (G1, G2, G3 <= 0)
    // -------------------------------------------------------------------------
    private static class HS367Ineq extends InequalityConstraint {

        HS367Ineq() {
            // All right-hand sides are 0
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        /**
         * Inequality constraints (G <= 0):
         *
         * G1 = 10 - (x1 + x2 + x3 + x4 + x5 + x6 + x7) <= 0
         * G2 = 5  - (x1 + x2 + x3 + x4)               <= 0
         * G3 = 5  - (x1 + x3 + x5 + x6^2 + x7^2)      <= 0
         */
        @Override
        public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);

            double g1 = 10.0 - (x1 + x2 + x3 + x4 + x5 + x6 + x7);
            double g2 = 5.0  - (x1 + x2 + x3 + x4);
            double g3 = 5.0  - (x1 + x3 + x5 + x6 * x6 + x7 * x7);

            return new ArrayRealVector(new double[] { g1, g2, g3 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // G1 gradient = [-1, -1, -1, -1, -1, -1, -1]
            for (int j = 0; j < DIM; j++) {
                J.setEntry(0, j, -1.0);
            }

            // G2 gradient = [-1, -1, -1, -1, 0, 0, 0]
            J.setEntry(1, 0, -1.0);
            J.setEntry(1, 1, -1.0);
            J.setEntry(1, 2, -1.0);
            J.setEntry(1, 3, -1.0);
            J.setEntry(1, 4,  0.0);
            J.setEntry(1, 5,  0.0);
            J.setEntry(1, 6,  0.0);

            // G3 gradient = [-1, 0, -1, 0, -1, -2*x6, -2*x7]
            J.setEntry(2, 0, -1.0);
            J.setEntry(2, 1,  0.0);
            J.setEntry(2, 2, -1.0);
            J.setEntry(2, 3,  0.0);
            J.setEntry(2, 4, -1.0);
            J.setEntry(2, 5, -2.0 * x6);
            J.setEntry(2, 6, -2.0 * x7);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS367_optimization() {

        // Initial guess: X(I) = 0.1
        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 0.1;
        }

        // Bounds: Xi >= 0, no upper bounds
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = 0.0;
            upper[i] = Double.POSITIVE_INFINITY;
        }
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS367Obj()),
                new HS367Eq(),     // 2 equality constraints
                new HS367Ineq(),   // 3 inequality constraints
                bounds
        );

        double f = sol.getValue();

        // LEX = .FALSE. → FEX is only a reference value, not exact:
        // rule: use FEX >= f (up to numerical tolerance)
        final double fExpected = -0.37412960e2; // -37.412960
        final double tolF = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(fExpected + tolF >= f,
                   "HS367: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
