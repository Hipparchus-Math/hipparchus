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

/*
 * HS232 (TP232)
 *
 * N    = 2
 * NILI = 3 (three linear inequalities)
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Objective (MODE=2):
 *   HV  = sqrt(3)
 *   f(x) = - (9 - (x1 - 3)^2) * x2^3 / (27 * HV)
 *
 * Fortran constraints (G(x) >= 0):
 *
 *   G1(x) = x1 / HV - x2
 *   G2(x) = x1 + HV * x2
 *   G3(x) = 6 - x1 - HV * x2
 *
 * Bounds:
 *   x1 >= 0
 *   x2 >= 0
 *
 * Reference solution (MODE=1):
 *   x*  = (3, sqrt(3))
 *   f*  = -1
 */
public class HS232Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 3;
    private static final int NUM_EQ   = 0;

    private static final double HV = FastMath.sqrt(3.0);

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS232Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            double term = 9.0 - (x1 - 3.0) * (x1 - 3.0);
            return -term * x2 * x2 * x2 / (27.0 * HV);
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // From Fortran:
            // GF(1)= 2*(x1-3)*x2^3 / (27*HV)
            // GF(2)= -(9-(x1-3)^2)*x2^2 / (9*HV)
            double df1 = 2.0 * (x1 - 3.0) * x2 * x2 * x2 / (27.0 * HV);
            double df2 = -(9.0 - (x1 - 3.0) * (x1 - 3.0)) * x2 * x2 / (9.0 * HV);

            return new ArrayRealVector(new double[]{df1, df2}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Start with zero Hessian; BFGS/SQP will update it.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints G(x) >= 0
    // -------------------------------------------------------------------------
    private static class HS232Ineq extends InequalityConstraint {

        HS232Ineq() {
            // RHS = 0 for all inequalities (G >= 0)
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            double G1 =  x1 / HV - x2;
            double G2 =  x1 + HV * x2;
            double G3 =  6.0 - x1 - HV * x2;

            return new ArrayRealVector(new double[]{G1, G2, G3}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // G1 = x1/HV - x2  → grad = [1/HV, -1]
            J.setEntry(0, 0, 1.0 / HV);
            J.setEntry(0, 1, -1.0);

            // G2 = x1 + HV*x2  → grad = [1, HV]
            J.setEntry(1, 0, 1.0);
            J.setEntry(1, 1, HV);

            // G3 = 6 - x1 - HV*x2 → grad = [-1, -HV]
            J.setEntry(2, 0, -1.0);
            J.setEntry(2, 1, -HV);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS232_optimization() {

        // Initial guess (MODE=1): X(1)=2, X(2)=0.5
        double[] x0 = new double[]{2.0, 0.5};

        // Bounds: x1 >= 0, x2 >= 0
        double[] lower = new double[]{0.0, 0.0};
        double[] upper = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS232Obj()),
                null,             // no equalities
                new HS232Ineq(),  // 3 linear inequalities
                bounds            // lower bounds x >= 0
        );

        double f = sol.getValue();

        final double fExpected = -1.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
