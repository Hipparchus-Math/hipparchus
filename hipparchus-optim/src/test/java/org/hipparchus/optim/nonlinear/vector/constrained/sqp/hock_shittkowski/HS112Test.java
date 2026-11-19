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

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS112Test {

    /** Constants C(i) from TP112. */
    private static final double[] C = {
        -6.089, -17.164, -34.054, -5.914, -24.721,
        -14.986, -24.100, -10.708, -26.662, -22.179
    };

    /** Reference optimum value f(x*). */
    private static final double F_REF = -0.47761086e2; // -47.761086

    /** Build A and b for the 3 linear equalities A·x = b. */
    private static RealMatrix buildA() {
        double[][] A = new double[3][10];
        // eq1: x1 + 2*x2 + 2*x3 + x6 + x10 = 2
        A[0][0] = 1;  A[0][1] = 2;  A[0][2] = 2;  A[0][5] = 1;  A[0][9] = 1;
        // eq2: x4 + 2*x5 + x6 + x7 = 1
        A[1][3] = 1;  A[1][4] = 2;  A[1][5] = 1;  A[1][6] = 1;
        // eq3: x3 + x7 + x8 + 2*x9 + x10 = 1
        A[2][2] = 1;  A[2][6] = 1;  A[2][7] = 1;  A[2][8] = 2;  A[2][9] = 1;
        return new Array2DRowRealMatrix(A, false);
    }
    private static RealVector buildB() {
        return new ArrayRealVector(new double[]{2.0, 1.0, 1.0}, false);
    }

    /** Objective f(x) = sum_i x_i * (C_i + log(x_i) - log(sum x)). Gradient provided. */
    static final class HS112Objective extends TwiceDifferentiableFunction {
        private static final double T_MIN = 1e-5;

        @Override public int dim() { return 10; }

        @Override public double value(RealVector x) {
            // Fortran fallback if sum too small or some xi < 0
            double sum = 0;
            for (int i = 0; i < 10; i++) sum += x.getEntry(i);
            boolean bad = (sum < T_MIN);
            for (int i = 0; i < 10 && !bad; i++) if (x.getEntry(i) < 0.0) bad = true;

            if (bad) {
                double s = 0;
                for (int i = 0; i < 10; i++) {
                    double xi = x.getEntry(i);
                    if (xi < 0.0) {
                        double d = xi - 5.0;
                        s += d * d;
                    }
                }
                return (s + 1.0e3 - 47.8);
            }

            final double logSum = Math.log(sum);
            double f = 0.0;
            for (int i = 0; i < 10; i++) {
                double xi = x.getEntry(i);
                // Bounds keep xi >= 1e-4, so log(xi) is safe in nominal path
                f += xi * (C[i] + Math.log(xi) - logSum);
            }
            return f;
        }

        @Override public RealVector gradient(RealVector x) {
            double sum = 0;
            for (int i = 0; i < 10; i++) sum += x.getEntry(i);
            boolean bad = (sum < T_MIN);
            for (int i = 0; i < 10 && !bad; i++) if (x.getEntry(i) < 0.0) bad = true;

            double[] g = new double[10];
            if (bad) {
                for (int i = 0; i < 10; i++) {
                    double xi = x.getEntry(i);
                    g[i] = (xi < 0.0) ? 2.0 * (xi - 5.0) : 0.0;
                }
                return new ArrayRealVector(g, false);
            }

            // grad_i = C_i + log(x_i) - log(sum x)
            final double logSum = Math.log(sum);
            for (int i = 0; i < 10; i++) {
                g[i] = C[i] + Math.log(x.getEntry(i)) - logSum;
            }
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            // Not required by the test suite
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

  
static final class HS112Eq extends EqualityConstraint {
    private final RealMatrix A;
    private final RealVector b;
    HS112Eq() {
        super(new ArrayRealVector(new double[] { 0.0, 0.0, 0.0 }));  // explicit zero RHS
        this.A = buildA();
        this.b = buildB();
    }

    @Override
    public int dim() { return 10; }

    @Override
    public RealVector value(RealVector x) {
        return A.operate(x).subtract(b);
    }

    @Override
    public RealMatrix jacobian(RealVector x) {
        return A;
    }
}

    /** Solve helper (adds requested debug printer). */
    static LagrangeSolution solve() {
        final double INF = Double.NEGATIVE_INFINITY;
        final double SUP = Double.POSITIVE_INFINITY;

        double[] start = new double[10];
        double[] lo    = new double[10];
        double[] up    = new double[10];
        for (int i = 0; i < 10; i++) {
            start[i] = 0.1;     // Fortran initial guess
            lo[i]    = 1.0e-4;  // XL(i) = 1e-4
            up[i]    = SUP;     // unbounded above
        }

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        return opt.optimize(
            new InitialGuess(start),
            new ObjectiveFunction(new HS112Objective()),
            new HS112Eq(),
            new SimpleBounds(lo, up)
        );
    }

    @Test
    public void testHS112() {
        LagrangeSolution sol = solve();
        double f = sol.getValue();
        HSProblemTestUtils.assertExpectedObjective(F_REF , sol);
    }
}
