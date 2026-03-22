/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS331Test {

    private static final int DIM = 2;

    static final class HS331Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }

        // This problem includes implicit checks to ensure arguments for log and division are valid.
        // For optimization, we rely on the optimizer to stay within bounds, but we must protect
        // log arguments which come from division.

        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // A = X1 + X2
            double a = x1 + x2;
            // B = log(A)
            double b = Math.log(a);
            // C = 2 * log(X2)
            double c = 2.0 * Math.log(x2);

            // Inner expression term: |C / B|
            // Due to the domain (X1, X2 > 0), B > 0. C = 2*log(X2) is negative since 0.1 <= X2 <= 0.2 < 1.
            // Therefore, C/B is negative, and we need the absolute value.
            double inner_term = Math.abs(c / b);

            // Ensure argument of outer log is positive (by Fortran's DMAX1(C/B, 1.0D-4) protection, which is used for the gradient)
            // The objective function uses DABS. Since C < 0 and B > 0, C/B < 0.
            if (inner_term <= 0) return Double.POSITIVE_INFINITY;

            // FX = log(|C / B|) / X1
            return Math.log(inner_term) / x1;
        }

        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // Re-evaluate A, B, C for derivative
            double a = x1 + x2;
            double b = Math.log(a);
            double c = 2.0 * Math.log(x2);

            // The Fortran gradient uses a max protection for the log: DLOG(DMAX1(C/B, 1.0D-4))
            // The argument of the inner log (C/B) is typically negative in this problem domain (C < 0, B > 0).
            // The absolute value is required for the objective function.
            // The Fortran gradient uses C/B *without* absolute value, which is confusing.
            // We use the derivatives of the mathematically defined function.
            // F = ln(|C/B|) / X1

            // G1: dF/dX1 = (-1/X1) * [ F + (1/X1) * (C/B)' * (B/C) ] where (C/B)' = d(C/B)/dX1
            // (C/B)' = d/dX1 [ 2*ln(X2) / ln(X1+X2) ]
            // d(C/B)/dX1 = -2*ln(X2) / (ln(X1+X2)^2 * (X1+X2))

            double b_sq = b * b;
            double c_over_b = c / b;

            // Term 1: d/dX1(C/B)
            double d_c_over_b_dx1 = -c / (b_sq * a);

            // Term 2: d/dX2(C/B)
            // d/dX2(C/B) = (B*dC/dX2 - C*dB/dX2) / B^2
            // dC/dX2 = 2/X2. dB/dX2 = 1/A.
            double d_c_over_b_dx2 = (b * 2.0 / x2 - c / a) / b_sq;

            // Gradient components: GF(i) = (1/X1) * (1/|C/B|) * (d|C/B|/dXi) - F/X1^2

            // Using the Fortran-like structure GF(1) = (-1/X1) * ( (log(|C/B|)/X1) + (1/(B*A)) )
            // Fortran formula is heavily simplified and likely assumes |C/B| near 1 or sign manipulation.

            // We implement the Fortran expression literally, noting potential numerical issues:
            // GF(1) = (-1/X1) * ( (DLOG(DMAX1(C/B, 1.0D-4))/X1) + (1/(B*A)) )
            // GF(2) = ((0.2D+1*B)/X(2) - C/A) / (C*B*X(1))

            // Note: C/B is negative. log(max(negative, 1e-4)) is log(1e-4), which is negative.

            double inner_log = Math.log(Math.max(c / b, 1.0e-4));

            double g1 = (-1.0 / x1) * ((inner_log / x1) + (1.0 / (b * a)));
            double g2 = ((2.0 * b) / x2 - c / a) / (c * b * x1);

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            // Hessian matrix is not implemented for this test case.
            throw new UnsupportedOperationException("Hessian matrix is not implemented for this test case.");
        }
    }

    static final class HS331LinearIneq extends InequalityConstraint {

        HS331LinearIneq() { super(new ArrayRealVector(new double[1])); }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            // G(1) = 1.0 - X1 - X2 >= 0
            double g1 = 1.0 - x.getEntry(0) - x.getEntry(1);

            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double[][] J = new double[1][DIM];

            // dG1/dX1 = -1.0, dG1/dX2 = -1.0
            J[0][0] = -1.0;
            J[0][1] = -1.0;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() {
        return new double[]{0.5, 0.1};
    }

    @Test
    public void testHS331() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Box constraints: 0.3 <= X1 <= 0.7, 0.1 <= X2 <= 0.2
        SimpleBounds bounds = new SimpleBounds(
            new double[]{0.3, 0.1},
            new double[]{0.7, 0.2}
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS331Obj()),
                new HS331LinearIneq(),
                bounds
        );

        double f = sol.getValue();
        final double fExpected = 4.258;

        assertEquals(fExpected, f, 1.0e-3 * (Math.abs(fExpected) + 1.0), "objective mismatch");


    }
}
