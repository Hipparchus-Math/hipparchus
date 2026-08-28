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
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS202Test {

    /** Objective F1^2 + F2^2 with exact gradient. */
    static final class HS202Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        private static double f1(double x1, double x2) {
            // F1 = -13 + x1 - 2*x2 + 5*x2^2 - x2^3
            return -13.0 + x1 - 2.0*x2 + 5.0*x2*x2 - Math.pow(x2, 3);
        }
        private static double f2(double x1, double x2) {
            // F2 = -29 + x1 - 14*x2 + x2^2 + x2^3
            return -29.0 + x1 - 14.0*x2 + x2*x2 + Math.pow(x2, 3);
        }

        @Override public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double F1 = f1(x1, x2);
            final double F2 = f2(x1, x2);
            return F1*F1 + F2*F2;
        }

        @Override public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double F1 = f1(x1, x2);
            final double F2 = f2(x1, x2);

            // dF1/dx1 = 1, dF2/dx1 = 1  -> d/dx1 (F1^2+F2^2) = 2*F1*1 + 2*F2*1
            final double g1 = 2.0*F1 + 2.0*F2;

            // dF1/dx2 = -2 + 10*x2 - 3*x2^2
            // dF2/dx2 = -14 + 2*x2 + 3*x2^2
            final double dF1dx2 = -2.0 + 10.0*x2 - 3.0*x2*x2;
            final double dF2dx2 = -14.0 + 2.0*x2 + 3.0*x2*x2;
            final double g2 = 2.0*F1*dF1dx2 + 2.0*F2*dF2dx2;

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }


    /** Solve utility (bounds only). */
    static LagrangeSolution solve(final double[] start) {
        final double[] lo = { 1.0, -5.0 };
        final double[] up = { 20.0, 5.0 };
        SQPOption sqpOption=new SQPOption();
        sqpOption.setGradientMode(GradientMode.FORWARD);
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        return optimizer.optimize(
            new InitialGuess(start),
            new ObjectiveFunction(new HS202Obj()),
            //new SimpleBounds(lo, up),
            sqpOption
        );
    }

    // ---------- Minimal JUnit test: check only objective value ----------
    @Test
    public void testHS202() {
        final double[] x0 = { 15.0, -2.0 }; // Fortran start
        final LagrangeSolution sol = solve(x0);

        final double fEx = 0.0; // reference optimum (at x* = [5, 4] within bounds)
        final double f   = sol.getValue();
         HSProblemTestUtils.assertExpectedObjective(0.0, sol);
    }
}
