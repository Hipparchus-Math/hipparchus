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
import org.junit.jupiter.api.Test;

import static java.lang.Math.*;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS107Test {

    // constants from the model
    private static final double V1 = 48.4 / 50.176;
    private static final double C  = V1 * sin(0.25);
    private static final double D  = V1 * cos(0.25);

    // -------- objective f, grad ----------
    private static final class HS107Objective extends TwiceDifferentiableFunction {
        @Override public int dim() { return 9; }

        @Override public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            // 3000*x1 + 1000*x1^3 + 2000*x2 + (2000/3)*x2^3
            return 3000.0 * x1 + 1000.0 * x1 * x1 * x1
                 + 2000.0 * x2 + (2000.0 / 3.0) * x2 * x2 * x2;
        }

        @Override public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double[] g = new double[9];
            g[0] = 3000.0 + 3000.0 * x1 * x1;        // d/dx1
            g[1] = 2000.0 + 2000.0 * x2 * x2;        // d/dx2
            // other partials are zero
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    // -------- six equalities h(x) = 0 with analytic Jacobian ----------
    private static final class HS107Eq extends EqualityConstraint {
        HS107Eq() { super(new ArrayRealVector(new double[6])); }
        @Override public int dim() { return 9; }

        @Override public RealVector value(RealVector x) {
            final double x1=x.getEntry(0), x2=x.getEntry(1), x3=x.getEntry(2), x4=x.getEntry(3),
                         x5=x.getEntry(4), x6=x.getEntry(5), x7=x.getEntry(6),
                         x8=x.getEntry(7), x9=x.getEntry(8);

            final double y1 = sin(x8);
            final double y2 = cos(x8);
            final double y3 = sin(x9);
            final double y4 = cos(x9);
            final double y5 = sin(x8 - x9);
            final double y6 = cos(x8 - x9);

            final double h1 = 0.4  - x1 + 2.0*C*x5*x5 - x5*x6*(D*y1 + C*y2) - x5*x7*(D*y3 + C*y4);
            final double h2 = 0.4  - x2 + 2.0*C*x6*x6 + x5*x6*(D*y1 - C*y2) + x6*x7*(D*y5 - C*y6);
            final double h3 = 0.8         + 2.0*C*x7*x7 + x5*x7*(D*y3 - C*y4) - x6*x7*(D*y5 + C*y6);
            final double h4 = 0.2  - x3 + 2.0*D*x5*x5 + x5*x6*(C*y1 - D*y2) + x5*x7*(C*y3 - D*y4);
            final double h5 = 0.2  - x4 + 2.0*D*x6*x6 - x5*x6*(C*y1 + D*y2) - x6*x7*(C*y5 + D*y6);
            final double h6 = -0.337      + 2.0*D*x7*x7 - x5*x7*(C*y3 + D*y4) + x6*x7*(C*y5 - D*y6);

            return new ArrayRealVector(new double[]{h1,h2,h3,h4,h5,h6}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            final double x1=x.getEntry(0), x2=x.getEntry(1), x3=x.getEntry(2), x4=x.getEntry(3),
                         x5=x.getEntry(4), x6=x.getEntry(5), x7=x.getEntry(6),
                         x8=x.getEntry(7), x9=x.getEntry(8);

            final double y1 = sin(x8);
            final double y2 = cos(x8);
            final double y3 = sin(x9);
            final double y4 = cos(x9);
            final double y5 = sin(x8 - x9);
            final double y6 = cos(x8 - x9);

            final double[][] J = new double[6][9];

            // h1 = 0.4 - x1 + 2C x5^2 - x5 x6 (D y1 + C y2) - x5 x7 (D y3 + C y4)
            J[0][0] = -1.0;
            J[0][4] = 4.0*C*x5 - x6*(D*y1 + C*y2) - x7*(D*y3 + C*y4);
            J[0][5] = -x5*(D*y1 + C*y2);
            J[0][6] = -x5*(D*y3 + C*y4);
            J[0][7] = -x5*x6*(D*y2 - C*y1);              // d/dx8 of -(D y1 + C y2)
            J[0][8] = -x5*x7*(D*y4 - C*y3);              // d/dx9 of -(D y3 + C y4)

            // h2 = 0.4 - x2 + 2C x6^2 + x5 x6 (D y1 - C y2) + x6 x7 (D y5 - C y6)
            final double v3  =  D*y6 + C*y5;              // d/dx8(y5)= y6 ; d/dx8(y6)= -y5
            final double v3m = -(D*y6 + C*y5);            // for d/dx9 of (D y5 - C y6)
            J[1][1] = -1.0;
            J[1][4] =  x6 * (D*y1 - C*y2);
            J[1][5] =  4.0*C*x6 + x5*(D*y1 - C*y2) + x7*(D*y5 - C*y6);
            J[1][6] =  x6 * (D*y5 - C*y6);
            J[1][7] =  x5*x6*(D*y2 + C*y1) + x6*x7 * v3;
            J[1][8] =  x6*x7 * v3m;

            // h3 = 0.8 + 2C x7^2 + x5 x7 (D y3 - C y4) - x6 x7 (D y5 + C y6)
            final double v7p =  D*y6 - C*y5;              // d/dx8 of (D y5 + C y6)
            final double v7m = -(D*y6 - C*y5);            // d/dx9 of (D y5 + C y6)
            J[2][4] =  x7 * (D*y3 - C*y4);
            J[2][5] = -x7 * (D*y5 + C*y6);
            J[2][6] =  4.0*C*x7 + x5*(D*y3 - C*y4) - x6*(D*y5 + C*y6);
            J[2][7] = -x6*x7 * v7p;
            J[2][8] =  x5*x7*(D*y4 + C*y3) - x6*x7 * v7m;

            // h4 = 0.2 - x3 + 2D x5^2 + x5 x6 (C y1 - D y2) + x5 x7 (C y3 - D y4)
            J[3][2] = -1.0;
            J[3][4] =  4.0*D*x5 + x6*(C*y1 - D*y2) + x7*(C*y3 - D*y4);
            J[3][5] =  x5 * (C*y1 - D*y2);
            J[3][6] =  x5 * (C*y3 - D*y4);
            J[3][7] =  x5*x6*(C*y2 + D*y1);
            J[3][8] =  x5*x7*(C*y4 + D*y3);

            // h5 = 0.2 - x4 + 2D x6^2 - x5 x6 (C y1 + D y2) - x6 x7 (C y5 + D y6)
            final double v12 = (C*y6 - D*y5) * x6;        // d/dx8 of -(x6*x7(Cy5+Dy6))
            J[4][3] = -1.0;
            J[4][4] = -x6 * (C*y1 + D*y2);
            J[4][5] =  4.0*D*x6 - x5*(C*y1 + D*y2) - x7*(C*y5 + D*y6);
            J[4][6] = -x6 * (C*y5 + D*y6);
            J[4][7] = -x5*x6*(C*y2 - D*y1) - x7 * v12;
            J[4][8] =  x7 * v12;

            // h6 = -0.337 + 2D x7^2 - x5 x7 (C y3 + D y4) + x6 x7 (C y5 - D y6)
            final double v15 = (C*y6 + D*y5) * x6 * x7;   // d/dx8 of +x6x7(Cy5-Dy6)
            J[5][4] = -x7 * (C*y3 + D*y4);
            J[5][5] =  x7 * (C*y5 - D*y6);
            J[5][6] =  4.0*D*x7 - x5*(C*y3 + D*y4) + x6*(C*y5 - D*y6);
            J[5][7] =  v15;
            J[5][8] = -x5*x7*(C*y4 - D*y3) - v15;

            return new Array2DRowRealMatrix(J, false);
        }
    }

    // -------- solve utility ----------
    private static LagrangeSolution solve() {
        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);

        // Start from the Fortran initial point (angles = 0)
        final double[] x0 = {
            0.8, 0.8, 0.2, 0.2, 1.0454, 1.0454, 1.0454, 0.0, 0.0
        };

        // Use large finite bounds instead of +/- infinity
        final double SUP = Double.POSITIVE_INFINITY; 
        final double INF = Double.NEGATIVE_INFINITY;
        final double[] lo = { 0.0, 0.0, INF, INF, 0.90909, 0.90909, 0.90909, INF, INF };
        final double[] up = { SUP,  SUP, SUP,  SUP, 1.0909,  1.0909,  1.0909,  SUP,  SUP };

        return optimizer.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS107Objective()),
                new HS107Eq(),
                new SimpleBounds(lo, up)
        );
    }

    // -------- minimal test: compare only objective (nonconvex tolerance) ----------
    @Test
    public void testHS107() {
        final double fEx = 0.505501180339e4;
        final LagrangeSolution sol = solve();
        final double f = sol.getValue();
        assertEquals(fEx, f, 1.0e-6 * (Math.abs(fEx) + 1.0), "objective mismatch");
    }
}
