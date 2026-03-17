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

import org.hipparchus.linear.*;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;

import static java.lang.Math.*;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS109Test {

    // Constants from the model
    private static final double A  = 50.176;
    private static final double RA = 1.0 / A;
    private static final double B  = sin(0.25);
    private static final double C  = cos(0.25);

    private static final double INF = Double.NEGATIVE_INFINITY;
    private static final double SUP = Double.POSITIVE_INFINITY;

    // Reference solution (XEX) and best-known objective (FEX)
    private static final double[] X_REF = {
            0.674888100445e3,  // x1
            0.113417039470e4,  // x2
            0.133569060261e0,  // x3
            -0.371152592466e0, // x4
            0.252e3,           // x5
            0.252e3,           // x6
            0.201464535316e3,  // x7
            0.426660777226e3,  // x8
            0.368494083867e3   // x9
    };
    private static final double F_REF = 0.536206927538e4;

    // ---------- Objective (value + gradient) ----------
    private static final class HS109Objective extends TwiceDifferentiableFunction {
        @Override public int dim() { return 9; }

        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return 3.0 * x1 + 1.0e-6 * x1*x1*x1
                 + 2.0 * x2 + 0.522074e-6 * x2*x2*x2;
        }

        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double[] g = new double[9];
            g[0] = 3.0 + 3.0e-6 * x1*x1;
            g[1] = 2.0 + 1.566222e-6 * x2*x2; // 3 * 0.522074e-6
            // others are zero
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    // ---------- Linear inequalities (2): g1,g2 >= 0 ----------
    // g1 = x4 - x3 + 0.55
    // g2 = x3 - x4 + 0.55
    private static final class HS109LinIneq extends InequalityConstraint {
        HS109LinIneq() { super(new ArrayRealVector(new double[2])); }
        @Override public int dim() { return 9; }

        @Override public RealVector value(RealVector x) {
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            return new ArrayRealVector(new double[]{
                    x4 - x3 + 0.55,
                    x3 - x4 + 0.55
            }, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double[][] J = new double[2][9];
            // g1 gradient: d/dx3=-1, d/dx4=+1
            J[0][2] = -1.0;
            J[0][3] = +1.0;
            // g2 gradient: d/dx3=+1, d/dx4=-1
            J[1][2] = +1.0;
            J[1][3] = -1.0;
            return new Array2DRowRealMatrix(J, false);
        }
    }

    // ---------- Nonlinear inequalities (2): g3,g4 >= 0 ----------
    // g3 = 2.25e6 - x1^2 - x8^2
    // g4 = 2.25e6 - x2^2 - x9^2
    private static final class HS109NonlinIneq extends InequalityConstraint {
        HS109NonlinIneq() { super(new ArrayRealVector(new double[2])); }
        @Override public int dim() { return 9; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x8 = x.getEntry(7);
            double x9 = x.getEntry(8);
            double g3 = 2.25e6 - x1*x1 - x8*x8;
            double g4 = 2.25e6 - x2*x2 - x9*x9;
            return new ArrayRealVector(new double[]{g3, g4}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x8 = x.getEntry(7);
            double x9 = x.getEntry(8);
            double[][] J = new double[2][9];
            // g3: d/dx1 = -2*x1, d/dx8 = -2*x8
            J[0][0] = -2.0 * x1;
            J[0][7] = -2.0 * x8;
            // g4: d/dx2 = -2*x2, d/dx9 = -2*x9
            J[1][1] = -2.0 * x2;
            J[1][8] = -2.0 * x9;
            return new Array2DRowRealMatrix(J, false);
        }
    }

    // ---------- Nonlinear equalities (6): g5..g10 = 0 ----------
    private static final class HS109Eq extends EqualityConstraint {
        HS109Eq() { super(new ArrayRealVector(new double[6])); }
        @Override public int dim() { return 9; }

        @Override public RealVector value(RealVector x) {
            double x1=x.getEntry(0), x2=x.getEntry(1), x3=x.getEntry(2),
                   x4=x.getEntry(3), x5=x.getEntry(4), x6=x.getEntry(5),
                   x7=x.getEntry(6), x8=x.getEntry(7), x9=x.getEntry(8);

            double y1 = sin(x8),           y2 = cos(x8);
            double y3 = sin(x9),           y4 = cos(x9);
            double y5 = sin(x8 - x9),      y6 = cos(x8 - x9);

            // Match Fortran formulas exactly (note the angles’ signs)
            double g5 = (x5*x6*sin(-x3 - 0.25)
                       + x5*x7*sin(-x4 - 0.25)
                       + 2.0*x5*x5*B)*RA + 400.0 - x1;

            double g6 = (x5*x6*sin( x3 - 0.25)
                       + x6*x7*sin( x3 - x4 - 0.25)
                       + 2.0*x6*x6*B)*RA + 400.0 - x2;

            double g7 = (x5*x7*sin( x4 - 0.25)
                       + x6*x7*sin( x4 - x3 - 0.25)
                       + 2.0*x7*x7*B)*RA + 881.779;

            double g8 = x8 + (x5*x6*cos(-x3 - 0.25)
                            + x5*x7*cos(-x4 - 0.25)
                            - 2.0*x5*x5*C)*RA + 0.7533e-3*x5*x5 - 200.0;

            double g9 = x9 + (x5*x6*cos( x3 - 0.25)
                            + x7*x6*cos( x3 - x4 - 0.25)
                            - 2.0*x6*x6*C)*RA + 0.7533e-3*x6*x6 - 200.0;

            double g10 = (x5*x7*cos( x4 - 0.25)
                        + x6*x7*cos( x4 - x3 - 0.25)
                        - 2.0*x7*x7*C)*RA + 0.7533e-3*x7*x7 - 22.938;

            return new ArrayRealVector(new double[]{g5,g6,g7,g8,g9,g10}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1=x.getEntry(0), x2=x.getEntry(1), x3=x.getEntry(2),
                   x4=x.getEntry(3), x5=x.getEntry(4), x6=x.getEntry(5),
                   x7=x.getEntry(6), x8=x.getEntry(7), x9=x.getEntry(8);

            double[][] J = new double[6][9];

            // g5 derivatives (Fortran lines “IF (.NOT.INDEX2(5)) ...”)
            {
                double V1 = sin(-x3 - 0.25);
                double V2 = sin(-x4 - 0.25);
                double V3 = x5 * RA;
                // d/dx1 = -1
                J[0][0] = -1.0;
                // d/dx3
                J[0][2] = -x6 * V3 * cos(-x3 - 0.25);
                // d/dx4
                J[0][3] = -x7 * V3 * cos(-x4 - 0.25);
                // d/dx5
                J[0][4] = (x6 * V1 + x7 * V2 + 4.0 * x5 * B) * RA;
                // d/dx6
                J[0][5] = V3 * V1;
                // d/dx7
                J[0][6] = V3 * V2;
            }

            // g6 derivatives
            {
                double HV1 = x3 - x4 - 0.25;
                double V3  = cos(HV1);
                double V4  = sin(x3 - 0.25);
                double V5  = x6 * RA;
                double V6  = sin(HV1);
                // d/dx2 = -1
                J[1][1] = -1.0;
                // d/dx3
                J[1][2] = x5 * V5 * cos(x3 - 0.25) + x7 * V5 * V3;
                // d/dx4
                J[1][3] = -x7 * V5 * V3;
                // d/dx5
                J[1][4] = V5 * V4;
                // d/dx6
                J[1][5] = (x5 * V4 + x7 * V6) * RA + 4.0 * V5 * B;
                // d/dx7
                J[1][6] = V5 * V6;
            }

            // g7 derivatives
            {
                double HV1 = x4 - x3 - 0.25;
                double V7  = x7 * RA;
                double V8  = cos(HV1);
                double V9  = sin(x4 - 0.25);
                double V10 = sin(HV1);
                // d/dx3
                J[2][2] = -x6 * V7 * V8;
                // d/dx4
                J[2][3] = x5 * V7 * cos(x4 - 0.25) + x6 * V7 * V8;
                // d/dx5
                J[2][4] = V7 * V9;
                // d/dx6
                J[2][5] = V7 * V10;
                // d/dx7
                J[2][6] = (x5 * V9 + x6 * V10) * RA + 4.0 * V7 * B;
            }

            // g8 derivatives
            {
                double V11 = x5 * RA;
                double V12 = cos(-x3 - 0.25) * RA;
                double V13 = cos(-x4 - 0.25) * RA;
                // d/dx8
                J[3][7] = 1.0;
                // d/dx3
                J[3][2] = x6 * V11 * sin(-x3 - 0.25);
                // d/dx4
                J[3][3] = x7 * V11 * sin(-x4 - 0.25);
                // d/dx5
                J[3][4] = x6 * V12 + x7 * V13 - 4.0 * V11 * C + 1.5066e-3 * x5;
                // d/dx6
                J[3][5] = x5 * V12;
                // d/dx7
                J[3][6] = x5 * V13;
            }

            // g9 derivatives
            {
                double HV1 = x3 - x4 - 0.25;
                double V14 = sin(HV1) * x6 * RA;
                double V15 = cos(x3 - 0.25) * RA;
                double V16 = cos(HV1) * RA;
                // d/dx9
                J[4][8] = 1.0;
                // d/dx3
                J[4][2] = -x5 * x6 * sin(x3 - 0.25) * RA - x7 * V14;
                // d/dx4
                J[4][3] = x7 * V14;
                // d/dx5
                J[4][4] = x6 * V15;
                // d/dx6
                J[4][5] = x5 * V15 + x7 * V16 - 4.0 * x6 * C * RA + 1.5066e-3 * x6;
                // d/dx7
                J[4][6] = x6 * V16;
            }

            // g10 derivatives
            {
                double HV1 = x4 - x3 - 0.25;
                double V17 = sin(HV1) * x6 * RA;
                double V18 = cos(x4 - 0.25) * RA;
                double V19 = cos(HV1) * RA;
                double V20 = x7 * RA;
                // d/dx3
                J[5][2] = x7 * V17;
                // d/dx4
                J[5][3] = -x5 * V20 * sin(x4 - 0.25) - x7 * V17;
                // d/dx5
                J[5][4] = x7 * V18;
                // d/dx6
                J[5][5] = x7 * V19;
                // d/dx7
                J[5][6] = x5 * V18 + x6 * V19 - 4.0 * V20 * C + 1.5066e-3 * x7;
            }

            return new Array2DRowRealMatrix(J, false);
        }
    }

    // ---------- Solve helper ----------
    private static LagrangeSolution solve() {
        // Bounds:
        // x1 >= 0, x2 >= 0
        // x3,x4 in [-0.55, +0.55]
        // x5,x6 in [196, 252]
        // x7 in [196, 252]
        // x8,x9 in [-400, 800]
        double[] lo = {0.0, 0.0, -0.55, -0.55, 196.0, 196.0, 196.0, -400.0, -400.0};
        double[] up = {SUP, SUP,  +0.55,  +0.55, 252.0, 252.0, 252.0,  800.0,  800.0};

        // Initial guess (as in the Fortran initialization)
        double[] x0 = {0, 0, 0, 0, 250, 250, 200, 0, 0};

        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        return optimizer.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS109Objective()),
                new HS109LinIneq(),
                new HS109NonlinIneq(),
                new HS109Eq(),
                new SimpleBounds(lo, up)
        );
    }

    // ---------- Minimal test: compare objective near reference ----------
    @Test
    public void testHS109() {
        LagrangeSolution sol = solve();
        double f = sol.getValue();

        // Objective comparison (looser tol due to nonconvexity & mixed eq/ineq)
        assertEquals(F_REF, f, 1e-6 * (F_REF + 1.0), "objective mismatch");
    }
}
