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
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;


public class HS101toHS103Test {

    // -------- Problem selector (0->HS101, 1->HS102, 2->HS103) --------
    static final class Problem {
        final int m;         // 0,1,2
        final double a;      // {-0.25, 0.125, 0.5}
        final double[] xEx;  // reference x*
        final double fEx;    // reference f(x*)

        Problem(int m) {
            this.m = m;
            this.a = new double[]{-0.25, 0.125, 0.5}[m];
            if (m == 0) {
                xEx = new double[]{2.85615855584, 0.610823030755, 2.15081256203, 4.71287370945, 0.999487540961, 1.34750750498, 0.0316527664991};
                fEx = 0.180976476556e4;
            } else if (m == 1) {
                xEx = new double[]{3.89625319099, 0.809358760118, 2.66438599373, 4.30091287458, 0.853554935267, 1.09528744459, 0.0273104596581};
                fEx = 0.911880571336e3;
            } else {
                xEx = new double[]{4.39410451026, 0.854468738817, 2.84323031380, 3.39997866779, 0.722926133025, 0.870406381840, 0.0246388263302};
                fEx = 0.543667958424e3;
            }
        }
    }

    // -------- Helpers --------
    private static double powAbs(double x, double p) { return Math.pow(Math.abs(x), p); }

    // -------- Objective f, grad --------
    static final class HSObjective extends TwiceDifferentiableFunction {
        private final Problem prob;
        HSObjective(Problem p) { this.prob = p; }
        @Override public int dim() { return 7; }

        @Override public double value(RealVector x) {
            for (int k = 0; k < 7; k++) if (x.getEntry(k) < 1e-8) {
                double sum = 0;
                for (int i = 0; i < 7; i++) {
                    double d = x.getEntry(i) - 5.0;
                    sum += d * d;
                }
                final double[] fmin = {1.8e3, 9.1e2, 5.4e2};
                return sum + 1.0e3 + fmin[prob.m];
            }
            double x1=x.getEntry(0), x2=x.getEntry(1), x3=x.getEntry(2),
                   x4=x.getEntry(3), x5=x.getEntry(4), x6=x.getEntry(5),
                   x7=x.getEntry(6);

            double term1 = 10.0 * x1 * x4*x4 * Math.pow(x7, prob.a) / (x2 * Math.pow(x6, 3));
            double term2 = 15.0 * x3 * x4 / (x1 * x2*x2 * x5 * Math.pow(x7, 0.5));
            double term3 = 20.0 * x2 * x6 / (x1*x1 * x4 * x5*x5);
            double term4 = 25.0 * x1*x1 * x2*x2 * Math.pow(x5, 0.5) * x7 / (x3 * x6*x6);
            return term1 + term2 + term3 + term4;
        }

        @Override public RealVector gradient(RealVector x) {
            for (int k = 0; k < 7; k++) if (x.getEntry(k) < 1e-8) {
                double[] g = new double[7];
                for (int i = 0; i < 7; i++) g[i] = 2.0 * (x.getEntry(i) - 5.0);
                return new ArrayRealVector(g, false);
            }
            final double a = prob.a;
            double x1=x.getEntry(0), x2=x.getEntry(1), x3=x.getEntry(2),
                   x4=x.getEntry(3), x5=x.getEntry(4), x6=x.getEntry(5),
                   x7=x.getEntry(6);

            double V1  = 10.0 * x4*x4;
            double V2  = Math.pow(x7, a);
            double V3  = x2 * Math.pow(x6, 3);
            double V4  = 15.0 * x3 * x4;
            double V5  = x1 * x2*x2 * x5 * Math.pow(x7, 0.5);
            double V6  = 20.0 * x2 * x6;
            double V7  = x1*x1 * x4 * x5*x5;          // x5^2
            double V8  = 25.0 * x1 * x2 * Math.pow(x5, 0.5) * x7;
            double V9  = x3 * x6*x6;
            double V10 = 12.5 * x1*x1 * x2*x2 * x7;
            double V11 = Math.pow(x5, 0.5);

            double[] g = new double[7];
            g[0] =  V1*V2/V3 - V4/(x1*V5) - 2.0*V6/(x1*V7) + 2.0*x2*V8/V9;
            g[1] = -V1*x1*V2/(x2*V3) - 2.0*V4/(x2*V5) + 20.0*x6/V7 + 2.0*x1*V8/V9;
            g[2] =  15.0*x4/V5 - x1*x2*V8/(x3*V9);
            g[3] =  20.0*x1*x4*V2/V3 + 15.0*x3/V5 - V6/(x4*V7);
            g[4] = -V4/(x5*V5) - 2.0*V6/(x5*V7) + V10/(V9*V11);
            g[5] = -3.0*V1*x1*V2/(x6*V3) + 20.0*x2/V7 - 4.0*V10*V11/(x6*V9);
            g[6] =  a*V1*x1*Math.pow(x7, a-1.0)/V3 - 0.5*V4/(V5*x7) + V8*x1*x2/(x7*V9);
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    // -------- Six nonlinear inequalities g(x) >= 0 (with Jacobian) --------
    static final class HSIneq extends InequalityConstraint {
        private final Problem prob;
        HSIneq(Problem p) { super(new ArrayRealVector(new double[6])); this.prob = p; }
        @Override public int dim() { return 7; }

        @Override public RealVector value(RealVector x) {
            for (int k = 0; k < 7; k++) if (x.getEntry(k) < 1e-8) {
                return new ArrayRealVector(new double[6]); // zeroes
            }
            double x1=x.getEntry(0), x2=x.getEntry(1), x3=x.getEntry(2),
                   x4=x.getEntry(3), x5=x.getEntry(4), x6=x.getEntry(5),
                   x7=x.getEntry(6);
            final double a = prob.a;

            double g1 = 1.0
                    - 0.5*powAbs(x1,0.5)*x7/(x3*Math.pow(x6,2))
                    - 0.7*Math.pow(x1,3)*x2*x6*powAbs(x7,0.5)/(x3*x3)
                    - 0.2*x3*powAbs(x6,2.0/3.0)*powAbs(x7,0.25)/(x2*powAbs(x4,0.5));

            double g2 = 1.0
                    - 1.3*x2*x6/(powAbs(x1,0.5)*x3*x5)
                    - 0.8*x3*x6*x6/(x4*x5)
                    - 3.1*powAbs(x2,0.5)*powAbs(x6,1.0/3.0)/(x1*x4*x4*x5);

            double g3 = 1.0
                    - 2.0*x1*x5*powAbs(x7,1.0/3.0)/(powAbs(x3,1.5)*x6)
                    - 0.1*x2*x5/(powAbs(x3*x7,0.5)*x6)
                    - x2*powAbs(x3,0.5)*x5/x1
                    - 0.65*x3*x5*x7/(x2*x2*x6);

            double g4 = 1.0
                    - 0.2*x2*powAbs(x5,0.5)*powAbs(x7,1.0/3.0)/(x1*x1*x4)
                    - 0.3*powAbs(x1,0.5)*x2*x2*x3*powAbs(x4,1.0/3.0)*powAbs(x7,0.25)/powAbs(x5,2.0/3.0)
                    - 0.4*x3*x5*powAbs(x7,0.75)/(x1*x1*x1*x2*x2)
                    - 0.5*x4*powAbs(x7,0.5)/(x3*x3);

            double term1 = 10.0 * x1 * x4*x4 * Math.pow(Math.abs(x7), a) / (x2 * Math.pow(x6, 3));
            double term2 = 15.0 * x3 * x4 / (x1 * x2*x2 * x5 * Math.pow(Math.abs(x7), 0.5));
            double term3 = 20.0 * x2 * x6 / (x1*x1 * x4 * x5*x5);
            double term4 = 25.0 * x1*x1 * x2*x2 * Math.pow(Math.abs(x5), 0.5) * x7 / (x3 * x6*x6);

            double g5 = term1 + term2 + term3 + term4 - 100.0;
            double g6 = -(term1 + term2 + term3 + term4) + 3.0e3;

            return new ArrayRealVector(new double[]{g1,g2,g3,g4,g5,g6}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            for (int k = 0; k < 7; k++) if (x.getEntry(k) < 1e-8) {
                return new Array2DRowRealMatrix(6,7);
            }
            double x1=x.getEntry(0), x2=x.getEntry(1), x3=x.getEntry(2),
                   x4=x.getEntry(3), x5=x.getEntry(4), x6=x.getEntry(5),
                   x7=x.getEntry(6);
            final double a = prob.a;

            double[][] J = new double[6][7];

            // g1
            double V1 = powAbs(x1,0.5);
            double V2 = x1*x1*x1;
            double V4 = x3*x3;
            double V6 = powAbs(x4,0.5);
            double V7 = x6*x6;
            double V8 = powAbs(x6,2.0/3.0);
            double V9 = powAbs(x7,0.5);
            double V10 = powAbs(x7,0.25);
            J[0][0] = -0.25*x7/(V1*x3*V7) - 2.1*x1*x1*x2*x6*V9/V4;
            J[0][1] = -0.7*V2*x6*V9/V4 + 0.2*x3*V8*V10/(x2*x2*V6);
            J[0][2] =  0.5*V1*x7/(V4*V7) + 1.4*V2*x2*x6*V9/(x3*V4) - 0.2*V8*V10/(x2*V6);
            J[0][3] =  0.1*x3*V8*V10/(x2*x4*V6);
            J[0][5] =  V1*x7/(x3*V7*x6) - 0.7*V2*x2*V9/V4 - (0.4/3.0)*x3*V10/(x2*V6*Math.pow(Math.abs(x6),1.0/3.0));
            J[0][6] = -0.5*V1/(x3*V7) - 0.35*V2*x2*x6/(V4*V9) - 0.05*x3*V8/(x2*V6*V9*V10);

            // g2
            double V11 = powAbs(x1,0.5);
            double V12 = powAbs(x2,0.5);
            double V13 = x4*x4;
            double V14 = x5*x5;
            double V15 = powAbs(x6,1.0/3.0);
            double V16 = x6*x6;
            J[1][0] =  0.65*x2*x6/(x1*V11*x3*x5) + 3.1*V12*V15/(x1*x1*V13*x5);
            J[1][1] = -1.3*x6/(V11*x3*x5) - 1.55*V15/(x1*V12*V13*x5);
            J[1][2] =  1.3*x2*x6/(V11*x3*x3*x5) - 0.8*V16/(x4*x5);
            J[1][3] =  0.8*x3*V16/(V13*x5) + 6.2*V12*V15/(x1*V13*x4*x5);
            J[1][4] =  1.3*x2*x6/(V11*x3*V14) + 0.8*x3*V16/(x4*V14) + 3.1*V12*V15/(x1*V13*V14);
            J[1][5] = -1.3*x2/(V11*x3*x5) - 1.6*x3*x6/(x4*x5) - (3.1/3.0)*V12/(x1*V13*x5*Math.pow(V15,2));

            // g3
            double V17 = x2*x2;
            double V18 = powAbs(x3,0.5);
            double V19 = V18*x3;
            double V20 = x6*x6;
            double V21 = powAbs(x7,1.0/3.0);
            double V22 = powAbs(x7,0.5);
            J[2][0] = -2.0*x5*V21/(V19*x6) + x2*V18*x5/(x1*x1);
            J[2][1] = -V18*x5/x1 + 1.3*x3*x5*x7/(V17*x2*x6) - 0.1*x5/(V18*V22*x6);
            J[2][2] =  3.0*x1*x5*V21/(x3*V19*x6) + 0.05*x2*x5/(V19*x6*V22) - 0.5*x2*x5/(x1*V18) - 0.65*x5*x7/(V17*x6);
            J[2][4] = -2.0*x1*V21/(V19*x6) - 0.1*x2/(V18*x6*V22) - x2*V18/x1 - 0.65*x3*x7/(V17*x6);
            J[2][5] =  2.0*x1*x5*V21/(V19*V20) + 0.1*x2*x5/(V18*V20*V22) + 0.65*x3*x5*x7/(V17*V20);
            J[2][6] = -2.0/3.0 * x1*x5/(V19*x6*Math.pow(V21,2)) + 0.05*x2*x5/(V18*x6*V22*x7) - 0.65*x3*x5/(V17*x6);

            // g4
            double V23 = powAbs(x1,0.5);
            double V24 = x1*x1;
            double V25 = V24*x1;
            double V26 = x2*x2;
            double V27 = x3*x3;
            double V28 = powAbs(x4,1.0/3.0);
            double V29 = powAbs(x5,2.0/3.0);
            double V30 = powAbs(x5,0.5);
            double V31 = powAbs(x7,0.25);
            double V32 = V31*V31;
            double V33 = V31*V32;
            double V34 = powAbs(x7,1.0/3.0);
            J[3][0] =  0.4*x2*V30*V34/(V25*x4) - 0.15*V26*x3*V28*V31/(V23*V29) + 1.2*x3*x5*V33/(V24*V24*V26);
            J[3][1] = -0.2*V30*V34/(V24*x4) - 0.6*V23*x2*x3*V28*V31/V29 + 0.8*x3*x5*V33/(V25*V26*x2);
            J[3][2] = -0.3*V23*V26*V28*V31/V29 - 0.4*x5*V33/(V25*V26) + x4*V32/(V27*x3);
            J[3][3] =  0.2*x2*V30*V34/(V24*x4*x4) - 0.1*V23*V26*x3*V31/(Math.pow(V28,2)*V29) - 0.5*V32/V27;
            J[3][4] = -0.1*x2*V34/(V24*x4*V30) + 0.2*V23*V26*x3*V28*V31/(x5*V29) - 0.4*x3*V33/(V25*V26);
            J[3][6] = -(0.2/3.0)*x2*V30/(V24*x4*Math.pow(V34,2)) - 0.075*V23*V26*x3*V28/(V29*V33) - 0.3*x3*x5/(V25*V26*V31) - 0.25*x4/(V27*V32);

            // g5 / g6 gradients (±)
            double V35 = 10.0 * x4*x4;
            double V36 = Math.pow(Math.abs(x7), a);
            double V37 = x2 * Math.pow(x6, 3);
            double V38 = 15.0 * x3 * x4;
            double V39 = x1 * x2*x2 * x5 * Math.pow(Math.abs(x7), 0.5);
            double V40 = 20.0 * x2 * x6;
            double V41 = x1*x1 * x4 * x5*x5;
            double V42 = 25.0 * x1 * x2 * Math.pow(Math.abs(x5), 0.5) * x7;
            double V43 = x3 * x6*x6;
            double V44 = 12.5 * x1*x1 * x2*x2 * x7;
            double V45 = Math.pow(Math.abs(x5), 0.5);

            double[] GV = new double[7];
            GV[0] =  V35*V36/V37 - V38/(x1*V39) - 2.0*V40/(x1*V41) + 2.0*x2*V42/V43;
            GV[1] = -V35*x1*V36/(x2*V37) - 2.0*V38/(x2*V39) + 20.0*x6/V41 + 2.0*x1*V42/V43;
            GV[2] =  15.0*x4/V39 - x1*x2*V42/(x3*V43);
            GV[3] =  20.0*x1*x4*V36/V37 + 15.0*x3/V39 - V40/(x4*V41);
            GV[4] = -V38/(x5*V39) - 2.0*V40/(x5*V41) + V44/(V43*V45);
            GV[5] = -3.0*V35*x1*V36/(x6*V37) + 20.0*x2/V41 - 4.0*V44*V45/(x6*V43);
            GV[6] =  a*V35*x1*Math.pow(Math.abs(x7), a-1.0)/V37 - 0.5*V38/(V39*x7) + V42*x1*x2/(x7*V43);

            System.arraycopy(GV, 0, J[4], 0, 7);
            for (int i = 0; i < 7; i++) J[5][i] = -GV[i];

            return new Array2DRowRealMatrix(J, false);
        }
    }

    // -------- Solve utility --------
    static LagrangeSolution solve(Problem p, double[] start, double[] lo, double[] up) {
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        return optimizer.optimize(
                new InitialGuess(start),
                new ObjectiveFunction(new HSObjective(p)),
                new HSIneq(p),
                new SimpleBounds(lo, up)
        );
    }

    // -------- JUnit tests kept INSIDE THE SAME FILE --------
    
        private void runCase(int m) {
            Problem prob = new Problem(m);
            double[] x0 = {6,6,6,6,6,6,6};
            double[] lo = {0.1,0.1,0.1,0.1,0.1,0.1,0.01};
            double[] up = {10,10,10,10,10,10,10};

            LagrangeSolution sol = solve(prob, x0, lo, up);

           
            // Objective ~ reference (nonconvex tolerance)
            double f = sol.getValue();
            assertEquals(prob.fEx, f, 1.e-6*(prob.fEx+1.0), "objective mismatch");

           
        }
        @Test 
        public void testHS101() { runCase(0); }
        
        @Test 
        public void testHS102() { runCase(1); }
        
        @Test 
        public void testHS103() { runCase(2); }
    
}
