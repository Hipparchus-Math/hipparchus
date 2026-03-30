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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/** TP88..TP92 combined: strict 1-based port of the nonlinear inequality; bounds separate. */
public class HS088toHS092 {

    private static final double ASSERT_TOL = 1e-5;

    /** Minimal 1-based array wrapper to mirror Fortran indexing. */
    static final class F1Array {
        private final double[] a; // index 0 unused
        F1Array(int n) { a = new double[n + 1]; }
        double get(int i) { return a[i]; }
        void set(int i, double v) { a[i] = v; }
        int len() { return a.length - 1; }
        static double powMinus1(int k) { return ((k & 1) == 0) ? 1.0 : -1.0; }
    }

    /** Exact 1:1 port of DOUBLE PRECISION FUNCTION GLEICH(P). */
    static double GLEICH(double p) {
        final double EPS = 1.0e-5;
        double y = p + 1.0, f, a;
        for (int it = 0; it < 10000; it++) { // safety cap
            f = y - p - FastMath.atan(1.0 / y);
            if (FastMath.abs(f) <= EPS) break;
            a = y * y + 1.0;
            a = (a + 1.0) / a;
            y = y - f / a;
        }
        final double EPS2 = EPS * EPS;
        return (y > EPS2) ? y : EPS2;
    }

    /** f(x) = sum x_i^2. */
    private static final class Obj extends TwiceDifferentiableFunction {
        private final int n;
        Obj(int n) { this.n = n; }
        @Override public int dim() { return n; }
        @Override public double value(RealVector x) {
            double s = 0.0;
            for (int i = 0; i < n; i++) s += x.getEntry(i) * x.getEntry(i);
            return s;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    /** g(x) = 1e-4 - V1(x) - 2/15  (feasible if g(x) >= 0). Everything 1-based like Fortran. */
    static final class HS88to92Ineq extends InequalityConstraint {
        private static final int M = 30;
        private static final double INTKO = 2.0 / 15.0;
        private static final double PI = FastMath.atan(1.0) * 4.0;

        private final int N;
        private final F1Array MUE = new F1Array(M);
        private final F1Array A   = new F1Array(M);
        private final F1Array DCOSKO = new F1Array(M);

        HS88to92Ineq(int N) {
            super(new ArrayRealVector(new double[]{0.0}));
            this.N = N;
            // MODE=1 precomputation: DO I=1,30
            for (int I = 1; I <= M; I++) {
                double Z = PI * (I - 1);
                double mu = GLEICH(Z);
                MUE.set(I, mu);
                double s = FastMath.sin(mu), c = FastMath.cos(mu);
                DCOSKO.set(I, (s / mu - c) / (mu * mu));
                A.set(I, 2.0 * s / (mu + s * c));
            }
        }

        @Override public RealVector value(RealVector x0) {
            // x1(1..N)
            F1Array x = new F1Array(N);
            for (int i = 1; i <= N; i++) x.set(i, x0.getEntry(i - 1));

            // T(N)=x(N)^2; T(N-I)=x(N-I)^2 + T(N-I+1)
            F1Array T = new F1Array(N);
            T.set(N, sq(x.get(N)));
            for (int I = 1; I <= N - 1; I++) {
                int idx = N - I;
                T.set(idx, sq(x.get(idx)) + T.get(idx + 1));
            }

            double V1 = 0.0;
            for (int J = 1; J <= M; J++) {
                double W  = MUE.get(J);
                double V3 = -W * W;

                double rho = F1Array.powMinus1(N);
                for (int I = 1; I <= N - 1; I++) {
                    int tIndex = N + 1 - I;
                    double a1 = V3 * T.get(tIndex);
                    double ep1 = (a1 > -100.0) ? FastMath.exp(a1) : 0.0;
                    rho += F1Array.powMinus1(N - I) * 2.0 * ep1;
                }
                {
                    double a1 = V3 * T.get(1);
                    double ep1 = (a1 > -100.0) ? FastMath.exp(a1) : 0.0;
                    rho = (rho + ep1) / V3;
                }

                double term = W * FastMath.sin(W) * rho - 2.0 * DCOSKO.get(J);
                V1 += -V3 * A.get(J) * rho * term;
            }

            double g = 1.0e-4 - V1 - INTKO; // Fortran G(1)
            return new ArrayRealVector(new double[]{ g });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return N; }

        private static double sq(double v) { return v * v; }
    }

    // ---------------------------- Tests (TP88..TP92) ----------------------------

    @Test 
    public void testHS088() { // N=2
        final int n = 2;
        final double[] xex = {
                .107431872940D+01,
                -0.456613707247D+00
        };
        final double fex = 0.136265680997e+01;

        runCaseAndCheck(n, fex, xex);
    }

    @Test 
    public void testHS089() { // N=3
        final int n = 3;
        final double[] xex = {
                0.107431872754e+01,
               -0.456613706239e+00,
                0.300836097604e-10
        };
        final double fex = 0.136265680508e+01;

        runCaseAndCheck(n, fex, xex);
    }

    @Test 
    public void testHS090() { // N=4
        final int n = 4;
        final double[] xex = {
                0.708479399007e+00,
                0.237919269592e-04,
                0.807599939006e+00,
               -0.456613723294e+00
        };
        final double fex = 0.136265681317e+01;

        runCaseAndCheck(n, fex, xex);
    }

    @Test 
    public void testHS091() { // N=5
        final int n = 5;
        final double[] xex = {
                0.701892928031e+00,
                0.221084326516e-11,
                0.813330836201e+00,
                0.456613707134e+00,
                0.899937588382e-11
        };
        final double fex = 0.136265680910e+01;

        runCaseAndCheck(n, fex, xex);
    }

    @Test 
    public void testHS092() { // N=6
        final int n = 6;
        final double[] xex = {
                0.494144465323e+00,
               -0.103530473697e-04,
                0.614950839550e+00,
               -0.242186612731e-05,
                0.729258528936e+00,
               -0.456613099133e+00
        };
        final double fex = 0.136265681213e+01;

        runCaseAndCheck(n, fex, xex);
    }

    // ---------------------------- Common runner ----------------------------

    private static void runCaseAndCheck(int n, double expectedFex, double[] xex) {
        // start: 0.5, -0.5, 0.5, -0.5, ...
        double[] x0 = new double[n];
        for (int i = 0; i < n; i++) x0[i] = (i % 2 == 0) ? 0.5 : -0.5;

        // bounds: x1 in [0.1, 10], others in [-10, 10]
        double[] lo = new double[n], hi = new double[n];
        for (int i = 0; i < n; i++) { lo[i] = -10.0; hi[i] = 10.0; }
        lo[0] = 0.1;

        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        

        // Run optimizer and verify objective quality against the reference value.
        LagrangeSolution sol = optimizer.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new Obj(n)),
                new HS88to92Ineq(n),
                new SimpleBounds(lo, hi)
        );

        HSProblemTestUtils.assertExpectedObjective(expectedFex, sol);
    }
}