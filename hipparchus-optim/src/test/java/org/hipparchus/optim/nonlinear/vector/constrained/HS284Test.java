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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS284 — 15 variables, 10 nonlinear inequality constraints.
 * Minimize f(x) = - sum(C_i * x_i)  subject to  g_i(x) = B_i - sum_j A_{ij} x_j^2 >= 0.
 */
public class HS284Test {

    /** Problem data: A is 10x15 exactly as in the Fortran DATA (rows = constraints i, cols = vars j). */
    static final class Data {
        // C(15)
        final double[] C = {20, 40, 400, 20, 80, 20, 40, 140, 380, 280, 80, 40, 140, 40, 120};

        // B(10)
        final double[] B = {385, 470, 560, 565, 645, 430, 485, 455, 390, 460};

        // A(10,15) — rows are constraints 1..10, columns are variables 1..15.
        // (Attenzione: numeri espansi senza scorciatoie Fortran; nessuna trasposizione.)
        final double[][] A = {
    {100.0, 100.0,  10.0,   5.0,  10.0,   0.0,   0.0,  25.0,   0.0,  10.0,  55.0,   5.0,  45.0,  20.0,  0.0},
    { 90.0, 100.0,  10.0,  35.0,  20.0,   5.0,   0.0,  35.0,  55.0,  25.0,  20.0,   0.0,  40.0,  25.0, 10.0},
    { 70.0,  50.0,   0.0,  55.0,  25.0, 100.0,  40.0,  50.0,   0.0,  30.0,  60.0,  10.0,  30.0,   0.0, 40.0},
    { 50.0,   0.0,   0.0,  65.0,  35.0, 100.0,  35.0,  60.0,   0.0,  15.0,   0.0,  75.0,  35.0,  30.0, 65.0},
    { 50.0,  10.0,  70.0,  60.0,  45.0,  45.0,   0.0,  35.0,  65.0,   5.0,  75.0, 100.0,  75.0,  10.0,  0.0},
    { 40.0,   0.0,  50.0,  95.0,  50.0,  35.0,  10.0,  60.0,   0.0,  45.0,  15.0,  20.0,   0.0,   5.0,  5.0},
    { 30.0,  60.0,  30.0,  90.0,   0.0,  30.0,   5.0,  25.0,   0.0,  70.0,  20.0,  25.0,  70.0,  15.0, 15.0},
    { 20.0,  30.0,  40.0,  25.0,  40.0,  25.0,  15.0,  10.0,  80.0,  20.0,  30.0,  30.0,   5.0,  65.0, 20.0},
    { 10.0,  70.0,  10.0,  35.0,  25.0,  65.0,   0.0,  30.0,   0.0,   0.0,  25.0,   0.0,  15.0,  50.0, 55.0},
    {  5.0,  10.0, 100.0,   5.0,  20.0,   5.0,  10.0,  35.0,  95.0,  70.0,  20.0,  10.0,  35.0,  10.0, 30.0}
};
    }

    /** Objective f(x) = - sum(C_i * x_i). */
    static final class HS284Objective extends TwiceDifferentiableFunction {
        private final double[] C;
        HS284Objective(double[] C) { this.C = C; }
        @Override public int dim() { return 15; }
        @Override public double value(RealVector x) {
            double s = 0.0;
            for (int i = 0; i < 15; i++) s += C[i] * x.getEntry(i);
            return -s;
        }
        @Override public RealVector gradient(RealVector x) {
            double[] g = new double[15];
            for (int i = 0; i < 15; i++) g[i] = -C[i];
            return new ArrayRealVector(g, false);
        }
        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    /** Ten nonlinear inequalities g(x) >= 0 with Jacobian (no transpose). */
    static final class HS284Ineq extends InequalityConstraint {
        private final double[][] A; // A[i][j]
        private final double[] B;
        HS284Ineq(double[][] A, double[] B) {
            super(new ArrayRealVector(new double[10])); // 10 constraints
            this.A = A;
            this.B = B;
        }
        @Override public int dim() { return 15; }

        @Override public RealVector value(RealVector x) {
            double[] g = new double[10];
            for (int i = 0; i < 10; i++) {
                double s = 0.0;
                for (int j = 0; j < 15; j++) {
                    final double v = x.getEntry(j);
                    s += A[i][j] * v * v;
                }
                g[i] = B[i] - s; // >= 0
            }
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double[][] J = new double[10][15];
            for (int i = 0; i < 10; i++) {
                for (int j = 0; j < 15; j++) {
                    J[i][j] = -2.0 * A[i][j] * x.getEntry(j);
                }
            }
            return new Array2DRowRealMatrix(J, false);
        }
    }

    private static LagrangeSolution solve(double[] x0) {
        Data d = new Data();
        SQPOptimizerS2 opt = new SQPOptimizerS2();
        opt.setDebugPrinter(System.out::println); // keep debug on
        SQPOption sqpOption=new SQPOption();

        return opt.optimize(
                sqpOption,
            new InitialGuess(x0),
            new ObjectiveFunction(new HS284Objective(d.C)),
            new HS284Ineq(d.A, d.B) // no bounds at all
        );
    }

    @Test
    public void testHS284() {
        Data d = new Data();

        // Start at zero (feasible: g = B >= 0)
        double[] x0 = new double[15];

        LagrangeSolution sol = solve(x0);

        // Expected: FEX = -sum(C)
        double fEx = 0.0;
        for (double v : d.C) fEx -= v;

        double f = sol.getValue();
        assertEquals(fEx, f, 1.0e-6 * (Math.abs(fEx) + 1.0), "objective mismatch");
    }
}
