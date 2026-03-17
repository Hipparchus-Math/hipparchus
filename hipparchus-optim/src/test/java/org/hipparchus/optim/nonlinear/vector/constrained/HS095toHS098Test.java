/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
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

public class HS095toHS098Test {

    /** Objective: f(x) = c·x (linear). */
    private static class TP95_98_Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 6; }
        @Override public double value(RealVector x) {
            final double x1 = x.getEntry(0), x2 = x.getEntry(1), x3 = x.getEntry(2),
                         x4 = x.getEntry(3), x5 = x.getEntry(4), x6 = x.getEntry(5);
            return 4.3 * x1 + 31.8 * x2 + 63.3 * x3 + 15.8 * x4 + 68.5 * x5 + 4.7 * x6;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    /**
     * Inequalities g(x) ≥ 0. This class internally flips the original forms
     * so that the optimizer always sees ≥ 0 constraints.
     */
    private static class TP95_98_Ineq extends InequalityConstraint {
        private final double B1, B2, B3, B4;

        TP95_98_Ineq(double B1, double B2, double B3, double B4) {
            // 4 inequalities
            super(new ArrayRealVector(new double[] {0, 0, 0, 0}));
            this.B1 = B1; this.B2 = B2; this.B3 = B3; this.B4 = B4;
        }

        @Override public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0), x2 = x.getEntry(1), x3 = x.getEntry(2),
                         x4 = x.getEntry(3), x5 = x.getEntry(4), x6 = x.getEntry(5);

            // Original forms were written as <= 0. We return -G so they become ≥ 0.
            final double G1f =
                    17.1*x1 + 38.2*x2 + 204.2*x3 + 212.3*x4 + 623.4*x5 + 1495.5*x6
                  - 169.0*x1*x3 - 3580.0*x3*x5 - 3810.0*x4*x5
                  - 18500.0*x4*x6 - 24300.0*x5*x6 - B1;

            final double G2f =
                    17.9*x1 + 36.8*x2 + 113.9*x3 + 169.7*x4 + 337.8*x5 + 1385.2*x6
                  - 139.0*x1*x3 - 2450.0*x4*x5 - 16600.0*x4*x6 - 17200.0*x5*x6 - B2;

            final double G3f =
                    -273.0*x2 - 70.0*x4 - 819.0*x5 + 26000.0*x4*x5 - B3;

            final double G4f =
                    159.9*x1 - 311.0*x2 + 587.0*x4 + 391.0*x5 + 2198.0*x6
                  - 14000.0*x1*x6 - B4;

            return new ArrayRealVector(new double[] { G1f, G2f, G3f, G4f });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 6; }
    }

    private static SimpleBounds bounds() {
        // 0 ≤ x ≤ XU
        final double[] lo = {0, 0, 0, 0, 0, 0};
        final double[] up = {0.31, 0.046, 0.068, 0.042, 0.028, 0.0134};
        return new SimpleBounds(lo, up);
    }

    private static InitialGuess guess() {
        // The original setup started from zeros.
        return new InitialGuess(new double[] {0, 0, 0, 0, 0, 0});
    }

    @Test
    public void testHS095() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        // KN1 = 1 → B1..B4 as below
        final double B1 = 4.97,  B2 = -1.88,  B3 = -29.08,  B4 = -78.02;

        LagrangeSolution sol = optimizer.optimize(
                guess(),
                new ObjectiveFunction(new TP95_98_Obj()),
                new TP95_98_Ineq(B1, B2, B3, B4),
                bounds()
        );

        // Expected f*
        final double fExpected = 0.0156195144282;
        HSProblemTestUtils.assertExpectedObjective(fExpected, sol);
    }

    @Test
    public void testHS096() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        // KN1 = 2
        final double B1 = 4.97,  B2 = -1.88,  B3 = -69.08,  B4 = -118.02;

        LagrangeSolution sol = optimizer.optimize(
                guess(),
                new ObjectiveFunction(new TP95_98_Obj()),
                new TP95_98_Ineq(B1, B2, B3, B4),
                bounds()
        );

        final double fExpected = 0.0156195134384;
        HSProblemTestUtils.assertExpectedObjective(fExpected, sol);
    }

    @Test
    public void testHS097() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        // KN1 = 3
        final double B1 = 32.97,  B2 = 25.12,  B3 = -29.08,  B4 = -78.02;

        LagrangeSolution sol = optimizer.optimize(
                guess(),
                new ObjectiveFunction(new TP95_98_Obj()),
                new TP95_98_Ineq(B1, B2, B3, B4),
                bounds()
        );

        final double fExpected = 3.1358089;
        HSProblemTestUtils.assertExpectedObjective(fExpected, sol);
    }

    @Test
    public void testHS098() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        // KN1 = 4
        final double B1 = 32.97,  B2 = 25.12,  B3 = -124.08,  B4 = -173.02;

        LagrangeSolution sol = optimizer.optimize(
                guess(),
                new ObjectiveFunction(new TP95_98_Obj()),
                new TP95_98_Ineq(B1, B2, B3, B4),
                bounds()
        );

        final double fExpected = 3.1358089;
        HSProblemTestUtils.assertExpectedObjective(fExpected, sol);
    }

    // Utility: quick check of f(x) at a given point (can help during debugging)
    @SuppressWarnings("unused")
    private static double evalObjective(double[] x) {
        return 4.3*x[0] + 31.8*x[1] + 63.3*x[2] + 15.8*x[3] + 68.5*x[4] + 4.7*x[5];
    }

    // Optional: known candidate points (can be used locally if needed)
    @SuppressWarnings("unused")
    private static double[] hs97_refPoint() {
        return new double[] { 0.268564912352, 0.0, 0.0, 0.0, 0.028, 0.0134 };
    }

    @SuppressWarnings("unused")
    private static double[] hs98_refPoint() {
        return new double[] { 0.268564912323, 0.0, 0.0, 0.0, 0.028, 0.0134 };
    }
}
