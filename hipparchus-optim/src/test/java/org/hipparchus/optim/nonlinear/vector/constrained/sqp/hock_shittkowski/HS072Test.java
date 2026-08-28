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
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS TP72 (Schittkowski). Minimize 1 + sum(x_i), 2 NL inequalities. */
public class HS072Test {

    // Bounds: XL(i)=1e-3; XU(i)=1e5*(5-i) for i=1..4
    private static final double[] LB = { 1e-3, 1e-3, 1e-3, 1e-3 };
    private static final double[] UB = { 4e5, 3e5, 2e5, 1e5 };

    /** f(x) = 1 + x1 + x2 + x3 + x4. */
    private static class TP72Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 4; }
        @Override public double value(RealVector X) {
            return 1.0 + X.getEntry(0) + X.getEntry(1) + X.getEntry(2) + X.getEntry(3);
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /**
     * Inequalities in forma g(x) ≥ 0 (come da tua convenzione).
     * g1 = -4/x1 - 2.25/x2 - 1/x3 - 0.25/x4 + 0.0401 ≥ 0
     * g2 = -0.16/x1 - 0.36/x2 - 0.64*(1/x3 + 1/x4) + 0.010085 ≥ 0
     */
    private static class TP72Ineq extends InequalityConstraint {
        TP72Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public int dim() { return 4; }
        @Override public RealVector value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x2 = X.getEntry(1);
            final double x3 = X.getEntry(2);
            final double x4 = X.getEntry(3);

            final double g1 = -4.0/x1 - 2.25/x2 - 1.0/x3 - 0.25/x4 + 0.0401;
            final double g2 = -0.16/x1 - 0.36/x2 - 0.64*(1.0/x3 + 1.0/x4) + 0.010085;

            return new ArrayRealVector(new double[]{ g1, g2 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS072() {
        final InitialGuess guess = new InitialGuess(new double[]{ 1.0, 1.0, 1.0, 1.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new TP72Obj()),
                new TP72Ineq(),
                bounds
        );

        // FEX = 0.72767936D+3
        HSProblemTestUtils.assertExpectedObjective(727.67936, sol);
    }
}
