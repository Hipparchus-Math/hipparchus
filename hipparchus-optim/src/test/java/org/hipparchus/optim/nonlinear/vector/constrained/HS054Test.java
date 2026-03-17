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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

public class HS054Test {

    private static final class HS054Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 6; }

        @Override
        public double value(final RealVector x) {
            final double v1 = x.getEntry(0) - 1.0e4;
            final double v2 = x.getEntry(1) - 1.0;
            final double v3 = x.getEntry(2) - 2.0e6;
            final double v4 = x.getEntry(3) - 10.0;
            final double v5 = x.getEntry(4) - 1.0e-3;
            final double v6 = x.getEntry(5) - 1.0e8;
            final double q = (1.5625e-8 * v1 * v1 + 5.0e-5 * v1 * v2 + v2 * v2) / 0.96 +
                             v3 * v3 / 4.9e13 + 4.0e-4 * v4 * v4 + 4.0e2 * v5 * v5 + 4.0e-18 * v6 * v6;
            return -FastMath.exp(-0.5 * q);
        }

        @Override public RealVector gradient(final RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(final RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static final class HS054Eq extends EqualityConstraint {
        HS054Eq() { super(new ArrayRealVector(new double[] {0.0})); }

        @Override
        public RealVector value(final RealVector x) {
            return new ArrayRealVector(new double[] {x.getEntry(0) + 4.0e3 * x.getEntry(1) - 1.76e4}, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            return new Array2DRowRealMatrix(new double[][] {{1.0, 4.0e3, 0.0, 0.0, 0.0, 0.0}}, false);
        }

        @Override public int dim() { return 6; }
    }

    @Test
    void testHS054() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final LagrangeSolution sol = optimizer.optimize(
                new InitialGuess(new double[] {6.0e3, 1.5, 4.0e6, 2.0, 3.0e-3, 5.0e7}),
                new ObjectiveFunction(new HS054Obj()),
                new HS054Eq(),
                new SimpleBounds(
                        new double[] {0.0, -10.0, 0.0, 0.0, 0.0, 0.0},
                        new double[] {2.0e4, 10.0, 1.0e7, 20.0, 1.0, 2.0e8}
                )
        );

        HSProblemTestUtils.assertExpectedObjective(-FastMath.exp(-27.0 / 280.0), sol);
    }
}
