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

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

public class HS025Test {

    private static final class HS025Obj extends TwiceDifferentiableFunction {
        @Override
        public int dim() {
            return 3;
        }

        @Override
        public double value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);

            double fx = 0.0;
            for (int i = 1; i <= 99; i++) {
                final double u = 25.0 + FastMath.pow(-50.0 * FastMath.log(0.01 * i), 2.0 / 3.0);
                final double v = u - x2;
                if (v < 0.0) {
                    final double p1 = x1 - 5.0;
                    final double p2 = x2 - 5.0;
                    final double p3 = x3 - 5.0;
                    return p1 * p1 + p2 * p2 + p3 * p3;
                }
                final double b = FastMath.exp(-FastMath.pow(v, x3) / x1);
                final double a = b - 0.01 * i;
                fx += a * a;
            }
            return fx;
        }

        @Override
        public RealVector gradient(final RealVector x) {
            throw new UnsupportedOperationException();
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            throw new UnsupportedOperationException();
        }
    }

    @Test
    void testHS025() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final LagrangeSolution sol = optimizer.optimize(
                new InitialGuess(new double[] {100.0, 12.5, 3.0}),
                new ObjectiveFunction(new HS025Obj()),
                new SimpleBounds(
                        new double[] {0.1, 1.0e-5, 1.0e-5},
                        new double[] {100.0, 25.6, 5.0}
                )
        );

        HSProblemTestUtils.assertExpectedObjective(0.0, sol);
    }
}
