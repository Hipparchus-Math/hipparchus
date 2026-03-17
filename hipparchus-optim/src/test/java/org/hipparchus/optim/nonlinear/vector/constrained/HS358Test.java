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

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS358Test {

    private static final class HS358Obj extends TwiceDifferentiableFunction {
        private static final double[] Y = {
            0.844,0.908,0.932,0.936,0.925,0.908,0.881,0.850,0.818,0.784,0.751,
            0.718,0.685,0.658,0.628,0.603,0.580,0.558,0.538,0.522,0.506,0.490,
            0.478,0.467,0.457,0.448,0.438,0.431,0.424,0.420,0.414,0.411,0.406
        };

        @Override public int dim() { return 5; }

        @Override
        public double value(final RealVector x) {
            double fx = 0.0;
            for (int i = 0; i < 33; i++) {
                final double ti = i * 0.1 + 2.0;
                final double fi = Y[i] - (x.getEntry(0) + x.getEntry(1) * FastMath.exp(-x.getEntry(3) * ti) +
                                          x.getEntry(2) * FastMath.exp(-x.getEntry(4) * ti));
                fx += fi * fi;
            }
            return fx;
        }

        @Override public RealVector gradient(final RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(final RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    void testHS358() {
        final HS358Obj obj = new HS358Obj();
        final double[] xEx = {0.3754, 1.9358, -1.4647, 0.01287, 0.02212};
        assertEquals(0.546e-4, obj.value(new org.hipparchus.linear.ArrayRealVector(xEx, false)), 5.0e-7);

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final LagrangeSolution sol = optimizer.optimize(
                new InitialGuess(new double[] {0.5, 1.5, -1.0, 0.01, 0.02}),
                new ObjectiveFunction(obj),
                new SimpleBounds(
                        new double[] {-0.5, 1.5, -2.0, 0.001, 0.001},
                        new double[] { 0.5, 2.5, -1.0, 0.1,   0.1  }
                )
        );

        HSProblemTestUtils.assertBetterObjective(0.546e-4, sol);
    }
}
