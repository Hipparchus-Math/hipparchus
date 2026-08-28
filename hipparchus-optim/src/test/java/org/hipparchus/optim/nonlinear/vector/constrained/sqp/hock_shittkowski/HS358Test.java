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
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS358Test {

    // TP358 MODE=1 bounds
    private static final double[] LB = {-0.5, 1.5, -2.0, 0.001, 0.001};
    private static final double[] UB = { 0.5, 2.5, -1.0, 0.1,   0.1  };

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
            final double x1 = FastMath.max(LB[0], FastMath.min(UB[0], x.getEntry(0)));
            final double x2 = FastMath.max(LB[1], FastMath.min(UB[1], x.getEntry(1)));
            final double x3 = FastMath.max(LB[2], FastMath.min(UB[2], x.getEntry(2)));
            final double x4 = FastMath.max(LB[3], FastMath.min(UB[3], x.getEntry(3)));
            final double x5 = FastMath.max(LB[4], FastMath.min(UB[4], x.getEntry(4)));
            for (int i = 0; i < 33; i++) {
                final double ti = i * 0.1 + 2.0;
                final double fi = Y[i] - (x1 + x2 * FastMath.exp(-x4 * ti) +
                                          x3 * FastMath.exp(-x5 * ti));
                fx += fi * fi;
            }
            return fx;
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final double[] g = new double[5];
            final double x1 = FastMath.max(LB[0], FastMath.min(UB[0], x.getEntry(0)));
            final double x2 = FastMath.max(LB[1], FastMath.min(UB[1], x.getEntry(1)));
            final double x3 = FastMath.max(LB[2], FastMath.min(UB[2], x.getEntry(2)));
            final double x4 = FastMath.max(LB[3], FastMath.min(UB[3], x.getEntry(3)));
            final double x5 = FastMath.max(LB[4], FastMath.min(UB[4], x.getEntry(4)));

            for (int i = 0; i < 33; i++) {
                final double ti = i * 0.1 + 2.0;
                final double exp4 = FastMath.exp(-x4 * ti);
                final double exp5 = FastMath.exp(-x5 * ti);
                final double fi = Y[i] - (x1 + x2 * exp4 + x3 * exp5);
                final double[] dfi = {
                        -1.0,
                        -exp4,
                        -exp5,
                        x2 * exp4 * ti,
                        x3 * exp5 * ti
                };
                for (int j = 0; j < 5; j++) {
                    g[j] += 2.0 * fi * dfi[j];
                }
            }
            return new ArrayRealVector(g, false);
        }
        @Override public RealMatrix hessian(final RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    void testHS358() {
        final HS358Obj obj = new HS358Obj();
        

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final LagrangeSolution sol = optimizer.optimize(
                HSProblemTestUtils.newExternalOption(),
                new InitialGuess(new double[] {0.5, 1.5, -1.0, 0.01, 0.02}),
                new ObjectiveFunction(obj),
                new SimpleBounds(LB, UB)
        );

        HSProblemTestUtils.assertBetterObjective(0.546e-4, sol);
    }
}