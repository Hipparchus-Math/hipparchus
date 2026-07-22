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
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

/**
 * PROB.FOR has no TP200 entry in this tree (it starts at TP201 in that range).
 * We keep HS200 aligned with the first available neighboring quadratic template.
 */
public class HS200Test {

    static final class HS200Objective extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        @Override
        public double value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            return 4.0 * (x1 - 5.0) * (x1 - 5.0) + (x2 - 6.0) * (x2 - 6.0);
        }

        @Override
        public RealVector gradient(final RealVector x) {
            return new ArrayRealVector(new double[] {
                8.0 * (x.getEntry(0) - 5.0),
                2.0 * (x.getEntry(1) - 6.0)
            }, false);
        }

        @Override public RealMatrix hessian(final RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    void testHS200() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final LagrangeSolution sol = optimizer.optimize(
                new InitialGuess(new double[] {8.0, 9.0}),
                new ObjectiveFunction(new HS200Objective()),
                new SimpleBounds(new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
                                 new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY})
        );

        HSProblemTestUtils.assertExpectedObjective(0.0, sol);
    }
}

