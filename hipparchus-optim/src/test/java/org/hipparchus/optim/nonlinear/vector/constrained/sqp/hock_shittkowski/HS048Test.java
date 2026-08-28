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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

public class HS048Test {

    private static final class HS048Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 5; }

        @Override
        public double value(final RealVector x) {
            return (x.getEntry(0) - 1.0) * (x.getEntry(0) - 1.0) +
                   (x.getEntry(1) - x.getEntry(2)) * (x.getEntry(1) - x.getEntry(2)) +
                   (x.getEntry(3) - x.getEntry(4)) * (x.getEntry(3) - x.getEntry(4));
        }

        @Override public RealVector gradient(final RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(final RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static final class HS048Eq extends EqualityConstraint {
        HS048Eq() { super(new ArrayRealVector(new double[] {0.0, 0.0})); }

        @Override
        public RealVector value(final RealVector x) {
            return new ArrayRealVector(new double[] {
                x.getEntry(0) + x.getEntry(1) + x.getEntry(2) + x.getEntry(3) + x.getEntry(4) - 5.0,
                x.getEntry(2) - 2.0 * (x.getEntry(3) + x.getEntry(4)) + 3.0
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            return new Array2DRowRealMatrix(new double[][] {
                {1.0, 1.0, 1.0, 1.0, 1.0},
                {0.0, 0.0, 1.0, -2.0, -2.0}
            }, false);
        }

        @Override public int dim() { return 5; }
    }

    @Test
    void testHS048() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final LagrangeSolution sol = optimizer.optimize(
                new InitialGuess(new double[] {3.0, 5.0, -3.0, 2.0, -2.0}),
                new ObjectiveFunction(new HS048Obj()),
                new HS048Eq()
        );

        HSProblemTestUtils.assertExpectedObjective(0.0, sol);
    }
}
