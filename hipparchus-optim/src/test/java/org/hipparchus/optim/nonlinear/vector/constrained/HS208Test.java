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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/** HS208: Rosenbrock with C = 100 (unconstrained). */
public class HS208Test {

    /** Objective f(x) = 100*(x2 - x1^2)^2 + (1 - x1)^2, with analytic gradient. */
    static final class HS208Objective extends TwiceDifferentiableFunction {
        private static final double C = 100.0;

        @Override public int dim() { return 2; }

        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return C * Math.pow(x2 - x1 * x1, 2) + Math.pow(1.0 - x1, 2);
        }

        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double g1 = -4.0 * C * (x2 - x1 * x1) * x1 - 2.0 * (1.0 - x1);
            double g2 =  2.0 * C * (x2 - x1 * x1);
            return new ArrayRealVector(new double[] { g1, g2 }, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            return new Array2DRowRealMatrix(2, 2); // not used by the solver
        }
    }

    @Test
    public void testHS208() {
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);

        double[] start = { -1.2, 1.0 };
        LagrangeSolution sol = optimizer.optimize(
                new InitialGuess(start),
                new ObjectiveFunction(new HS208Objective())
        );

        double expected = 0.0; // known optimum at (1,1)
        assertEquals(expected, sol.getValue(),
                     1.0e-6 * (Math.abs(expected) + 1.0),
                     "objective mismatch");
    }
}

