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
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS205Test {

    /** Beale-type least-squares (same model as HS203, different start). */
    static final class HS205Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        @Override public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double f1 = 1.5   - x1 * (1.0 - Math.pow(x2, 1));
            final double f2 = 2.25  - x1 * (1.0 - Math.pow(x2, 2));
            final double f3 = 2.625 - x1 * (1.0 - Math.pow(x2, 3));
            return f1*f1 + f2*f2 + f3*f3;
        }

        @Override public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double f1 = 1.5   - x1 * (1.0 - Math.pow(x2, 1));
            final double f2 = 2.25  - x1 * (1.0 - Math.pow(x2, 2));
            final double f3 = 2.625 - x1 * (1.0 - Math.pow(x2, 3));

            final double df1dx1 = (x2 - 1.0);
            final double df1dx2 = x1;

            final double df2dx1 = (x2*x2 - 1.0);
            final double df2dx2 = 2.0 * x1 * x2;

            final double df3dx1 = (x2*x2*x2 - 1.0);
            final double df3dx2 = 3.0 * x1 * x2 * x2;

            final double g1 = 2.0 * (f1*df1dx1 + f2*df2dx1 + f3*df3dx1);
            final double g2 = 2.0 * (f1*df1dx2 + f2*df2dx2 + f3*df3dx2);

            return new ArrayRealVector(new double[]{ g1, g2 }, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            return new Array2DRowRealMatrix(2, 2); // not used
        }
    }

    private static LagrangeSolution solve() {
        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        final double[] x0 = { 0.0, 0.0 }; // Fortran start
        return opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS205Obj())
        );
    }

    @Test
    public void testHS205() {
        final LagrangeSolution sol = solve();
        final double expected = 0.0; // fEx
        final double f = sol.getValue();
        assertEquals(expected, f, 1.0e-6 * (Math.abs(expected) + 1.0), "objective mismatch at optimum");
    }
}
