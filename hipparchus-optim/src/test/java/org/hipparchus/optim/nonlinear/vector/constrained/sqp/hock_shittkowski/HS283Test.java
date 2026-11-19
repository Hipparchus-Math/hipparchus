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
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/** HS283 (n=10): f(x) = ( Σ i^3 (x_i-1)^2 )^3, minimum at x = 1. */
public class HS283Test {

    /** Objective and gradient. */
    static final class HS283Objective extends TwiceDifferentiableFunction {
        @Override public int dim() { return 10; }

        @Override public double value(RealVector x) {
            double s = 0.0;
            for (int i = 0; i < 10; i++) {
                double wi = (1.0*i + 1.0)*(1.0*i + 1.0)*(1.0*i + 1.0);         
                double di = x.getEntry(i) - 1.0;
                s += wi * di * di;
            }
            return s * s * s; // s^3
        }

        @Override public RealVector gradient(RealVector x) {
            // s = Σ w_i (x_i-1)^2 ; g_i = 6 w_i (x_i-1) s^2
            double s = 0.0;
            double[] w = new double[10];
            double[] d = new double[10];
            for (int i = 0; i < 10; i++) {
                w[i] = Math.pow(i + 1, 3);
                d[i] = x.getEntry(i) - 1.0;
                s += w[i] * d[i] * d[i];
            }
            double coeff = 6.0 * s * s;
            double[] g = new double[10];
            for (int i = 0; i < 10; i++) {
                g[i] = coeff * w[i] * d[i];
            }
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            // Not required by the solver.
            return new Array2DRowRealMatrix(10, 10);
        }
    }

    private LagrangeSolution solve() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final double[] start = new double[10]; // all zeros as in the Fortran setup

      

        return optimizer.optimize(
                new InitialGuess(start),
                new ObjectiveFunction(new HS283Objective())
                
        );
    }

    @Test
    public void testHS283() {
        LagrangeSolution sol = solve();
        double f = sol.getValue();
        double fEx = 0.0;
        assertEquals(fEx, f, 1.0e-6 * (Math.abs(fEx) + 1.0), "objective mismatch at optimum");
    }
}

