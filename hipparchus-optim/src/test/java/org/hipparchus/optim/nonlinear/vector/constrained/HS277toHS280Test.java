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

/**
 * HS277–HS280 (TP277–TP280).
 *
 * <p>Fortran model summary:
 * <ul>
 *   <li>Objective: f(x) = Σ_i c_i x_i, with c_i = Σ_j 1 / (i + j - 1).</li>
 *   <li>Inequalities: g_i(x) = Σ_j (x_j - 1) / (i + j - 1) >= 0, i = 1..N.</li>
 *   <li>Lower bounds x_i >= 0 are declared, but they are redundant for this model.</li>
 *   <li>Reference solution x* = 1, f* = Σ_i c_i.</li>
 * </ul>
 */
public class HS277toHS280Test {

    private static RealMatrix hilbert(final int n) {
        final double[][] h = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                h[i][j] = 1.0 / (i + j + 1.0);
            }
        }
        return new Array2DRowRealMatrix(h, false);
    }

    private static RealVector rowSums(final RealMatrix h) {
        final int n = h.getRowDimension();
        final double[] c = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += h.getEntry(i, j);
            }
            c[i] = sum;
        }
        return new ArrayRealVector(c, false);
    }

    private static final class HS277Objective extends TwiceDifferentiableFunction {
        private final RealVector c;

        HS277Objective(final RealVector c) {
            this.c = c;
        }

        @Override
        public int dim() {
            return c.getDimension();
        }

        @Override
        public double value(final RealVector x) {
            return c.dotProduct(x);
        }

        @Override
        public RealVector gradient(final RealVector x) {
            return c.copy();
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            return new Array2DRowRealMatrix(c.getDimension(), c.getDimension());
        }
    }

    private static final class HS277Ineq extends InequalityConstraint {
        private final int n;

        HS277Ineq(final int n) {
            super(new ArrayRealVector(n));
            this.n = n;
        }

        @Override
        public int dim() {
            return n;
        }

        @Override
        public RealVector value(final RealVector x) {
            final double[] g = new double[n];
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    sum += (x.getEntry(j) - 1.0) / (i + j + 1.0);
                }
                g[i] = sum;
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            final double[][] j = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < n; k++) {
                    j[i][k] = 1.0 / (i + k + 1.0);
                }
            }
            return new Array2DRowRealMatrix(j, false);
        }
    }

    private void runCase(final int n) {
        final RealMatrix h = hilbert(n);
        final RealVector c = rowSums(h);
        final double expected = c.dotProduct(new ArrayRealVector(n, 1.0));

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);

        final double[] start = new double[n]; // Fortran starts from 0.

        final LagrangeSolution sol = optimizer.optimize(
                option,
                new InitialGuess(start),
                new ObjectiveFunction(new HS277Objective(c)),
                null,
                new HS277Ineq(n),
                null);

        assertEquals(expected, sol.getValue(), 1.0e-6 * (Math.abs(expected) + 1.0), "objective mismatch");
    }

    @Test
    public void testHS277() {
        runCase(4);
    }

    @Test
    public void testHS278() {
        runCase(6);
    }

    @Test
    public void testHS279() {
        runCase(8);
    }

    @Test
    public void testHS280() {
        runCase(10);
    }
}