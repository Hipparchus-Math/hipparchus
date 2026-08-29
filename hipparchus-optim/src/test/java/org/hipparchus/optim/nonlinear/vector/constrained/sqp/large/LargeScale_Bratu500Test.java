/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.large;


import java.util.Locale;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Reduced Bratu benchmark with 500 variables and 498 nonlinear equalities.
 */
public class LargeScale_Bratu500Test {

    /** Number of variables. */
    private static final int DIMENSION = 500;

    /** Number of nonlinear equality constraints. */
    private static final int CONSTRAINTS = DIMENSION - 2;

    /** Bratu parameter. */
    private static final double LAMBDA = 1.0;

    /** Mesh spacing used by the original Rust implementation. */
    private static final double H = 1.0 / (DIMENSION + 1.0);

    /** Squared mesh spacing. */
    private static final double H2 = H * H;

    /** Maximum number of iterations. */
    private static final int MAX_ITERATIONS = 3000;

    /** Residual tolerance. */
    private static final double RESIDUAL_TOLERANCE = 1.0e-6;

    /** Constant zero objective. */
    private static final class BratuObjective
        extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIMENSION;
        }

        @Override
        public double value(final RealVector x) {
            return 0.0;
        }

        @Override
        public RealVector gradient(final RealVector x) {
            return new ArrayRealVector(DIMENSION);
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            return MatrixUtils.createRealMatrix(DIMENSION, DIMENSION);
        }
    }

    /** Bratu finite-difference equations. */
    private static final class BratuEquality
        extends EqualityConstraint {

        BratuEquality() {
            super(new ArrayRealVector(CONSTRAINTS));
        }

        @Override
        public int dim() {
            return DIMENSION;
        }

        @Override
        public RealVector value(final RealVector x) {

            final RealVector residual =
                    new ArrayRealVector(CONSTRAINTS);

            for (int j = 0; j < CONSTRAINTS; ++j) {
                final int i = j + 1;

                residual.setEntry(
                        j,
                        (-x.getEntry(i - 1) +
                         2.0 * x.getEntry(i) -
                         x.getEntry(i + 1)) / H2 -
                        LAMBDA * FastMath.exp(x.getEntry(i)));
            }

            return residual;
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {

            final RealMatrix jacobian =
                    MatrixUtils.createRealMatrix(
                            CONSTRAINTS,
                            DIMENSION);

            final double offDiagonal = -1.0 / H2;

            for (int j = 0; j < CONSTRAINTS; ++j) {
                final int i = j + 1;

                jacobian.setEntry(j, i - 1, offDiagonal);

                jacobian.setEntry(
                        j,
                        i,
                        2.0 / H2 -
                        LAMBDA * FastMath.exp(x.getEntry(i)));

                jacobian.setEntry(j, i + 1, offDiagonal);
            }

            return jacobian;
        }
    }

    /** Original zero starting point. */
    private static double[] initialPoint() {
        return new double[DIMENSION];
    }

    /**
     * Bounds reproducing the Rust benchmark:
     * the two boundary variables are fixed to zero.
     */
    private static SimpleBounds bounds() {

        final double[] lower = new double[DIMENSION];
        final double[] upper = new double[DIMENSION];

        for (int i = 0; i < DIMENSION; ++i) {
            lower[i] = Double.NEGATIVE_INFINITY;
            upper[i] = Double.POSITIVE_INFINITY;
        }

        lower[0] = 0.0;
        upper[0] = 0.0;

        lower[DIMENSION - 1] = 0.0;
        upper[DIMENSION - 1] = 0.0;

        return new SimpleBounds(lower, upper);
    }

    @Test
    public void testBratu500() {

        final BratuObjective objective =
                new BratuObjective();

        final BratuEquality equality =
                new BratuEquality();

        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        option.setMaxIteration(1000);
        final SQPOptimizerS2 optimizer =
                LargeScaleProblemTestUtils.newOptimizer();

        final long startTime = System.nanoTime();

        final LagrangeSolution solution = optimizer.optimize(
                new MaxIter(MAX_ITERATIONS),
                new InitialGuess(initialPoint()),
                new ObjectiveFunction(objective),
                equality,
                bounds(),
                option);

        final double elapsedSeconds =
                (System.nanoTime() - startTime) * 1.0e-9;

        final RealVector point = solution.getX();
        final double residualNorm =
                equality.value(point).getLInfNorm();

        System.out.printf(
                Locale.US,
                "%nBratu 500%n" +
                "----------------------------%n" +
                "variables          : %d%n" +
                "equalities         : %d%n" +
                "objective          : %.16e%n" +
                "constraint LInf    : %.16e%n" +
                "time               : %.6f s%n",
                DIMENSION,
                CONSTRAINTS,
                solution.getValue(),
                residualNorm,
                elapsedSeconds);

        assertTrue(
                Double.isFinite(residualNorm),
                "The Bratu residual is not finite.");

        assertTrue(
                residualNorm <= RESIDUAL_TOLERANCE,
                "Bratu 500 did not reach feasibility. Residual = " +
                residualNorm);
    }
}
