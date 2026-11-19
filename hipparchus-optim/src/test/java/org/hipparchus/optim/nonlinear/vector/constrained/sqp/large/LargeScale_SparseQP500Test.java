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
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

/**
 * SparseQP benchmark with 500 variables and 500 cyclic inequalities.
 *
 * <p>The mathematical problem is:</p>
 *
 * <pre>
 * minimize 0.5 x' Q x - sum(x_i)
 *
 * subject to
 *
 *     x_j + x_(j+1) + x_(j+2) <= 2.5
 *     0 <= x_i <= 10
 *
 * where Q has diagonal 4 and first off-diagonal -1.
 * </pre>
 */
public class LargeScale_SparseQP500Test {

    /** Number of variables and inequalities. */
    private static final int DIMENSION = 500;

    /** Maximum iterations. */
    private static final int MAX_ITERATIONS = 3000;

    /** Independently computed optimal objective. */
    private static final double EXPECTED_OBJECTIVE =
            -124.81698729810776;

    /** Objective tolerance. */
    private static final double OBJECTIVE_TOLERANCE = 1.0e-6;

    /** Feasibility tolerance. */
    private static final double FEASIBILITY_TOLERANCE = 1.0e-8;

    /** Convex quadratic objective. */
    private static final class SparseQPObjective
        extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIMENSION;
        }

        @Override
        public double value(final RealVector x) {

            double value = 0.0;

            for (int i = 0; i < DIMENSION; ++i) {
                final double xi = x.getEntry(i);

                value += 2.0 * xi * xi;
                value -= xi;

                if (i < DIMENSION - 1) {
                    value -= xi * x.getEntry(i + 1);
                }
            }

            return value;
        }

        @Override
        public RealVector gradient(final RealVector x) {

            final RealVector gradient =
                    new ArrayRealVector(DIMENSION);

            for (int i = 0; i < DIMENSION; ++i) {

                double entry =
                        4.0 * x.getEntry(i) - 1.0;

                if (i > 0) {
                    entry -= x.getEntry(i - 1);
                }

                if (i < DIMENSION - 1) {
                    entry -= x.getEntry(i + 1);
                }

                gradient.setEntry(i, entry);
            }

            return gradient;
        }

        @Override
        public RealMatrix hessian(final RealVector x) {

            final RealMatrix hessian =
                    MatrixUtils.createRealMatrix(
                            DIMENSION,
                            DIMENSION);

            for (int i = 0; i < DIMENSION; ++i) {
                hessian.setEntry(i, i, 4.0);

                if (i < DIMENSION - 1) {
                    hessian.setEntry(i, i + 1, -1.0);
                    hessian.setEntry(i + 1, i, -1.0);
                }
            }

            return hessian;
        }
    }

    /**
     * Cyclic inequalities converted to Hipparchus convention:
     *
     * <pre>
     * 2.5 - x_j - x_(j+1) - x_(j+2) >= 0.
     * </pre>
     */
    private static final class SparseQPInequality
        extends InequalityConstraint {

        SparseQPInequality() {
            super(new ArrayRealVector(DIMENSION));
        }

        @Override
        public int dim() {
            return DIMENSION;
        }

        @Override
        public RealVector value(final RealVector x) {

            final RealVector residual =
                    new ArrayRealVector(DIMENSION);

            for (int j = 0; j < DIMENSION; ++j) {

                final int j1 = (j + 1) % DIMENSION;
                final int j2 = (j + 2) % DIMENSION;

                residual.setEntry(
                        j,
                        2.5 -
                        x.getEntry(j) -
                        x.getEntry(j1) -
                        x.getEntry(j2));
            }

            return residual;
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {

            final RealMatrix jacobian =
                    MatrixUtils.createRealMatrix(
                            DIMENSION,
                            DIMENSION);

            for (int j = 0; j < DIMENSION; ++j) {
                jacobian.setEntry(j, j, -1.0);
                jacobian.setEntry(
                        j,
                        (j + 1) % DIMENSION,
                        -1.0);
                jacobian.setEntry(
                        j,
                        (j + 2) % DIMENSION,
                        -1.0);
            }

            return jacobian;
        }
    }

    /** Original initial point, x_i = 0.5. */
    private static double[] initialPoint() {

        final double[] initial =
                new double[DIMENSION];

        for (int i = 0; i < DIMENSION; ++i) {
            initial[i] = 0.5;
        }

        return initial;
    }

    /** Original bounds, 0 <= x_i <= 10. */
    private static SimpleBounds bounds() {

        final double[] lower =
                new double[DIMENSION];

        final double[] upper =
                new double[DIMENSION];

        for (int i = 0; i < DIMENSION; ++i) {
            lower[i] = 0.0;
            upper[i] = 10.0;
        }

        return new SimpleBounds(lower, upper);
    }

    @Test
    public void testSparseQP500() {

        final SparseQPObjective objective =
                new SparseQPObjective();

        final SparseQPInequality inequality =
                new SparseQPInequality();

        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);

        // just to improve test coverage…
        option.setGradientEps(SQPOption.DEFAULT_GRAD_EPS);
        option.setEpsQP(SQPOption.DEFAULT_EPSILON_QP);
        option.setAlphaMin(SQPOption.DEFAULT_ALPHA_MIN);
        option.setHistory(0);
        option.setMaxIteration(SQPOption.DEFAULT_MAX_ITERATION);
        option.setConvergenceFunction(null);

        final SQPOptimizerS2 optimizer =
                LargeScaleProblemTestUtils.newOptimizer();

        final long startTime = System.nanoTime();

        final LagrangeSolution solution = optimizer.optimize(
                new MaxIter(MAX_ITERATIONS),
                new InitialGuess(initialPoint()),
                new ObjectiveFunction(objective),
                inequality,
                bounds(),
                option);

        final double elapsedSeconds =
                (System.nanoTime() - startTime) * 1.0e-9;

        final double objectiveValue =
                solution.getValue();

        final RealVector inequalityValue =
                inequality.value(solution.getX());

        double maximumViolation = 0.0;

        for (int i = 0; i < DIMENSION; ++i) {
            maximumViolation =
                    Math.max(
                            maximumViolation,
                            -inequalityValue.getEntry(i));
        }

        System.out.printf(
                Locale.US,
                "%nSparseQP 500%n" +
                "----------------------------%n" +
                "variables          : %d%n" +
                "inequalities       : %d%n" +
                "objective          : %.16e%n" +
                "expected objective : %.16e%n" +
                "maximum violation  : %.16e%n" +
                "time               : %.6f s%n",
                DIMENSION,
                DIMENSION,
                objectiveValue,
                EXPECTED_OBJECTIVE,
                maximumViolation,
                elapsedSeconds);

        LargeScaleProblemTestUtils.assertExpectedObjective(EXPECTED_OBJECTIVE,solution);
    }
}
