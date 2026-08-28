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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Linear-quadratic optimal-control benchmark with 999 variables.
 */
public class LargeScale_OptimalControl999Test {

    /** Number of time intervals. */
    private static final int TIME_STEPS = 499;

    /** Number of state variables y_0,...,y_T. */
    private static final int STATE_COUNT = TIME_STEPS + 1;

    /** Total variables: T+1 states and T controls. */
    private static final int DIMENSION =
            2 * TIME_STEPS + 1;

    /** Number of equality constraints. */
    private static final int CONSTRAINTS =
            TIME_STEPS + 1;

    /** Integration step. */
    private static final double H =
            1.0 / TIME_STEPS;

    /** Control regularization. */
    private static final double ALPHA = 0.01;

    /** Maximum iterations. */
    private static final int MAX_ITERATIONS = 3000;

    /** Independently computed optimal objective. */
    private static final double EXPECTED_OBJECTIVE =
            0.11824299608478428;

    /** Objective tolerance. */
    private static final double OBJECTIVE_TOLERANCE = 1.0e-6;

    /** Equality residual tolerance. */
    private static final double FEASIBILITY_TOLERANCE = 1.0e-8;

    /** Quadratic tracking objective. */
    private static final class OptimalControlObjective
        extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIMENSION;
        }

        @Override
        public double value(final RealVector x) {

            double value = 0.0;

            for (int i = 0; i <= TIME_STEPS; ++i) {
                final double difference =
                        x.getEntry(i) - 1.0;

                value += H * difference * difference;
            }

            for (int i = 0; i < TIME_STEPS; ++i) {
                final double control =
                        x.getEntry(STATE_COUNT + i);

                value +=
                        ALPHA * H *
                        control * control;
            }

            return value;
        }

        @Override
        public RealVector gradient(final RealVector x) {

            final RealVector gradient =
                    new ArrayRealVector(DIMENSION);

            for (int i = 0; i <= TIME_STEPS; ++i) {
                gradient.setEntry(
                        i,
                        2.0 * H *
                        (x.getEntry(i) - 1.0));
            }

            for (int i = 0; i < TIME_STEPS; ++i) {
                final int index =
                        STATE_COUNT + i;

                gradient.setEntry(
                        index,
                        2.0 * ALPHA * H *
                        x.getEntry(index));
            }

            return gradient;
        }

        @Override
        public RealMatrix hessian(final RealVector x) {

            final double[] diagonal =
                    new double[DIMENSION];

            for (int i = 0; i <= TIME_STEPS; ++i) {
                diagonal[i] = 2.0 * H;
            }

            for (int i = 0; i < TIME_STEPS; ++i) {
                diagonal[STATE_COUNT + i] =
                        2.0 * ALPHA * H;
            }

            return MatrixUtils.createRealDiagonalMatrix(
                    diagonal);
        }
    }

    /** Discretized state equations. */
    private static final class OptimalControlEquality
        extends EqualityConstraint {

        OptimalControlEquality() {
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

            residual.setEntry(0, x.getEntry(0));

            for (int i = 0; i < TIME_STEPS; ++i) {
                residual.setEntry(
                        i + 1,
                        x.getEntry(i + 1) -
                        (1.0 - H) * x.getEntry(i) -
                        H * x.getEntry(STATE_COUNT + i));
            }

            return residual;
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {

            final RealMatrix jacobian =
                    MatrixUtils.createRealMatrix(
                            CONSTRAINTS,
                            DIMENSION);

            jacobian.setEntry(0, 0, 1.0);

            for (int i = 0; i < TIME_STEPS; ++i) {

                final int row = i + 1;

                jacobian.setEntry(
                        row,
                        i,
                        -(1.0 - H));

                jacobian.setEntry(
                        row,
                        i + 1,
                        1.0);

                jacobian.setEntry(
                        row,
                        STATE_COUNT + i,
                        -H);
            }

            return jacobian;
        }
    }

    /** Original zero initial point. */
    private static double[] initialPoint() {
        return new double[DIMENSION];
    }

    @Test
    public void testOptimalControl999() {

        final OptimalControlObjective objective =
                new OptimalControlObjective();

        final OptimalControlEquality equality =
                new OptimalControlEquality();

        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);

        final SQPOptimizerS2 optimizer =
                LargeScaleProblemTestUtils.newOptimizer();

        final long startTime = System.nanoTime();

        final LagrangeSolution solution = optimizer.optimize(
                new MaxIter(MAX_ITERATIONS),
                new InitialGuess(initialPoint()),
                new ObjectiveFunction(objective),
                equality,
                option);

        final double elapsedSeconds =
                (System.nanoTime() - startTime) * 1.0e-9;

        final double objectiveValue =
                solution.getValue();

        final double residualNorm =
                equality.value(solution.getX())
                        .getLInfNorm();

        System.out.printf(
                Locale.US,
                "%nOptimal Control 999%n" +
                "----------------------------%n" +
                "variables          : %d%n" +
                "equalities         : %d%n" +
                "objective          : %.16e%n" +
                "expected objective : %.16e%n" +
                "constraint LInf    : %.16e%n" +
                "time               : %.6f s%n",
                DIMENSION,
                CONSTRAINTS,
                objectiveValue,
                EXPECTED_OBJECTIVE,
                residualNorm,
                elapsedSeconds);

       LargeScaleProblemTestUtils.assertExpectedObjective(EXPECTED_OBJECTIVE,solution);
    }
}