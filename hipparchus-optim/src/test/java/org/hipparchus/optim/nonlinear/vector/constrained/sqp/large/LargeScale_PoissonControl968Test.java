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
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Poisson optimal-control benchmark on a 22 by 22 grid.
 *
 * <p>The problem contains 484 state variables, 484 control variables
 * and 484 linear equality constraints.</p>
 */
public class LargeScale_PoissonControl968Test {

    /** Grid size in each spatial direction. */
    private static final int GRID = 22;

    /** Number of grid points. */
    private static final int GRID_SIZE = GRID * GRID;

    /** State plus control variables. */
    private static final int DIMENSION = 2 * GRID_SIZE;

    /** Number of equality constraints. */
    private static final int CONSTRAINTS = GRID_SIZE;

    /** Mesh spacing. */
    private static final double H =
            1.0 / (GRID + 1.0);

    /** Squared mesh spacing. */
    private static final double H2 = H * H;

    /** Control regularization. */
    private static final double ALPHA = 0.01;

    /** Maximum number of iterations. */
    private static final int MAX_ITERATIONS = 3000;

    /** Independently computed reference objective. */
    private static final double EXPECTED_OBJECTIVE =
            0.09940761183819619;

    /** Objective tolerance. */
    private static final double OBJECTIVE_TOLERANCE = 1.0e-6;

    /** Equality residual tolerance. */
    private static final double FEASIBILITY_TOLERANCE = 1.0e-7;

    /** Return the state-variable index. */
    private static int stateIndex(final int i,
                                  final int j) {
        return i + j * GRID;
    }

    /** Return the control-variable index. */
    private static int controlIndex(final int i,
                                    final int j) {
        return GRID_SIZE + i + j * GRID;
    }

    /** Desired state at a grid point. */
    private static double desiredState(final int i,
                                       final int j) {

        final double x =
                (i + 1.0) * H;

        final double y =
                (j + 1.0) * H;

        return FastMath.sin(FastMath.PI * x) *
               FastMath.sin(FastMath.PI * y);
    }

    /** Quadratic state-tracking objective. */
    private static final class PoissonObjective
        extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIMENSION;
        }

        @Override
        public double value(final RealVector x) {

            double value = 0.0;

            for (int j = 0; j < GRID; ++j) {
                for (int i = 0; i < GRID; ++i) {

                    final double state =
                            x.getEntry(stateIndex(i, j));

                    final double difference =
                            state - desiredState(i, j);

                    value +=
                            0.5 * H2 *
                            difference * difference;

                    final double control =
                            x.getEntry(controlIndex(i, j));

                    value +=
                            0.5 * ALPHA * H2 *
                            control * control;
                }
            }

            return value;
        }

        @Override
        public RealVector gradient(final RealVector x) {

            final RealVector gradient =
                    new ArrayRealVector(DIMENSION);

            for (int j = 0; j < GRID; ++j) {
                for (int i = 0; i < GRID; ++i) {

                    final int state =
                            stateIndex(i, j);

                    final int control =
                            controlIndex(i, j);

                    gradient.setEntry(
                            state,
                            H2 *
                            (x.getEntry(state) -
                             desiredState(i, j)));

                    gradient.setEntry(
                            control,
                            ALPHA * H2 *
                            x.getEntry(control));
                }
            }

            return gradient;
        }

        @Override
        public RealMatrix hessian(final RealVector x) {

            final double[] diagonal =
                    new double[DIMENSION];

            for (int i = 0; i < GRID_SIZE; ++i) {
                diagonal[i] = H2;
                diagonal[GRID_SIZE + i] =
                        ALPHA * H2;
            }

            return MatrixUtils.createRealDiagonalMatrix(
                    diagonal);
        }
    }

    /** Discrete Poisson equations. */
    private static final class PoissonEquality
        extends EqualityConstraint {

        PoissonEquality() {
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

            for (int j = 0; j < GRID; ++j) {
                for (int i = 0; i < GRID; ++i) {

                    final int row =
                            j * GRID + i;

                    double laplacian =
                            4.0 *
                            x.getEntry(stateIndex(i, j));

                    if (i > 0) {
                        laplacian -=
                                x.getEntry(
                                        stateIndex(i - 1, j));
                    }

                    if (i < GRID - 1) {
                        laplacian -=
                                x.getEntry(
                                        stateIndex(i + 1, j));
                    }

                    if (j > 0) {
                        laplacian -=
                                x.getEntry(
                                        stateIndex(i, j - 1));
                    }

                    if (j < GRID - 1) {
                        laplacian -=
                                x.getEntry(
                                        stateIndex(i, j + 1));
                    }

                    residual.setEntry(
                            row,
                            laplacian / H2 -
                            x.getEntry(controlIndex(i, j)));
                }
            }

            return residual;
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {

            final RealMatrix jacobian =
                    MatrixUtils.createRealMatrix(
                            CONSTRAINTS,
                            DIMENSION);

            final double center = 4.0 / H2;
            final double neighbor = -1.0 / H2;

            for (int j = 0; j < GRID; ++j) {
                for (int i = 0; i < GRID; ++i) {

                    final int row =
                            j * GRID + i;

                    jacobian.setEntry(
                            row,
                            stateIndex(i, j),
                            center);

                    if (i > 0) {
                        jacobian.setEntry(
                                row,
                                stateIndex(i - 1, j),
                                neighbor);
                    }

                    if (i < GRID - 1) {
                        jacobian.setEntry(
                                row,
                                stateIndex(i + 1, j),
                                neighbor);
                    }

                    if (j > 0) {
                        jacobian.setEntry(
                                row,
                                stateIndex(i, j - 1),
                                neighbor);
                    }

                    if (j < GRID - 1) {
                        jacobian.setEntry(
                                row,
                                stateIndex(i, j + 1),
                                neighbor);
                    }

                    jacobian.setEntry(
                            row,
                            controlIndex(i, j),
                            -1.0);
                }
            }

            return jacobian;
        }
    }

    /** Original zero initial point. */
    private static double[] initialPoint() {
        return new double[DIMENSION];
    }

    @Test
    public void testPoissonControl968() {

        final PoissonObjective objective =
                new PoissonObjective();

        final PoissonEquality equality =
                new PoissonEquality();

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
                "%nPoisson Control 968%n" +
                "----------------------------%n" +
                "grid               : %d x %d%n" +
                "variables          : %d%n" +
                "equalities         : %d%n" +
                "objective          : %.16e%n" +
                "expected objective : %.16e%n" +
                "constraint LInf    : %.16e%n" +
                "time               : %.6f s%n",
                GRID,
                GRID,
                DIMENSION,
                CONSTRAINTS,
                objectiveValue,
                EXPECTED_OBJECTIVE,
                residualNorm,
                elapsedSeconds);

         LargeScaleProblemTestUtils.assertExpectedObjective(EXPECTED_OBJECTIVE,solution);
    }
}