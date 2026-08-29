/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.CUTEst;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.OpenMapRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

/**
 * CUTEst/COPS benchmark {@code ELEC}.
 *
 * <p>The problem places 200 electrons on the unit sphere and minimizes
 * their Coulomb potential. The exact CUTEst dimensions are retained:
 * 600 variables and 200 nonlinear equality constraints.</p>
 */
@Disabled
public class ELECTest {

    /** Number of electrons. */
    private static final int POINTS = 200;

    /** Three Cartesian coordinates for every electron. */
    private static final int N = 3 * POINTS;

    /** One unit-sphere equality for every electron. */
    private static final int M = POINTS;

    /** Standard CUTEst reference objective for NP = 200. */
    private static final double EXPECTED_OBJECTIVE = 1.84389e+04;

    /**
     * Coulomb potential:
     *
     * <pre>
     * f(p) = sum(i=0..POINTS-2) sum(j=i+1..POINTS-1)
     *        1 / ||p_i - p_j||_2.
     * </pre>
     */
    private static final class Objective extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            double value = 0.0;

            for (int i = 0; i < POINTS - 1; ++i) {
                final double xi = point.getEntry(xIndex(i));
                final double yi = point.getEntry(yIndex(i));
                final double zi = point.getEntry(zIndex(i));

                for (int j = i + 1; j < POINTS; ++j) {
                    final double dx = xi - point.getEntry(xIndex(j));
                    final double dy = yi - point.getEntry(yIndex(j));
                    final double dz = zi - point.getEntry(zIndex(j));

                    final double distanceSquared =
                            dx * dx + dy * dy + dz * dz;

                    value += 1.0 / FastMath.sqrt(distanceSquared);
                }
            }

            return value;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            final RealVector gradient = new ArrayRealVector(N);

            for (int i = 0; i < POINTS - 1; ++i) {
                final double xi = point.getEntry(xIndex(i));
                final double yi = point.getEntry(yIndex(i));
                final double zi = point.getEntry(zIndex(i));

                for (int j = i + 1; j < POINTS; ++j) {
                    final double dx = xi - point.getEntry(xIndex(j));
                    final double dy = yi - point.getEntry(yIndex(j));
                    final double dz = zi - point.getEntry(zIndex(j));

                    final double distanceSquared =
                            dx * dx + dy * dy + dz * dz;

                    final double inverseDistance =
                            1.0 / FastMath.sqrt(distanceSquared);

                    final double inverseDistanceCubed =
                            inverseDistance / distanceSquared;

                    final double gx = -dx * inverseDistanceCubed;
                    final double gy = -dy * inverseDistanceCubed;
                    final double gz = -dz * inverseDistanceCubed;

                    gradient.addToEntry(xIndex(i), gx);
                    gradient.addToEntry(yIndex(i), gy);
                    gradient.addToEntry(zIndex(i), gz);

                    gradient.addToEntry(xIndex(j), -gx);
                    gradient.addToEntry(yIndex(j), -gy);
                    gradient.addToEntry(zIndex(j), -gz);
                }
            }

            return gradient;
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            final RealMatrix hessian = new OpenMapRealMatrix(N, N);
            final double[] difference = new double[3];

            for (int i = 0; i < POINTS - 1; ++i) {
                for (int j = i + 1; j < POINTS; ++j) {
                    difference[0] =
                            point.getEntry(xIndex(i)) -
                            point.getEntry(xIndex(j));

                    difference[1] =
                            point.getEntry(yIndex(i)) -
                            point.getEntry(yIndex(j));

                    difference[2] =
                            point.getEntry(zIndex(i)) -
                            point.getEntry(zIndex(j));

                    final double distanceSquared =
                            difference[0] * difference[0] +
                            difference[1] * difference[1] +
                            difference[2] * difference[2];

                    final double inverseDistance =
                            1.0 / FastMath.sqrt(distanceSquared);

                    final double inverseDistanceCubed =
                            inverseDistance / distanceSquared;

                    final double inverseDistanceFifth =
                            inverseDistanceCubed / distanceSquared;

                    for (int row = 0; row < 3; ++row) {
                        final int rowI = coordinateIndex(row, i);
                        final int rowJ = coordinateIndex(row, j);

                        for (int column = 0; column < 3; ++column) {
                            final int columnI =
                                    coordinateIndex(column, i);

                            final int columnJ =
                                    coordinateIndex(column, j);

                            final double blockEntry =
                                    3.0 *
                                    difference[row] *
                                    difference[column] *
                                    inverseDistanceFifth -
                                    (row == column ?
                                     inverseDistanceCubed :
                                     0.0);

                            hessian.addToEntry(
                                    rowI,
                                    columnI,
                                    blockEntry);

                            hessian.addToEntry(
                                    rowJ,
                                    columnJ,
                                    blockEntry);

                            hessian.addToEntry(
                                    rowI,
                                    columnJ,
                                    -blockEntry);

                            hessian.addToEntry(
                                    rowJ,
                                    columnI,
                                    -blockEntry);
                        }
                    }
                }
            }

            return hessian;
        }
    }

    /**
     * Unit-sphere equalities:
     *
     * <pre>
     * x_i^2 + y_i^2 + z_i^2 - 1 = 0,
     * i = 0,...,199.
     * </pre>
     */
    private static final class SphereEqualities
            extends EqualityConstraint {

        SphereEqualities() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {
            final RealVector constraints =
                    new ArrayRealVector(M);

            for (int i = 0; i < POINTS; ++i) {
                final double x = point.getEntry(xIndex(i));
                final double y = point.getEntry(yIndex(i));
                final double z = point.getEntry(zIndex(i));

                constraints.setEntry(
                        i,
                        x * x + y * y + z * z - 1.0);
            }

            return constraints;
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {
            final RealMatrix jacobian =
                    new OpenMapRealMatrix(M, N);

            for (int i = 0; i < POINTS; ++i) {
                jacobian.setEntry(
                        i,
                        xIndex(i),
                        2.0 * point.getEntry(xIndex(i)));

                jacobian.setEntry(
                        i,
                        yIndex(i),
                        2.0 * point.getEntry(yIndex(i)));

                jacobian.setEntry(
                        i,
                        zIndex(i),
                        2.0 * point.getEntry(zIndex(i)));
            }

            return jacobian;
        }
    }

    /**
     * Deterministic quasi-uniform feasible starting point from ELEC.SIF.
     */
    private static double[] initialPoint() {
        final double[] initial = new double[N];

        for (int i = 0; i < POINTS; ++i) {
            final double theta =
                    2.0 * FastMath.PI * (i + 1.0) / POINTS;

            final double phi =
                    FastMath.PI * i / POINTS;

            final double sinPhi =
                    FastMath.sin(phi);

            initial[xIndex(i)] =
                    FastMath.cos(theta) * sinPhi;

            initial[yIndex(i)] =
                    FastMath.sin(theta) * sinPhi;

            initial[zIndex(i)] =
                    FastMath.cos(phi);
        }

        return initial;
    }

    private static int xIndex(final int point) {
        return 3 * point;
    }

    private static int yIndex(final int point) {
        return 3 * point + 1;
    }

    private static int zIndex(final int point) {
        return 3 * point + 2;
    }

    private static int coordinateIndex(final int coordinate,
                                       final int point) {
        switch (coordinate) {
            case 0:
                return xIndex(point);
            case 1:
                return yIndex(point);
            default:
                return zIndex(point);
        }
    }

    @Test
    public void testELEC() {
        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();
        option.setMaxIteration(5000);
        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(initialPoint()),
                        new ObjectiveFunction(new Objective()),
                        new SphereEqualities(),
                        SimpleBounds.unbounded(N),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}
