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

import java.util.Arrays;

import org.hipparchus.linear.ArrayRealVector;
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
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code CATENARY}, original 166-beam instance.
 *
 * <p>This is the deliberately erroneous catenary model. The third squared
 * term of beam {@code i} is {@code (X(i) - Z(i - 1))^2}, exactly as written
 * in {@code CATENARY.SIF}. The corrected model is {@code CATENA}.</p>
 *
 * <p>This class keeps all 501 variables declared by the original SIF:
 * three coordinates for each of the 167 joints. Four variables are fixed
 * through bounds: {@code X(0)}, {@code Y(0)}, {@code Z(0)} and
 * {@code X(166)}.</p>
 *
 * <p>All derivatives are evaluated using forward finite differences.</p>
 */
@Disabled
public class CATENARYTest {

    /** Number of beams, equal to the SIF parameter N+1. */
    private static final int BEAMS = 166;

    /** Number of joints, including both endpoints. */
    private static final int JOINTS = BEAMS + 1;

    /** Three coordinates for every joint. */
    private static final int N = 3 * JOINTS;

    /** One equality constraint for every beam. */
    private static final int M = BEAMS;

    /** Gravitational acceleration. */
    private static final double GRAVITY = 9.81;

    /** Total chain mass. */
    private static final double TOTAL_MASS = 500.0;

    /** Beam length. */
    private static final double BEAM_LENGTH = 1.0;

    /** Horizontal shortening factor. */
    private static final double FRACTION = 0.6;

    /** Fixed final horizontal coordinate. */
    private static final double FINAL_X =
            BEAM_LENGTH * BEAMS * FRACTION;

    /** Mass times gravity for one beam. */
    private static final double MG =
            TOTAL_MASS * GRAVITY / BEAMS;

    /** Half endpoint coefficient in the trapezoidal objective. */
    private static final double HALF_MG =
            0.5 * MG;

    /** Objective reference printed in CATENARY.SIF for 166 beams. */
    private static final double EXPECTED_OBJECTIVE =
            -348403.164505;

    /** Linear gravitational-potential objective. */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {

            double objective =
                    HALF_MG * point.getEntry(yIndex(0));

            for (int joint = 1; joint < BEAMS; ++joint) {
                objective += MG * point.getEntry(yIndex(joint));
            }

            objective +=
                    HALF_MG * point.getEntry(yIndex(BEAMS));

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Gradient is evaluated by finite differences.");
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Hessian is not required by the finite-difference option.");
        }
    }

    /** Beam-length equalities copied literally from CATENARY.SIF. */
    private static final class BeamEqualities
            extends EqualityConstraint {

        BeamEqualities() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {

            final double[] constraints = new double[M];

            for (int beam = 1; beam <= BEAMS; ++beam) {

                final double dx =
                        point.getEntry(xIndex(beam)) -
                        point.getEntry(xIndex(beam - 1));

                final double dy =
                        point.getEntry(yIndex(beam)) -
                        point.getEntry(yIndex(beam - 1));

                /* Intentional CATENARY.SIF error. */
                final double erroneousDz =
                        point.getEntry(xIndex(beam)) -
                        point.getEntry(zIndex(beam - 1));

                constraints[beam - 1] =
                        dx * dx +
                        dy * dy +
                        erroneousDz * erroneousDz -
                        BEAM_LENGTH * BEAM_LENGTH;
            }

            return new ArrayRealVector(constraints, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Jacobian is evaluated by finite differences.");
        }
    }

    /** SIF starting point after enforcing the fixed-variable bounds. */
    private static double[] initialPoint() {

        final double[] initial = new double[N];
        final double tmp = FINAL_X / BEAMS;

        for (int joint = 1; joint <= BEAMS; ++joint) {
            initial[xIndex(joint)] = tmp;
        }

        initial[xIndex(BEAMS)] = FINAL_X;
        return initial;
    }

    /** Free default bounds plus the four fixed variables from CATENARY.SIF. */
    private static SimpleBounds bounds() {

        final double[] lower = new double[N];
        final double[] upper = new double[N];

        Arrays.fill(lower, Double.NEGATIVE_INFINITY);
        Arrays.fill(upper, Double.POSITIVE_INFINITY);

        fix(lower, upper, xIndex(0), 0.0);
        fix(lower, upper, yIndex(0), 0.0);
        fix(lower, upper, zIndex(0), 0.0);
        fix(lower, upper, xIndex(BEAMS), FINAL_X);

        return new SimpleBounds(lower, upper);
    }

    /** Fix one variable by assigning identical lower and upper bounds. */
    private static void fix(final double[] lower,
                            final double[] upper,
                            final int index,
                            final double value) {
        lower[index] = value;
        upper[index] = value;
    }

    private static int xIndex(final int joint) {
        return 3 * joint;
    }

    private static int yIndex(final int joint) {
        return 3 * joint + 1;
    }

    private static int zIndex(final int joint) {
        return 3 * joint + 2;
    }

    @Test
    public void testCATENARY() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(10000),
                        new InitialGuess(initialPoint()),
                        new ObjectiveFunction(new Objective()),
                        new BeamEqualities(),
                        bounds(),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}