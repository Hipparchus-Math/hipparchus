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
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code CRESC100}.
 *
 * <p>The problem finds a minimum-area crescent containing 100 prescribed
 * planar points. The SIF parametrization uses two circle centers and radii
 * constructed from six variables so that the crescent area remains
 * computable throughout the bounded domain.</p>
 *
 * <p>The problem has six variables, 200 nonlinear inequality constraints and
 * no equality constraints. Objective gradients, constraint Jacobians and the
 * Hessian are deliberately evaluated through the finite-difference/BFGS
 * configuration used by the SQPOptimizerS2 CUTEst test suite.</p>
 */
public class CRESC100Test {

    /** Number of variables. */
    private static final int N = 6;

    /** Number of data points. */
    private static final int POINTS = 100;

    /** Two inequalities are generated for every point. */
    private static final int M = 2 * POINTS;

    /** Reference objective reported for CRESC100. */
    private static final double EXPECTED_OBJECTIVE =
            0.5676027;

    /** Variable indices in the original SIF order. */
    private static final int V1 = 0;
    private static final int W1 = 1;
    private static final int D  = 2;
    private static final int A  = 3;
    private static final int T  = 4;
    private static final int R  = 5;

    /** Official SIF starting point. */
    private static final double[] START = {
        -40.0,
          5.0,
          1.0,
          2.0,
          1.5,
          0.75
    };

    /** X coordinates of the 100 points. */
    private static final double[] X = {
        0.544, 0.714, 0.594, 0.474, 0.470, 0.241, 0.503, 0.854,
        0.438, 0.294, 0.479, 0.413, 0.722, 0.358, 0.836, 0.648,
        0.267, 0.362, 0.232, 0.667, 0.476, 0.868, 0.603, 0.788,
        0.745, 0.627, 0.394, 0.220, 0.548, 0.446, 0.463, 0.541,
        0.631, 0.735, 0.674, 0.620, 0.257, 0.659, 0.475, 0.708,
        0.545, 0.586, 0.796, 0.525, 0.498, 0.790, 0.632, 0.699,
        0.758, 0.656, 0.652, 0.801, 0.730, 0.689, 0.527, 0.696,
        0.911, 0.805, 0.764, 0.331, 0.078, 0.656, 0.480, 0.503,
        0.412, 0.338, 0.920, 0.548, 0.826, 0.071, 0.635, 0.591,
        0.489, 0.565, 0.791, 0.725, 0.091, 0.401, 0.226, 0.518,
        0.651, 0.738, 0.509, 0.833, 0.669, 0.121, 0.809, 0.477,
        0.267, 0.412, 0.663, 0.830, 0.324, 0.225, 0.870, 0.343,
        0.849, 0.655, 0.558, 0.830
    };

    /** Y coordinates of the 100 points. */
    private static final double[] Y = {
        0.492, 0.505, 0.094, 0.398, 0.804, 0.955, 0.097, 0.359,
        0.477, 0.887, 0.372, 0.911, 0.082, 0.788, 0.225, 0.424,
        0.932, 0.697, 0.969, 0.527, 0.735, 0.141, 0.550, 0.399,
        0.365, 0.762, 0.779, 0.869, 0.625, 0.778, 0.762, 0.068,
        0.670, 0.151, 0.490, 0.565, 0.862, 0.010, 0.585, 0.475,
        0.371, 0.085, 0.309, 0.573, 0.530, 0.217, 0.169, 0.024,
        0.436, 0.662, 0.251, 0.511, 0.611, 0.277, 0.612, 0.074,
        0.052, 0.075, 0.067, 0.759, 0.983, 0.340, 0.232, 0.140,
        0.683, 0.624, 0.139, 0.474, 0.443, 0.969, 0.728, 0.570,
        0.427, 0.105, 0.554, 0.174, 0.951, 0.790, 0.968, 0.663,
        0.166, 0.146, 0.280, 0.176, 0.044, 0.921, 0.458, 0.436,
        0.840, 0.630, 0.648, 0.086, 0.862, 0.853, 0.331, 0.629,
        0.521, 0.714, 0.005, 0.526
    };

    /**
     * Crescent-area objective translated directly from the SIF {@code SC}
     * nonlinear element.
     */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {

            final double d =
                    point.getEntry(D);

            final double a =
                    point.getEntry(A);

            final double r =
                    point.getEntry(R);

            /*
             * These are exactly the SIF temporaries A, B and D used by
             * the SC element. The temporary names are expanded here to
             * avoid confusion with the optimization variable A.
             */
            final double firstRadiusTerm =
                    d + r;

            final double secondRadiusTerm =
                    a * d + r;

            final double centerDistance =
                    a * d;

            final double firstRadiusSquared =
                    firstRadiusTerm *
                    firstRadiusTerm;

            final double secondRadiusSquared =
                    secondRadiusTerm *
                    secondRadiusTerm;

            final double distanceSquared =
                    centerDistance *
                    centerDistance;

            final double firstCosine =
                    (distanceSquared +
                     secondRadiusSquared -
                     firstRadiusSquared) /
                    (2.0 *
                     secondRadiusTerm *
                     centerDistance);

            final double secondCosine =
                    -(distanceSquared -
                      secondRadiusSquared +
                      firstRadiusSquared) /
                    (2.0 *
                     firstRadiusTerm *
                     centerDistance);

            final double firstAngle =
                    FastMath.acos(firstCosine);

            final double secondAngle =
                    FastMath.acos(secondCosine);

            return firstRadiusSquared *
                   secondAngle -
                   secondRadiusSquared *
                   firstAngle +
                   firstRadiusTerm *
                   centerDistance *
                   FastMath.sin(secondAngle);
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

    /**
     * Point-containment inequalities.
     *
     * <p>The first 100 components impose that every point lies inside the
     * second circle. The final 100 components impose that every point lies
     * outside the first circle. Hipparchus inequalities are returned in
     * nonnegative form.</p>
     */
    private static final class CrescentInequalities
            extends InequalityConstraint {

        CrescentInequalities() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {

            final double v1 =
                    point.getEntry(V1);

            final double w1 =
                    point.getEntry(W1);

            final double d =
                    point.getEntry(D);

            final double a =
                    point.getEntry(A);

            final double t =
                    point.getEntry(T);

            final double r =
                    point.getEntry(R);

            final double ad =
                    a * d;

            final double cosT =
                    FastMath.cos(t);

            final double sinT =
                    FastMath.sin(t);

            final double secondCenterX =
                    v1 +
                    ad * cosT;

            final double secondCenterY =
                    w1 +
                    ad * sinT;

            final double secondRadius =
                    d + r;

            final double firstRadius =
                    ad + r;

            final double secondRadiusSquared =
                    secondRadius *
                    secondRadius;

            final double firstRadiusSquared =
                    firstRadius *
                    firstRadius;

            final double[] constraints =
                    new double[M];

            for (int i = 0;
                 i < POINTS;
                 ++i) {

                final double secondDx =
                        secondCenterX -
                        X[i];

                final double secondDy =
                        secondCenterY -
                        Y[i];

                /*
                 * SIF group IS2(i):
                 * distance-to-circle-2 squared minus radius squared <= 0.
                 * It is negated here to obtain the Hipparchus >= 0 form.
                 */
                constraints[i] =
                        secondRadiusSquared -
                        secondDx * secondDx -
                        secondDy * secondDy;

                final double firstDx =
                        v1 -
                        X[i];

                final double firstDy =
                        w1 -
                        Y[i];

                /*
                 * SIF group OS1(i):
                 * distance-to-circle-1 squared minus radius squared >= 0.
                 */
                constraints[POINTS + i] =
                        firstDx * firstDx +
                        firstDy * firstDy -
                        firstRadiusSquared;
            }

            return new ArrayRealVector(
                    constraints,
                    false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Jacobian is evaluated by finite differences.");
        }
    }

    /**
     * Original SIF variable bounds.
     *
     * <p>Although the descriptive comment states {@code d <= 1}, the actual
     * BOUNDS section specifies only {@code D >= 1e-8}. This test follows the
     * executable SIF definition exactly.</p>
     */
    private static SimpleBounds bounds() {

        final double[] lower =
                new double[N];

        final double[] upper =
                new double[N];

        Arrays.fill(
                lower,
                Double.NEGATIVE_INFINITY);

        Arrays.fill(
                upper,
                Double.POSITIVE_INFINITY);

        lower[D] =
                1.0e-8;

        lower[A] =
                1.0;

        lower[T] =
                0.0;

        upper[T] =
                6.2831852;

        lower[R] =
                0.39;

        return new SimpleBounds(
                lower,
                upper);
    }

    @Test
    public void testCRESC100() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(10000),
                        new InitialGuess(START.clone()),
                        new ObjectiveFunction(new Objective()),
                        new CrescentInequalities(),
                        bounds(),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}