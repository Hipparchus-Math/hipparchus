/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
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
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code PFIT1}.
 *
 * <p>The problem fits a model containing a pole by solving three nonlinear
 * equations in the three variables A, R and H. The SIF solution is
 * {@code (1, 3, 2)}.</p>
 *
 * <p>This is a square nonlinear-equation problem: the number of equality
 * constraints is equal to the number of variables.</p>
 */
public class PFIT1Test {

    /** Number of variables. */
    private static final int N = 3;

    /** Number of equality constraints. */
    private static final int M = 3;

    /** First target constant from the SIF. */
    private static final double CF = -8.0;

    /** Second target constant from the SIF. */
    private static final double CG = -18.66666666;

    /** Third target constant from the SIF. */
    private static final double CH = -23.11111111;

    /** Expected constant objective. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    /** Official starting point: A = 1, R = 0, H = 1. */
    private static final double[] START = {
        1.0, 0.0, 1.0
    };

    /** Lower bounds; only H is bounded below. */
    private static final double[] LOWER = {
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        -0.5
    };

    /** Upper bounds. */
    private static final double[] UPPER = {
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY
    };

    /**
     * Constant zero objective used by the nonlinear-equation formulation.
     */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            return 0.0;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            return new ArrayRealVector(N);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            return new OpenMapRealMatrix(N, N);
        }
    }

    /**
     * Three nonlinear equality constraints from the PFIT1 SIF.
     */
    private static final class Equalities
            extends EqualityConstraint {

        Equalities() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {
            return new ArrayRealVector(
                    residuals(point),
                    false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {

            return residualJacobian(point);
        }
    }

    /**
     * Evaluate the three PFIT1 equality residuals.
     *
     * @param point current point, ordered as A, R, H
     * @return residual array
     */
    private static double[] residuals(final RealVector point) {

        final double a = point.getEntry(0);
        final double r = point.getEntry(1);
        final double h = point.getEntry(2);

        final double y = 1.0 + h;

        final double power0 = FastMath.pow(y, -a);
        final double power1 = FastMath.pow(y, -(a + 1.0));
        final double power2 = FastMath.pow(y, -(a + 2.0));

        final double ea =
                a * (a + 1.0) * r * h * h;
        final double eb =
                a * r * h * (1.0 - power1);
        final double ec =
                a * r * h;
        final double ed =
                r * (1.0 - power0);
        final double ee =
                ea * (1.0 - power2);

        return new double[] {
            -0.5 * ea + ec - ed - CF,
            -ea + eb - CG,
            -ee - CH
        };
    }

    /**
     * Evaluate the equality-constraint Jacobian.
     *
     * @param point current point, ordered as A, R, H
     * @return constraint Jacobian
     */
    private static RealMatrix residualJacobian(final RealVector point) {

        final double a = point.getEntry(0);
        final double r = point.getEntry(1);
        final double h = point.getEntry(2);

        final double y = 1.0 + h;
        final double logY = FastMath.log(y);

        final double power0 = FastMath.pow(y, -a);
        final double power1 = FastMath.pow(y, -(a + 1.0));
        final double power2 = FastMath.pow(y, -(a + 2.0));

        final double b0 = 1.0 - power0;
        final double b1 = 1.0 - power1;
        final double b2 = 1.0 - power2;

        final double ea =
                a * (a + 1.0) * r * h * h;

        final double[] dEa = {
            (2.0 * a + 1.0) * r * h * h,
            a * (a + 1.0) * h * h,
            2.0 * a * (a + 1.0) * r * h
        };

        final double[] dEc = {
            r * h,
            a * h,
            a * r
        };

        final double[] dEd = {
            r * logY * power0,
            b0,
            r * a * power0 / y
        };

        final double[] dEb = {
            r * h * b1 +
                a * r * h * logY * power1,
            a * h * b1,
            a * r * b1 +
                a * r * h * (a + 1.0) * power1 / y
        };

        final double[] dEe = {
            dEa[0] * b2 +
                ea * logY * power2,
            dEa[1] * b2,
            dEa[2] * b2 +
                ea * (a + 2.0) * power2 / y
        };

        final RealMatrix jacobian =
                new OpenMapRealMatrix(M, N);

        for (int column = 0; column < N; ++column) {
            jacobian.setEntry(
                    0,
                    column,
                    -0.5 * dEa[column] +
                    dEc[column] -
                    dEd[column]);

            jacobian.setEntry(
                    1,
                    column,
                    -dEa[column] +
                    dEb[column]);

            jacobian.setEntry(
                    2,
                    column,
                    -dEe[column]);
        }

        return jacobian;
    }

    @Test
    public void testPFIT1() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(START),
                        new ObjectiveFunction(new Objective()),
                        new Equalities(),
                        new SimpleBounds(LOWER, UPPER),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}