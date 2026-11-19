/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class SQPOptimizerS2Test {

    /** Tolerance on the solution of the small convex test problems. */
    private static final double TOL = 1.0e-6;

    @Test
    public void testSQPProblemWithAllConstraints() {
        // min x² + y² + z²  s.t.  x + y = 1, x + z ≥ 2, 0 ≤ x, y, z: solution (1, 0, 1)
        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        final BoundedConstraint box = new LinearBoundedConstraint(MatrixUtils.createRealIdentityMatrix(3),
                                                                  new ArrayRealVector(3, 0.0),
                                                                  new ArrayRealVector(3, 1.0));
        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        final LagrangeSolution sol = optimizer.optimize(option, new DummyProblem(true, true, true), box);

        Assertions.assertSame(box, optimizer.getBoxConstraint());
        Assertions.assertNotNull(optimizer.getEqConstraint());
        Assertions.assertNotNull(optimizer.getIqConstraint());
        Assertions.assertNotNull(optimizer.getSimpleBounds());
        Assertions.assertArrayEquals(new double[] { 1.0, 0.0, 1.0 }, sol.getX().toArray(), TOL);
        Assertions.assertEquals(2.0, sol.getValue(), TOL);
    }

    @Test
    public void testSQPProblemWithoutConstraints() {
        // min x² + y² + z², unbounded: the initial guess is the origin, which is the solution
        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        final LagrangeSolution sol = optimizer.optimize(option, new DummyProblem(false, false, false));

        Assertions.assertNull(optimizer.getEqConstraint());
        Assertions.assertNull(optimizer.getIqConstraint());
        Assertions.assertNull(optimizer.getSimpleBounds());
        Assertions.assertArrayEquals(new double[] { 0.0, 0.0, 0.0 }, sol.getX().toArray(), TOL);
        Assertions.assertEquals(0.0, sol.getValue(), TOL);
    }

    @Test
    public void testInitialGuessFromPartialBounds() {
        // no start point: the initial guess is built from bounds that are
        // two-sided, upper only, lower only and absent; the minimum is inside the bounds
        // (the default eps = 1e-7 leaves x about 4e-5 away from the minimum, so it is tightened)
        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        option.setEps(1.0e-12);
        final SimpleBounds bounds = new SimpleBounds(new double[] { 0.0,  Double.NEGATIVE_INFINITY, -5.0, Double.NEGATIVE_INFINITY },
                                                     new double[] { 10.0, 5.0, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY });
        final LagrangeSolution sol = new SQPOptimizerS2().optimize(option, bounds,
                                                                   new ObjectiveFunction(shiftedQuadratic(1.0, 2.0, 3.0, 4.0)));

        Assertions.assertArrayEquals(new double[] { 1.0, 2.0, 3.0, 4.0 }, sol.getX().toArray(), TOL);
        Assertions.assertEquals(0.0, sol.getValue(), TOL);
    }

    @Test
    public void testFixedVariableForwardGradient() {
        doTestFixedVariable(GradientMode.FORWARD);
    }

    @Test
    public void testFixedVariableCentralGradient() {
        doTestFixedVariable(GradientMode.CENTRAL);
    }

    private void doTestFixedVariable(final GradientMode gradientMode) {
        // min (x₀-1)² + (x₁-3)² + (x₂-1)²  s.t.  x₁ + x₂ = 5, x₂ ≥ 0, x₀ = 2 (fixed by bounds), x₁ ≥ 3.8
        // solution (2, 3.8, 1.2), objective 1 + 0.64 + 0.04
        // x₀ cannot be perturbed at all, and x₁ starts on its lower bound so it can only be
        // perturbed upwards
        final SQPOption option = new SQPOption();
        option.setGradientMode(gradientMode);
        final SimpleBounds bounds = new SimpleBounds(new double[] { 2.0,  3.8, -10.0 },
                                                     new double[] { 2.0, 10.0,  10.0 });
        final LinearEqualityConstraint eq =
                new LinearEqualityConstraint(new double[][] { { 0.0, 1.0, 1.0 } }, new double[] { 5.0 });
        final LinearInequalityConstraint iq =
                new LinearInequalityConstraint(new double[][] { { 0.0, 0.0, 1.0 } }, new double[] { 0.0 });
        final LagrangeSolution sol = new SQPOptimizerS2().optimize(option, bounds, eq, iq,
                                                                   new InitialGuess(new double[] { 2.0, 3.8, 0.5 }),
                                                                   new ObjectiveFunction(shiftedQuadratic(1.0, 3.0, 1.0)));

        Assertions.assertArrayEquals(new double[] { 2.0, 3.8, 1.2 }, sol.getX().toArray(), TOL);
        Assertions.assertEquals(1.68, sol.getValue(), TOL);
    }

    @Test
    public void testCustomConvergenceCriterion() {
        final int[] calls = new int[1];
        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        option.setConvergenceFunction((alpha, dxNorm, dxHdx, complSlack, kkt, viol, funDiff, xNorm, funNorm) -> {
            ++calls[0];
            return true;
        });
        new SQPOptimizerS2().optimize(option,
                                      new InitialGuess(new double[] { 5.0, -5.0 }),
                                      new ObjectiveFunction(shiftedQuadratic(1.0, 2.0)));

        // the optimizer stopped at the first convergence check
        Assertions.assertEquals(1, calls[0]);
    }

    /** f(x) = Σ (xᵢ - cᵢ)² = ½ xᵀ (2I) x - 2cᵀx + cᵀc. */
    private static QuadraticFunction shiftedQuadratic(final double... center) {
        final RealVector c = new ArrayRealVector(center);
        return new QuadraticFunction(MatrixUtils.createRealIdentityMatrix(center.length).scalarMultiply(2.0),
                                     c.mapMultiply(-2.0), c.dotProduct(c));
    }

}
