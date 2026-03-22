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
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS395TestWithSQPProblemInterface {

    /*
     * SQPProblem implementation for Hock–Schittkowski problem 395.
     *
     * Objective:
     *   f(x) = Σ_{i=1}^{50} i * (x_i^2 + x_i^4)
     *
     * Equality constraint:
     *   h(x) = Σ_{i=1}^{50} x_i^2 - 1 = 0
     *
     * No inequality constraints, no box bounds.
     */
    private static class HS395Problem implements SQPProblem, OptimizationData {

        /** Dimension of the decision vector. */
        private static final int N = 50;

        @Override
        public int getdim() {
            return N;
        }

        // -------- Initial guess / bounds flags --------

        @Override
        public boolean hasInitialGuess() {
            // In this test we provide the initial guess explicitly via InitialGuess,
            // so we return false here.
            return false;
        }

        @Override
        public double[] getInitialGuess() {
            // Not used because hasInitialGuess() == false.
            // Could be implemented if the optimizer prefers retrieving
            // the guess from the problem itself.
            double[] start = new double[N];
            for (int i = 0; i < N; i++) {
                start[i] = 2.0;
            }
            return start;
        }

        @Override
        public boolean hasBounds() {
            return false;
        }

        @Override
        public double[] getBoxConstraintLB() {
            return null;
        }

        @Override
        public double[] getBoxConstraintUB() {
            return null;
        }

        // -------- Objective function --------

        @Override
        public double getObjectiveEvaluation(final RealVector x) {
            double sum = 0.0;
            for (int i = 0; i < x.getDimension(); i++) {
                final double xi  = x.getEntry(i);
                final double idx = i + 1; // 1-based index
                sum += idx * (xi * xi + FastMath.pow(xi, 4));
            }
            return sum;
        }

        @Override
        public RealVector getObjectiveGradient(final RealVector x) {
            // As in the original test, no analytical gradient is provided.
            // The optimizer is expected to approximate gradients numerically
            // if needed.
            throw new UnsupportedOperationException("Analytical gradient not provided for HS395.");
        }

        // -------- Equality constraints --------

        @Override
        public boolean hasEquality() {
            return true;
        }

        @Override
        public RealVector getEqCostraintEvaluation(final RealVector x) {
            double sum = 0.0;
            for (int i = 0; i < x.getDimension(); i++) {
                sum += FastMath.pow(x.getEntry(i), 2);
            }
            // h(x) = Σ x_i^2 - 1 = 0
            return new ArrayRealVector(new double[] { sum - 1.0 });
        }

        @Override
        public RealMatrix getEqCostraintJacobian(final RealVector x) {
            // As in the original HS395Eq, we do not provide the Jacobian.
            // The optimizer may approximate it numerically.
            throw new UnsupportedOperationException("Analytical Jacobian for equality constraints not provided.");
        }

        @Override
        public RealVector getEqCostraintLB() {
            // Equality constraints are modeled as h(x) = 0; lower bound is 0.
            return new ArrayRealVector(new double[] { 0.0 });
        }

        // -------- Inequality constraints (none for HS395) --------

        @Override
        public boolean hasInequality() {
            return false;
        }

        @Override
        public RealVector getIneqConstraintEvaluation(final RealVector x) {
            return null;
        }

        @Override
        public RealMatrix getIneqCostraintJacobian(final RealVector x) {
            return null;
        }

        @Override
        public RealVector getIneqCostraintLB() {
            return null;
        }
    }

    @Test
    public void testHS395() {
        // Explicit initial guess (as in the original test).
        double[] start = new double[50];
        for (int i = 0; i < 50; i++) {
            start[i] = 2.0;
        }

        InitialGuess guess = new InitialGuess(start);
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }

        // Expected optimal objective value.
        double expected = 1.9166668;

        // Use the new SQPProblem-based interface.
        HS395Problem problem = new HS395Problem();

        LagrangeSolution sol = optimizer.optimize(
                problem

        );

        assertEquals(expected, sol.getValue(), 1e-6);
    }
}
