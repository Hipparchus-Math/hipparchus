/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.analysis.solvers;

import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;

/**
 * Implements <a href="http://mathworld.wolfram.com/NewtonsMethod.html">
 * Newton's Method</a> for finding zeros of real univariate differentiable
 * functions.
 *
 */
public class NewtonRaphsonSolver extends AbstractUnivariateDifferentiableSolver {
    /** Default absolute accuracy. */
    private static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    /**
     * Construct a solver.
     */
    public NewtonRaphsonSolver() {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }
    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy Absolute accuracy.
     */
    public NewtonRaphsonSolver(double absoluteAccuracy) {
        super(absoluteAccuracy);
    }

    /**
     * Find a zero near the midpoint of {@code min} and {@code max}.
     *
     * @param f Function to solve.
     * @param min Lower bound for the interval.
     * @param max Upper bound for the interval.
     * @param maxEval Maximum number of evaluations.
     * @return the value where the function is zero.
     * @throws org.hipparchus.exception.MathIllegalStateException
     * if the maximum evaluation count is exceeded.
     * @throws org.hipparchus.exception.MathIllegalArgumentException
     * if {@code min >= max}.
     */
    @Override
    public double solve(int maxEval, final UnivariateDifferentiableFunction f,
                        final double min, final double max)
        throws MathIllegalStateException {
        return super.solve(maxEval, f, UnivariateSolverUtils.midpoint(min, max));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected double doSolve()
        throws MathIllegalStateException {
        final double startValue = getStartValue();
        final double absoluteAccuracy = getAbsoluteAccuracy();

        double x0 = startValue;
        double x1;
        while (true) {
            final DerivativeStructure y0 = computeObjectiveValueAndDerivative(x0);
            x1 = x0 - (y0.getValue() / y0.getPartialDerivative(1));
            if (FastMath.abs(x1 - x0) <= absoluteAccuracy) {
                return x1;
            }

            x0 = x1;
        }
    }
}
