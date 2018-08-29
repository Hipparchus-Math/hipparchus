/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
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

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test case for {@link RegulaFalsiSolver Regula Falsi} solver.
 *
 */
public final class RegulaFalsiSolverTest extends BaseSecantSolverAbstractTest {
    /** {@inheritDoc} */
    @Override
    protected UnivariateSolver getSolver() {
        UnivariateSolver solver = new RegulaFalsiSolver();
        checktype(solver, BaseSecantSolver.Method.REGULA_FALSI);
        return solver;
    }

    /** {@inheritDoc} */
    @Override
    protected int[] getQuinticEvalCounts() {
        // While the Regula Falsi method guarantees convergence, convergence
        // may be extremely slow. The last test case does not converge within
        // even a million iterations. As such, it was disabled.
        return new int[] {3, 7, 8, 19, 18, 11, 67, 55, 288, 151, -1};
    }

    @Test(expected=MathIllegalStateException.class)
    public void testIssue631() {
        final UnivariateFunction f = new UnivariateFunction() {
                /** {@inheritDoc} */
                public double value(double x) {
                    return FastMath.exp(x) - FastMath.pow(Math.PI, 3.0);
                }
            };

        final UnivariateSolver solver = new RegulaFalsiSolver();
        final double root = solver.solve(3624, f, 1, 10);
        Assert.assertEquals(3.4341896575482003, root, 1e-15);
    }
}
