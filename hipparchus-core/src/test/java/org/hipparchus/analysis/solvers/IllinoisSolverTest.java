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

import org.junit.Test;

/**
 * Test case for {@link IllinoisSolver Illinois} solver.
 *
 */
public final class IllinoisSolverTest extends BaseSecantSolverAbstractTest {
    /** {@inheritDoc} */
    @Override
    protected UnivariateSolver getSolver() {
        UnivariateSolver solver = new IllinoisSolver();
        checktype(solver, BaseSecantSolver.Method.ILLINOIS);
        return solver;
    }

    /** {@inheritDoc} */
    @Override
    protected int[] getQuinticEvalCounts() {
        return new int[] {3, 7, 9, 10, 10, 10, 12, 12, 14, 15, 20};
    }

    @Test
    public void testGitHubIssue44() {
        checktype(new IllinoisSolver(1.0e-6, 1.0e-14, 1.0e-15), BaseSecantSolver.Method.ILLINOIS);
    }

}
