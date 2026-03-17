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

import static org.junit.jupiter.api.Assertions.assertEquals;

import static org.junit.jupiter.api.Assertions.assertTrue;

import org.hipparchus.util.FastMath;

final class HSProblemTestUtils {

    private static final boolean DEBUG = Boolean.getBoolean("hipparchus.debug.sqp");

    private HSProblemTestUtils() {
        // utility class
    }

    static SQPOptimizerS2 newOptimizer() {
        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (DEBUG) {
            optimizer.setDebugPrinter(System.out::println);
        }
        return optimizer;
    }

    static void assertExpectedObjective(final double fExpected, final LagrangeSolution sol) {
        final double tol = 1.0e-2 * (FastMath.abs(fExpected) + 1.0);
        assertEquals(fExpected, sol.getValue(), tol);
    }
    
    static void assertBetterObjective(final double fExpected, final LagrangeSolution sol) {
       final double tol = 1.0e-4 * (FastMath.abs(fExpected) + 1.0); 
       assertTrue(sol.getValue() < fExpected+tol,"equal or better Solution");
    }

}
