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
package org.hipparchus.optim;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class SimplePointCheckerTest {
    @Test
    void testIterationCheckPrecondition() {
        assertThrows(MathIllegalArgumentException.class, () -> new SimplePointChecker<PointValuePair>(1e-1, 1e-2, 0));
    }

    @Test
    void testIterationCheck() {
        final int max = 10;
        final SimplePointChecker<PointValuePair> checker
            = new SimplePointChecker<PointValuePair>(1e-1, 1e-2, max);
        assertTrue(checker.converged(max, null, null));
        assertTrue(checker.converged(max + 1, null, null));
    }

    @Test
    void testIterationCheckDisabled() {
        final SimplePointChecker<PointValuePair> checker
            = new SimplePointChecker<PointValuePair>(1e-8, 1e-8);

        final PointValuePair a = new PointValuePair(new double[] { 1d }, 1d);
        final PointValuePair b = new PointValuePair(new double[] { 10d }, 10d);

        assertFalse(checker.converged(-1, a, b));
        assertFalse(checker.converged(0, a, b));
        assertFalse(checker.converged(1000000, a, b));

        assertTrue(checker.converged(-1, a, a));
        assertTrue(checker.converged(-1, b, b));
    }
}
