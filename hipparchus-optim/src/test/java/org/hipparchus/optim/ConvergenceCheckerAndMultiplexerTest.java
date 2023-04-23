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
package org.hipparchus.optim;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class ConvergenceCheckerAndMultiplexerTest {

    @Test
    public void testFalseFalse() {
        Assert.assertFalse(new ConvergenceCheckerAndMultiplexer<>(buildCheckers(false, false)).converged(0, null, null));
    }

    @Test
    public void testFalseTrue() {
        Assert.assertFalse(new ConvergenceCheckerAndMultiplexer<>(buildCheckers(false, true)).converged(0, null, null));
    }

    @Test
    public void testTrueFalse() {
        Assert.assertFalse(new ConvergenceCheckerAndMultiplexer<>(buildCheckers(true, false)).converged(0, null, null));
    }

    @Test
    public void testTrueTrue() {
        Assert.assertTrue(new ConvergenceCheckerAndMultiplexer<>(buildCheckers(true, true)).converged(0, null, null));
    }

    private List<ConvergenceChecker<Object>> buildCheckers(final boolean result1, final boolean result2) {
        final List<ConvergenceChecker<Object>> checkers = new ArrayList<>();
        checkers.add((iteration, previous, current) -> result1);
        checkers.add((iteration, previous, current) -> result2);
        return checkers;
    }

}
