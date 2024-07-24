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
package org.hipparchus.analysis.integration.gauss;

import org.hipparchus.util.Pair;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test for {@link AbstractRuleFactory}.
 *
 */
class RuleFactoryTest {
    /**
     * Tests that a given rule rule will be computed and added once to the cache
     * whatever the number of times this rule is called concurrently.
     */
    @Test
    void testConcurrentCreation() throws InterruptedException,
        ExecutionException {
        // Number of times the same rule will be called.
        final int numTasks = 20;

        final ThreadPoolExecutor exec
            = new ThreadPoolExecutor(3, numTasks, 1, TimeUnit.SECONDS,
                                     new ArrayBlockingQueue<Runnable>(2));

        final List<Future<Pair<double[], double[]>>> results
            = new ArrayList<Future<Pair<double[], double[]>>>();
        for (int i = 0; i < numTasks; i++) {
            results.add(exec.submit(new RuleBuilder()));
        }

        // Ensure that all computations have completed.
        for (Future<Pair<double[], double[]>> f : results) {
            f.get();
        }

        // Assertion would fail if "getRuleInternal" were not "synchronized".
        final int n = RuleBuilder.getNumberOfCalls();
        assertEquals(1, n, "Rule computation was called " + n + " times");
    }

    private static class RuleBuilder implements Callable<Pair<double[], double[]>> {
        private static final DummyRuleFactory factory = new DummyRuleFactory();

        public Pair<double[], double[]> call() {
            final int dummy = 2; // Always request the same rule.
            return factory.getRule(dummy);
        }

        public static int getNumberOfCalls() {
            return factory.getNumberOfCalls();
        }
    }

    private static class DummyRuleFactory extends AbstractRuleFactory {
        /** Rule computations counter. */
        private static AtomicInteger nCalls = new AtomicInteger();

        @Override
        protected Pair<double[], double[]> computeRule(int order) {
            // Tracks whether this computation has been called more than once.
            nCalls.getAndIncrement();

            assertDoesNotThrow(() -> {
                // Sleep to simulate computation time.
                Thread.sleep(20);
            }, "Unexpected interruption");

            // Dummy rule (but contents must exist).
            final double[] p = new double[order];
            final double[] w = new double[order];
            for (int i = 0; i < order; i++) {
                p[i] = Double.valueOf(i);
                w[i] = Double.valueOf(i);
            }
            return new Pair<>(p, w);
        }

        public int getNumberOfCalls() {
            return nCalls.get();
        }

    }

}

