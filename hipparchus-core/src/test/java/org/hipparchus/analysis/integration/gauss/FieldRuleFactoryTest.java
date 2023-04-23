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

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.Pair;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for {@link FieldAbstractRuleFactory}.
 *
 */
public class FieldRuleFactoryTest {
    /**
     * Tests that a given rule rule will be computed and added once to the cache
     * whatever the number of times this rule is called concurrently.
     */
    @Test
        public void testConcurrentCreation() throws InterruptedException,
                                                    ExecutionException {
        // Number of times the same rule will be called.
        final int numTasks = 20;

        final ThreadPoolExecutor exec
            = new ThreadPoolExecutor(3, numTasks, 1, TimeUnit.SECONDS,
                                     new ArrayBlockingQueue<Runnable>(2));

        final List<Future<Pair<Binary64[], Binary64[]>>> results
            = new ArrayList<Future<Pair<Binary64[], Binary64[]>>>();
        for (int i = 0; i < numTasks; i++) {
            results.add(exec.submit(new RuleBuilder()));
        }

        // Ensure that all computations have completed.
        for (Future<Pair<Binary64[], Binary64[]>> f : results) {
            f.get();
        }

        // Assertion would fail if "getRuleInternal" were not "synchronized".
        final int n = RuleBuilder.getNumberOfCalls();
        Assert.assertEquals("Rule computation was called " + n + " times", 1, n);
    }

    private static class RuleBuilder implements Callable<Pair<Binary64[], Binary64[]>> {
        private static final DummyRuleFactory factory = new DummyRuleFactory();

        public Pair<Binary64[], Binary64[]> call() {
            final int dummy = 2; // Always request the same rule.
            return factory.getRule(dummy);
        }

        public static int getNumberOfCalls() {
            return factory.getNumberOfCalls();
        }
    }

    private static class DummyRuleFactory extends FieldAbstractRuleFactory<Binary64> {
        /** Rule computations counter. */
        private static AtomicInteger nCalls = new AtomicInteger();

        
        DummyRuleFactory() {
            super(Binary64Field.getInstance());
        }

        @Override
        protected Pair<Binary64[], Binary64[]> computeRule(int order) {
            // Tracks whether this computation has been called more than once.
            nCalls.getAndIncrement();

            try {
                // Sleep to simulate computation time.
                Thread.sleep(20);
            } catch (InterruptedException e) {
                Assert.fail("Unexpected interruption");
            }

            // Dummy rule (but contents must exist).
            final Binary64[] p = new Binary64[order];
            final Binary64[] w = new Binary64[order];
            for (int i = 0; i < order; i++) {
                p[i] = new Binary64(i);
                w[i] = new Binary64(i);
            }
            return new Pair<>(p, w);
        }

        public int getNumberOfCalls() {
            return nCalls.get();
        }

    }

}
