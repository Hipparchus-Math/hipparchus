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
package org.hipparchus.stat.descriptive;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.stat.descriptive.moment.SecondMoment;
import org.hipparchus.util.FastMath;
import org.junit.Test;

/**
 * Test cases for {@link StorelessUnivariateStatistic} classes.
 */
public abstract class StorelessUnivariateStatisticAbstractTest
    extends UnivariateStatisticAbstractTest {

    /** Small sample arrays */
    protected double[][] smallSamples = {{}, {1}, {1,2}, {1,2,3}, {1,2,3,4}};

    /** Return a new instance of the statistic */
    @Override
    public abstract StorelessUnivariateStatistic getUnivariateStatistic();

    /**Expected value for the testArray defined in UnivariateStatisticAbstractTest */
    @Override
    public abstract double expectedValue();

    /**
     * Verifies that increment() and incrementAll work properly.
     */
    @Test
    public void testIncrementation() {

        StorelessUnivariateStatistic statistic = getUnivariateStatistic();

        // Add testArray one value at a time and check result
        for (int i = 0; i < testArray.length; i++) {
            statistic.increment(testArray[i]);
        }

        assertEquals(expectedValue(), statistic.getResult(), getTolerance());
        assertEquals(testArray.length, statistic.getN());

        statistic.clear();

        // Add testArray all at once and check again
        statistic.incrementAll(testArray);
        assertEquals(expectedValue(), statistic.getResult(), getTolerance());
        assertEquals(testArray.length, statistic.getN());

        statistic.clear();

        // Cleared
        checkClearValue(statistic);
        assertEquals(0, statistic.getN());
    }

    protected void checkClearValue(StorelessUnivariateStatistic statistic){
        assertTrue(Double.isNaN(statistic.getResult()));
    }

    @Test
    public void testSerialization() {
        StorelessUnivariateStatistic statistic = getUnivariateStatistic();

        UnitTestUtils.checkSerializedEquality(statistic);
        statistic.clear();

        for (int i = 0; i < testArray.length; i++) {
            statistic.increment(testArray[i]);
            if(i % 5 == 0) {
                statistic = (StorelessUnivariateStatistic)UnitTestUtils.serializeAndRecover(statistic);
            }
        }

        UnitTestUtils.checkSerializedEquality(statistic);
        assertEquals(expectedValue(), statistic.getResult(), getTolerance());

        statistic.clear();
        checkClearValue(statistic);
    }

    @Test
    public void testEqualsAndHashCode() {
        StorelessUnivariateStatistic statistic = getUnivariateStatistic();
        StorelessUnivariateStatistic statistic2 = null;

        assertTrue("non-null, compared to null", !statistic.equals(statistic2));
        assertTrue("reflexive, non-null", statistic.equals(statistic));

        int emptyHash = statistic.hashCode();
        statistic2 = getUnivariateStatistic();
        assertTrue("empty stats should be equal", statistic.equals(statistic2));
        assertEquals("empty stats should have the same hashcode",
                     emptyHash, statistic2.hashCode());

        statistic.increment(1d);
        assertTrue("reflexive, non-empty", statistic.equals(statistic));
        assertTrue("non-empty, compared to empty", !statistic.equals(statistic2));
        assertTrue("non-empty, compared to empty", !statistic2.equals(statistic));
        assertTrue("non-empty stat should have different hashcode from empty stat",
                   statistic.hashCode() != emptyHash);

        statistic2.increment(1d);
        assertTrue("stats with same data should be equal", statistic.equals(statistic2));
        assertEquals("stats with same data should have the same hashcode",
                     statistic.hashCode(), statistic2.hashCode());

        statistic.increment(Double.POSITIVE_INFINITY);
        assertTrue("stats with different n's should not be equal", !statistic2.equals(statistic));
        assertTrue("stats with different n's should have different hashcodes",
                   statistic.hashCode() != statistic2.hashCode());

        statistic2.increment(Double.POSITIVE_INFINITY);
        assertTrue("stats with same data should be equal", statistic.equals(statistic2));
        assertEquals("stats with same data should have the same hashcode",
                     statistic.hashCode(), statistic2.hashCode());

        statistic.clear();
        statistic2.clear();
        assertTrue("cleared stats should be equal", statistic.equals(statistic2));
        assertEquals("cleared stats should have thashcode of empty stat",
                     emptyHash, statistic2.hashCode());
        assertEquals("cleared stats should have thashcode of empty stat",
                     emptyHash, statistic.hashCode());

    }

    @Test
    public void testMomentSmallSamples() {
        UnivariateStatistic stat = getUnivariateStatistic();
        if (stat instanceof SecondMoment) {
            SecondMoment moment = (SecondMoment) getUnivariateStatistic();
            assertTrue(Double.isNaN(moment.getResult()));
            moment.increment(1d);
            assertEquals(0d, moment.getResult(), 0);
        }
    }

    /**
     * Make sure that evaluate(double[]) and inrementAll(double[]),
     * getResult() give same results.
     */
    @Test
    public void testConsistency() {
        StorelessUnivariateStatistic stat = getUnivariateStatistic();
        stat.incrementAll(testArray);
        assertEquals(stat.getResult(), stat.evaluate(testArray), getTolerance());
        for (int i = 0; i < smallSamples.length; i++) {
            stat.clear();
            for (int j =0; j < smallSamples[i].length; j++) {
                stat.increment(smallSamples[i][j]);
            }
            UnitTestUtils.assertEquals(stat.getResult(), stat.evaluate(smallSamples[i]), getTolerance());
        }
    }

    /**
     * Verifies that copied statistics remain equal to originals when
     * incremented the same way.
     */
    @Test
    public void testCopyConsistency() {

        StorelessUnivariateStatistic master = getUnivariateStatistic();

        StorelessUnivariateStatistic replica = null;

        // Randomly select a portion of testArray to load first
        long index = FastMath.round((FastMath.random()) * testArray.length);

        // Put first half in master and copy master to replica
        master.incrementAll(testArray, 0, (int) index);
        replica = master.copy();

        // Check same
        assertTrue(replica.equals(master));
        assertTrue(master.equals(replica));

        // Now add second part to both and check again
        master.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        replica.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        assertTrue(replica.equals(master));
        assertTrue(master.equals(replica));
    }

    @Test
    public void testSerial() {
        StorelessUnivariateStatistic s = getUnivariateStatistic();
        assertEquals(s, UnitTestUtils.serializeAndRecover(s));
    }

    /**
     * Make sure that evaluate(double[]) does not alter the internal state.
     */
    @Test
    public void testEvaluateInternalState() {
        StorelessUnivariateStatistic stat = getUnivariateStatistic();
        stat.evaluate(testArray);
        assertEquals(0, stat.getN());

        stat.incrementAll(testArray);

        StorelessUnivariateStatistic savedStatistic = stat.copy();

        assertNotEquals(stat.getResult(), stat.evaluate(testArray, 0, 5), getTolerance());

        assertEquals(savedStatistic.getResult(), stat.getResult(), 0.0);
        assertEquals(savedStatistic.getN(), stat.getN());
    }

    /**
     * Test that the aggregate operation is consistent with individual increment.
     */
    @Test
    @SuppressWarnings("unchecked")
    public <T> void testAggregate() {
        // Union of both statistics.
        StorelessUnivariateStatistic statU = getUnivariateStatistic();
        if (!(statU instanceof AggregatableStatistic<?>)) {
            return;
        }

        // Aggregated statistic.
        AggregatableStatistic<T> aggregated = (AggregatableStatistic<T>) getUnivariateStatistic();
        StorelessUnivariateStatistic statAgg = (StorelessUnivariateStatistic) aggregated;

        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(100);
        // Create Set A
        StorelessUnivariateStatistic statA = getUnivariateStatistic();
        for (int i = 0; i < 10; i++) {
            final double val = randomDataGenerator.nextGaussian();
            statA.increment(val);
            statU.increment(val);
        }

        aggregated.aggregate((T) statA);
        assertEquals(statA.getN(), statAgg.getN(), getTolerance());
        assertEquals(statA.getResult(), statAgg.getResult(), getTolerance());

        // Create Set B
        StorelessUnivariateStatistic statB = getUnivariateStatistic();
        for (int i = 0; i < 4; i++) {
            final double val = randomDataGenerator.nextGaussian();
            statB.increment(val);
            statU.increment(val);
        }

        aggregated.aggregate((T) statB);
        assertEquals(statU.getN(), statAgg.getN(), getTolerance());
        assertEquals(statU.getResult(), statAgg.getResult(), getTolerance());
    }

    @Test
    public void testConsume() {
        StorelessUnivariateStatistic stat = getUnivariateStatistic();

        Arrays.stream(testArray)
              .forEach(stat);

        assertEquals(expectedValue(), stat.getResult(), getTolerance());

        StorelessUnivariateStatistic stat2 = getUnivariateStatistic();
        stat2.incrementAll(testArray);
        assertEquals(stat2.getResult(), stat.getResult(), getTolerance());
    }
}
