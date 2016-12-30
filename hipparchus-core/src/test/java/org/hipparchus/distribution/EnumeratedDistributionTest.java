/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.Pair;
import org.junit.Test;

public class EnumeratedDistributionTest {

    @Test
    public void testCheckAndNormalizeBadArguments() {
        double[] bad = new double[] {-1, 0, 1, 1};
        try {
            EnumeratedDistribution.checkAndNormalize(bad);
            fail("Expecting IAE - negative probability");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        bad = new double[] {0, Double.NaN, 1, 1};
        try {
            EnumeratedDistribution.checkAndNormalize(bad);
            fail("Expecting IAE - NaN probability");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        bad = new double[] {0, Double.POSITIVE_INFINITY, 1, 1};
        try {
            EnumeratedDistribution.checkAndNormalize(bad);
            fail("Expecting IAE - infinite probability");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        bad = new double[] {0, 0, 0, 0};
        try {
            EnumeratedDistribution.checkAndNormalize(bad);
            fail("Expecting IAE - no positive probabilities");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        bad = new double[] {};
        try {
            EnumeratedDistribution.checkAndNormalize(bad);
            fail("Expecting IAE - empty probability array");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        bad = null;
        try {
            EnumeratedDistribution.checkAndNormalize(bad);
            fail("Expecting IAE - empty probability array");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testCheckAndNormalize() {
        double[] p = new double[] {0, 2, 2, 1};
        double[] normP = EnumeratedDistribution.checkAndNormalize(p);
        assertEquals(0, normP[0], 0);
        assertEquals(0.4, normP[1], 0);
        assertEquals(0.4, normP[2], 0);
        assertEquals(0.2, normP[3], 0);
        p = new double[] {0.2, 0.2, 0.4, 0.2};
        assertTrue(MathArrays.equals(p, EnumeratedDistribution.checkAndNormalize(p)));
    }

    @Test
    public void testNullValues() {
        final List<Pair<String, Double>> pmf = new ArrayList<>();
        pmf.add(new Pair<>("a", 0.5));
        pmf.add(new Pair<>(null, 0.5));
        final EnumeratedDistribution<String> dist = new EnumeratedDistribution<>(pmf);
        assertEquals(0.5, dist.probability(null), 0);
        assertEquals(0.5, dist.probability("a"), 0);
        assertEquals(0, dist.probability("b"), 0);
    }

    @Test
    public void testRepeatedValues() {
        final List<Pair<String, Double>> pmf = new ArrayList<>();
        pmf.add(new Pair<>("a", 0.5));
        pmf.add(new Pair<>("a", 0.5));
        pmf.add(new Pair<>("b", 0.0));
        final EnumeratedDistribution<String> dist = new EnumeratedDistribution<>(pmf);
        assertEquals(0, dist.probability(null), 0);
        assertEquals(1, dist.probability("a"), 0);
        assertEquals(0, dist.probability("b"), 0);
        assertEquals(0, dist.probability("c"), 0);
    }
}
