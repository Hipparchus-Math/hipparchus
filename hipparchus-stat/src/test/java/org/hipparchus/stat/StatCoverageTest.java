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
package org.hipparchus.stat;

import org.hipparchus.stat.descriptive.moment.SemiVariance;
import org.junit.jupiter.api.Test;
import java.util.Arrays;
import java.util.List;
import static org.junit.jupiter.api.Assertions.*;

public class StatCoverageTest {

    @Test
    public void testLongFrequencyCoverage() {
        LongFrequency freq = new LongFrequency();
        freq.addValue(1);
        freq.incrementValue(2, 3);
        
        assertEquals(1, freq.getCount(1));
        assertEquals(3, freq.getCount(2));
        assertEquals(4, freq.getSumFreq());
        
        assertEquals(0.25, freq.getPct(1));
        assertEquals(0.25, freq.getCumPct(1));
        assertEquals(1, freq.getCumFreq(1));
        
        assertEquals(1.0, freq.getCumPct(2));
        assertEquals(4, freq.getCumFreq(2));
    }

    @Test
    public void testFrequencyEqualsHashCode() {
        Frequency<String> f1 = new Frequency<>();
        Frequency<String> f2 = new Frequency<>();
        
        assertEquals(f1, f2);
        assertEquals(f1.hashCode(), f2.hashCode());
        
        f1.addValue("A");
        assertNotEquals(f1, f2);
        
        f2.addValue("A");
        assertEquals(f1, f2);
        
        f1.addValue("B");
        assertNotEquals(f1, f2);
        
        assertNotEquals(f1, "not a frequency");
    }

    @Test
    public void testFrequencyMerge() {
        Frequency<String> f1 = new Frequency<>();
        f1.addValue("A");
        Frequency<String> f2 = new Frequency<>();
        f2.addValue("B");
        
        f1.merge(f2);
        assertEquals(1, f1.getCount("A"));
        assertEquals(1, f1.getCount("B"));
        
        Frequency<String> f3 = new Frequency<>();
        f3.addValue("C");
        f1.merge(Arrays.asList(f3));
        assertEquals(1, f1.getCount("C"));
    }

    @Test
    public void testStatUtilsDifference() {
        double[] s1 = { 1.0, 2.0, 3.0 };
        double[] s2 = { 0.5, 1.5, 2.5 };
        
        assertEquals(1.5, StatUtils.sumDifference(s1, s2), 1e-10);
        assertEquals(0.5, StatUtils.meanDifference(s1, s2), 1e-10);
        assertEquals(0.0, StatUtils.varianceDifference(s1, s2, 0.5), 1e-10);
    }
    
    @Test
    public void testStatUtilsMode() {
        double[] sample = { 1.0, 2.0, 2.0, 3.0, 3.0 };
        double[] mode = StatUtils.mode(sample);
        assertEquals(2, mode.length);
        assertEquals(2.0, mode[0]);
        assertEquals(3.0, mode[1]);
        
        assertEquals(0, StatUtils.mode(new double[0]).length);
    }

    @Test
    public void testSemiVarianceCoverage() {
        SemiVariance sv = new SemiVariance();
        SemiVariance copy = sv.copy();
        assertEquals(sv.isBiasCorrected(), copy.isBiasCorrected());
        assertEquals(sv.getVarianceDirection(), copy.getVarianceDirection());
        
        SemiVariance sv2 = new SemiVariance(SemiVariance.Direction.UPSIDE);
        assertEquals(SemiVariance.Direction.UPSIDE, sv2.getVarianceDirection());
        
        SemiVariance sv3 = new SemiVariance(false, SemiVariance.Direction.DOWNSIDE);
        assertFalse(sv3.isBiasCorrected());
        
        SemiVariance sv4 = new SemiVariance(sv3);
        assertFalse(sv4.isBiasCorrected());
        
        SemiVariance sv5 = sv.withBiasCorrected(false);
        assertFalse(sv5.isBiasCorrected());
    }
}
