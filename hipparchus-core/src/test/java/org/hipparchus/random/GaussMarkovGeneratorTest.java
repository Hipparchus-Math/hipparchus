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
package org.hipparchus.random;

import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class GaussMarkovGeneratorTest {

    @Test
    void testExpectation() {
        final double          tau             = 3600.0;
        final double          stationarySigma = 0.2;
        final RandomGenerator random          = new Well1024a(0xb8005f29892534a8L);
        GaussMarkovGenerator  gm              = new GaussMarkovGenerator(tau,
                                                                         stationarySigma,
                                                                         random);
        assertEquals(tau,             gm.getTau(),             1.0e-15);
        assertEquals(stationarySigma, gm.getStationarySigma(), 1.0e-15);

        double sum = 0;
        int count = 0;
        final double deltaT = 0.1;
        for (double t = 0; t < 1000000; t += deltaT) {
            sum += gm.next(deltaT);

            ++count;
        }
        assertEquals(0.0, sum / count, 0.014);
    }

    @Test
    void testVariance() {
        final double          tau             = 3600.0;
        final double          stationarySigma = 0.2;
        final RandomGenerator random          = new Well1024a(0x09efbd4e87e7791eL);
        GaussMarkovGenerator gm = new GaussMarkovGenerator(tau,
                                                           stationarySigma,
                                                           random);

        // here, we already assume expectation is 0
        double sum2 = 0;
        int count = 0;
        final double deltaT = 0.1;
        for (double t = 0; t < 1000000; t += deltaT) {
            final double v = gm.next(tau);
            ++count;
            sum2 += v * v;
        }

        assertEquals(gm.getStationarySigma(), FastMath.sqrt(sum2 / count), 3.7e-5);

    }

}
