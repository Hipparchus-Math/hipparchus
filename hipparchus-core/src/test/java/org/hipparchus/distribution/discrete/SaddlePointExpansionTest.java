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
package org.hipparchus.distribution.discrete;

import org.junit.Assert;
import org.junit.Test;

public class SaddlePointExpansionTest {

    @Test
    public void testSmallInteger() {
        for (int n = 4; n < 16; ++n) {
             Assert.assertEquals(alternateStirlingErrorImplementation(n, 8),
                                SaddlePointExpansion.getStirlingError(n),
                                1.0e-10);
        }
    }

    @Test
    public void testSmallNonInteger() {
        for (double z = 3.75; z < 14.8; z += 1.0) {
            Assert.assertEquals(alternateStirlingErrorImplementation(z, 12),
                                SaddlePointExpansion.getStirlingError(z),
                                1.0e-10);
        }
    }

    @Test
    public void testLargeValues() {
        for (double z = 15.25; z < 25.5; z += 1.0) {
            Assert.assertEquals(alternateStirlingErrorImplementation(z, 21),
                                SaddlePointExpansion.getStirlingError(z),
                                1.0e-15);
        }
    }

    @Test
    public void testSpecialValues() {
        Assert.assertEquals(0.0, SaddlePointExpansion.logBinomialProbability(0, 0, 0.6, 0.4), 1.0e-15);
        Assert.assertEquals(Double.NEGATIVE_INFINITY, SaddlePointExpansion.logBinomialProbability(1, 0, 0.6, 0.4), 1.0e-15);
    }

    private double alternateStirlingErrorImplementation(final double z, final int kMax) {
        // numerators of coefficients in Stirling's expansion for log(Γ(z)) from https://oeis.org/A046968
        final double[] A046968 = {
            1.0, -1.0, 1.0, -1.0, 1.0, -691.0, 1.0,
            -3617.0, 43867.0, -174611.0, 77683.0, -236364091.0, 657931.0, -3392780147.0,
            1723168255201.0, -7709321041217.0, 151628697551.0, -26315271553053477373.0,
            154210205991661.0, -261082718496449122051.0, 1520097643918070802691.0
        };
        // denominators of coefficients in Stirling's expansion for log(Γ(z)) from https://oeis.org/A046969
        final double[] A046969 = {
            12.0, 360.0, 1260.0, 1680.0, 1188.0, 360360.0, 156.0,
            122400.0, 244188.0, 125400.0, 5796.0, 1506960.0, 300.0, 93960.0,
            2492028.0, 505920.0, 396.0, 2418179400.0, 444.0, 21106800.0, 3109932.0
        };
        double zk   = z;
        double expansion = 0.0;
        for (int i = 0; i < kMax; ++i) {
            expansion += A046968[i] / (zk * A046969[i]);
            zk *= z * z;
        }
        return expansion;
    }
}
