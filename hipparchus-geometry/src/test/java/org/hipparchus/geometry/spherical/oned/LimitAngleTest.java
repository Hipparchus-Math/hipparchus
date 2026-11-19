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
package org.hipparchus.geometry.spherical.oned;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathUtils;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class LimitAngleTest {

    @Test
    void testReversedLimit() {
        for (int k = -2; k < 3; ++k) {
            LimitAngle l  = new LimitAngle(new S1Point(1.0 + k * MathUtils.TWO_PI), false, 1.0e-10);
            assertEquals(l.getLocation().getAlpha(), l.getReverse().getLocation().getAlpha(), 1.0e-10);
            assertEquals(l.getTolerance(), l.getReverse().getTolerance(), 1.0e-10);
            assertTrue(l.sameOrientationAs(l));
            assertFalse(l.sameOrientationAs(l.getReverse()));
            assertEquals(MathUtils.TWO_PI, l.wholeSpace().getSize(), 1.0e-10);
            assertEquals(MathUtils.TWO_PI, l.getReverse().wholeSpace().getSize(), 1.0e-10);
        }
    }

    @Test
    void testTooSmallTolerance() {
        assertThrows(MathIllegalArgumentException.class, () -> new LimitAngle(new S1Point(1.0), false, 0.9 * Sphere1D.SMALLEST_TOLERANCE));
    }

}
