/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus;

import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public abstract class RealFieldElementAbstractTest<T extends RealFieldElement<T>>
    extends CalculusFieldElementAbstractTest<T> {

    @Test
    public void testRemainderField() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(FastMath.IEEEremainder(x, y), build(x).remainder(build(y)));
            }
        }
    }

    @Test
    public void testRemainderDouble() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3.2; y < 3.2; y += 0.25) {
                checkRelative(FastMath.IEEEremainder(x, y), build(x).remainder(y));
            }
        }
    }

    @Test
    public void testAbs() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.abs(x), build(x).abs());
        }
    }

    @Test
    public void testCeil() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.ceil(x), build(x).ceil());
        }
    }

    @Test
    public void testFloor() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.floor(x), build(x).floor());
        }
    }

    @Test
    public void testRint() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.rint(x), build(x).rint());
        }
    }

    @Test
    public void testRound() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            Assert.assertEquals(FastMath.round(x), build(x).round());
        }
    }

    @Test
    public void testSignum() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.signum(x), build(x).signum());
        }
    }

    @Test
    public void testCopySignField() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(FastMath.copySign(x, y), build(x).copySign(build(y)));
            }
        }
    }

    @Test
    public void testCopySignDouble() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(FastMath.copySign(x, y), build(x).copySign(y));
            }
        }
    }

}
