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
package org.hipparchus.complex;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;


/**
 * Unit tests for the {@link RootsOfUnity} class.
 *
 */
class RootsOfUnityTest {

    @Test
    void testMathIllegalState1() {
        assertThrows(MathIllegalStateException.class, () -> {
            final RootsOfUnity roots = new RootsOfUnity();
            roots.getReal(0);
        });
    }

    @Test
    void testMathIllegalState2() {
        assertThrows(MathIllegalStateException.class, () -> {
            final RootsOfUnity roots = new RootsOfUnity();
            roots.getImaginary(0);
        });
    }

    @Test
    void testMathIllegalState3() {
        assertThrows(MathIllegalStateException.class, () -> {
            final RootsOfUnity roots = new RootsOfUnity();
            roots.isCounterClockWise();
        });
    }

    @Test
    void testZeroNumberOfRoots() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final RootsOfUnity roots = new RootsOfUnity();
            roots.computeRoots(0);
        });
    }

    @Test
    void testGetNumberOfRoots() {
        final RootsOfUnity roots = new RootsOfUnity();
        assertEquals(0, roots.getNumberOfRoots(), "");
        roots.computeRoots(5);
        assertEquals(5, roots.getNumberOfRoots(), "");
        /*
         * Testing -5 right after 5 is important, as the roots in this case are
         * not recomputed.
         */
        roots.computeRoots(-5);
        assertEquals(5, roots.getNumberOfRoots(), "");
        roots.computeRoots(6);
        assertEquals(6, roots.getNumberOfRoots(), "");
    }

    @Test
    public void testErrorConditions() {
        final RootsOfUnity roots = new RootsOfUnity();
        roots.computeRoots(5);
        try {
            roots.getReal(-2);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_ROOT_OF_UNITY_INDEX,
                         e.getSpecifier());
            assertEquals(-2, e.getParts()[0]);
            assertEquals( 0, e.getParts()[1]);
            assertEquals( 4, e.getParts()[2]);
        }
        try {
            roots.getImaginary(-2);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_ROOT_OF_UNITY_INDEX,
                         e.getSpecifier());
            assertEquals(-2, e.getParts()[0]);
            assertEquals( 0, e.getParts()[1]);
            assertEquals( 4, e.getParts()[2]);
        }
        try {
            roots.getReal(5);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_ROOT_OF_UNITY_INDEX,
                         e.getSpecifier());
            assertEquals( 5, e.getParts()[0]);
            assertEquals( 0, e.getParts()[1]);
            assertEquals( 4, e.getParts()[2]);
        }
        try {
            roots.getImaginary(5);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_ROOT_OF_UNITY_INDEX,
                         e.getSpecifier());
            assertEquals( 5, e.getParts()[0]);
            assertEquals( 0, e.getParts()[1]);
            assertEquals( 4, e.getParts()[2]);
        }
    }

    @Test
    void testComputeRoots() {
        final RootsOfUnity roots = new RootsOfUnity();
        for (int n = -10; n < 11; n++) {
            /*
             * Testing -n right after n is important, as the roots in this case
             * are not recomputed.
             */
            if (n != 0) {
                roots.computeRoots(n);
                doTestComputeRoots(roots);
                roots.computeRoots(-n);
                doTestComputeRoots(roots);
            }
        }
    }

    private void doTestComputeRoots(final RootsOfUnity roots) {
        final int n = roots.isCounterClockWise() ? roots.getNumberOfRoots() :
            -roots.getNumberOfRoots();
        final double tol = 10 * Math.ulp(1.0);
        for (int k = 0; k < n; k++) {
            final double t = 2.0 * FastMath.PI * k / n;
            final String msg = String.format("n = %d, k = %d", n, k);
            assertEquals(FastMath.cos(t), roots.getReal(k), tol, msg);
            assertEquals(FastMath.sin(t), roots.getImaginary(k), tol, msg);
        }
    }
}
