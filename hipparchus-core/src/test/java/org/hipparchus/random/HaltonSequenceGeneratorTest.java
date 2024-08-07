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
package org.hipparchus.random;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

class HaltonSequenceGeneratorTest {

    private double[][] referenceValues = {
            { 0.0,    0.0,    0.0  },
            { 0.5,    0.6667, 0.6  },
            { 0.25,   0.3333, 0.2  },
            { 0.75,   0.2223, 0.8  },
            { 0.125,  0.8888, 0.4  },
            { 0.625,  0.5555, 0.12 },
            { 0.375,  0.1111, 0.72 },
            { 0.875,  0.7777, 0.32 },
            { 0.0625, 0.4444, 0.92 },
            { 0.5625, 0.0740, 0.52 }
    };

    private double[][] referenceValuesUnscrambled = {
            { 0.0,    0.0    },
            { 0.5,    0.3333 },
            { 0.25,   0.6666 },
            { 0.75,   0.1111 },
            { 0.125,  0.4444 },
            { 0.625,  0.7777 },
            { 0.375,  0.2222 },
            { 0.875,  0.5555 },
            { 0.0625, 0.8888 },
            { 0.5625, 0.0370 }
    };

    private HaltonSequenceGenerator generator;

    @BeforeEach
    void setUp() {
        generator = new HaltonSequenceGenerator(3);
    }

    @Test
    void test3DReference() {
        for (int i = 0; i < referenceValues.length; i++) {
            double[] result = generator.nextVector();
            assertArrayEquals(referenceValues[i], result, 1e-3);
            assertEquals(i + 1, generator.getNextIndex());
        }
    }

    @Test
    void test2DUnscrambledReference() {
        generator = new HaltonSequenceGenerator(2, new int[] {2, 3}, null);
        for (int i = 0; i < referenceValuesUnscrambled.length; i++) {
            double[] result = generator.nextVector();
            assertArrayEquals(referenceValuesUnscrambled[i], result, 1e-3);
            assertEquals(i + 1, generator.getNextIndex());
        }
    }

    @Test
    void testConstructor() {
        try {
            new HaltonSequenceGenerator(0);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        try {
            new HaltonSequenceGenerator(41);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
    }

    @Test
    void testConstructor2() throws Exception{
        try {
            new HaltonSequenceGenerator(2, new int[] { 1 }, null);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        try {
            new HaltonSequenceGenerator(2, null, null);
            fail("an exception should have been thrown");
        } catch (NullArgumentException e) {
            // expected
        }

        try {
            new HaltonSequenceGenerator(2, new int[] { 1, 1 }, new int[] { 1 });
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
    }

    @Test
    void testSkip() {
        double[] result = generator.skipTo(5);
        assertArrayEquals(referenceValues[5], result, 1e-3);
        assertEquals(6, generator.getNextIndex());

        for (int i = 6; i < referenceValues.length; i++) {
            result = generator.nextVector();
            assertArrayEquals(referenceValues[i], result, 1e-3);
            assertEquals(i + 1, generator.getNextIndex());
        }
    }

}
