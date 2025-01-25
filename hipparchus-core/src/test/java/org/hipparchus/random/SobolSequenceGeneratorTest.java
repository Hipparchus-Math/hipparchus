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
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.InputStream;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

class SobolSequenceGeneratorTest {

    private final double[][] referenceValues = {
            { 0.0, 0.0, 0.0 },
            { 0.5, 0.5, 0.5 },
            { 0.75, 0.25, 0.25 },
            { 0.25, 0.75, 0.75 },
            { 0.375, 0.375, 0.625 },
            { 0.875, 0.875, 0.125 },
            { 0.625, 0.125, 0.875 },
            { 0.125, 0.625, 0.375 },
            { 0.1875, 0.3125, 0.9375 },
            { 0.6875, 0.8125, 0.4375 }
    };

    private SobolSequenceGenerator generator;

    @BeforeEach
    void setUp() {
        generator = new SobolSequenceGenerator(3);
    }

    @Test
    void test3DReference() {
        for (int i = 0; i < referenceValues.length; i++) {
            double[] result = generator.nextVector();
            assertArrayEquals(referenceValues[i], result, 1e-6);
            assertEquals(i + 1, generator.getNextIndex());
        }
    }

    @Test
    void testConstructor() {
        try {
            new SobolSequenceGenerator(0);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        try {
            new SobolSequenceGenerator(21202);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
    }

    @Test
    void testConstructor2() throws Exception{
        try {
            final String RESOURCE_NAME = "/assets/org/hipparchus/random/new-joe-kuo-6.21201";
            final InputStream is = getClass().getResourceAsStream(RESOURCE_NAME);
            new SobolSequenceGenerator(21202, is);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        try {
            new SobolSequenceGenerator(21202);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
    }

    @Test
    void testSkip() {
        double[] result = generator.skipTo(5);
        assertArrayEquals(referenceValues[5], result, 1e-6);
        assertEquals(6, generator.getNextIndex());

        for (int i = 6; i < referenceValues.length; i++) {
            result = generator.nextVector();
            assertArrayEquals(referenceValues[i], result, 1e-6);
            assertEquals(i + 1, generator.getNextIndex());
        }
    }

}
