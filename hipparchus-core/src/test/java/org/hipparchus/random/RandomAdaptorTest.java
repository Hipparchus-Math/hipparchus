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

import org.hipparchus.exception.NullArgumentException;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test cases for the RandomAdaptor class.
 */
class RandomAdaptorTest {

    @Test
    void testAdaptor() {
        ConstantGenerator generator = new ConstantGenerator();
        Random random = RandomAdaptor.of(generator);
        checkConstant(random);
        RandomAdaptor randomAdaptor = new RandomAdaptor(generator);
        checkConstant(randomAdaptor);
    }

    private void checkConstant(Random random) {
        byte[] bytes = new byte[] {0};
        random.nextBytes(bytes);
        assertEquals(0, bytes[0]);
        assertFalse(random.nextBoolean());
        assertEquals(0, random.nextDouble(), 0);
        assertEquals(0, random.nextFloat(), 0);
        assertEquals(0, random.nextGaussian(), 0);
        assertEquals(0, random.nextInt());
        assertEquals(0, random.nextInt(1));
        assertEquals(0, random.nextLong());
        random.setSeed(100);
        assertEquals(0, random.nextDouble(), 0);
    }

    @SuppressWarnings("unused")
    @Test
    void testNullGenerator(){
        assertThrows(NullArgumentException.class, () -> {
            RandomGenerator nullGenerator = null;
            Random random = new RandomAdaptor(nullGenerator);
        });
    }

    @SuppressWarnings("unused")
    @Test
    void testNullGenerator2(){
        assertThrows(NullArgumentException.class, () -> {
            RandomGenerator nullGenerator = null;
            Random random = RandomAdaptor.of(nullGenerator);
        });
    }

    /*
     * "Constant" generator to test Adaptor delegation.
     * "Powered by Eclipse ;-)"
     *
     */
    public static class ConstantGenerator implements RandomGenerator {

        private final double value;

        public ConstantGenerator() {
            value = 0;
        }

        public ConstantGenerator(double value) {
            this.value = value;
        }

        @Override
        public boolean nextBoolean() {
            return false;
        }

        @Override
        public void nextBytes(byte[] bytes) {
        }

        @Override
        public void nextBytes(byte[] bytes, int offset, int len) {
        }

        @Override
        public double nextDouble() {
            return value;
        }

        @Override
        public float nextFloat() {
            return (float) value;
        }

        @Override
        public double nextGaussian() {
            return value;
        }

        @Override
        public int nextInt() {
            return (int) value;
        }

        @Override
        public int nextInt(int n) {
            return (int) value;
        }

        @Override
        public long nextLong() {
            return (int) value;
        }

        @Override
        public long nextLong(long n) {
            return (int) value;
        }

        @Override
        public void setSeed(int seed) {
        }

        @Override
        public void setSeed(int[] seed) {
        }

        @Override
        public void setSeed(long seed) {
        }

    }
}
