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
package org.hipparchus.random;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Base class for RandomGenerator tests.
 * <p>
 * Tests RandomGenerator methods directly and also executes RandomDataTest
 * test cases against a RandomDataImpl created using the provided generator.
 * <p>
 * RandomGenerator test classes should extend this class, implementing
 * makeGenerator() to provide a concrete generator to test. The generator
 * returned by makeGenerator should be seeded with a fixed seed.
 */
public abstract class RandomGeneratorAbstractTest extends RandomDataGeneratorTest {

    /** RandomGenerator under test */
    protected RandomGenerator generator;

    /**
     * Override this method in subclasses to provide a concrete generator to test.
     * Return a generator seeded with a fixed seed.
     */
    protected abstract RandomGenerator makeGenerator();

    /**
     * Initialize generator and randomData instance in superclass.
     */
    public RandomGeneratorAbstractTest() {
        generator = makeGenerator();
        randomData = RandomDataGenerator.of(generator);
    }

    /**
     * Set a fixed seed for the tests
     */
    @Before
    public void setUp() {
        generator = makeGenerator();
    }

    /**
     * Tests uniformity of nextInt(int) distribution by generating 1000
     * samples for each of 10 test values and for each sample performing
     * a chi-square test of homogeneity of the observed distribution with
     * the expected uniform distribution.  Tests are performed at the .01
     * level and an average failure rate higher than 2% (i.e. more than 20
     * null hypothesis rejections) causes the test case to fail.
     *
     * All random values are generated using the generator instance used by
     * other tests and the generator is not reseeded, so this is a fixed seed
     * test.
     */
    @Test
    public void testNextIntDirect() {
        // Set up test values - end of the array filled randomly
        int[] testValues = new int[] {4, 10, 12, 32, 100, 10000, 0, 0, 0, 0};
        for (int i = 6; i < 10; i++) {
            final int val = generator.nextInt();
            testValues[i] = val < 0 ? -val : val + 1;
        }

        final int numTests = 1000;
        for (int i = 0; i < testValues.length; i++) {
            final int n = testValues[i];
            // Set up bins
            int[] binUpperBounds;
            if (n < 32) {
                binUpperBounds = new int[n];
                for (int k = 0; k < n; k++) {
                    binUpperBounds[k] = k;
                }
            } else {
                binUpperBounds = new int[10];
                final int step = n / 10;
                for (int k = 0; k < 9; k++) {
                    binUpperBounds[k] = (k + 1) * step;
                }
                binUpperBounds[9] = n - 1;
            }
            // Run the tests
            int numFailures = 0;
            final int binCount = binUpperBounds.length;
            final long[] observed = new long[binCount];
            final double[] expected = new double[binCount];
            expected[0] = binUpperBounds[0] == 0 ? (double) smallSampleSize / (double) n :
                (double) ((binUpperBounds[0] + 1) * smallSampleSize) / (double) n;
            for (int k = 1; k < binCount; k++) {
                expected[k] = (double) smallSampleSize *
                (double) (binUpperBounds[k] - binUpperBounds[k - 1]) / n;
            }
            for (int j = 0; j < numTests; j++) {
                Arrays.fill(observed, 0);
                for (int k = 0; k < smallSampleSize; k++) {
                    final int value = generator.nextInt(n);
                    assertTrue("nextInt range",(value >= 0) && (value < n));
                    for (int l = 0; l < binCount; l++) {
                        if (binUpperBounds[l] >= value) {
                            observed[l]++;
                            break;
                        }
                    }
                }
                if (UnitTestUtils.chiSquareTest(expected, observed) < 0.01) {
                    numFailures++;
                }
            }
            if ((double) numFailures / (double) numTests > 0.02) {
                fail("Too many failures for n = " + n +
                     " " + numFailures + " out of " + numTests + " tests failed.");
            }
        }
    }

    @Test
    public void testNextLongDirect() {
        long q1 = Long.MAX_VALUE/4;
        long q2 = 2 *  q1;
        long q3 = 3 * q1;

        long[] observed = new long[4];
        long val = 0;
        for (int i = 0; i < smallSampleSize; i++) {
            val = generator.nextLong();
            val = val < 0 ? -val : val;
            if (val < q1) {
               observed[0] = ++observed[0];
            } else if (val < q2) {
               observed[1] = ++observed[1];
            } else if (val < q3) {
               observed[2] = ++observed[2];
            } else {
               observed[3] = ++observed[3];
            }
        }

        /* Use ChiSquare dist with df = 4-1 = 3, alpha = .001
         * Change to 11.34 for alpha = .01
         */
        assertTrue("chi-square test -- will fail about 1 in 1000 times",
                   UnitTestUtils.chiSquare(expected,observed) < 16.27);
    }

    @Test
    public void testNextBooleanDirect() {
        long halfSampleSize = smallSampleSize / 2;
        double[] expected = {halfSampleSize, halfSampleSize};
        long[] observed = new long[2];
        for (int i=0; i<smallSampleSize; i++) {
            if (generator.nextBoolean()) {
                observed[0]++;
            } else {
                observed[1]++;
            }
        }
        /* Use ChiSquare dist with df = 2-1 = 1, alpha = .001
         * Change to 6.635 for alpha = .01
         */
        assertTrue("chi-square test -- will fail about 1 in 1000 times",
                   UnitTestUtils.chiSquare(expected,observed) < 10.828);
    }

    @Test
    public void testNextFloatDirect() {
        long[] observed = new long[4];
        for (int i=0; i<smallSampleSize; i++) {
            float val = generator.nextFloat();
            if (val < 0.25) {
                observed[0] = ++observed[0];
            } else if (val < 0.5) {
                observed[1] = ++observed[1];
            } else if (val < 0.75) {
                observed[2] = ++observed[2];
            } else {
                observed[3] = ++observed[3];
            }
        }

        /* Use ChiSquare dist with df = 4-1 = 3, alpha = .001
         * Change to 11.34 for alpha = .01
         */
        assertTrue("chi-square test -- will fail about 1 in 1000 times",
                   UnitTestUtils.chiSquare(expected,observed) < 16.27);
    }

    @Test
    public void testNextDouble() {
        final double[] sample = new double[10000];
        final double[] expected = new double[100];
        final long[] observed = new long[100];
        Arrays.fill(expected, 100d);
        for (int i = 0; i < sample.length; i++) {
            sample[i] = generator.nextDouble();
            int j = 0;
            while (sample[i] < j / 100) {
                j++;
            }
            observed[j] = observed[j]++;
        }
        UnitTestUtils.assertChiSquareAccept(expected, observed, 0.01);
    }


    @Test(expected=MathIllegalArgumentException.class)
    public void testNextIntPrecondition1() {
        generator.nextInt(-1);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNextIntPrecondition2() {
        generator.nextInt(0);
    }

    @Test
    public void testNextInt2() {
        int walk = 0;
        final int N = 10000;
        for (int k = 0; k < N; ++k) {
           if (generator.nextInt() >= 0) {
               ++walk;
           } else {
               --walk;
           }
        }
        assertTrue("Walked too far astray: " + walk + "\nNote: This " +
                   "test will fail randomly about 1 in 100 times.",
                   FastMath.abs(walk) < FastMath.sqrt(N) * 2.576);
    }

    @Test
    public void testNextLong2() {
        int walk = 0;
        final int N = 1000;
        for (int k = 0; k < N; ++k) {
           if (generator.nextLong() >= 0) {
               ++walk;
           } else {
               --walk;
           }
        }
        assertTrue("Walked too far astray: " + walk + "\nNote: This " +
                   "test will fail randomly about 1 in 100 times.",
                   FastMath.abs(walk) < FastMath.sqrt(N) * 2.576);
    }

    @Test
    public void testNextBoolean2() {
        int walk = 0;
        final int N = 10000;
        for (int k = 0; k < N; ++k) {
           if (generator.nextBoolean()) {
               ++walk;
           } else {
               --walk;
           }
        }
        assertTrue("Walked too far astray: " + walk + "\nNote: This " +
                   "test will fail randomly about 1 in 100 times.",
                   FastMath.abs(walk) < FastMath.sqrt(N) * 2.576);
    }

    @Test
    public void testNextBytes() {
        long[] count = new long[256];
        byte[] bytes = new byte[10];
        double[] expected = new double[256];
        final int sampleSize = 100000;

        for (int i = 0; i < 256; i++) {
            expected[i] = (double) sampleSize / 265f;
        }

        for (int k = 0; k < sampleSize; ++k) {
           generator.nextBytes(bytes);
           for (byte b : bytes) {
               ++count[b + 128];
           }
        }

        UnitTestUtils.assertChiSquareAccept(expected, count, 0.001);
    }

    // MATH-1300
    @Test
    public void testNextBytesChunks() {
        final int[] chunkSizes = { 4, 8, 12, 16 };
        final int[] chunks = { 1, 2, 3, 4, 5 };
        for (int chunkSize : chunkSizes) {
            for (int numChunks : chunks) {
                checkNextBytesChunks(chunkSize, numChunks);
            }
        }
    }

    @Test
    public void testSeeding() {
        // makeGenerator initializes with fixed seed
        RandomGenerator gen = makeGenerator();
        RandomGenerator gen1 = makeGenerator();
        checkSameSequence(gen, gen1);
        // reseed, but recreate the second one
        // verifies MATH-723
        gen.setSeed(100);
        gen1 = makeGenerator();
        gen1.setSeed(100);
        checkSameSequence(gen, gen1);
    }

    private void checkSameSequence(RandomGenerator gen1, RandomGenerator gen2) {
        final int len = 11;  // Needs to be an odd number to check MATH-723
        final double[][] values = new double[2][len];
        for (int i = 0; i < len; i++) {
            values[0][i] = gen1.nextDouble();
        }
        for (int i = 0; i < len; i++) {
            values[1][i] = gen2.nextDouble();
        }
        assertTrue(Arrays.equals(values[0], values[1]));
        for (int i = 0; i < len; i++) {
            values[0][i] = gen1.nextFloat();
        }
        for (int i = 0; i < len; i++) {
            values[1][i] = gen2.nextFloat();
        }
        assertTrue(Arrays.equals(values[0], values[1]));
        for (int i = 0; i < len; i++) {
            values[0][i] = gen1.nextInt();
        }
        for (int i = 0; i < len; i++) {
            values[1][i] = gen2.nextInt();
        }
        assertTrue(Arrays.equals(values[0], values[1]));
        for (int i = 0; i < len; i++) {
            values[0][i] = gen1.nextLong();
        }
        for (int i = 0; i < len; i++) {
            values[1][i] = gen2.nextLong();
        }
        assertTrue(Arrays.equals(values[0], values[1]));
        for (int i = 0; i < len; i++) {
            values[0][i] = gen1.nextInt(len);
        }
        for (int i = 0; i < len; i++) {
            values[1][i] = gen2.nextInt(len);
        }
        assertTrue(Arrays.equals(values[0], values[1]));
        for (int i = 0; i < len; i++) {
            values[0][i] = gen1.nextBoolean() ? 1 : 0;
        }
        for (int i = 0; i < len; i++) {
            values[1][i] = gen2.nextBoolean() ? 1 : 0;
        }
        assertTrue(Arrays.equals(values[0], values[1]));
        for (int i = 0; i < len; i++) {
            values[0][i] = gen1.nextGaussian();
        }
        for (int i = 0; i < len; i++) {
            values[1][i] = gen2.nextGaussian();
        }
        assertTrue(Arrays.equals(values[0], values[1]));
    }

    // MATH-1300
    private void checkNextBytesChunks(int chunkSize,
                                      int numChunks) {
        final RandomGenerator rg = makeGenerator();
        final long seed = 1234567L;

        final byte[] b1 = new byte[chunkSize * numChunks];
        final byte[] b2 = new byte[chunkSize];

        // Generate the chunks in a single call.
        rg.setSeed(seed);
        rg.nextBytes(b1);

        // Reset.
        rg.setSeed(seed);
        // Generate the chunks in consecutive calls.
        for (int i = 0; i < numChunks; i++) {
            rg.nextBytes(b2);
        }

        // Store last 128 bytes chunk of b1 into b3.
        final byte[] b3 = new byte[chunkSize];
        System.arraycopy(b1, b1.length - b3.length, b3, 0, b3.length);

        // Sequence of calls must be the same.
        Assert.assertArrayEquals("chunkSize=" + chunkSize + " numChunks=" + numChunks,
                                 b2, b3);
    }

    @Override
    public void testNextZipf() {
        // Skip this test for the individual generators
    }
}
