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
package org.hipparchus.util;

import java.io.Serializable;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.RandomGenerator;

/**
 * A helper class containing various {@link PivotingStrategy} implementations.
 */
class PivotingStrategies {

    /** Utility class, prevent instantiation. */
    private PivotingStrategies() {}

    /** Returns the central pivoting strategy. */
    static PivotingStrategy getCentralPivotingStrategy() {
        return CentralPivotingStrategy.INSTANCE;
    }

    /** Returns the median of 3 pivoting strategy. */
    static PivotingStrategy getMedianOf3PivotingStrategy() {
        return MedianOf3PivotingStrategy.INSTANCE;
    }

    /** Returns a new random pivoting strategy using the given random generator. */
    static PivotingStrategy getRandomPivotingStrategy(RandomGenerator random) {
        return new RandomPivotingStrategy(random);
    }

    /**
     * A mid point strategy based on the average of begin and end indices.
     */
    static class CentralPivotingStrategy implements PivotingStrategy, Serializable {

        /** Serializable UID. */
        private static final long serialVersionUID = 20140713L;
        /** The singleton instance. */
        private static final PivotingStrategy INSTANCE = new CentralPivotingStrategy();

        /** Singleton class, prevent instantiation. */
        private CentralPivotingStrategy() {}

        /**
         * {@inheritDoc}
         * This in particular picks a average of begin and end indices
         * @return The index corresponding to a simple average of
         * the first and the last element indices of the array slice
         * @throws MathIllegalArgumentException when indices exceeds range
         */
        @Override
        public int pivotIndex(final double[] work, final int begin, final int end)
            throws MathIllegalArgumentException {
            MathArrays.verifyValues(work, begin, end-begin);
            return begin + (end - begin)/2;
        }

        private Object readResolve() {
            return INSTANCE;
        }

    }

    /**
     * Classic median of 3 strategy given begin and end indices.
     */
    static class MedianOf3PivotingStrategy implements PivotingStrategy, Serializable {

        /** Serializable UID. */
        private static final long serialVersionUID = 20140713L;
        /** The singleton instance. */
        private static final PivotingStrategy INSTANCE = new MedianOf3PivotingStrategy();

        /** Singleton class, prevent instantiation. */
        private MedianOf3PivotingStrategy() {}

        /**{@inheritDoc}
         * This in specific makes use of median of 3 pivoting.
         * @return The index corresponding to a pivot chosen between the
         * first, middle and the last indices of the array slice
         * @throws MathIllegalArgumentException when indices exceeds range
         */
        @Override
        public int pivotIndex(final double[] work, final int begin, final int end)
            throws MathIllegalArgumentException {
            MathArrays.verifyValues(work, begin, end-begin);
            final int inclusiveEnd = end - 1;
            final int middle = begin + (inclusiveEnd - begin) / 2;
            final double wBegin = work[begin];
            final double wMiddle = work[middle];
            final double wEnd = work[inclusiveEnd];

            if (wBegin < wMiddle) {
                if (wMiddle < wEnd) {
                    return middle;
                } else {
                    return wBegin < wEnd ? inclusiveEnd : begin;
                }
            } else {
                if (wBegin < wEnd) {
                    return begin;
                } else {
                    return wMiddle < wEnd ? inclusiveEnd : middle;
                }
            }
        }

        private Object readResolve() {
            return INSTANCE;
        }

    }

    /**
     * A strategy of selecting random index between begin and end indices.
     */
    static class RandomPivotingStrategy implements PivotingStrategy, Serializable {

        /** Serializable UID. */
        private static final long serialVersionUID = 20140713L;

        /** Random generator to use for selecting pivot. */
        private final RandomGenerator random;

        /** Simple constructor.
         * @param random random generator to use for selecting pivot
         */
        private RandomPivotingStrategy(final RandomGenerator random) {
            this.random = random;
        }

        /**
         * {@inheritDoc}
         * A uniform random pivot selection between begin and end indices
         * @return The index corresponding to a random uniformly selected
         * value between first and the last indices of the array slice
         * @throws MathIllegalArgumentException when indices exceeds range
         */
        @Override
        public int pivotIndex(final double[] work, final int begin, final int end)
            throws MathIllegalArgumentException {
            MathArrays.verifyValues(work, begin, end-begin);
            return begin + random.nextInt(end - begin - 1);
        }

    }

}
