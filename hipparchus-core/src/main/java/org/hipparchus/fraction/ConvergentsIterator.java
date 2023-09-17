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
package org.hipparchus.fraction;

import java.util.Iterator;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Predicate;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Pair;
import org.hipparchus.util.Precision;

/**
 * Generator for convergents.
 */
class ConvergentsIterator {

    /** Unused constructor.
     */
    private ConvergentsIterator() {
    } // static use only

    /** Container for one convergent step. */
    static class ConvergenceStep {
        /** Numerator of previous convergent. */
        private final long   p0;

        /** Denominator of previous convergent. */
        private final long   q0;

        /** Numerator of current convergent. */
        private final long   p1;

        /** Denominator of current convergent. */
        private final long   q1;

        /** Remainder of current convergent. */
        private final double r1;

        private ConvergenceStep(final long p0, final long q0, final long p1, final long q1, final double r1) {
            this.p0 = p1;
            this.q0 = q1;
            final long a1 = (long) FastMath.floor(r1);
            try {
                this.p1 = FastMath.addExact(Math.multiplyExact(a1, p1), p0);
                this.q1 = FastMath.addExact(Math.multiplyExact(a1, q1), q0);
                this.r1 = 1.0 / (r1 - a1);
            } catch (ArithmeticException e) { // unlike the name implies FastMath's multiplyExact() is slower
                throw new MathIllegalStateException(e, LocalizedCoreFormats.FRACTION_CONVERSION_OVERFLOW, r1, p1, q1);
            }
        }

        /** Builder from a double value.
         * @param value value to approximate
         * @return first step in approximation
         */
        public static ConvergenceStep start(double value) {
            return new ConvergenceStep(0, 1, 1, 0, value);
        }

        /** Compute next step in convergence.
         * @return next convergence step
         */
        public ConvergenceStep next() {
            return new ConvergenceStep(p0, q0, p1, q1, r1);
        }

        /** Get the numerator of current convergent.
         * @return numerator of current convergent
         */
        public long getNumerator() {
            return p1;
        }

        /** Get the denominator of current convergent.
         * @return denominator of current convergent
         */
        public long getDenominator() {
            return q1;
        }

        /** Compute double value of current convergent.
         * @return double value of current convergent
         */
        public double getFractionValue() {
            return getNumerator() / (double) getDenominator();
        }

        /** Convert convergent to string representation.
         * @return string representation of convergent
         */
        @Override
        public String toString() {
            return getNumerator() + "/" + getDenominator();
        }

    }

    /**
     * Returns the last element of the series of convergent-steps to approximate the
     * given value.
     * <p>
     * The series terminates either at the first step that satisfies the given
     * {@code convergenceTest} or after at most {@code maxConvergents} elements. The
     * returned Pair consists of that terminal step and a {@link Boolean} that
     * indicates if it satisfies the given convergence tests. If the returned pair's
     * value is {@code false} the element at position {@code maxConvergents} was
     * examined but failed to satisfy the {@code convergenceTest}.
     *
     * @param value           value to approximate
     * @param maxConvergents  maximum number of convergents to examine
     * @param convergenceTests the test if the series has converged at a step
     * @return the pair of last element of the series of convergents and a boolean
     *         indicating if that element satisfies the specified convergent test
     */
    static Pair<ConvergenceStep, Boolean> convergent(double value, int maxConvergents,
            Predicate<ConvergenceStep> convergenceTests) {
        ConvergenceStep step = ConvergenceStep.start(value);
        for (int i = 1; i < maxConvergents; i++) { // start performs first iteration
            if (convergenceTests.test(step)) {
                return Pair.create(step, Boolean.TRUE);
            }
            step = step.next();
        }
        return Pair.create(step, convergenceTests.test(step));
    }

    /**
     * Generate a {@link Stream stream} of {@code ConvergenceStep convergent-steps}
     * from a real number.
     *
     * @param value           value to approximate
     * @param maxConvergents maximum number of convergent steps.
     * @return stream of {@link ConvergenceStep convergent-steps} approximating
     *         {@code value}
     */
    static Stream<ConvergenceStep> convergents(final double value, final int maxConvergents) {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(new Iterator<ConvergenceStep>() {

            /** Next convergent. */
            private ConvergenceStep next = ConvergenceStep.start(value);

            /** {@inheritDoc} */
            @Override
            public boolean hasNext() {
                return next != null;
            }

            /** {@inheritDoc} */
            @Override
            public ConvergenceStep next() {
                final ConvergenceStep ret = next;
                if (Precision.equals(ret.getFractionValue(), value, 1)) {
                    next = null; // stop if precision has been reached
                } else {
                    next = next.next();
                }
                return ret;
            }

        }, Spliterator.DISTINCT | Spliterator.NONNULL | Spliterator.IMMUTABLE | Spliterator.ORDERED),
        false).
        limit(maxConvergents);
    }
}
