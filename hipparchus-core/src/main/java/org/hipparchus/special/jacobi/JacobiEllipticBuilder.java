/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.special.jacobi;

/** Builder for algorithms compmuting Jacobi elliptic functions.
 * The Jacobi elliptic functions are related to elliptic integrals.
 * We use hare the notations from <a
 * href="https://en.wikipedia.org/wiki/Abramowitz_and_Stegun">Abramowitz and
 * Stegun</a> (Ch. 16) with a parameter {@code m}. The notations from
 * <a href="https://dlmf.nist.gov/22">Digital Library of Mathematical Functions (Ch. 22)</a>
 * are different as they use modulus {@code k} instead of parameter {@code m},
 * with {@code k² = m}.
 * @since 2.0
 */
public class JacobiEllipticBuilder {

    /** Threshold near 0 for using specialized algorithm. */
    private static final double NEAR_ZERO = 1.0e-9;

    /** Threshold near 1 for using specialized algorithm. */
    private static final double NEAR_ONE = 1.0 - NEAR_ZERO;

    /** Private constructor for utility class.
     */
    private JacobiEllipticBuilder() {
        // nothing to do
    }

    /** Build an algorithm for computing Jacobi elliptic functions.
     * @param m parameter of the Jacobi elliptic function
     * @return selected algorithm
     */
    public static JacobiElliptic build(final double m) {
        if (m < 0) {
            return new NegativeParameter(m);
        } else if (m > 1) {
            return new BigParameter(m);
        } else if (m < NEAR_ZERO) {
            return new NearZeroParameter(m);
        } else if (m > NEAR_ONE) {
            return new NearOneParameter(m);
        } else {
            return new BoundedParameter(m);
        }
    }

}
