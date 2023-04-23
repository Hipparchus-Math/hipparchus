/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.special.elliptic.jacobi;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.FieldComplex;

/** Builder for algorithms compmuting Jacobi elliptic functions.
 * <p>
 * The Jacobi elliptic functions are related to elliptic integrals.
 * </p>
 * <p>
 * There are different conventions to interpret the arguments of
 * Jacobi elliptic functions. The first argument may be  the amplitude φ,
 * but is more often the variable u (with sn(u) = sin(φ) and cn(u) = cos(φ)).
 * The second argument  is either the modulus k or the parameter m with m = k².
 * In Hipparchus, we adopted the convention to use u and m.
 * </p>
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

    /** Build an algorithm for computing Jacobi elliptic functions.
     * @param m parameter of the Jacobi elliptic function
     * @param <T> type of the field elements
     * @return selected algorithm
     */
    public static <T extends CalculusFieldElement<T>> FieldJacobiElliptic<T> build(final T m) {
        if (m.getReal() < 0) {
            return new FieldNegativeParameter<>(m);
        } else if (m.getReal() > 1) {
            return new FieldBigParameter<>(m);
        } else if (m.getReal() < NEAR_ZERO) {
            return new FieldNearZeroParameter<>(m);
        } else if (m.getReal() > NEAR_ONE) {
            return new FieldNearOneParameter<>(m);
        } else {
            return new FieldBoundedParameter<>(m);
        }
    }

    /** Build an algorithm for computing Jacobi elliptic functions.
     * @param m parameter of the Jacobi elliptic function
     * @return selected algorithm
     */
    public static FieldJacobiElliptic<Complex> build(final Complex m) {
        return new ComplexParameter(m);
    }

    /** Build an algorithm for computing Jacobi elliptic functions.
     * @param m parameter of the Jacobi elliptic function
     * @param <T> type of the field elements
     * @return selected algorithm
     */
    public static <T extends CalculusFieldElement<T>> FieldJacobiElliptic<FieldComplex<T>> build(final FieldComplex<T> m) {
        return new FieldComplexParameter<>(m);
    }

}
