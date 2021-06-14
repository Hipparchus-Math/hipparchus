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
package org.hipparchus.special.elliptic.jacobi;

import org.hipparchus.CalculusFieldElement;

/** Builder for algorithms compmuting Jacobi elliptic functions.
 * <p>
 * The Jacobi elliptic functions are related to elliptic integrals.
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
     * <p>
     * Beware that {@link org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral
     * Legendre elliptic integrals} are defined in terms of elliptic modulus {@code k} whereas
     * {@link JacobiEllipticBuilder#build(double) Jacobi elliptic functions} (which
     * are their inverses) are defined in terms of parameter {@code m} and {@link
     * JacobiTheta#JacobiTheta(double) Jacobi theta functions} are defined
     * in terms of the {@link
     * org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral#nome(double)
     * nome q}. All are related as {@code k² = m} and the nome can be computed from ratios of complete
     * elliptic integrals.
     * </p>
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
     * <p>
     * Beware that {@link org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral
     * Legendre elliptic integrals} are defined in terms of elliptic modulus {@code k} whereas
     * {@link JacobiEllipticBuilder#build(CalculusFieldElement) Jacobi elliptic functions} (which
     * are their inverses) are defined in terms of parameter {@code m} and {@link
     * FieldJacobiTheta#FieldJacobiTheta(CalculusFieldElement) Jacobi theta functions} are defined
     * in terms of the {@link
     * org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral#nome(CalculusFieldElement)
     * nome q}. All are related as {@code k² = m} and the nome can be computed from ratios of complete
     * elliptic integrals.
     * </p>
     * @param m parameter of the Jacobi elliptic function
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

}
