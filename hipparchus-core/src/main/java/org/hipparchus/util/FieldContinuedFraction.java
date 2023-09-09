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
package org.hipparchus.util;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;

/**
 * Provides a generic means to evaluate continued fractions.  Subclasses simply
 * provided the a and b coefficients to evaluate the continued fraction.
 * <p>
 * References:
 * <ul>
 * <li><a href="http://mathworld.wolfram.com/ContinuedFraction.html">
 * Continued Fraction</a></li>
 * </ul>
 *
 */
public abstract class FieldContinuedFraction {
    /** Maximum allowed numerical error. */
    private static final double DEFAULT_EPSILON = 10e-9;

    /**
     * Default constructor.
     */
    protected FieldContinuedFraction() {
        super();
    }

    /**
     * Access the n-th a coefficient of the continued fraction.  Since a can be
     * a function of the evaluation point, x, that is passed in as well.
     * @param n the coefficient index to retrieve.
     * @param x the evaluation point.
     * @param <T> type of the field elements.
     * @return the n-th a coefficient.
     */
    public abstract <T extends CalculusFieldElement<T>> T getA(int n, T x);

    /**
     * Access the n-th b coefficient of the continued fraction.  Since b can be
     * a function of the evaluation point, x, that is passed in as well.
     * @param n the coefficient index to retrieve.
     * @param x the evaluation point.
     * @param <T> type of the field elements.
     * @return the n-th b coefficient.
     */
    public abstract <T extends CalculusFieldElement<T>> T getB(int n, T x);

    /**
     * Evaluates the continued fraction at the value x.
     * @param x the evaluation point.
     * @param <T> type of the field elements.
     * @return the value of the continued fraction evaluated at x.
     * @throws MathIllegalStateException if the algorithm fails to converge.
     */
    public <T extends CalculusFieldElement<T>> T evaluate(T x) throws MathIllegalStateException {
        return evaluate(x, DEFAULT_EPSILON, Integer.MAX_VALUE);
    }

    /**
     * Evaluates the continued fraction at the value x.
     * @param x the evaluation point.
     * @param epsilon maximum error allowed.
     * @param <T> type of the field elements.
     * @return the value of the continued fraction evaluated at x.
     * @throws MathIllegalStateException if the algorithm fails to converge.
     */
    public <T extends CalculusFieldElement<T>> T evaluate(T x, double epsilon) throws MathIllegalStateException {
        return evaluate(x, epsilon, Integer.MAX_VALUE);
    }

    /**
     * Evaluates the continued fraction at the value x.
     * @param x the evaluation point.
     * @param maxIterations maximum number of convergents
     * @param <T> type of the field elements.
     * @return the value of the continued fraction evaluated at x.
     * @throws MathIllegalStateException if the algorithm fails to converge.
     * @throws MathIllegalStateException if maximal number of iterations is reached
     */
    public <T extends CalculusFieldElement<T>> T evaluate(T x, int maxIterations)
        throws MathIllegalStateException {
        return evaluate(x, DEFAULT_EPSILON, maxIterations);
    }

    /**
     * Evaluates the continued fraction at the value x.
     * <p>
     * The implementation of this method is based on the modified Lentz algorithm as described
     * on page 18 ff. in:
     * </p>
     * <ul>
     *   <li>
     *   I. J. Thompson,  A. R. Barnett. "Coulomb and Bessel Functions of Complex Arguments and Order."
     *   <a target="_blank" href="http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf">
     *   http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf</a>
     *   </li>
     * </ul>
     * <p>
     * <b>Note:</b> the implementation uses the terms a<sub>i</sub> and b<sub>i</sub> as defined in
     * <a href="http://mathworld.wolfram.com/ContinuedFraction.html">Continued Fraction @ MathWorld</a>.
     * </p>
     *
     * @param x the evaluation point.
     * @param epsilon maximum error allowed.
     * @param maxIterations maximum number of convergents
     * @param <T> type of the field elements.
     * @return the value of the continued fraction evaluated at x.
     * @throws MathIllegalStateException if the algorithm fails to converge.
     * @throws MathIllegalStateException if maximal number of iterations is reached
     */
    public <T extends CalculusFieldElement<T>> T evaluate(T x, double epsilon, int maxIterations)
        throws MathIllegalStateException {
        final T zero = x.getField().getZero();
        final T one  = x.getField().getOne();

        final double small      = 1e-50;
        final T      smallField = one.multiply(small);

        T hPrev = getA(0, x);

        // use the value of small as epsilon criteria for zero checks
        if (Precision.equals(hPrev.getReal(), 0.0, small)) {
            hPrev = one.multiply(small);
        }

        int n     = 1;
        T   dPrev = zero;
        T   cPrev = hPrev;
        T   hN    = hPrev;

        while (n < maxIterations) {
            final T a = getA(n, x);
            final T b = getB(n, x);

            T dN = a.add(b.multiply(dPrev));
            if (Precision.equals(dN.getReal(), 0.0, small)) {
                dN = smallField;
            }
            T cN = a.add(b.divide(cPrev));
            if (Precision.equals(cN.getReal(), 0.0, small)) {
                cN = smallField;
            }

            dN = dN.reciprocal();
            final T deltaN = cN.multiply(dN);
            hN = hPrev.multiply(deltaN);

            if (hN.isInfinite()) {
                throw new MathIllegalStateException(LocalizedCoreFormats.CONTINUED_FRACTION_INFINITY_DIVERGENCE, x);
            }
            if (hN.isNaN()) {
                throw new MathIllegalStateException(LocalizedCoreFormats.CONTINUED_FRACTION_NAN_DIVERGENCE, x);
            }

            if (deltaN.subtract(1.0).abs().getReal() < epsilon) {
                break;
            }

            dPrev = dN;
            cPrev = cN;
            hPrev = hN;
            n++;
        }

        if (n >= maxIterations) {
            throw new MathIllegalStateException(LocalizedCoreFormats.NON_CONVERGENT_CONTINUED_FRACTION,
                                                maxIterations, x);
        }

        return hN;
    }

}
