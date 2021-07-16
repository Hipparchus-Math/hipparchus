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
package org.hipparchus.special.elliptic.legendre;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.special.elliptic.carlson.CarlsonEllipticIntegral;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/** Complete and incomplete elliptic integrals in Legendre form.
 * <p>
 * The elliptic integrals are related to Jacobi elliptic functions.
 * </p>
 * <p>
 * Beware that {@link org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral
 * Legendre elliptic integrals} are defined in terms of elliptic modulus {@code k} whereas
 * {@link JacobiEllipticBuilder#build(double) Jacobi elliptic functions} (which
 * are their inverses) are defined in terms of parameter {@code m} and {@link
 * org.hipparchus.special.elliptic.jacobi.JacobiTheta#JacobiTheta(double) Jacobi theta
 * functions} are defined in terms of the {@link
 * org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral#nome(double)
 * nome q}. All are related as {@code k² = m} and the nome can be computed from ratios of complete
 * elliptic integrals.
 * </p>
 * @since 2.0
 */
public class LegendreEllipticIntegral {

    /** Private constructor for a utility class.
     */
    private LegendreEllipticIntegral() {
        // nothing to do
    }

    /** Get the nome q.
     * @param k elliptic modulus
     * @return nome q
     */
    public static double nome(final double k) {
        if (k < 1.0e-8) {
            // first terms of infinite series in Abramowitz and Stegun 17.3.21
            final double m16 = k * k * 0.0625;
            return m16 * (1 + 8 * m16);
        } else {
            return FastMath.exp(-FastMath.PI * bigKPrime(k) / bigK(k));
        }
    }

    /** Get the nome q.
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return nome q
     */
    public static <T extends CalculusFieldElement<T>> T nome(final T k) {
        final T one = k.getField().getOne();
        if (k.norm() < 1.0e7 * one.ulp().getReal()) {
            // first terms of infinite series in Abramowitz and Stegun 17.3.21
            final T m16 = k.multiply(k).multiply(0.0625);
            return m16.multiply(m16.multiply(8).add(1));
        } else {
            return FastMath.exp(bigKPrime(k).divide(bigK(k)).multiply(one.getPi().negate()));
        }
    }

    /** Get the complete elliptic integral of the first kind K(k).
     * <p>
     * The complete elliptic integral of the first kind K(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}}
     * \]
     * it corresponds to the real quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @return complete elliptic integral of the first kind K(k)
     * @see #bigKPrime(double)
     * @see #bigF(double, double)
     */
    public static double bigK(final double k) {
        final double m = k * k;
        if (m < 1.0e-8) {
            // first terms of infinite series in Abramowitz and Stegun 17.3.11
            return (1 + 0.25 * m) * MathUtils.SEMI_PI;
        } else {
            return CarlsonEllipticIntegral.rF(0, 1.0 - m, 1);
        }
    }

    /** Get the complete elliptic integral of the first kind K(k).
     * <p>
     * The complete elliptic integral of the first kind K(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}}
     * \]
     * it corresponds to the real quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the first kind K(k)
     * @see #bigKPrime(CalculusFieldElement)
     * @see #bigF(CalculusFieldElement, CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigK(final T k) {
        final T zero = k.getField().getZero();
        final T one  = k.getField().getOne();
        final T m    = k.multiply(k);
        if (m.norm() < 1.0e7 * one.ulp().getReal()) {

            // first terms of infinite series in Abramowitz and Stegun 17.3.11
            return one.add(m.multiply(0.25)).multiply(zero.getPi().multiply(0.5));

        } else {
            return CarlsonEllipticIntegral.rF(zero, one.subtract(m), one);
        }
    }

    /** Get the complete elliptic integral of the first kind K(k).
     * <p>
     * The complete elliptic integral of the first kind K(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}}
     * \]
     * it corresponds to the real quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @return complete elliptic integral of the first kind K(k)
     * @see #bigKPrime(Complex)
     * @see #bigF(Complex, Complex)
     */
    public static Complex bigK(final Complex k) {
        final Complex m = k.multiply(k);
        if (m.norm() < 1.0e-8) {
            // first terms of infinite series in Abramowitz and Stegun 17.3.11
            return Complex.ONE.add(m.multiply(0.25)).multiply(MathUtils.SEMI_PI);
        } else {
            return CarlsonEllipticIntegral.rF(Complex.ZERO, Complex.ONE.subtract(m), Complex.ONE);
        }
    }

    /** Get the complete elliptic integral of the first kind K(k).
     * <p>
     * The complete elliptic integral of the first kind K(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}}
     * \]
     * it corresponds to the real quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the first kind K(k)
     * @see #bigKPrime(FieldComplex)
     * @see #bigF(FieldComplex, FieldComplex)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigK(final FieldComplex<T> k) {
        final FieldComplex<T> zero = k.getField().getZero();
        final FieldComplex<T> one  = k.getField().getOne();
        final FieldComplex<T> m    = k.multiply(k);
        if (m.norm() < 1.0e7 * one.ulp().getReal()) {

            // first terms of infinite series in Abramowitz and Stegun 17.3.11
            return one.add(m.multiply(0.25)).multiply(zero.getPi().multiply(0.5));

        } else {
            return CarlsonEllipticIntegral.rF(zero, one.subtract(m), one);
        }
    }

    /** Get the complete elliptic integral of the first kind K'(k).
     * <p>
     * The complete elliptic integral of the first kind K'(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-(1-k^2) \sin^2\theta}}
     * \]
     * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @return complete elliptic integral of the first kind K'(k)
     * @see #bigK(double)
     */
    public static double bigKPrime(final double k) {
        final double m = k * k;
        return CarlsonEllipticIntegral.rF(0, m, 1);
    }

    /** Get the complete elliptic integral of the first kind K'(k).
     * <p>
     * The complete elliptic integral of the first kind K'(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-(1-k^2) \sin^2\theta}}
     * \]
     * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the first kind K'(k)
     * @see #bigK(CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigKPrime(final T k) {
        final T zero = k.getField().getZero();
        final T one  = k.getField().getOne();
        final T m    = k.multiply(k);
        return CarlsonEllipticIntegral.rF(zero, m, one);
    }

    /** Get the complete elliptic integral of the first kind K'(k).
     * <p>
     * The complete elliptic integral of the first kind K'(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-(1-k^2) \sin^2\theta}}
     * \]
     * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @return complete elliptic integral of the first kind K'(k)
     * @see #bigK(Complex)
     */
    public static Complex bigKPrime(final Complex k) {
        final Complex m = k.multiply(k);
        return CarlsonEllipticIntegral.rF(Complex.ZERO, m, Complex.ONE);
    }

    /** Get the complete elliptic integral of the first kind K'(k).
     * <p>
     * The complete elliptic integral of the first kind K'(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-(1-k^2) \sin^2\theta}}
     * \]
     * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the first kind K'(k)
     * @see #bigK(FieldComplex)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigKPrime(final FieldComplex<T> k) {
        final FieldComplex<T> zero = k.getField().getZero();
        final FieldComplex<T> one  = k.getField().getOne();
        final FieldComplex<T> m    = k.multiply(k);
        return CarlsonEllipticIntegral.rF(zero, m, one);
    }

    /** Get the complete elliptic integral of the second kind E(k).
     * <p>
     * The complete elliptic integral of the second kind E(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \sqrt{1-k^2 \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @return complete elliptic integral of the second kind E(k)
     * @see #bigE(double, double)
     */
    public static double bigE(final double k) {
        return CarlsonEllipticIntegral.rG(0, 1 - k * k, 1) * 2;
    }

    /** Get the complete elliptic integral of the second kind E(k).
     * <p>
     * The complete elliptic integral of the second kind E(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \sqrt{1-k^2 \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the second kind E(k)
     * @see #bigE(CalculusFieldElement, CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigE(final T k) {
        final T zero = k.getField().getZero();
        final T one  = k.getField().getOne();
        return CarlsonEllipticIntegral.rG(zero, one.subtract(k.multiply(k)), one).multiply(2);
    }

    /** Get the complete elliptic integral of the second kind E(k).
     * <p>
     * The complete elliptic integral of the second kind E(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \sqrt{1-k^2 \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @return complete elliptic integral of the second kind E(k)
     * @see #bigE(Complex, Complex)
     */
    public static Complex bigE(final Complex k) {
        return CarlsonEllipticIntegral.rG(Complex.ZERO,
                                          Complex.ONE.subtract(k.multiply(k)),
                                          Complex.ONE).multiply(2);
    }

    /** Get the complete elliptic integral of the second kind E(k).
     * <p>
     * The complete elliptic integral of the second kind E(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \sqrt{1-k^2 \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the second kind E(k)
     * @see #bigE(FieldComplex, FieldComplex)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigE(final FieldComplex<T> k) {
        final FieldComplex<T> zero = k.getField().getZero();
        final FieldComplex<T> one  = k.getField().getOne();
        return CarlsonEllipticIntegral.rG(zero, one.subtract(k.multiply(k)), one).multiply(2);
    }

    /** Get the complete elliptic integral D(k) = [K(k) - E(k)]/k².
     * <p>
     * The complete elliptic integral D(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{\sin^2\theta}{\sqrt{1-k^2 \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @return complete elliptic integral D(k)
     * @see #bigD(double, double)
     */
    public static double bigD(final double k) {
        return CarlsonEllipticIntegral.rD(0, 1 - k * k, 1) / 3;
    }

    /** Get the complete elliptic integral D(k) = [K(k) - E(k)]/k².
     * <p>
     * The complete elliptic integral D(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{\sin^2\theta}{\sqrt{1-k^2 \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral D(k)
     * @see #bigD(CalculusFieldElement, CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigD(final T k) {
        final T zero = k.getField().getZero();
        final T one  = k.getField().getOne();
        return CarlsonEllipticIntegral.rD(zero, one.subtract(k.multiply(k)), one).divide(3);
    }

    /** Get the complete elliptic integral D(k) = [K(k) - E(k)]/k².
     * <p>
     * The complete elliptic integral D(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{\sin^2\theta}{\sqrt{1-k^2 \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @return complete elliptic integral D(k)
     * @see #bigD(Complex, Complex)
     */
    public static Complex bigD(final Complex k) {
        return CarlsonEllipticIntegral.rD(Complex.ZERO, Complex.ONE.subtract(k.multiply(k)), Complex.ONE).divide(3);
    }

    /** Get the complete elliptic integral D(k) = [K(k) - E(k)]/k².
     * <p>
     * The complete elliptic integral D(k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{\sin^2\theta}{\sqrt{1-k^2 \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral D(k)
     * @see #bigD(FieldComplex, FieldComplex)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigD(final FieldComplex<T> k) {
        final FieldComplex<T> zero = k.getField().getZero();
        final FieldComplex<T> one  = k.getField().getOne();
        return CarlsonEllipticIntegral.rD(zero, one.subtract(k.multiply(k)), one).divide(3);
    }

    /** Get the complete elliptic integral of the third kind Π(α², k).
     * <p>
     * The complete elliptic integral of the third kind Π(α², k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}(1-\alpha^2 \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param alpha2 α² parameter (already squared)
     * @param k elliptic modulus
     * @return complete elliptic integral of the third kind Π(α², k)
     * @see #bigPi(double, double, double)
     */
    public static double bigPi(final double alpha2, final double k) {
        final double kPrime2 = 1 - k * k;
        return CarlsonEllipticIntegral.rF(0, kPrime2, 1) +
               CarlsonEllipticIntegral.rJ(0, kPrime2, 1, 1 - alpha2) * alpha2 / 3;
    }

    /** Get the complete elliptic integral of the third kind Π(α², k).
     * <p>
     * The complete elliptic integral of the third kind Π(α², k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}(1-\alpha^2 \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param alpha2 α² parameter (already squared)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the third kind Π(α², k)
     * @see #bigPi(CalculusFieldElement, CalculusFieldElement, CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigPi(final T alpha2, final T k) {
        final T zero    = k.getField().getZero();
        final T one     = k.getField().getOne();
        final T kPrime2 = one.subtract(k.multiply(k));
        return CarlsonEllipticIntegral.rF(zero, kPrime2, one).
               add(CarlsonEllipticIntegral.rJ(zero, kPrime2, one, one.subtract(alpha2)).multiply(alpha2).divide(3));
    }

    /** Get the complete elliptic integral of the third kind Π(α², k).
     * <p>
     * The complete elliptic integral of the third kind Π(α², k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}(1-\alpha^2 \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param alpha2 α² parameter (already squared)
     * @param k elliptic modulus
     * @return complete elliptic integral of the third kind Π(α², k)
     * @see #bigPi(Complex, Complex, Complex)
     */
    public static Complex bigPi(final Complex alpha2, final Complex k) {
        final Complex kPrime2 = Complex.ONE.subtract(k.multiply(k));
        return CarlsonEllipticIntegral.rF(Complex.ZERO, kPrime2, Complex.ONE).
               add(CarlsonEllipticIntegral.rJ(Complex.ZERO, kPrime2, Complex.ONE, Complex.ONE.subtract(alpha2)).multiply(alpha2).divide(3));
    }

    /** Get the complete elliptic integral of the third kind Π(α², k).
     * <p>
     * The complete elliptic integral of the third kind Π(α², k) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}(1-\alpha^2 \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param alpha2 α² parameter (already squared)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the third kind Π(α², k)
     * @see #bigPi(FieldComplex, FieldComplex, FieldComplex)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigPi(final FieldComplex<T> alpha2, final FieldComplex<T> k) {
        final FieldComplex<T> zero = k.getField().getZero();
        final FieldComplex<T> one  = k.getField().getOne();
        final FieldComplex<T> kPrime2 = one.subtract(k.multiply(k));
        return CarlsonEllipticIntegral.rF(zero, kPrime2, one).
               add(CarlsonEllipticIntegral.rJ(zero, kPrime2, one, one.subtract(alpha2)).multiply(alpha2).divide(3));
    }

    /** Get the incomplete elliptic integral of the first kind F(Φ, k).
     * <p>
     * The incomplete elliptic integral of the first kind F(Φ, k) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @return incomplete elliptic integral of the first kind K(k)
     * @see #bigK(double)
     */
    public static double bigF(final double phi, final double k) {
        final double csc  = 1.0 / FastMath.sin(phi);
        final double c    = csc * csc;
        final double cM1  = c - 1.0;
        final double cMk2 = c - k * k;
        return CarlsonEllipticIntegral.rF(cM1, cMk2, c);
    }

    /** Get the incomplete elliptic integral of the first kind F(Φ, k).
     * <p>
     * The incomplete elliptic integral of the first kind F(Φ, k) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the first kind K(k)
     * @see #bigK(CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigF(final T phi, final T k) {
        final T one  = k.getField().getOne();
        final T csc  = FastMath.sin(phi).reciprocal();
        final T c    = csc.multiply(csc);
        final T cM1  = c.subtract(one);
        final T cMk2 = c.subtract(k.multiply(k));
        return CarlsonEllipticIntegral.rF(cM1, cMk2, c);
    }

    /** Get the incomplete elliptic integral of the first kind F(Φ, k).
     * <p>
     * The incomplete elliptic integral of the first kind F(Φ, k) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @return incomplete elliptic integral of the first kind K(k)
     * @see #bigK(Complex)
     */
    public static Complex bigF(final Complex phi, final Complex k) {
        final Complex csc  = FastMath.sin(phi).reciprocal();
        final Complex c    = csc.multiply(csc);
        final Complex cM1  = c.subtract(Complex.ONE);
        final Complex cMk2 = c.subtract(k.multiply(k));
        return CarlsonEllipticIntegral.rF(cM1, cMk2, c);
    }

    /** Get the incomplete elliptic integral of the first kind F(Φ, k).
     * <p>
     * The incomplete elliptic integral of the first kind F(Φ, k) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the first kind K(k)
     * @see #bigK(CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigF(final FieldComplex<T> phi, final FieldComplex<T> k) {
        final FieldComplex<T> one  = k.getField().getOne();
        final FieldComplex<T> csc  = FastMath.sin(phi).reciprocal();
        final FieldComplex<T> c    = csc.multiply(csc);
        final FieldComplex<T> cM1  = c.subtract(one);
        final FieldComplex<T> cMk2 = c.subtract(k.multiply(k));
        return CarlsonEllipticIntegral.rF(cM1, cMk2, c);
    }

    /** Get the incomplete elliptic integral of the second kind E(Φ, k).
     * <p>
     * The incomplete elliptic integral of the second kind E(Φ, k) is
     * \[
     *    \int_0^{\phi} \sqrt{1-k^2 \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @return incomplete elliptic integral of the second kind E(Φ, k)
     * @see #bigE(double)
     */
    public static double bigE(final double phi, final double k) {
        final double csc  = 1.0 / FastMath.sin(phi);
        final double c    = csc * csc;
        final double k2   = k * k;
        final double cM1  = c - 1.0;
        final double cMk2 = c - k2;
        return CarlsonEllipticIntegral.rF(cM1, cMk2, c) -
               CarlsonEllipticIntegral.rD(cM1, cMk2, c) * (k2 / 3);
    }

    /** Get the incomplete elliptic integral of the second kind E(Φ, k).
     * <p>
     * The incomplete elliptic integral of the second kind E(Φ, k) is
     * \[
     *    \int_0^{\phi} \sqrt{1-k^2 \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the second kind E(Φ, k)
     * @see #bigE(CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigE(final T phi, final T k) {
        final T one  = k.getField().getOne();
        final T csc  = FastMath.sin(phi).reciprocal();
        final T c    = csc.multiply(csc);
        final T k2   = k.multiply(k);
        final T cM1  = c.subtract(one);
        final T cMk2 = c.subtract(k2);
        return CarlsonEllipticIntegral.rF(cM1, cMk2, c).
               subtract(CarlsonEllipticIntegral.rD(cM1, cMk2, c).multiply(k2.divide(3)));
    }

    /** Get the incomplete elliptic integral of the second kind E(Φ, k).
     * <p>
     * The incomplete elliptic integral of the second kind E(Φ, k) is
     * \[
     *    \int_0^{\phi} \sqrt{1-k^2 \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @return incomplete elliptic integral of the second kind E(Φ, k)
     * @see #bigE(Complex)
     */
    public static Complex bigE(final Complex phi, final Complex k) {
        final Complex csc  = FastMath.sin(phi).reciprocal();
        final Complex c    = csc.multiply(csc);
        final Complex k2   = k.multiply(k);
        final Complex cM1  = c.subtract(Complex.ONE);
        final Complex cMk2 = c.subtract(k2);
        return CarlsonEllipticIntegral.rF(cM1, cMk2, c).
               subtract(CarlsonEllipticIntegral.rD(cM1, cMk2, c).multiply(k2.divide(3)));
    }

    /** Get the incomplete elliptic integral of the second kind E(Φ, k).
     * <p>
     * The incomplete elliptic integral of the second kind E(Φ, k) is
     * \[
     *    \int_0^{\phi} \sqrt{1-k^2 \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the second kind E(Φ, k)
     * @see #bigE(FieldComplex)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigE(final FieldComplex<T> phi, final FieldComplex<T> k) {
        final FieldComplex<T> one  = k.getField().getOne();
        final FieldComplex<T> csc  = FastMath.sin(phi).reciprocal();
        final FieldComplex<T> c    = csc.multiply(csc);
        final FieldComplex<T> k2   = k.multiply(k);
        final FieldComplex<T> cM1  = c.subtract(one);
        final FieldComplex<T> cMk2 = c.subtract(k2);
        return CarlsonEllipticIntegral.rF(cM1, cMk2, c).
               subtract(CarlsonEllipticIntegral.rD(cM1, cMk2, c).multiply(k2.divide(3)));
    }

    /** Get the incomplete elliptic integral D(Φ, k) = [F(Φ, k) - E(Φ, k)]/k².
     * <p>
     * The incomplete elliptic integral D(Φ, k) is
     * \[
     *    \int_0^{\phi} \frac{\sin^2\theta}{\sqrt{1-k^2 \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @return incomplete elliptic integral D(Φ, k)
     * @see #bigD(double)
     */
    public static double bigD(final double phi, final double k) {
        final double csc  = 1.0 / FastMath.sin(phi);
        final double c    = csc * csc;
        final double k2   = k * k;
        final double cM1  = c - 1.0;
        final double cMk2 = c - k2;
        return CarlsonEllipticIntegral.rD(cM1, cMk2, 1) / 3;
    }

    /** Get the incomplete elliptic integral D(Φ, k) = [F(Φ, k) - E(Φ, k)]/k².
     * <p>
     * The incomplete elliptic integral D(Φ, k) is
     * \[
     *    \int_0^{\phi} \frac{\sin^2\theta}{\sqrt{1-k^2 \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral D(Φ, k)
     * @see #bigD(CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigD(final T phi, final T k) {
        final T one  = k.getField().getOne();
        final T csc  = FastMath.sin(phi).reciprocal();
        final T c    = csc.multiply(csc);
        final T k2   = k.multiply(k);
        final T cM1  = c.subtract(one);
        final T cMk2 = c.subtract(k2);
        return CarlsonEllipticIntegral.rD(cM1, cMk2, one).divide(3);
    }

    /** Get the incomplete elliptic integral D(Φ, k) = [F(Φ, k) - E(Φ, k)]/k².
     * <p>
     * The incomplete elliptic integral D(Φ, k) is
     * \[
     *    \int_0^{\phi} \frac{\sin^2\theta}{\sqrt{1-k^2 \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @return incomplete elliptic integral D(Φ, k)
     * @see #bigD(Complex)
     */
    public static Complex bigD(final Complex phi, final Complex k) {
        final Complex csc  = FastMath.sin(phi).reciprocal();
        final Complex c    = csc.multiply(csc);
        final Complex k2   = k.multiply(k);
        final Complex cM1  = c.subtract(Complex.ONE);
        final Complex cMk2 = c.subtract(k2);
        return CarlsonEllipticIntegral.rD(cM1, cMk2, Complex.ONE).divide(3);
    }

    /** Get the incomplete elliptic integral D(Φ, k) = [F(Φ, k) - E(Φ, k)]/k².
     * <p>
     * The incomplete elliptic integral D(Φ, k) is
     * \[
     *    \int_0^{\phi} \frac{\sin^2\theta}{\sqrt{1-k^2 \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral D(Φ, k)
     * @see #bigD(CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigD(final FieldComplex<T> phi, final FieldComplex<T> k) {
        final FieldComplex<T> one  = k.getField().getOne();
        final FieldComplex<T> csc  = FastMath.sin(phi).reciprocal();
        final FieldComplex<T> c    = csc.multiply(csc);
        final FieldComplex<T> k2   = k.multiply(k);
        final FieldComplex<T> cM1  = c.subtract(one);
        final FieldComplex<T> cMk2 = c.subtract(k2);
        return CarlsonEllipticIntegral.rD(cM1, cMk2, one).divide(3);
    }

    /** Get the incomplete elliptic integral of the third kind Π(Φ, α², k).
     * <p>
     * The incomplete elliptic integral of the third kind Π(Φ, α², k) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}(1-\alpha^2 \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param alpha2 α² parameter (already squared)
     * @param k elliptic modulus
     * @return incomplete elliptic integral of the third kind Π(Φ, α², k)
     * @see #bigPi(double, double)
     */
    public static double bigPi(final double phi, final double alpha2, final double k) {
        final double csc  = 1.0 / FastMath.sin(phi);
        final double c    = csc * csc;
        final double k2   = k * k;
        final double cM1  = c - 1.0;
        final double cMk2 = c - k2;
        final double cMa2 = c - alpha2;
        return bigF(phi, k) +
               CarlsonEllipticIntegral.rJ(cM1, cMk2, c, cMa2) * alpha2 / 3;
    }

    /** Get the incomplete elliptic integral of the third kind Π(Φ, α², k).
     * <p>
     * The incomplete elliptic integral of the third kind Π(Φ, α², k) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}(1-\alpha^2 \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param alpha2 α² parameter (already squared)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the third kind Π(Φ, α², k)
     * @see #bigPi(CalculusFieldElement, CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigPi(final T phi, final T alpha2, final T k) {
        final T one  = k.getField().getOne();
        final T csc  = FastMath.sin(phi).reciprocal();
        final T c    = csc.multiply(csc);
        final T k2   = k.multiply(k);
        final T cM1  = c.subtract(one);
        final T cMk2 = c.subtract(k2);
        final T cMa2 = c.subtract(alpha2);
        return bigF(phi, k).
               add(CarlsonEllipticIntegral.rJ(cM1, cMk2, c, cMa2).multiply(alpha2).divide(3));
    }

    /** Get the incomplete elliptic integral of the third kind Π(Φ, α², k).
     * <p>
     * The incomplete elliptic integral of the third kind Π(Φ, α², k) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}(1-\alpha^2 \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param alpha2 α² parameter (already squared)
     * @param k elliptic modulus
     * @return incomplete elliptic integral of the third kind Π(Φ, α², k)
     * @see #bigPi(Complex, Complex)
     */
    public static Complex bigPi(final Complex phi, final Complex alpha2, final Complex k) {
        final Complex one  = k.getField().getOne();
        final Complex csc  = FastMath.sin(phi).reciprocal();
        final Complex c    = csc.multiply(csc);
        final Complex k2   = k.multiply(k);
        final Complex cM1  = c.subtract(one);
        final Complex cMk2 = c.subtract(k2);
        final Complex cMa2 = c.subtract(alpha2);
        return bigF(phi, k).
               add(CarlsonEllipticIntegral.rJ(cM1, cMk2, c, cMa2).multiply(alpha2).divide(3));
    }

    /** Get the incomplete elliptic integral of the third kind Π(Φ, α², k).
     * <p>
     * The incomplete elliptic integral of the third kind Π(Φ, α², k) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-k^2 \sin^2\theta}(1-\alpha^2 \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param alpha2 α² parameter (already squared)
     * @param k elliptic modulus
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the third kind Π(Φ, α², k)
     * @see #bigPi(FieldComplex, FieldComplex)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigPi(final FieldComplex<T> phi,
                                                                            final FieldComplex<T> alpha2,
                                                                            final FieldComplex<T> k) {
        final FieldComplex<T> one  = k.getField().getOne();
        final FieldComplex<T> csc  = FastMath.sin(phi).reciprocal();
        final FieldComplex<T> c    = csc.multiply(csc);
        final FieldComplex<T> k2   = k.multiply(k);
        final FieldComplex<T> cM1  = c.subtract(one);
        final FieldComplex<T> cMk2 = c.subtract(k2);
        final FieldComplex<T> cMa2 = c.subtract(alpha2);
        return bigF(phi, k).
               add(CarlsonEllipticIntegral.rJ(cM1, cMk2, c, cMa2).multiply(alpha2).divide(3));
    }

}
