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
package org.hipparchus.special.elliptic.legendre;

import java.util.function.DoubleFunction;
import java.util.function.Function;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexUnivariateIntegrator;
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.complex.FieldComplexUnivariateIntegrator;
import org.hipparchus.special.elliptic.carlson.CarlsonEllipticIntegral;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/** Complete and incomplete elliptic integrals in Legendre form.
 * <p>
 * The elliptic integrals are related to Jacobi elliptic functions.
 * </p>
 * <p>
 * <em>
 * Beware that when computing elliptic integrals in the complex plane,
 * many issues arise due to branch cuts. See the
 * <a href="https://www.hipparchus.org/hipparchus-core/special.html#Elliptic_functions_and_integrals">user guide</a>
 * for a thorough explanation.
 * </em>
 * </p>
 * <p>
 * There are different conventions to interpret the arguments of
 * Legendre elliptic integrals. In mathematical texts, these conventions show
 * up using the separator between arguments. So for example for the incomplete
 * integral of the first kind F we have:
 * </p>
 * <ul>
 *   <li>F(φ, k): the first argument φ is an angle and the second argument k
 *       is the elliptic modulus: this is the trigonometric form of the integral</li>
 *   <li>F(φ; m): the first argument φ is an angle and the second argument m=k²
 *       is the parameter: this is also a trigonometric form of the integral</li>
 *   <li>F(x|m): the first argument x=sin(φ) is not an angle anymore and the
 *       second argument m=k² is the parameter: this is the Legendre form</li>
 *   <li>F(φ\α): the first argument φ is an angle and the second argument α is the
 *       modular angle</li>
 * </ul>
 * <p>
 * As we have no separator in a method call, we have to adopt one convention
 * and stick to it. In Hipparchus, we adopted the Legendre form (i.e. F(x|m),
 * with x=sin(φ) and m=k². These conventions are consistent with Wolfram Alpha
 * functions EllipticF, EllipticE, ElliptiPI…
 * </p>
 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
 * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
 * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
 * @since 2.0
 */
public class LegendreEllipticIntegral {

    /** Private constructor for a utility class.
     */
    private LegendreEllipticIntegral() {
        // nothing to do
    }

    /** Get the nome q.
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return nome q
     */
    public static double nome(final double m) {
        if (m < 1.0e-16) {
            // first terms of infinite series in Abramowitz and Stegun 17.3.21
            final double m16 = m * 0.0625;
            return m16 * (1 + 8 * m16);
        } else {
            return FastMath.exp(-FastMath.PI * bigKPrime(m) / bigK(m));
        }
    }

    /** Get the nome q.
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return nome q
     */
    public static <T extends CalculusFieldElement<T>> T nome(final T m) {
        final T one = m.getField().getOne();
        if (m.norm() < 100 * one.ulp().getReal()) {
            // first terms of infinite series in Abramowitz and Stegun 17.3.21
            final T m16 = m.multiply(0.0625);
            return m16.multiply(m16.multiply(8).add(1));
        } else {
            return FastMath.exp(bigKPrime(m).divide(bigK(m)).multiply(one.getPi().negate()));
        }
    }

    /** Get the complete elliptic integral of the first kind K(m).
     * <p>
     * The complete elliptic integral of the first kind K(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * it corresponds to the real quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral of the first kind K(m)
     * @see #bigKPrime(double)
     * @see #bigF(double, double)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static double bigK(final double m) {
        if (m < 1.0e-8) {
            // first terms of infinite series in Abramowitz and Stegun 17.3.11
            return (1 + 0.25 * m) * MathUtils.SEMI_PI;
        } else {
            return CarlsonEllipticIntegral.rF(0, 1.0 - m, 1);
        }
    }

    /** Get the complete elliptic integral of the first kind K(m).
     * <p>
     * The complete elliptic integral of the first kind K(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * it corresponds to the real quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the first kind K(m)
     * @see #bigKPrime(CalculusFieldElement)
     * @see #bigF(CalculusFieldElement, CalculusFieldElement)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> T bigK(final T m) {
        final T zero = m.getField().getZero();
        final T one  = m.getField().getOne();
        if (m.norm() < 1.0e7 * one.ulp().getReal()) {

            // first terms of infinite series in Abramowitz and Stegun 17.3.11
            return one.add(m.multiply(0.25)).multiply(zero.getPi().multiply(0.5));

        } else {
            return CarlsonEllipticIntegral.rF(zero, one.subtract(m), one);
        }
    }

    /** Get the complete elliptic integral of the first kind K(m).
     * <p>
     * The complete elliptic integral of the first kind K(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * it corresponds to the real quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral of the first kind K(m)
     * @see #bigKPrime(Complex)
     * @see #bigF(Complex, Complex)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigK(final Complex m) {
        if (m.norm() < 1.0e-8) {
            // first terms of infinite series in Abramowitz and Stegun 17.3.11
            return Complex.ONE.add(m.multiply(0.25)).multiply(MathUtils.SEMI_PI);
        } else {
            return CarlsonEllipticIntegral.rF(Complex.ZERO, Complex.ONE.subtract(m), Complex.ONE);
        }
    }

    /** Get the complete elliptic integral of the first kind K(m).
     * <p>
     * The complete elliptic integral of the first kind K(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * it corresponds to the real quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the first kind K(m)
     * @see #bigKPrime(FieldComplex)
     * @see #bigF(FieldComplex, FieldComplex)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigK(final FieldComplex<T> m) {
        final FieldComplex<T> zero = m.getField().getZero();
        final FieldComplex<T> one  = m.getField().getOne();
        if (m.norm() < 1.0e7 * one.ulp().getReal()) {

            // first terms of infinite series in Abramowitz and Stegun 17.3.11
            return one.add(m.multiply(0.25)).multiply(zero.getPi().multiply(0.5));

        } else {
            return CarlsonEllipticIntegral.rF(zero, one.subtract(m), one);
        }
    }

    /** Get the complete elliptic integral of the first kind K'(m).
     * <p>
     * The complete elliptic integral of the first kind K'(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-(1-m) \sin^2\theta}}
     * \]
     * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral of the first kind K'(m)
     * @see #bigK(double)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static double bigKPrime(final double m) {
        return CarlsonEllipticIntegral.rF(0, m, 1);
    }

    /** Get the complete elliptic integral of the first kind K'(m).
     * <p>
     * The complete elliptic integral of the first kind K'(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-(1-m) \sin^2\theta}}
     * \]
     * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the first kind K'(m)
     * @see #bigK(CalculusFieldElement)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> T bigKPrime(final T m) {
        final T zero = m.getField().getZero();
        final T one  = m.getField().getOne();
        return CarlsonEllipticIntegral.rF(zero, m, one);
    }

    /** Get the complete elliptic integral of the first kind K'(m).
     * <p>
     * The complete elliptic integral of the first kind K'(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-(1-m) \sin^2\theta}}
     * \]
     * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral of the first kind K'(m)
     * @see #bigK(Complex)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigKPrime(final Complex m) {
        return CarlsonEllipticIntegral.rF(Complex.ZERO, m, Complex.ONE);
    }

    /** Get the complete elliptic integral of the first kind K'(m).
     * <p>
     * The complete elliptic integral of the first kind K'(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-(1-m) \sin^2\theta}}
     * \]
     * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the first kind K'(m)
     * @see #bigK(FieldComplex)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigKPrime(final FieldComplex<T> m) {
        final FieldComplex<T> zero = m.getField().getZero();
        final FieldComplex<T> one  = m.getField().getOne();
        return CarlsonEllipticIntegral.rF(zero, m, one);
    }

    /** Get the complete elliptic integral of the second kind E(m).
     * <p>
     * The complete elliptic integral of the second kind E(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral of the second kind E(m)
     * @see #bigE(double, double)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static double bigE(final double m) {
        return CarlsonEllipticIntegral.rG(0, 1 - m, 1) * 2;
    }

    /** Get the complete elliptic integral of the second kind E(m).
     * <p>
     * The complete elliptic integral of the second kind E(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the second kind E(m)
     * @see #bigE(CalculusFieldElement, CalculusFieldElement)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> T bigE(final T m) {
        final T zero = m.getField().getZero();
        final T one  = m.getField().getOne();
        return CarlsonEllipticIntegral.rG(zero, one.subtract(m), one).multiply(2);
    }

    /** Get the complete elliptic integral of the second kind E(m).
     * <p>
     * The complete elliptic integral of the second kind E(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral of the second kind E(m)
     * @see #bigE(Complex, Complex)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigE(final Complex m) {
        return CarlsonEllipticIntegral.rG(Complex.ZERO,
                                          Complex.ONE.subtract(m),
                                          Complex.ONE).multiply(2);
    }

    /** Get the complete elliptic integral of the second kind E(m).
     * <p>
     * The complete elliptic integral of the second kind E(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the second kind E(m)
     * @see #bigE(FieldComplex, FieldComplex)
     * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigE(final FieldComplex<T> m) {
        final FieldComplex<T> zero = m.getField().getZero();
        final FieldComplex<T> one  = m.getField().getOne();
        return CarlsonEllipticIntegral.rG(zero, one.subtract(m), one).multiply(2);
    }

    /** Get the complete elliptic integral D(m) = [K(m) - E(m)]/m.
     * <p>
     * The complete elliptic integral D(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{\sin^2\theta}{\sqrt{1-m \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral D(m)
     * @see #bigD(double, double)
     */
    public static double bigD(final double m) {
        return CarlsonEllipticIntegral.rD(0, 1 - m, 1) / 3;
    }

    /** Get the complete elliptic integral D(m) = [K(m) - E(m)]/m.
     * <p>
     * The complete elliptic integral D(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{\sin^2\theta}{\sqrt{1-m \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral D(m)
     * @see #bigD(CalculusFieldElement, CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigD(final T m) {
        final T zero = m.getField().getZero();
        final T one  = m.getField().getOne();
        return CarlsonEllipticIntegral.rD(zero, one.subtract(m), one).divide(3);
    }

    /** Get the complete elliptic integral D(m) = [K(m) - E(m)]/m.
     * <p>
     * The complete elliptic integral D(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{\sin^2\theta}{\sqrt{1-m \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral D(m)
     * @see #bigD(Complex, Complex)
     */
    public static Complex bigD(final Complex m) {
        return CarlsonEllipticIntegral.rD(Complex.ZERO, Complex.ONE.subtract(m), Complex.ONE).divide(3);
    }

    /** Get the complete elliptic integral D(m) = [K(m) - E(m)]/m.
     * <p>
     * The complete elliptic integral D(m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{\sin^2\theta}{\sqrt{1-m \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral D(m)
     * @see #bigD(FieldComplex, FieldComplex)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigD(final FieldComplex<T> m) {
        final FieldComplex<T> zero = m.getField().getZero();
        final FieldComplex<T> one  = m.getField().getOne();
        return CarlsonEllipticIntegral.rD(zero, one.subtract(m), one).divide(3);
    }

    /** Get the complete elliptic integral of the third kind Π(n, m).
     * <p>
     * The complete elliptic integral of the third kind Π(n, m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param n elliptic characteristic
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral of the third kind Π(n, m)
     * @see #bigPi(double, double, double)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static double bigPi(final double n, final double m) {
        final double kPrime2 = 1 - m;
        final double delta   = n * (m - n) * (n - 1);
        return CarlsonEllipticIntegral.rF(0, kPrime2, 1) +
               CarlsonEllipticIntegral.rJ(0, kPrime2, 1, 1 - n, delta) * n / 3;
    }

    /** Get the complete elliptic integral of the third kind Π(n, m).
     * <p>
     * The complete elliptic integral of the third kind Π(n, m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param n elliptic characteristic
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the third kind Π(n, m)
     * @see #bigPi(CalculusFieldElement, CalculusFieldElement, CalculusFieldElement)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> T bigPi(final T n, final T m) {
        final T zero    = m.getField().getZero();
        final T one     = m.getField().getOne();
        final T kPrime2 = one.subtract(m);
        final T delta   = n.multiply(m.subtract(n)).multiply(n.subtract(1));
        return CarlsonEllipticIntegral.rF(zero, kPrime2, one).
               add(CarlsonEllipticIntegral.rJ(zero, kPrime2, one, one.subtract(n), delta).multiply(n).divide(3));
    }

    /** Get the complete elliptic integral of the third kind Π(n, m).
     * <p>
     * The complete elliptic integral of the third kind Π(n, m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param n elliptic characteristic
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return complete elliptic integral of the third kind Π(n, m)
     * @see #bigPi(Complex, Complex, Complex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigPi(final Complex n, final Complex m) {
        final Complex kPrime2 = Complex.ONE.subtract(m);
        final Complex delta   = n.multiply(m.subtract(n)).multiply(n.subtract(1));
        return CarlsonEllipticIntegral.rF(Complex.ZERO, kPrime2, Complex.ONE).
               add(CarlsonEllipticIntegral.rJ(Complex.ZERO, kPrime2, Complex.ONE, Complex.ONE.subtract(n), delta).multiply(n).divide(3));
    }

    /** Get the complete elliptic integral of the third kind Π(n, m).
     * <p>
     * The complete elliptic integral of the third kind Π(n, m) is
     * \[
     *    \int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param n elliptic characteristic
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return complete elliptic integral of the third kind Π(n, m)
     * @see #bigPi(FieldComplex, FieldComplex, FieldComplex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigPi(final FieldComplex<T> n, final FieldComplex<T> m) {
        final FieldComplex<T> zero    = m.getField().getZero();
        final FieldComplex<T> one     = m.getField().getOne();
        final FieldComplex<T> kPrime2 = one.subtract(m);
        final FieldComplex<T> delta   = n.multiply(m.subtract(n)).multiply(n.subtract(1));
        return CarlsonEllipticIntegral.rF(zero, kPrime2, one).
               add(CarlsonEllipticIntegral.rJ(zero, kPrime2, one, one.subtract(n), delta).multiply(n).divide(3));
    }

    /** Get the incomplete elliptic integral of the first kind F(φ, m).
     * <p>
     * The incomplete elliptic integral of the first kind F(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return incomplete elliptic integral of the first kind F(φ, m)
     * @see #bigK(double)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static double bigF(final double phi, final double m) {

        // argument reduction
        final DoubleArgumentReduction ar = new DoubleArgumentReduction(phi, m, n -> bigK(n));

        // integrate part between 0 and π/2
        final double cM1 = ar.csc2 - 1.0;
        final double cMm = ar.csc2 - m;
        final double incomplete =  CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2);

        // combine complete and incomplete parts
        return ar.negate ? ar.complete - incomplete : ar.complete + incomplete;

    }

    /** Get the incomplete elliptic integral of the first kind F(φ, m).
     * <p>
     * The incomplete elliptic integral of the first kind F(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the first kind F(φ, m)
     * @see #bigK(CalculusFieldElement)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> T bigF(final T phi, final T m) {

        // argument reduction
        final FieldArgumentReduction<T> ar = new FieldArgumentReduction<>(phi, m, n -> bigK(n));

        // integrate part between 0 and π/2
        final T cM1        = ar.csc2.subtract(1);
        final T cMm        = ar.csc2.subtract(m);
        final T incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2);

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the first kind F(φ, m).
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the first kind F(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return incomplete elliptic integral of the first kind F(φ, m)
     * @see #bigK(Complex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigF(final Complex phi, final Complex m) {

        // argument reduction
        final FieldArgumentReduction<Complex> ar = new FieldArgumentReduction<>(phi, m, n -> bigK(n));

        // integrate part between 0 and π/2
        final Complex cM1        = ar.csc2.subtract(1);
        final Complex cMm        = ar.csc2.subtract(m);
        final Complex incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2);

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the first kind F(φ, m) using numerical integration.
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the first kind F(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on numerical integration.
     * If integration path comes too close to a pole of the integrand, then integration will fail
     * with a {@link org.hipparchus.exception.MathIllegalStateException MathIllegalStateException}
     * even for very large {@code maxEval}. This is normal behavior.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param integrator integrator to use
     * @param maxEval maximum number of evaluations (real and imaginary
     * parts are evaluated separately, so up to twice this number may be used)
     * @return incomplete elliptic integral of the first kind F(φ, m)
     * @see #bigK(Complex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigF(final Complex phi, final Complex m,
                               final ComplexUnivariateIntegrator integrator, final int maxEval) {
        return integrator.integrate(maxEval, new First<>(m), phi.getField().getZero(), phi);
    }

    /** Get the incomplete elliptic integral of the first kind F(φ, m).
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the first kind F(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the first kind F(φ, m)
     * @see #bigK(CalculusFieldElement)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigF(final FieldComplex<T> phi, final FieldComplex<T> m) {

        // argument reduction
        final FieldArgumentReduction<FieldComplex<T>> ar = new FieldArgumentReduction<>(phi, m, n -> bigK(n));

        // integrate part between 0 and π/2
        final FieldComplex<T> cM1        = ar.csc2.subtract(1);
        final FieldComplex<T> cMm        = ar.csc2.subtract(m);
        final FieldComplex<T> incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2);

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the first kind F(φ, m).
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the first kind F(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on numerical integration.
     * If integration path comes too close to a pole of the integrand, then integration will fail
     * with a {@link org.hipparchus.exception.MathIllegalStateException MathIllegalStateException}
     * even for very large {@code maxEval}. This is normal behavior.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param integrator integrator to use
     * @param maxEval maximum number of evaluations (real and imaginary
     * parts are evaluated separately, so up to twice this number may be used)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the first kind F(φ, m)
     * @see #bigK(CalculusFieldElement)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigF(final FieldComplex<T> phi, final FieldComplex<T> m,
                                                                           final FieldComplexUnivariateIntegrator<T> integrator,
                                                                           final int maxEval) {
        return integrator.integrate(maxEval, new First<>(m), phi.getField().getZero(), phi);
    }

    /** Get the incomplete elliptic integral of the second kind E(φ, m).
     * <p>
     * The incomplete elliptic integral of the second kind E(φ, m) is
     * \[
     *    \int_0^{\phi} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return incomplete elliptic integral of the second kind E(φ, m)
     * @see #bigE(double)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static double bigE(final double phi, final double m) {

        // argument reduction
        final DoubleArgumentReduction ar = new DoubleArgumentReduction(phi, m, n -> bigE(n));

        // integrate part between 0 and π/2
        final double cM1        = ar.csc2 - 1.0;
        final double cMm        = ar.csc2 - m;
        final double incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2) -
                                  CarlsonEllipticIntegral.rD(cM1, cMm, ar.csc2) * (m / 3);

        // combine complete and incomplete parts
        return ar.negate ? ar.complete - incomplete : ar.complete + incomplete;

    }

    /** Get the incomplete elliptic integral of the second kind E(φ, m).
     * <p>
     * The incomplete elliptic integral of the second kind E(φ, m) is
     * \[
     *    \int_0^{\phi} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the second kind E(φ, m)
     * @see #bigE(CalculusFieldElement)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> T bigE(final T phi, final T m) {

        // argument reduction
        final FieldArgumentReduction<T> ar = new FieldArgumentReduction<>(phi, m, n -> bigE(n));

        // integrate part between 0 and π/2
        final T cM1        = ar.csc2.subtract(1);
        final T cMm        = ar.csc2.subtract(m);
        final T incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2).
                             subtract(CarlsonEllipticIntegral.rD(cM1, cMm, ar.csc2).multiply(m.divide(3)));

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the second kind E(φ, m).
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the second kind E(φ, m) is
     * \[
     *    \int_0^{\phi} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return incomplete elliptic integral of the second kind E(φ, m)
     * @see #bigE(Complex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigE(final Complex phi, final Complex m) {

        // argument reduction
        final FieldArgumentReduction<Complex> ar = new FieldArgumentReduction<>(phi, m, n -> bigE(n));

        // integrate part between 0 and π/2
        final Complex cM1        = ar.csc2.subtract(1);
        final Complex cMm        = ar.csc2.subtract(m);
        final Complex incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2).
                                   subtract(CarlsonEllipticIntegral.rD(cM1, cMm, ar.csc2).multiply(m.divide(3)));

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the second kind E(φ, m) using numerical integration.
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the second kind E(φ, m) is
     * \[
     *    \int_0^{\phi} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on numerical integration.
     * If integration path comes too close to a pole of the integrand, then integration will fail
     * with a {@link org.hipparchus.exception.MathIllegalStateException MathIllegalStateException}
     * even for very large {@code maxEval}. This is normal behavior.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param integrator integrator to use
     * @param maxEval maximum number of evaluations (real and imaginary
     * parts are evaluated separately, so up to twice this number may be used)
     * @return incomplete elliptic integral of the second kind E(φ, m)
     * @see #bigE(Complex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigE(final Complex phi, final Complex m,
                               final ComplexUnivariateIntegrator integrator, final int maxEval) {
        return integrator.integrate(maxEval, new Second<>(m), phi.getField().getZero(), phi);
    }

    /** Get the incomplete elliptic integral of the second kind E(φ, m).
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the second kind E(φ, m) is
     * \[
     *    \int_0^{\phi} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the second kind E(φ, m)
     * @see #bigE(FieldComplex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigE(final FieldComplex<T> phi, final FieldComplex<T> m) {

        // argument reduction
        final FieldArgumentReduction<FieldComplex<T>> ar = new FieldArgumentReduction<>(phi, m, n -> bigE(n));

        // integrate part between 0 and π/2
        final FieldComplex<T> cM1        = ar.csc2.subtract(1);
        final FieldComplex<T> cMm        = ar.csc2.subtract(m);
        final FieldComplex<T> incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2).
                                           subtract(CarlsonEllipticIntegral.rD(cM1, cMm, ar.csc2).multiply(m.divide(3)));

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the second kind E(φ, m).
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the second kind E(φ, m) is
     * \[
     *    \int_0^{\phi} \sqrt{1-m \sin^2\theta} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on numerical integration.
     * If integration path comes too close to a pole of the integrand, then integration will fail
     * with a {@link org.hipparchus.exception.MathIllegalStateException MathIllegalStateException}
     * even for very large {@code maxEval}. This is normal behavior.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param integrator integrator to use
     * @param maxEval maximum number of evaluations (real and imaginary
     * parts are evaluated separately, so up to twice this number may be used)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the second kind E(φ, m)
     * @see #bigE(FieldComplex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigE(final FieldComplex<T> phi, final FieldComplex<T> m,
                                                                           final FieldComplexUnivariateIntegrator<T> integrator,
                                                                           final int maxEval) {
        return integrator.integrate(maxEval, new Second<>(m), phi.getField().getZero(), phi);
    }

    /** Get the incomplete elliptic integral D(φ, m) = [F(φ, m) - E(φ, m)]/m.
     * <p>
     * The incomplete elliptic integral D(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{\sin^2\theta}{\sqrt{1-m \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return incomplete elliptic integral D(φ, m)
     * @see #bigD(double)
     */
    public static double bigD(final double phi, final double m) {

        // argument reduction
        final DoubleArgumentReduction ar = new DoubleArgumentReduction(phi, m, n -> bigD(n));

        // integrate part between 0 and π/2
        final double cM1        = ar.csc2 - 1.0;
        final double cMm        = ar.csc2 - m;
        final double incomplete = CarlsonEllipticIntegral.rD(cM1, cMm, ar.csc2) / 3;

        // combine complete and incomplete parts
        return ar.negate ? ar.complete - incomplete : ar.complete + incomplete;

    }

    /** Get the incomplete elliptic integral D(φ, m) = [F(φ, m) - E(φ, m)]/m.
     * <p>
     * The incomplete elliptic integral D(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{\sin^2\theta}{\sqrt{1-m \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral D(φ, m)
     * @see #bigD(CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> T bigD(final T phi, final T m) {

        // argument reduction
        final FieldArgumentReduction<T> ar = new FieldArgumentReduction<>(phi, m, n -> bigD(n));

        // integrate part between 0 and π/2
        final T cM1        = ar.csc2.subtract(1);
        final T cMm        = ar.csc2.subtract(m);
        final T incomplete = CarlsonEllipticIntegral.rD(cM1, cMm, ar.csc2).divide(3);

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral D(φ, m) = [F(φ, m) - E(φ, m)]/m.
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral D(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{\sin^2\theta}{\sqrt{1-m \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return incomplete elliptic integral D(φ, m)
     * @see #bigD(Complex)
     */
    public static Complex bigD(final Complex phi, final Complex m) {

        // argument reduction
        final FieldArgumentReduction<Complex> ar = new FieldArgumentReduction<>(phi, m, n -> bigD(n));

        // integrate part between 0 and π/2
        final Complex cM1        = ar.csc2.subtract(1);
        final Complex cMm        = ar.csc2.subtract(m);
        final Complex incomplete = CarlsonEllipticIntegral.rD(cM1, cMm, ar.csc2).divide(3);

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral D(φ, m) = [F(φ, m) - E(φ, m)]/m.
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral D(φ, m) is
     * \[
     *    \int_0^{\phi} \frac{\sin^2\theta}{\sqrt{1-m \sin^2\theta}} d\theta
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral D(φ, m)
     * @see #bigD(CalculusFieldElement)
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigD(final FieldComplex<T> phi, final FieldComplex<T> m) {

        // argument reduction
        final FieldArgumentReduction<FieldComplex<T>> ar = new FieldArgumentReduction<>(phi, m, n -> bigD(n));

        // integrate part between 0 and π/2
        final FieldComplex<T> cM1        = ar.csc2.subtract(1);
        final FieldComplex<T> cMm        = ar.csc2.subtract(m);
        final FieldComplex<T> incomplete = CarlsonEllipticIntegral.rD(cM1, cMm, ar.csc2).divide(3);

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the third kind Π(n, φ, m).
     * <p>
     * The incomplete elliptic integral of the third kind Π(n, φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param n elliptic characteristic
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return incomplete elliptic integral of the third kind Π(n, φ, m)
     * @see #bigPi(double, double)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static double bigPi(final double n, final double phi, final double m) {

        // argument reduction
        final DoubleArgumentReduction ar = new DoubleArgumentReduction(phi, m, parameter -> bigPi(n, parameter));

        // integrate part between 0 and π/2
        final double cM1        = ar.csc2 - 1.0;
        final double cMm        = ar.csc2 - m;
        final double cMn        = ar.csc2 - n;
        final double delta      = n * (m - n) * (n - 1);
        final double incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2) +
                                  CarlsonEllipticIntegral.rJ(cM1, cMm, ar.csc2, cMn, delta) * n / 3;

        // combine complete and incomplete parts
        return ar.negate ? ar.complete - incomplete : ar.complete + incomplete;

    }

    /** Get the incomplete elliptic integral of the third kind Π(n, φ, m).
     * <p>
     * The incomplete elliptic integral of the third kind Π(n, φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param n elliptic characteristic
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the third kind Π(n, φ, m)
     * @see #bigPi(CalculusFieldElement, CalculusFieldElement)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> T bigPi(final T n, final T phi, final T m) {

        // argument reduction
        final FieldArgumentReduction<T> ar = new FieldArgumentReduction<>(phi, m, parameter -> bigPi(n, parameter));

        // integrate part between 0 and π/2
        final T cM1        = ar.csc2.subtract(1);
        final T cMm        = ar.csc2.subtract(m);
        final T cMn        = ar.csc2.subtract(n);
        final T delta      = n.multiply(m.subtract(n)).multiply(n.subtract(1));
        final T incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2).
                             add(CarlsonEllipticIntegral.rJ(cM1, cMm, ar.csc2, cMn, delta).multiply(n).divide(3));

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the third kind Π(n, φ, m).
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the third kind Π(n, φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param n elliptic characteristic
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @return incomplete elliptic integral of the third kind Π(n, φ, m)
     * @see #bigPi(Complex, Complex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigPi(final Complex n, final Complex phi, final Complex m) {

        // argument reduction
        final FieldArgumentReduction<Complex> ar = new FieldArgumentReduction<>(phi, m, parameter -> bigPi(n, parameter));

        // integrate part between 0 and π/2
        final Complex cM1        = ar.csc2.subtract(1);
        final Complex cMm        = ar.csc2.subtract(m);
        final Complex cMn        = ar.csc2.subtract(n);
        final Complex delta      = n.multiply(m.subtract(n)).multiply(n.subtract(1));
        final Complex incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2).
                                   add(CarlsonEllipticIntegral.rJ(cM1, cMm, ar.csc2, cMn, delta).multiply(n).divide(3));

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the third kind Π(n, φ, m) using numerical integration.
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the third kind Π(n, φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on numerical integration.
     * If integration path comes too close to a pole of the integrand, then integration will fail
     * with a {@link org.hipparchus.exception.MathIllegalStateException MathIllegalStateException}
     * even for very large {@code maxEval}. This is normal behavior.
     * </p>
     * @param n elliptic characteristic
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param integrator integrator to use
     * @param maxEval maximum number of evaluations (real and imaginary
     * @return incomplete elliptic integral of the third kind Π(n, φ, m)
     * @see #bigPi(Complex, Complex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static Complex bigPi(final Complex n, final Complex phi, final Complex m,
                                final ComplexUnivariateIntegrator integrator, final int maxEval) {
         return integrator.integrate(maxEval, new Third<>(n, m), phi.getField().getZero(), phi);
    }

    /** Get the incomplete elliptic integral of the third kind Π(n, φ, m).
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the third kind Π(n, φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on {@link CarlsonEllipticIntegral
     * Carlson elliptic integrals}.
     * </p>
     * @param n elliptic characteristic
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the third kind Π(n, φ, m)
     * @see #bigPi(FieldComplex, FieldComplex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigPi(final FieldComplex<T> n,
                                                                            final FieldComplex<T> phi,
                                                                            final FieldComplex<T> m) {

        // argument reduction
        final FieldArgumentReduction<FieldComplex<T>> ar = new FieldArgumentReduction<>(phi, m, parameter -> bigPi(n, parameter));

        // integrate part between 0 and π/2
        final FieldComplex<T> cM1        = ar.csc2.subtract(1);
        final FieldComplex<T> cMm        = ar.csc2.subtract(m);
        final FieldComplex<T> cMn        = ar.csc2.subtract(n);
        final FieldComplex<T> delta      = n.multiply(m.subtract(n)).multiply(n.subtract(1));
        final FieldComplex<T> incomplete = CarlsonEllipticIntegral.rF(cM1, cMm, ar.csc2).
                                           add(CarlsonEllipticIntegral.rJ(cM1, cMm, ar.csc2, cMn, delta).multiply(n).divide(3));

        // combine complete and incomplete parts
        return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);

    }

    /** Get the incomplete elliptic integral of the third kind Π(n, φ, m).
     * <p>
     * <em>
     * BEWARE! Elliptic integrals for complex numbers in the incomplete case
     * are considered experimental for now, they have known issues.
     * </em>
     * </p>
     * <p>
     * The incomplete elliptic integral of the third kind Π(n, φ, m) is
     * \[
     *    \int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}
     * \]
     * </p>
     * <p>
     * The algorithm for evaluating the functions is based on numerical integration.
     * If integration path comes too close to a pole of the integrand, then integration will fail
     * with a {@link org.hipparchus.exception.MathIllegalStateException MathIllegalStateException}
     * even for very large {@code maxEval}. This is normal behavior.
     * </p>
     * @param n elliptic characteristic
     * @param phi amplitude (i.e. upper bound of the integral)
     * @param m parameter (m=k² where k is the elliptic modulus)
     * @param integrator integrator to use
     * @param maxEval maximum number of evaluations (real and imaginary
     * parts are evaluated separately, so up to twice this number may be used)
     * @param <T> the type of the field elements
     * @return incomplete elliptic integral of the third kind Π(n, φ, m)
     * @see #bigPi(FieldComplex, FieldComplex)
     * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
     * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> bigPi(final FieldComplex<T> n,
                                                                            final FieldComplex<T> phi,
                                                                            final FieldComplex<T> m,
                                                                            final FieldComplexUnivariateIntegrator<T> integrator,
                                                                            final int maxEval) {
        return integrator.integrate(maxEval, new Third<>(n, m), phi.getField().getZero(), phi);
    }

    /** Argument reduction for an incomplete integral. */
    private static class DoubleArgumentReduction {

        /** Complete part. */
        private final double complete;

        /** Squared cosecant of the Jacobi amplitude. */
        private final double csc2;

        /** Indicator for negated Jacobi amplitude. */
        private boolean negate;

        /** Simple constructor.
         * @param phi amplitude (i.e. upper bound of the integral)
         * @param m parameter (m=k² where k is the elliptic modulus)
         * @param integral provider for complete integral
         */
        DoubleArgumentReduction(final double phi, final double m, final DoubleFunction<Double> integral) {
            final double sin = FastMath.sin(phi);
            final int    p   = (int) FastMath.rint(phi / FastMath.PI);
            complete         = p == 0 ? 0 : integral.apply(m) * 2 * p;
            negate           = sin < 0 ^ (p & 0x1) == 1;
            csc2             = 1.0 / (sin * sin);
        }

    }

    /** Argument reduction for an incomplete integral.
     * @param <T> type fo the field elements
     */
    private static class FieldArgumentReduction<T extends CalculusFieldElement<T>> {

        /** Complete part. */
        private final T complete;

        /** Squared cosecant of the Jacobi amplitude. */
        private final T csc2;

        /** Indicator for negated Jacobi amplitude. */
        private boolean negate;

        /** Simple constructor.
         * @param phi amplitude (i.e. upper bound of the integral)
         * @param m parameter (m=k² where k is the elliptic modulus)
         * @param integral provider for complete integral
         */
        FieldArgumentReduction(final T phi, final T m, final Function<T, T> integral) {
            final T   sin = FastMath.sin(phi);
            final int p   = (int) FastMath.rint(phi.getReal() / FastMath.PI);
            complete      = p == 0 ? phi.getField().getZero() : integral.apply(m).multiply(2 * p);
            negate        = sin.getReal() < 0 ^ (p & 0x1) == 1;
            csc2          = sin.multiply(sin).reciprocal();
        }

    }

    /** Integrand for elliptic integrals of the first kind.
     * @param <T> type of the field elements
     */
    private static class First<T extends CalculusFieldElement<T>> implements CalculusFieldUnivariateFunction<T> {

        /** Parameter. */
        private final T m;

        /** Simple constructor.
         * @param m parameter (m=k² where k is the elliptic modulus)
         */
        First(final T m) {
            this.m = m;
        }

        /** {@inheritDoc} */
        @Override
        public T value(final T theta) {
            final T sin  = theta.sin();
            final T sin2 = sin.multiply(sin);
            return sin2.multiply(m).negate().add(1).sqrt().reciprocal();
        }

    }

    /** Integrand for elliptic integrals of the second kind.
     * @param <T> type of the field elements
     */
    private static class Second<T extends CalculusFieldElement<T>> implements CalculusFieldUnivariateFunction<T> {

        /** Parameter. */
        private final T m;

        /** Simple constructor.
         * @param m parameter (m=k² where k is the elliptic modulus)
         */
        Second(final T m) {
            this.m = m;
        }

        /** {@inheritDoc} */
        @Override
        public T value(final T theta) {
            final T sin = theta.sin();
            final T sin2 = sin.multiply(sin);
            return sin2.multiply(m).negate().add(1).sqrt();
        }

    }

    /** Integrand for elliptic integrals of the third kind.
     * @param <T> type of the field elements
     */
    private static class Third<T extends CalculusFieldElement<T>> implements CalculusFieldUnivariateFunction<T> {

        /** Elliptic characteristic. */
        private final T n;

        /** Parameter. */
        private final T m;

        /** Simple constructor.
         * @param n elliptic characteristic
         * @param m parameter (m=k² where k is the elliptic modulus)
         */
        Third(final T n, final T m) {
            this.n = n;
            this.m = m;
        }

        /** {@inheritDoc} */
        @Override
        public T value(final T theta) {
            final T sin  = theta.sin();
            final T sin2 = sin.multiply(sin);
            final T d1   = sin2.multiply(m).negate().add(1).sqrt();
            final T da   = sin2.multiply(n).negate().add(1);
            return d1.multiply(da).reciprocal();
        }

    }
}
