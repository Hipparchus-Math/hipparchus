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
package org.hipparchus.special.elliptic;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;

/** Elliptic integrals in Carlson symmetric form.
 * <p>
 * This utility class computes the various symmetric elliptic
 * integrals defined as:
 * \[
 *   \left\{\begin{align}
 *   R_F(x,y,z)   &= \frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)}\\
 *   R_J(x,y,z,p) &= \frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)(t+p)}\\
 *   R_G(x,y,z)   &= \frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
                     \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t\\
 *   R_D(x,y,z)   &= R_J(x,y,z,z)\\
 *   R_C(x,y)     &= R_F(x,y,y)
 *   \end{align}\right.
 * \]
 * </p>
 * <p>
 * where
 * \[
 *   s(t) = \sqrt{t+x}\sqrt{t+y}\sqrt{t+z}
 * \]
 * </p>
 * <p>
 * The algorithms used are based on the duplication method as described in
 * B. C. Carlson 1995 paper "Numerical computation of real or complex
 * elliptic integrals", with the improvements described in the appendix
 * of B. C. Carlson and James FitzSimons 2000 paper "Reduction theorems
 * for elliptic integrands with the square root of two quadratic factors".
 * They are also described in <a href="https://dlmf.nist.gov/19.36#i">section 19.36(i)</a>
 * of Digital Library of Mathematical Functions.
 * </p>
 * @since 2.0
 */
public class CarlsonEllipticIntegral {

    /** Max number of iterations of the duplication method. */
    private static final int DUPLICATION_MAX = 100;

    /** Constant term in R<sub>C</sub> polynomial. */
    private static final double RC_0 = 80080;

    /** Coefficient of s² in R<sub>C</sub> polynomial. */
    private static final double RC_2 = 24024;

    /** Coefficient of s³ in R<sub>C</sub> polynomial. */
    private static final double RC_3 = 11440;

    /** Coefficient of s⁴ in R<sub>C</sub> polynomial. */
    private static final double RC_4 = 30030;

    /** Coefficient of s⁵ in R<sub>C</sub> polynomial. */
    private static final double RC_5 = 32760;

    /** Coefficient of s⁶ in R<sub>C</sub> polynomial. */
    private static final double RC_6 = 61215;

    /** Coefficient of s⁷ in R<sub>C</sub> polynomial. */
    private static final double RC_7 = 90090;

    /** Denominator in R<sub>C</sub> polynomial. */
    private static final double RC_DENOMINATOR = 80080;

    /** Constant term in R<sub>F</sub> polynomial. */
    private static final double RF_1 = 240240;

    /** Coefficient of E₂ in R<sub>F</sub> polynomial. */
    private static final double RF_E2 = -24024;

    /** Coefficient of E₃ in R<sub>F</sub> polynomial. */
    private static final double RF_E3 = 17160;

    /** Coefficient of E₂² in R<sub>F</sub> polynomial. */
    private static final double RF_E2_E2 = 10010;

    /** Coefficient of E₂E₃ in R<sub>F</sub> polynomial. */
    private static final double RF_E2_E3 = -16380;

    /** Coefficient of E₃² in R<sub>F</sub> polynomial. */
    private static final double RF_E3_E3 = 6930;

    /** Coefficient of E₂³ in R<sub>F</sub> polynomial. */
    private static final double RF_E2_E2_E2 = -5775;

    /** Denominator in R<sub>F</sub> polynomial. */
    private static final double RF_DENOMINATOR = 240240;

    /** Constant term in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_1 = 4084080;

    /** Coefficient of E₂ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E2 = -875160;

    /** Coefficient of E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E3 = 680680;

    /** Coefficient of E₂² in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E2_E2 = 417690;

    /** Coefficient of E₄ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E4 = -556920;

    /** Coefficient of E₂E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E2_E3 = -706860;

    /** Coefficient of E₅ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E5 = 471240;

    /** Coefficient of E₂³ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E2_E2_E2 = -255255;

    /** Coefficient of E₃² in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E3_E3 = 306306;

    /** Coefficient of E₂E₄ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E2_E4 = 612612;

    /** Coefficient of E₂²E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E2_E2_E3 = 675675;

    /** Coefficient of E₃E₄+E₂E₅ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_E3_E4_P_E2_E5 = -540540;

    /** Denominator in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double RJD_DENOMINATOR = 4084080;

    /** Private constructor for a utility class.
     */
    private CarlsonEllipticIntegral() {
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     */
    public static Complex rC(final Complex x, final Complex y) {
        if (y.getImaginaryPart() == 0 && y.getRealPart() < 0) {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            final Complex xMy = x.subtract(y);
            return FastMath.sqrt(x.divide(xMy)).multiply(computeRc(xMy, y.negate()));
        } else {
            return computeRc(x, y);
        }
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param <T> type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rC(final FieldComplex<T> x, final FieldComplex<T> y) {
        if (y.getImaginaryPart().isZero() && y.getRealPart().getReal() < 0) {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            final FieldComplex<T> xMy = x.subtract(y);
            return FastMath.sqrt(x.divide(xMy)).multiply(computeRc(xMy, y.negate()));
        } else {
            return computeRc(x, y);
        }
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>, with restrictions on {@code y}.
     * <p>
     * Here {@code y} must not be real and negative.
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     */
    private static <T extends CalculusFieldElement<T>> T computeRc(final T x, final T y) {

        // compute tolerance
        final T      a0 = x.add(y.multiply(2)).divide(3.0);
        final double r = FastMath.ulp(1.0 / FastMath.sqrt(a0).norm());
        final double q = a0.subtract(x).norm() / FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(3 * r)));

        // duplication iterations
        T xM   = x;
        T yM   = y;
        T aM   = a0;
        double fourM = 1.0;
        for (int i = 0; i < DUPLICATION_MAX; ++i) {

            if (q < fourM * aM.norm()) {
                // convergence has been reached

                // compute the single polynomial independent variable
                final T s = y.subtract(a0).divide(aM.multiply(fourM));

                // evaluate integral using equation 2.13 in Carlson[1995]
                final T poly = s.multiply(RC_7).
                               add(RC_6).multiply(s).
                               add(RC_5).multiply(s).
                               add(RC_4).multiply(s).
                               add(RC_3).multiply(s).
                               add(RC_2).multiply(s).
                               multiply(s).
                               add(RC_0).
                               divide(RC_DENOMINATOR);
                return poly.divide(FastMath.sqrt(aM));

            }

            // apply duplication once more (we know that Complex.sqrt() returns the root with nonnegative real part)
            final T sqrtXM  = xM.sqrt();
            final T sqrtYM  = yM.sqrt();
            final T lambdaM = sqrtXM.multiply(sqrtYM).multiply(2).add(yM);

            // update symmetric integral variables and their mean
            xM = xM.add(lambdaM).multiply(0.25);
            yM = yM.add(lambdaM).multiply(0.25);
            aM = aM.add(lambdaM).multiply(0.25);

            fourM *= 4;

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     */
    public static Complex rF(final Complex x, final Complex y, final Complex z) {
        return computeRf(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rF(final FieldComplex<T> x, final FieldComplex<T> y, final FieldComplex<T> z) {
        return computeRf(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     */
    private static <T extends CalculusFieldElement<T>> T computeRf(final T x, final T y, final T z) {

        // compute tolerance
        final T      a0 = x.add(y).add(z).divide(3.0);
        final double r = FastMath.ulp(1.0 / FastMath.sqrt(a0).norm());
        final double q = FastMath.max(FastMath.max(a0.subtract(x).norm(),
                                                   a0.subtract(y).norm()),
                                      a0.subtract(z).norm()) /
                         FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(3 * r)));

        // duplication iterations
        T xM   = x;
        T yM   = y;
        T zM   = z;
        T aM   = a0;
        double fourM = 1.0;
        for (int i = 0; i < DUPLICATION_MAX; ++i) {

            if (q < fourM * aM.norm()) {
                // convergence has been reached

                // compute symmetric differences
                final T inv  = aM.multiply(fourM).reciprocal();
                final T bigX = a0.subtract(x).multiply(inv);
                final T bigY = a0.subtract(y).multiply(inv);
                final T bigZ = bigX.add(bigY).negate();

                // compute elementary symmetric functions (we already know e1 = 0 by construction)
                final T e2  = bigX.multiply(bigY).subtract(bigZ.multiply(bigZ));
                final T e3  = bigX.multiply(bigY).multiply(bigZ);

                final T e2e2   = e2.multiply(e2);
                final T e2e3   = e2.multiply(e3);
                final T e3e3   = e3.multiply(e3);
                final T e2e2e2 = e2e2.multiply(e2);

                // evaluate integral using equation 19.36.1 in DLMF
                // (which add more terms than equation 2.7 in Carlson[1995])
                final T poly = e2e2e2.multiply(RF_E2_E2_E2).
                               add(e3e3.multiply(RF_E3_E3)).
                               add(e2e3.multiply(RF_E2_E3)).
                               add(e2e2.multiply(RF_E2_E2)).
                               add(e3.multiply(RF_E3)).
                               add(e2.multiply(RF_E2)).
                               add(RF_1).
                               divide(RF_DENOMINATOR);
                return poly.divide(FastMath.sqrt(aM));

            }

            // apply duplication once more (we know that Complex.sqrt() returns the root with nonnegative real part)
            final T sqrtXM  = xM.sqrt();
            final T sqrtYM  = yM.sqrt();
            final T sqrtZM  = zM.sqrt();
            final T lambdaM = sqrtXM.multiply(sqrtYM.add(sqrtZM)).add(sqrtYM.multiply(sqrtZM));

            // update symmetric integral variables and their mean
            xM = xM.add(lambdaM).multiply(0.25);
            yM = yM.add(lambdaM).multiply(0.25);
            zM = zM.add(lambdaM).multiply(0.25);
            aM = aM.add(lambdaM).multiply(0.25);

            fourM *= 4;

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     */
    public static Complex rJ(final Complex x, final Complex y, final Complex z, final Complex p) {
        return computeRj(x, y, z, p);
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param <T> type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rJ(final FieldComplex<T> x, final FieldComplex<T> y,
                                                                         final FieldComplex<T> z, final FieldComplex<T> p) {
        return computeRj(x, y, z, p);
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     */
    private static <T extends CalculusFieldElement<T>> T computeRj(final T x, final T y, final T z, final T p) {

        // compute tolerance
        final T a0     = x.add(y).add(z).add(p.multiply(2)).divide(5.0);
        final T delta  = p.subtract(x).multiply(p.subtract(y)).multiply(p.subtract(z));
        final double r = FastMath.ulp(1.0 / FastMath.sqrt(a0).norm());
        final double q = FastMath.max(FastMath.max(a0.subtract(x).norm(),
                                                   a0.subtract(y).norm()),
                                      FastMath.max(a0.subtract(z).norm(),
                                                   a0.subtract(p).norm())) /
                         FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(0.25 * r)));

        // duplication iterations
        T xM     = x;
        T yM     = y;
        T zM     = z;
        T pM     = p;

        T sM     = null;
        T aM     = a0;
        double fourM   = 1.0;
        for (int i = 0; i < DUPLICATION_MAX; ++i) {

            if (i > 0 && q < fourM * aM.norm()) {
                // convergence has been reached

                // compute symmetric differences
                final T inv     = aM.multiply(fourM).reciprocal();
                final T bigX    = a0.subtract(x).multiply(inv);
                final T bigY    = a0.subtract(y).multiply(inv);
                final T bigZ    = a0.subtract(z).multiply(inv);
                final T bigP    = bigX.add(bigY).add(bigZ).multiply(-0.5);
                final T bigP2   = bigP.multiply(bigP);

                // compute elementary symmetric functions (we already know e1 = 0 by construction)
                final T xyz = bigX.multiply(bigY).multiply(bigZ);
                final T e2  = bigX.multiply(bigY.add(bigZ)).add(bigY.multiply(bigZ)).
                                    subtract(bigP.multiply(bigP).multiply(3));
                final T e3  = xyz.add(bigP.multiply(2).multiply(e2.add(bigP2.multiply(2))));
                final T e4  = xyz.multiply(2).add(bigP.multiply(e2.add(bigP2.multiply(3)))).multiply(bigP);
                final T e5  = xyz.multiply(bigP2);

                final T e2e2   = e2.multiply(e2);
                final T e2e3   = e2.multiply(e3);
                final T e2e4   = e2.multiply(e4);
                final T e2e5   = e2.multiply(e5);
                final T e3e3   = e3.multiply(e3);
                final T e3e4   = e3.multiply(e4);
                final T e2e2e2 = e2e2.multiply(e2);
                final T e2e2e3 = e2e2.multiply(e3);

                // evaluate integral using equation 19.36.1 in DLMF
                // (which add more terms than equation 2.7 in Carlson[1995])
                final T poly = e3e4.add(e2e5).multiply(RJD_E3_E4_P_E2_E5).
                                     add(e2e2e3.multiply(RJD_E2_E2_E3)).
                                     add(e2e4.multiply(RJD_E2_E4)).
                                     add(e3e3.multiply(RJD_E3_E3)).
                                     add(e2e2e2.multiply(RJD_E2_E2_E2)).
                                     add(e5.multiply(RJD_E5)).
                                     add(e2e3.multiply(RJD_E2_E3)).
                                     add(e4.multiply(RJD_E4)).
                                     add(e2e2.multiply(RJD_E2_E2)).
                                     add(e3.multiply(RJD_E3)).
                                     add(e2.multiply(RJD_E2)).
                                     add(RJD_1).
                                     divide(RJD_DENOMINATOR);
                final T polyTerm = poly.divide(aM.multiply(FastMath.sqrt(aM)).multiply(fourM));

                // compute a single R_C term
                final T rcTerm = computeRc(a0.getField().getOne(), delta.divide(sM.multiply(sM).multiply(fourM)).add(1)).
                                           multiply(3).divide(sM);

                return polyTerm.add(rcTerm);

            }

            // the following equations are from appendix of Carlson[2000],
            // which improves on Carlson[1995]
            T sqrtXM = xM.sqrt();
            T sqrtYM = yM.sqrt();
            T sqrtZM = zM.sqrt();
            T sqrtPM = pM.sqrt();

            final T dM =          sqrtPM.add(sqrtXM).
                               multiply(sqrtPM.add(sqrtYM)).
                               multiply(sqrtPM.add(sqrtZM));
            if (i == 0) {
                sM = dM.multiply(0.5);
            } else {
                // equation A.3 in Carlson[2000]
                final T rM = sM.multiply(delta.divide(sM.multiply(sM).multiply(fourM)).add(1.0).sqrt().add(1.0));
                sM = dM.multiply(rM).subtract(delta.divide(fourM * fourM)).
                     divide(dM.add(rM.divide(fourM)).multiply(2));
            }

            // apply duplication once more (we know that T.sqrt() returns the root with nonnegative real part)
            final T lambdaM = sqrtXM.multiply(sqrtYM.add(sqrtZM)).add(sqrtYM.multiply(sqrtZM));

            // update symmetric integral variables and their mean
            xM = xM.add(lambdaM).multiply(0.25);
            yM = yM.add(lambdaM).multiply(0.25);
            zM = zM.add(lambdaM).multiply(0.25);
            pM = pM.add(lambdaM).multiply(0.25);
            aM = aM.add(lambdaM).multiply(0.25);

            fourM *= 4;

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     */
    public static Complex rD(final Complex x, final Complex y, final Complex z) {
        return computeRd(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rD(final FieldComplex<T> x, final FieldComplex<T> y,
                                                                         final FieldComplex<T> z) {
        return computeRd(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     */
    private static <T extends CalculusFieldElement<T>> T computeRd(final T x, final T y, final T z) {

        // compute tolerance
        final T a0     = x.add(y).add(z.multiply(3)).divide(5.0);
        final double r = FastMath.ulp(1.0 / FastMath.sqrt(a0).norm());
        final double q = FastMath.max(FastMath.max(a0.subtract(x).norm(),
                                                   a0.subtract(y).norm()),
                                      a0.subtract(z).norm()) /
                         FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(0.25 * r)));

        // duplication iterations
        T xM   = x;
        T yM   = y;
        T zM   = z;
        T aM   = a0;
        T sum  = x.getField().getZero();
        double fourM = 1.0;
        for (int i = 0; i < DUPLICATION_MAX; ++i) {

            if (i > 0 && q < fourM * aM.norm()) {
                // convergence has been reached

                // compute symmetric differences
                final T inv   = aM.multiply(fourM).reciprocal();
                final T bigX  = a0.subtract(x).multiply(inv);
                final T bigY  = a0.subtract(y).multiply(inv);
                final T bigZ  = bigX.add(bigY).divide(-3);
                final T bigXY = bigX.multiply(bigY);
                final T bigZ2 = bigZ.multiply(bigZ);

                // compute elementary symmetric functions (we already know e1 = 0 by construction)
                final T e2  = bigXY.subtract(bigZ2.multiply(6));
                final T e3  = bigXY.multiply(3).subtract(bigZ2.multiply(8)).multiply(bigZ);
                final T e4  = bigXY.subtract(bigZ2).multiply(3).multiply(bigZ2);
                final T e5  = bigXY.multiply(bigZ2).multiply(bigZ);

                final T e2e2   = e2.multiply(e2);
                final T e2e3   = e2.multiply(e3);
                final T e2e4   = e2.multiply(e4);
                final T e2e5   = e2.multiply(e5);
                final T e3e3   = e3.multiply(e3);
                final T e3e4   = e3.multiply(e4);
                final T e2e2e2 = e2e2.multiply(e2);
                final T e2e2e3 = e2e2.multiply(e3);

                // evaluate integral using equation 19.36.1 in DLMF
                // (which add more terms than equation 2.7 in Carlson[1995])
                final T poly = e3e4.add(e2e5).multiply(RJD_E3_E4_P_E2_E5).
                                     add(e2e2e3.multiply(RJD_E2_E2_E3)).
                                     add(e2e4.multiply(RJD_E2_E4)).
                                     add(e3e3.multiply(RJD_E3_E3)).
                                     add(e2e2e2.multiply(RJD_E2_E2_E2)).
                                     add(e5.multiply(RJD_E5)).
                                     add(e2e3.multiply(RJD_E2_E3)).
                                     add(e4.multiply(RJD_E4)).
                                     add(e2e2.multiply(RJD_E2_E2)).
                                     add(e3.multiply(RJD_E3)).
                                     add(e2.multiply(RJD_E2)).
                                     add(RJD_1).
                                     divide(RJD_DENOMINATOR);
                final T polyTerm = poly.divide(aM.multiply(FastMath.sqrt(aM)).multiply(fourM));

                return polyTerm.add(sum.multiply(3));

            }

            // apply duplication once more (we know that Complex.sqrt() returns the root with nonnegative real part)
            final T sqrtXM  = xM.sqrt();
            final T sqrtYM  = yM.sqrt();
            final T sqrtZM  = zM.sqrt();
            final T lambdaM = sqrtXM.multiply(sqrtYM.add(sqrtZM)).add(sqrtYM.multiply(sqrtZM));

            sum = sum.add(zM.add(lambdaM).multiply(sqrtZM).multiply(fourM).reciprocal());

            // update symmetric integral variables and their mean
            xM  = xM.add(lambdaM).multiply(0.25);
            yM  = yM.add(lambdaM).multiply(0.25);
            zM  = zM.add(lambdaM).multiply(0.25);
            aM  = aM.add(lambdaM).multiply(0.25);

            fourM *= 4;

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     */
    public static Complex rG(final Complex x, final Complex y, final Complex z) {
        return computeRg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     * @param <T> type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rG(final FieldComplex<T> x,
                                                                         final FieldComplex<T> y,
                                                                         final FieldComplex<T> z) {
        return computeRg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * This corresponds to equation 19.21.11 in DLFM.
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     */
    private static <T extends CalculusFieldElement<T>> T computeRg(final T x, final T y, final T z) {
        return d(x, y, z).add(d(y, z, x)).add(d(z, x, y)).divide(6);
    }

    /** Compute one R<sub>D</sub> term from DLFM 19.21.11.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     */
    private static <T extends CalculusFieldElement<T>> T d(final T u, final T v, final T w) {
        return u.isZero() ? u : u.multiply(v.add(w)).multiply(computeRd(v, w, u));
    }

}
