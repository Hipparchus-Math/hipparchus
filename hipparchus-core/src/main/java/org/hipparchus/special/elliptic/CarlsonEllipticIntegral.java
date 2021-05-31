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

import org.hipparchus.complex.Complex;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;

/** Elliptic integrals in Carlson symmetric form.
 * <p>
 * This utility class computes the various symmetric elliptic
 * integrals defined as:
 * \[
 *   \begin{eqnarray}
 *   R_F(x,y,z)   = \frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)}\\
 *   R_J(x,y,z,p) = \frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)(t+p)}\\
 *   R_D(x,y,z)   = R_J(x,y,z,z)\\
 *   R_C(x,y)     = R_F(x,y,y)
 *   \end{eqnarray}
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

        // compute tolerance
        final Complex a0 = x.add(y.multiply(2)).divide(3.0);
        final double r = FastMath.ulp(1.0 / FastMath.sqrt(a0).norm());
        final double q = a0.subtract(x).norm() / FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(3 * r)));

        // duplication iterations
        Complex xM   = x;
        Complex yM   = y;
        Complex aM   = a0;
        double fourM = 1.0;
        for (int i = 0; i < DUPLICATION_MAX; ++i) {

            if (q < fourM * aM.norm()) {
                // convergence has been reached

                // compute the single polynomial independent variable
                final Complex s = y.subtract(a0).divide(aM.multiply(fourM));

                // evaluate integral using equation 2.13 in Carlson[1995]
                final Complex poly = s.multiply(RC_7).
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
            final Complex sqrtXM  = xM.sqrt();
            final Complex sqrtYM  = yM.sqrt();
            final Complex lambdaM = sqrtXM.multiply(sqrtYM).multiply(2).add(yM);

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

        // compute tolerance
        final Complex a0 = x.add(y).add(z).divide(3.0);
        final double r = FastMath.ulp(1.0 / FastMath.sqrt(a0).norm());
        final double q = FastMath.max(FastMath.max(a0.subtract(x).norm(),
                                                   a0.subtract(y).norm()),
                                      a0.subtract(z).norm()) /
                         FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(3 * r)));

        // duplication iterations
        Complex xM   = x;
        Complex yM   = y;
        Complex zM   = z;
        Complex aM   = a0;
        double fourM = 1.0;
        for (int i = 0; i < DUPLICATION_MAX; ++i) {

            if (q < fourM * aM.norm()) {
                // convergence has been reached

                // compute symmetric differences
                final Complex inv  = aM.multiply(fourM).reciprocal();
                final Complex bigX = a0.subtract(x).multiply(inv);
                final Complex bigY = a0.subtract(y).multiply(inv);
                final Complex bigZ = bigX.add(bigY).negate();

                // compute elementary symmetric functions (we already know e1 = 0 by construction)
                final Complex e2  = bigX.multiply(bigY).subtract(bigZ.multiply(bigZ));
                final Complex e3  = bigX.multiply(bigY).multiply(bigZ);

                final Complex e2e2   = e2.multiply(e2);
                final Complex e2e3   = e2.multiply(e3);
                final Complex e3e3   = e3.multiply(e3);
                final Complex e2e2e2 = e2e2.multiply(e2);

                // evaluate integral using equation 19.36.1 in DLMF
                // (which add more terms than equation 2.7 in Carlson[1995])
                final Complex poly = e2e2e2.multiply(RF_E2_E2_E2).
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
            final Complex sqrtXM  = xM.sqrt();
            final Complex sqrtYM  = yM.sqrt();
            final Complex sqrtZM  = zM.sqrt();
            final Complex lambdaM = sqrtXM.multiply(sqrtYM.add(sqrtZM)).add(sqrtYM.multiply(sqrtZM));

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
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)},\]
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     */
    public static Complex rJ(final Complex x, final Complex y, final Complex z, final Complex p) {

        // compute tolerance
        final Complex a0    = x.add(y).add(z).add(p.multiply(2)).divide(5.0);
        final Complex delta = p.subtract(x).multiply(p.subtract(y)).multiply(p.subtract(z));
        final double r      = FastMath.ulp(1.0 / FastMath.sqrt(a0).norm());
        final double q      = FastMath.max(FastMath.max(a0.subtract(x).norm(),
                                                        a0.subtract(y).norm()),
                                           FastMath.max(a0.subtract(z).norm(),
                                                        a0.subtract(p).norm())) /
                              FastMath.cbrt(FastMath.sqrt(0.25 * r));

        // duplication iterations
        Complex xM   = x;
        Complex yM   = y;
        Complex zM   = z;
        Complex pM   = p;
        Complex aM   = a0;
        double fourM = 1.0;
        for (int i = 0; i < DUPLICATION_MAX; ++i) {

            if (q < fourM * aM.norm()) {
                // convergence has been reached

                // compute symmetric differences
                final Complex inv  = aM.multiply(fourM).reciprocal();
                final Complex bigX = a0.subtract(x).multiply(inv);
                final Complex bigY = a0.subtract(y).multiply(inv);
                final Complex bigZ = bigX.add(bigY).negate();

                // compute elementary symmetric functions (we already know e1 = 0 by construction)
                final Complex e2  = bigX.multiply(bigY).subtract(bigZ.multiply(bigZ));
                final Complex e3  = bigX.multiply(bigY).multiply(bigZ);

                final Complex e2e2   = e2.multiply(e2);
                final Complex e2e3   = e2.multiply(e3);
                final Complex e3e3   = e3.multiply(e3);
                final Complex e2e2e2 = e2e2.multiply(e2);

                // evaluate integral using equation 19.36.1 in DLMF
                // (which add more terms than equation 2.7 in Carlson[1995])
                final Complex poly = e2e2e2.multiply(RF_E2_E2_E2).
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
            final Complex sqrtXM  = xM.sqrt();
            final Complex sqrtYM  = yM.sqrt();
            final Complex sqrtZM  = zM.sqrt();
            final Complex lambdaM = sqrtXM.multiply(sqrtYM.add(sqrtZM)).add(sqrtYM.multiply(sqrtZM));

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

}
