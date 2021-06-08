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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.special.elliptic.jacobi.JacobiTheta;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/** Elliptic integrals in Legendre form.
 * <p>
 * The elliptic integrals are related to Jacobi elliptic functions.
 * </p>
 * @since 2.0
 */
public class LegendreEllipticIntegral {

    /** Max number of iterations of the AGM scale. */
    private static final int N_MAX = 16;

    /** Elliptic modulus. */
    private final double k;

    /** Simple constructor.
     * <p>
     * Beware that {@link LegendreEllipticIntegral#LegendreEllipticIntegral(double)
     * Legendre elliptic integrals} are defined in terms of elliptic modulus {@code k}
     * whereas {@link JacobiEllipticBuilder#build(double) Jacobi elliptic functions}
     * (which are their inverses) are defined in terms of parameter {@code m} and {@link
     * JacobiTheta#JacobiTheta(double) Jacobi theta functions} are defined in terms of
     * the {@link LegendreEllipticIntegral#getNome() nome q}. All are related as {@code k² = m}
     * and the nome can be computed from ratios of complete elliptic integrals.
     * </p>
     * @param k elliptic modulus
     */
    public LegendreEllipticIntegral(final double k) {
        this.k = k;
    }

    /** Get the elliptic modulus.
     * @return elliptic modulus
     */
    public double getK() {
        return k;
    }

    /** Get the nome q.
     * @return nome q
     */
    public double getNome() {
        return FastMath.exp(-FastMath.PI * getBigKPrime() / getBigK());
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
     * The algorithm for evaluating the functions is based on arithmetic-geometric
     * mean. It is given in Abramowitz and Stegun, section 17.6.
     * </p>
     * @return complete elliptic integral of the first kind K(k)
     * @see #getBigKPrime()
     */
    public double getBigK() {
        return MathUtils.SEMI_PI / arithmeticGeometricMean(1, FastMath.sqrt(1.0 - k * k));
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
     * The algorithm for evaluating the functions is based on arithmetic-geometric
     * mean. It is given in Abramowitz and Stegun, section 17.6.
     * </p>
     * @return complete elliptic integral of the first kind K'(k)
     * @see #getBigK()
     */
    public double getBigKPrime() {
        return MathUtils.SEMI_PI / arithmeticGeometricMean(1, k);
    }

    /** Compute arithmetic-geometric mean.
     * @param a first term
     * @param b second term
     * @return arithmetic-geometric mean of a and b
     */
    private double arithmeticGeometricMean(final double a, final double b) {

        double aCur = a;
        double bCur = b;

        // iterate down
        for (int i = 1; i < N_MAX; ++i) {

            final double aPrev = aCur;
            final double bPrev = bCur;

            // arithmetic mean
            aCur = 0.5 * (aPrev + bPrev);

            // geometric mean
            bCur = FastMath.sqrt(aPrev * bPrev);

            // convergence (by the inequality of arithmetic and geometric means, this is non-negative)
            if (0.5 * (aPrev - bPrev) <= FastMath.ulp(aCur)) {
                // convergence has been reached
                return aCur;
            }

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

}
