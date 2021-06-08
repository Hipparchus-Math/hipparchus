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
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.special.elliptic.jacobi.FieldJacobiTheta;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/** Elliptic integrals in Legendre form.
 * <p>
 * The elliptic integrals are related to Jacobi elliptic functions.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
public class FieldLegendreEllipticIntegral<T extends CalculusFieldElement<T>> {

    /** Max number of iterations of the AGM scale. */
    private static final int N_MAX = 16;

    /** Elliptic modulus. */
    private final T k;

    /** Simple constructor.
     * <p>
     * Beware that {@link
     * FieldLegendreEllipticIntegral#FieldLegendreEllipticIntegral(CalculusFieldElement)
     * elliptic integrals} are defined in terms of elliptic modulus {@code k} whereas
     * {@link JacobiEllipticBuilder#build(CalculusFieldElement) Jacobi elliptic
     * functions} (which are their inverses) are defined in terms of parameter
     * {@code m} and {@link FieldJacobiTheta#FieldJacobiTheta(CalculusFieldElement) Jacobi theta
     * functions} are defined in terms of the {@link FieldLegendreEllipticIntegral#getNome()
     * nome q}. All are related as {@code k² = m} and the nome can be computed
     * from ratios of complete elliptic integrals.
     * </p>
     * @param k elliptic modulus
     */
    public FieldLegendreEllipticIntegral(final T k) {
        this.k = k;
    }

    /** Get the elliptic modulus.
     * @return elliptic modulus
     */
    public T getK() {
        return k;
    }

    /** Get the nome q.
     * @return nome q
     */
    public T getNome() {
        return FastMath.exp(getBigKPrime().divide(getBigK()).multiply(- FastMath.PI));
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
    public T getBigK() {
        final T one = k.getField().getOne();
        return arithmeticGeometricMean(one, FastMath.sqrt(one.subtract(k.multiply(k)))).
               reciprocal().multiply(MathUtils.SEMI_PI);
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
    public T getBigKPrime() {
        final T one = k.getField().getOne();
        return arithmeticGeometricMean(one, k).
               reciprocal().multiply(MathUtils.SEMI_PI);
    }

    /** Compute arithmetic-geometric mean.
     * @param a first term
     * @param b second term
     * @return arithmetic-geometric mean of a and b
     */
    private T arithmeticGeometricMean(final T a, final T b) {

        T aCur = a;
        T bCur = b;

        // iterate down
        for (int i = 1; i < N_MAX; ++i) {

            final T aPrev = aCur;
            final T bPrev = bCur;

            // arithmetic mean
            aCur = aPrev.add(bPrev).multiply(0.5);

            // geometric mean
            bCur = FastMath.sqrt(aPrev.multiply(bPrev));

            // convergence (by the inequality of arithmetic and geometric means, this is non-negative)
            if (aPrev.subtract(bPrev).multiply(0.5).getReal() <= FastMath.ulp(aCur).getReal()) {
                // convergence has been reached
                return aCur;
            }

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }


}
