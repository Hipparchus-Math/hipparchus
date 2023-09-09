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
package org.hipparchus.special.elliptic.carlson;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/** Duplication algorithm for Carlson symmetric forms.
 * <p>
 * The algorithms are described in B. C. Carlson 1995 paper
 * "Numerical computation of real or complex elliptic integrals", with
 * improvements described in the appendix of B. C. Carlson and James FitzSimons
 * 2000 paper "Reduction theorems for elliptic integrands with the square root
 * of two quadratic factors". They are also described in
 * <a href="https://dlmf.nist.gov/19.36#i">section 19.36(i)</a>
 * of Digital Library of Mathematical Functions.
 * </p>
 * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
 * @since 2.0
 */
abstract class FieldDuplication<T extends CalculusFieldElement<T>> {

    /** Max number of iterations. */
    private static final int M_MAX = 16;

    /** Symmetric variables of the integral, plus mean point. */
    private final T[] initialVA;

    /** Convergence criterion. */
    private final double q;

    /** Constructor.
     * @param v symmetric variables of the integral
     */
    @SafeVarargs
    FieldDuplication(final T... v) {

        final Field<T> field = v[0].getField();
        final int n = v.length;
        initialVA = MathArrays.buildArray(field, n + 1);
        System.arraycopy(v, 0, initialVA, 0, n);
        initialMeanPoint(initialVA);

        T max = field.getZero();
        final T a0 = initialVA[n];
        for (final T vi : v) {
            max = FastMath.max(max, a0.subtract(vi).abs());
        }
        this.q = convergenceCriterion(FastMath.ulp(field.getOne()), max).getReal();

    }

    /** Get the i<sup>th</sup> symmetric variable.
     * @param i index of the variable
     * @return i<sup>th</sup> symmetric variable
     */
    protected T getVi(final int i) {
        return initialVA[i];
    }

    /** Compute initial mean point.
     * <p>
     * The initial mean point is put as the last array element
     * </>
     * @param va symmetric variables of the integral (plus placeholder for initial mean point)
     */
    protected abstract void initialMeanPoint(T[] va);

    /** Compute convergence criterion.
     * @param r relative tolerance
     * @param max max(|a0-v[i]|)
     * @return convergence criterion
     */
    protected abstract T convergenceCriterion(T r, T max);

    /** Update reduced variables in place.
     * <ul>
     *  <li>vₘ₊₁|i] ← (vₘ[i] + λₘ) / 4</li>
     *  <li>aₘ₊₁ ← (aₘ + λₘ) / 4</li>
     * </ul>
     * @param m iteration index
     * @param vaM reduced variables and mean point (updated in place)
     * @param sqrtM square roots of reduced variables
     * @param fourM 4<sup>m</sup>
     */
    protected abstract void update(int m, T[] vaM, T[] sqrtM, double fourM);

    /** Evaluate integral.
     * @param va0 initial symmetric variables and mean point of the integral
     * @param aM reduced mean point
     * @param fourM 4<sup>m</sup>
     * @return convergence criterion
     */
    protected abstract T evaluate(T[] va0, T aM, double fourM);

    /** Compute Carlson elliptic integral.
     * @return Carlson elliptic integral
     */
    public T integral() {

        // duplication iterations
        final int n     = initialVA.length - 1;
        final T[] vaM   = initialVA.clone();
        final T[] sqrtM = MathArrays.buildArray(initialVA[0].getField(), n);
        double    fourM = 1.0;
        for (int m = 0; m < M_MAX; ++m) {

            if (m > 0 && q < fourM * vaM[n].norm()) {
                // convergence reached
                break;
            }

            // apply duplication once more
            // (we know that {Field}Complex.sqrt() returns the root with nonnegative real part)
            for (int i = 0; i < n; ++i) {
                sqrtM[i] = vaM[i].sqrt();
            }
            update(m, vaM, sqrtM, fourM);

            fourM *= 4;

        }

        return evaluate(initialVA, vaM[n], fourM);

    }

}
