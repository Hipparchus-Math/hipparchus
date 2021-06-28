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
package org.hipparchus.special.elliptic.carlson;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.util.FastMath;

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

    /** Symmetric variables of the integral. */
    private final T[] initialV;

    /** Mean point. */
    private final T initialA;

    /** Convergence criterion. */
    private final double q;

    /** Constructor.
     * @param v symmetric variables of the integral
     */
    @SafeVarargs
    FieldDuplication(final T... v) {

        this.initialV = v;
        this.initialA = initialMeanPoint(initialV);

        T max = initialA.getField().getZero();
        for (final T vi : v) {
            max = FastMath.max(max, initialA.subtract(vi).abs());
        }
        this.q = convergenceCriterion(FastMath.ulp(initialA.getField().getOne()), max).getReal();

    }

    /** Get the i<sup>th</sup> symmetric variable.
     * @param i index of the variable
     * @return i<sup>th</sup> symmetric variable
     */
    protected T getVi(final int i) {
        return initialV[i];
    }

    /** Compute initial mean point.
     * @param v symmetric variables of the integral
     * @return initial mean point
     */
    protected abstract T initialMeanPoint(T[] v);

    /** Compute convergence criterion.
     * @param r relative tolerance
     * @param max max(|a0-v[i]|)
     * @return convergence criterion
     */
    protected abstract T convergenceCriterion(T r, T max);

    /** Compute λₘ.
     * @param m iteration index
     * @param vM reduced variables
     * @param sqrtM square roots of reduced variables
     * @param fourM 4<sup>m</sup>
     * @return λₘ
     */
    protected abstract T lambda(int m, T[] vM, T[] sqrtM, double fourM);

    /** Evaluate integral.
     * @param v0 symmetric variables of the integral
     * @param a0 initial mean point
     * @param aM reduced mean point
     * @param fourM 4<sup>m</sup>
     * @return convergence criterion
     */
    protected abstract T evaluate(T[] v0, T a0, T aM, double fourM);

    /** Compute Carlson elliptic integral.
     * @return Carlson elliptic integral
     */
    public T integral() {

        // duplication iterations
        final T[] vM    = initialV.clone();
        final T[] sqrtM = initialV.clone();
        T         aM    = initialA;
        double fourM = 1.0;
        for (int m = 0; m < M_MAX; ++m) {

            if (m > 0 && q < fourM * aM.norm()) {
                // convergence reached
                break;
            }

            // apply duplication once more
            // (we know that {Field}Complex.sqrt() returns the root with nonnegative real part)
            for (int i = 0; i < vM.length; ++i) {
                sqrtM[i] = vM[i].sqrt();
            }
            final T lambdaN = lambda(m, vM, sqrtM, fourM);

            // update symmetric integral variables and their mean
            for (int i = 0; i < vM.length; ++i) {
                vM[i] = vM[i].add(lambdaN).multiply(0.25);
            }
            aM = aM.add(lambdaN).multiply(0.25);

            fourM *= 4;

        }

        return evaluate(initialV, initialA, aM, fourM);

    }

}
