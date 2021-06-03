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

/** Duplication algorithm for Carlson R<sub>F</sub> elliptic integral.
 * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
 * @since 2.0
 */
class RfDuplication<T extends CalculusFieldElement<T>> extends Duplication<T> {

    /** Max number of iterations in the AGM scale. */
    private static final int AGM_MAX = 16;

    /** Constant term in R<sub>F</sub> polynomial. */
    private static final double CONSTANT = 240240;

    /** Coefficient of E₂ in R<sub>F</sub> polynomial. */
    private static final double E2 = -24024;

    /** Coefficient of E₃ in R<sub>F</sub> polynomial. */
    private static final double E3 = 17160;

    /** Coefficient of E₂² in R<sub>F</sub> polynomial. */
    private static final double E2_E2 = 10010;

    /** Coefficient of E₂E₃ in R<sub>F</sub> polynomial. */
    private static final double E2_E3 = -16380;

    /** Coefficient of E₃² in R<sub>F</sub> polynomial. */
    private static final double E3_E3 = 6930;

    /** Coefficient of E₂³ in R<sub>F</sub> polynomial. */
    private static final double E2_E2_E2 = -5775;

    /** Denominator in R<sub>F</sub> polynomial. */
    private static final double DENOMINATOR = 240240;

    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     */
    RfDuplication(final T x, final T y, final T z) {
        super(x, y, z);
    }

    /** {@inheritDoc} */
    @Override
    protected T initialMeanPoint(T[] v) {
        return v[0].add(v[1]).add(v[2]).divide(3.0);
    }

    /** {@inheritDoc} */
    @Override
    protected T convergenceCriterion(final T r, final T max) {
        return max.divide(FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(r.multiply(3.0)))));
    }

    /** {@inheritDoc} */
    @Override
    protected T lambda(final int m, final T[] vM, final T[] sqrtM, final  double fourM) {
        return sqrtM[0].multiply(sqrtM[1].add(sqrtM[2])).add(sqrtM[1].multiply(sqrtM[2]));
    }

    /** {@inheritDoc} */
    @Override
    protected T evaluate(final T[] v0, final T a0, final T aM, final  double fourM) {

        // compute symmetric differences
        final T inv  = aM.multiply(fourM).reciprocal();
        final T bigX = a0.subtract(v0[0]).multiply(inv);
        final T bigY = a0.subtract(v0[1]).multiply(inv);
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
        final T poly = e2e2e2.multiply(E2_E2_E2).
                       add(e3e3.multiply(E3_E3)).
                       add(e2e3.multiply(E2_E3)).
                       add(e2e2.multiply(E2_E2)).
                       add(e3.multiply(E3)).
                       add(e2.multiply(E2)).
                       add(CONSTANT).
                       divide(DENOMINATOR);
        return poly.divide(FastMath.sqrt(aM));

    }

    /** {@inheritDoc} */
    @Override
    public T integral() {
        final T x = getVi(0);
        final T y = getVi(1);
        final T z = getVi(2);
        if (x.isZero()) {
            return completeIntegral(y, z);
        } else if (y.isZero()) {
            return completeIntegral(x, z);
        } else if (z.isZero()) {
            return completeIntegral(x, y);
        } else {
            return super.integral();
        }
    }

    /** Compute Carlson complete elliptic integral R<sub>F</sub>(u, v, 0).
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     */
    private T completeIntegral(final T x, final T y) {

        T xM = x.sqrt();
        T yM = y.sqrt();

        // iterate down
        for (int i = 1; i < AGM_MAX; ++i) {

            final T xM1 = xM;
            final T yM1 = yM;

            // arithmetic mean
            xM = xM1.add(yM1).multiply(0.5);

            // geometric mean
            yM = xM1.multiply(yM1).sqrt();

            // convergence (by the inequality of arithmetic and geometric means, this is non-negative)
            if (xM.subtract(yM).norm() <= FastMath.ulp(xM).getReal()) {
                // convergence has been reached
                return xM.add(yM).reciprocal().multiply(FastMath.PI);
            }

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

}
