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

/** Duplication algorithm for Carlson R<sub>J</sub> elliptic integral.
 * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
 * @since 2.0
 */
class RjFieldDuplication<T extends CalculusFieldElement<T>> extends FieldDuplication<T> {

    /** Constant term in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double CONSTANT = 4084080;

    /** Coefficient of E₂ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E2 = -875160;

    /** Coefficient of E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E3 = 680680;

    /** Coefficient of E₂² in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E2_E2 = 417690;

    /** Coefficient of E₄ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E4 = -556920;

    /** Coefficient of E₂E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E2_E3 = -706860;

    /** Coefficient of E₅ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E5 = 471240;

    /** Coefficient of E₂³ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E2_E2_E2 = -255255;

    /** Coefficient of E₃² in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E3_E3 = 306306;

    /** Coefficient of E₂E₄ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E2_E4 = 612612;

    /** Coefficient of E₂²E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E2_E2_E3 = 675675;

    /** Coefficient of E₃E₄+E₂E₅ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double E3_E4_P_E2_E5 = -540540;

    /** Denominator in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    private static final double DENOMINATOR = 4084080;

    /** Delta product. */
    private T delta;

    /** sₘ iteration parameter. */
    private T sM;

    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     */
    RjFieldDuplication(final T x, final T y, final T z, final T p) {
        super(x, y, z, p);
        delta = p.subtract(x).multiply(p.subtract(y)).multiply(p.subtract(z));
    }

    /** {@inheritDoc} */
    @Override
    protected T initialMeanPoint(T[] v) {
        return v[0].add(v[1]).add(v[2]).add(v[3].multiply(2)).divide(5.0);
    }

    /** {@inheritDoc} */
    @Override
    protected T convergenceCriterion(final T r, final T max) {
        return max.divide(FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(r.multiply(0.25)))));
    }

    /** {@inheritDoc} */
    @Override
    protected T lambda(final int m, final T[] vM, final T[] sqrtM, final  double fourM) {
        final T dM =          sqrtM[3].add(sqrtM[0]).
                     multiply(sqrtM[3].add(sqrtM[1])).
                     multiply(sqrtM[3].add(sqrtM[2]));
        if (m == 0) {
            sM = dM.multiply(0.5);
        } else {
            // equation A.3 in Carlson[2000]
            final T rM = sM.multiply(delta.divide(sM.multiply(sM).multiply(fourM)).add(1.0).sqrt().add(1.0));
            sM = dM.multiply(rM).subtract(delta.divide(fourM * fourM)).
                 divide(dM.add(rM.divide(fourM)).multiply(2));
        }

        return sqrtM[0].multiply(sqrtM[1].add(sqrtM[2])).add(sqrtM[1].multiply(sqrtM[2]));

    }

    /** {@inheritDoc} */
    @Override
    protected T evaluate(final T[] v0, final T a0, final T aM, final  double fourM) {

        // compute symmetric differences
        final T inv    = aM.multiply(fourM).reciprocal();
        final T bigX   = a0.subtract(v0[0]).multiply(inv);
        final T bigY   = a0.subtract(v0[1]).multiply(inv);
        final T bigZ   = a0.subtract(v0[2]).multiply(inv);
        final T bigP   = bigX.add(bigY).add(bigZ).multiply(-0.5);
        final T bigP2  = bigP.multiply(bigP);

        // compute elementary symmetric functions (we already know e1 = 0 by construction)
        final T xyz    = bigX.multiply(bigY).multiply(bigZ);
        final T e2     = bigX.multiply(bigY.add(bigZ)).add(bigY.multiply(bigZ)).
                               subtract(bigP.multiply(bigP).multiply(3));
        final T e3     = xyz.add(bigP.multiply(2).multiply(e2.add(bigP2.multiply(2))));
        final T e4     = xyz.multiply(2).add(bigP.multiply(e2.add(bigP2.multiply(3)))).multiply(bigP);
        final T e5     = xyz.multiply(bigP2);

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
        final T poly = e3e4.add(e2e5).multiply(E3_E4_P_E2_E5).
                       add(e2e2e3.multiply(E2_E2_E3)).
                       add(e2e4.multiply(E2_E4)).
                       add(e3e3.multiply(E3_E3)).
                       add(e2e2e2.multiply(E2_E2_E2)).
                       add(e5.multiply(E5)).
                       add(e2e3.multiply(E2_E3)).
                       add(e4.multiply(E4)).
                       add(e2e2.multiply(E2_E2)).
                       add(e3.multiply(E3)).
                       add(e2.multiply(E2)).
                       add(CONSTANT).
                       divide(DENOMINATOR);
        final T polyTerm = poly.divide(aM.multiply(FastMath.sqrt(aM)).multiply(fourM));

        // compute a single R_C term
        final T rcTerm = new RcFieldDuplication<>(a0.getField().getOne(), delta.divide(sM.multiply(sM).multiply(fourM)).add(1)).
                         integral().
                         multiply(3).divide(sM);

        return polyTerm.add(rcTerm);

    }

}
