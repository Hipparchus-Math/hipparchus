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
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.util.FastMath;

/** Duplication algorithm for Carlson R<sub>D</sub> elliptic integral.
 * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
 * @since 2.0
 */
class RdFieldDuplication<T extends CalculusFieldElement<T>> extends FieldDuplication<T> {

    /** Partial sum. */
    private T sum;

    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     */
    RdFieldDuplication(final T x, final T y, final T z) {
        super(x, y, z);
        sum = x.getField().getZero();
    }

    /** {@inheritDoc} */
    @Override
    protected void initialMeanPoint(final T[] va) {
        va[3] = va[0].add(va[1]).add(va[2].multiply(3)).divide(5.0);
    }

    /** {@inheritDoc} */
    @Override
    protected T convergenceCriterion(final T r, final T max) {
        return max.divide(FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(r.multiply(0.25)))));
    }

    /** {@inheritDoc} */
    @Override
    protected void update(final int m, final T[] vaM, final T[] sqrtM, final  double fourM) {

        // equation 2.29 in Carlson[1995]
        final T lambdaA = sqrtM[0].multiply(sqrtM[1]);
        final T lambdaB = sqrtM[0].multiply(sqrtM[2]);
        final T lambdaC = sqrtM[1].multiply(sqrtM[2]);

        // running sum in equation 2.34 in Carlson[1995]
        final T lambda = lambdaA.add(lambdaB).add(lambdaC);
        sum = sum.add(vaM[2].add(lambda).multiply(sqrtM[2]).multiply(fourM).reciprocal());

        // equations 2.29 and 2.30 in Carlson[1995]
        vaM[0] = vaM[0].linearCombination(0.25, vaM[0], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // xₘ
        vaM[1] = vaM[1].linearCombination(0.25, vaM[1], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // yₘ
        vaM[2] = vaM[2].linearCombination(0.25, vaM[2], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // zₘ
        vaM[3] = vaM[3].linearCombination(0.25, vaM[3], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // aₘ

    }

    /** {@inheritDoc} */
    @Override
    protected T evaluate(final T[] va0, final T aM, final  double fourM) {

        // compute symmetric differences
        final T inv   = aM.multiply(fourM).reciprocal();
        final T bigX  = va0[3].subtract(va0[0]).multiply(inv);
        final T bigY  = va0[3].subtract(va0[1]).multiply(inv);
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
        final T poly = e3e4.add(e2e5).multiply(RdRealDuplication.E3_E4_P_E2_E5).
                       add(e2e2e3.multiply(RdRealDuplication.E2_E2_E3)).
                       add(e2e4.multiply(RdRealDuplication.E2_E4)).
                       add(e3e3.multiply(RdRealDuplication.E3_E3)).
                       add(e2e2e2.multiply(RdRealDuplication.E2_E2_E2)).
                       add(e5.multiply(RdRealDuplication.E5)).
                       add(e2e3.multiply(RdRealDuplication.E2_E3)).
                       add(e4.multiply(RdRealDuplication.E4)).
                       add(e2e2.multiply(RdRealDuplication.E2_E2)).
                       add(e3.multiply(RdRealDuplication.E3)).
                       add(e2.multiply(RdRealDuplication.E2)).
                       add(RdRealDuplication.CONSTANT).
                       divide(RdRealDuplication.DENOMINATOR);
        final T polyTerm = poly.divide(aM.multiply(FastMath.sqrt(aM)).multiply(fourM));

        return polyTerm.add(sum.multiply(3));

    }

}
