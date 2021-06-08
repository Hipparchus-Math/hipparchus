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

/** Duplication algorithm for Carlson R<sub>C</sub> elliptic integral.
 * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
 * @since 2.0
 */
class RcFieldDuplication<T extends CalculusFieldElement<T>> extends FieldDuplication<T> {

    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     */
    RcFieldDuplication(final T x, final T y) {
        super(x, y);
    }

    /** {@inheritDoc} */
    @Override
    protected T initialMeanPoint(final T[] v) {
        return v[0].add(v[1].multiply(2)).divide(3.0);
    }

    /** {@inheritDoc} */
    @Override
    protected T convergenceCriterion(final T r, final T max) {
        return max.divide(FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(r.multiply(3.0)))));
    }

    /** {@inheritDoc} */
    @Override
    protected T lambda(final int m, final T[] vM, final T[] sqrtM, final  double fourM) {
        return sqrtM[0].multiply(sqrtM[1]).multiply(2).add(vM[1]);
    }

    /** {@inheritDoc} */
    @Override
    protected T evaluate(final T[] v0, final T a0, final T aM, final  double fourM) {

        // compute the single polynomial independent variable
        final T s = v0[1].subtract(a0).divide(aM.multiply(fourM));

        // evaluate integral using equation 2.13 in Carlson[1995]
        final T poly = s.multiply(RcRealDuplication.S7).
                       add(RcRealDuplication.S6).multiply(s).
                       add(RcRealDuplication.S5).multiply(s).
                       add(RcRealDuplication.S4).multiply(s).
                       add(RcRealDuplication.S3).multiply(s).
                       add(RcRealDuplication.S2).multiply(s).
                       multiply(s).
                       add(RcRealDuplication.S0).
                       divide(RcRealDuplication.DENOMINATOR);
        return poly.divide(FastMath.sqrt(aM));

    }

}
