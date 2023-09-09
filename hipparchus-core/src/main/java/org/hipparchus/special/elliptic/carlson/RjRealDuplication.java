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

import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/** Duplication algorithm for Carlson R<sub>J</sub> elliptic integral.
 * @since 2.0
 */
class RjRealDuplication extends RealDuplication {

    /** Delta product. */
    private double delta;

    /** sₘ iteration parameter. */
    private double sM;

    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     */
    RjRealDuplication(final double x, final double y, final double z, final double p, final double delta) {
        super(x, y, z, p);
        this.delta = delta;
    }

    /** {@inheritDoc} */
    @Override
    protected void initialMeanPoint(final double[] va) {
        va[4] = (va[0] + va[1] + va[2] + va[3] * 2) / 5.0;
    }

    /** {@inheritDoc} */
    @Override
    protected double convergenceCriterion(final double r, final double max) {
        return max / (FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(r * 0.25))));
    }

    /** {@inheritDoc} */
    @Override
    protected void update(final int m, final double[] vaM, final double[] sqrtM, final  double fourM) {
        final double dM = (sqrtM[3] + sqrtM[0]) * (sqrtM[3] + sqrtM[1]) * (sqrtM[3] + sqrtM[2]);
        if (m == 0) {
            sM = dM * 0.5;
        } else {
            // equation A.3 in Carlson[2000]
            final double rM = sM * (FastMath.sqrt(delta / (sM * sM * fourM) + 1.0) + 1.0);
            sM = (dM * rM - delta / (fourM * fourM)) / ((dM + rM / fourM) * 2);
        }

        // equation 2.19 in Carlson[1995]
        final double lambdaA = sqrtM[0] * sqrtM[1];
        final double lambdaB = sqrtM[0] * sqrtM[2];
        final double lambdaC = sqrtM[1] * sqrtM[2];

        // equations 2.19 and 2.20 in Carlson[1995]
        vaM[0] = MathArrays.linearCombination(0.25, vaM[0], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // xₘ
        vaM[1] = MathArrays.linearCombination(0.25, vaM[1], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // yₘ
        vaM[2] = MathArrays.linearCombination(0.25, vaM[2], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // zₘ
        vaM[3] = MathArrays.linearCombination(0.25, vaM[3], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // pₘ
        vaM[4] = MathArrays.linearCombination(0.25, vaM[4], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // aₘ

    }

    /** {@inheritDoc} */
    @Override
    protected double evaluate(final double[] va0, final double aM, final  double fourM) {

        // compute symmetric differences
        final double inv    = 1.0 / (aM * fourM);
        final double bigX   = (va0[4] - va0[0]) * inv;
        final double bigY   = (va0[4] - va0[1]) * inv;
        final double bigZ   = (va0[4] - va0[2]) * inv;
        final double bigP   = (bigX + bigY + bigZ) * -0.5;
        final double bigP2  = bigP * bigP;

        // compute elementary symmetric functions (we already know e1 = 0 by construction)
        final double xyz    = bigX * bigY * bigZ;
        final double e2     = bigX * (bigY + bigZ) + bigY * bigZ - bigP * bigP * 3;
        final double e3     = xyz + bigP * 2 * (e2 + bigP2 * 2);
        final double e4     = (xyz * 2 + bigP * (e2 + bigP2 * 3)) * bigP;
        final double e5     = xyz * bigP2;

        final double e2e2   = e2   * e2;
        final double e2e3   = e2   * e3;
        final double e2e4   = e2   * e4;
        final double e2e5   = e2   * e5;
        final double e3e3   = e3   * e3;
        final double e3e4   = e3   * e4;
        final double e2e2e2 = e2e2 * e2;
        final double e2e2e3 = e2e2 * e3;

        // evaluate integral using equation 19.36.1 in DLMF
        // (which add more terms than equation 2.7 in Carlson[1995])
        final double poly = ((e3e4 + e2e5) * RdRealDuplication.E3_E4_P_E2_E5 +
                              e2e2e3       * RdRealDuplication.E2_E2_E3 +
                              e2e4         * RdRealDuplication.E2_E4 +
                              e3e3         * RdRealDuplication.E3_E3 +
                              e2e2e2       * RdRealDuplication.E2_E2_E2 +
                              e5           * RdRealDuplication.E5 +
                              e2e3         * RdRealDuplication.E2_E3 +
                              e4           * RdRealDuplication.E4 +
                              e2e2         * RdRealDuplication.E2_E2 +
                              e3           * RdRealDuplication.E3 +
                              e2           * RdRealDuplication.E2 +
                              RdRealDuplication.CONSTANT) /
                             RdRealDuplication.DENOMINATOR;
        final double polyTerm = poly / (aM * FastMath.sqrt(aM) * fourM);

        // compute a single R_C term
        final double rcTerm = new RcRealDuplication(1.0, delta / (sM * sM * fourM) + 1.0).integral() * 3 / sM;

        return polyTerm + rcTerm;

    }

}
