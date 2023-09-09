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

/** Duplication algorithm for Carlson R<sub>F</sub> elliptic integral.
 * @since 2.0
 */
class RfRealDuplication extends RealDuplication {

    /** Max number of iterations in the AGM scale. */
    static final int AGM_MAX = 32;

    /** Constant term in R<sub>F</sub> polynomial. */
    static final double CONSTANT = 240240;

    /** Coefficient of E₂ in R<sub>F</sub> polynomial. */
    static final double E2 = -24024;

    /** Coefficient of E₃ in R<sub>F</sub> polynomial. */
    static final double E3 = 17160;

    /** Coefficient of E₂² in R<sub>F</sub> polynomial. */
    static final double E2_E2 = 10010;

    /** Coefficient of E₂E₃ in R<sub>F</sub> polynomial. */
    static final double E2_E3 = -16380;

    /** Coefficient of E₃² in R<sub>F</sub> polynomial. */
    static final double E3_E3 = 6930;

    /** Coefficient of E₂³ in R<sub>F</sub> polynomial. */
    static final double E2_E2_E2 = -5775;

    /** Denominator in R<sub>F</sub> polynomial. */
    static final double DENOMINATOR = 240240;

    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     */
    RfRealDuplication(final double x, final double y, final double z) {
        super(x, y, z);
    }

    /** {@inheritDoc} */
    @Override
    protected void initialMeanPoint(final double[] va) {
        va[3] = (va[0] + va[1] + va[2]) / 3.0;
    }

    /** {@inheritDoc} */
    @Override
    protected double convergenceCriterion(final double r, final double max) {
        return max / FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(r * 3.0)));
    }

    /** {@inheritDoc} */
    @Override
    protected void update(final int m, final double[] vaM, final double[] sqrtM, final  double fourM) {

        // equation 2.3 in Carlson[1995]
        final double lambdaA = sqrtM[0] * sqrtM[1];
        final double lambdaB = sqrtM[0] * sqrtM[2];
        final double lambdaC = sqrtM[1] * sqrtM[2];

        // equations 2.3 and 2.4 in Carlson[1995]
        vaM[0] = MathArrays.linearCombination(0.25, vaM[0], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // xₘ
        vaM[1] = MathArrays.linearCombination(0.25, vaM[1], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // yₘ
        vaM[2] = MathArrays.linearCombination(0.25, vaM[2], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // zₘ
        vaM[3] = MathArrays.linearCombination(0.25, vaM[3], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // aₘ

    }

    /** {@inheritDoc} */
    @Override
    protected double evaluate(final double[] va0, final double aM, final  double fourM) {

        // compute symmetric differences
        final double inv  = 1.0 / (aM * fourM);
        final double bigX = (va0[3] - va0[0]) * inv;
        final double bigY = (va0[3] - va0[1]) * inv;
        final double bigZ = -(bigX + bigY);

        // compute elementary symmetric functions (we already know e1 = 0 by construction)
        final double e2  = bigX * bigY - bigZ * bigZ;
        final double e3  = bigX * bigY * bigZ;

        final double e2e2   =   e2 * e2;
        final double e2e3   =   e2 * e3;
        final double e3e3   =   e3 * e3;
        final double e2e2e2 = e2e2 * e2;

        // evaluate integral using equation 19.36.1 in DLMF
        // (which add more terms than equation 2.7 in Carlson[1995])
        final double poly = (e2e2e2 * E2_E2_E2 +
                             e3e3   * E3_E3 +
                             e2e3   * E2_E3 +
                             e2e2   * E2_E2 +
                             e3     * E3 +
                             e2     * E2 +
                             CONSTANT) /
                        DENOMINATOR;
        return poly / FastMath.sqrt(aM);

    }

    /** {@inheritDoc} */
    @Override
    public double integral() {
        final double x = getVi(0);
        final double y = getVi(1);
        final double z = getVi(2);
        if (x == 0) {
            return completeIntegral(y, z);
        } else if (y == 0) {
            return completeIntegral(x, z);
        } else if (z == 0) {
            return completeIntegral(x, y);
        } else {
            return super.integral();
        }
    }

    /** Compute Carlson complete elliptic integral R<sub>F</sub>(u, v, 0).
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @return Carlson complete elliptic integral R<sub>F</sub>(u, v, 0)
     */
    private double completeIntegral(final double x, final double y) {

        double xM = FastMath.sqrt(x);
        double yM = FastMath.sqrt(y);

        // iterate down
        for (int i = 1; i < AGM_MAX; ++i) {

            final double xM1 = xM;
            final double yM1 = yM;

            // arithmetic mean
            xM = (xM1 + yM1) * 0.5;

            // geometric mean
            yM = FastMath.sqrt(xM1 * yM1);

            // convergence (by the inequality of arithmetic and geometric means, this is non-negative)
            if (FastMath.abs(xM - yM) <= 4 * FastMath.ulp(xM)) {
                // convergence has been reached
                break;
            }

        }

        return FastMath.PI / (xM + yM);

    }

}
