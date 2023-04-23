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

/** Duplication algorithm for Carlson R<sub>D</sub> elliptic integral.
 * @since 2.0
 */
class RdRealDuplication extends RealDuplication {

    /** Constant term in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double CONSTANT = 4084080;

    /** Coefficient of E₂ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E2 = -875160;

    /** Coefficient of E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E3 = 680680;

    /** Coefficient of E₂² in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E2_E2 = 417690;

    /** Coefficient of E₄ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E4 = -556920;

    /** Coefficient of E₂E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E2_E3 = -706860;

    /** Coefficient of E₅ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E5 = 471240;

    /** Coefficient of E₂³ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E2_E2_E2 = -255255;

    /** Coefficient of E₃² in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E3_E3 = 306306;

    /** Coefficient of E₂E₄ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E2_E4 = 612612;

    /** Coefficient of E₂²E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E2_E2_E3 = 675675;

    /** Coefficient of E₃E₄+E₂E₅ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double E3_E4_P_E2_E5 = -540540;

    /** Denominator in R<sub>J</sub> and R<sub>D</sub> polynomials. */
    static final double DENOMINATOR = 4084080;

    /** Partial sum. */
    private double sum;

    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     */
    RdRealDuplication(final double x, final double y, final double z) {
        super(x, y, z);
        sum = 0;
    }

    /** {@inheritDoc} */
    @Override
    protected void initialMeanPoint(final double[] va) {
        va[3] = (va[0] + va[1] + va[2] * 3.0) / 5.0;
    }

    /** {@inheritDoc} */
    @Override
    protected double convergenceCriterion(final double r, final double max) {
        return max / (FastMath.sqrt(FastMath.sqrt(FastMath.sqrt(r * 0.25))));
    }

    /** {@inheritDoc} */
    @Override
    protected void update(final int m, final double[] vaM, final double[] sqrtM, final  double fourM) {

        // equation 2.29 in Carlson[1995]
        final double lambdaA = sqrtM[0] * sqrtM[1];
        final double lambdaB = sqrtM[0] * sqrtM[2];
        final double lambdaC = sqrtM[1] * sqrtM[2];

        // running sum in equation 2.34 in Carlson[1995]
        final double lambda = lambdaA + lambdaB + lambdaC;
        sum += 1.0 / ((vaM[2] + lambda) * sqrtM[2] * fourM);

        // equations 2.29 and 2.30 in Carlson[1995]
        vaM[0] = MathArrays.linearCombination(0.25, vaM[0], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // xₘ
        vaM[1] = MathArrays.linearCombination(0.25, vaM[1], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // yₘ
        vaM[2] = MathArrays.linearCombination(0.25, vaM[2], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // zₘ
        vaM[3] = MathArrays.linearCombination(0.25, vaM[3], 0.25, lambdaA, 0.25, lambdaB, 0.25, lambdaC); // aₘ

    }

    /** {@inheritDoc} */
    @Override
    protected double evaluate(final double[] va0, final double aM, final  double fourM) {

        // compute symmetric differences
        final double inv   = 1.0 / (aM * fourM);
        final double bigX  = (va0[3] - va0[0]) * inv;
        final double bigY  = (va0[3] - va0[1]) * inv;
        final double bigZ  = (bigX + bigY) / -3;
        final double bigXY = bigX * bigY;
        final double bigZ2 = bigZ * bigZ;

        // compute elementary symmetric functions (we already know e1 = 0 by construction)
        final double e2  = bigXY - bigZ2 * 6;
        final double e3  = (bigXY * 3 - bigZ2 * 8) * bigZ;
        final double e4  = (bigXY - bigZ2) * 3 * bigZ2;
        final double e5  = bigXY * bigZ2 * bigZ;

        final double e2e2   =   e2 * e2;
        final double e2e3   =   e2 * e3;
        final double e2e4   =   e2 * e4;
        final double e2e5   =   e2 * e5;
        final double e3e3   =   e3 * e3;
        final double e3e4   =   e3 * e4;
        final double e2e2e2 = e2e2 * e2;
        final double e2e2e3 = e2e2 * e3;

        // evaluate integral using equation 19.36.1 in DLMF
        // (which add more terms than equation 2.7 in Carlson[1995])
        final double poly = ((e3e4 + e2e5) * E3_E4_P_E2_E5 +
                              e2e2e3 * E2_E2_E3 +
                              e2e4 * E2_E4 +
                              e3e3 * E3_E3 +
                              e2e2e2 * E2_E2_E2 +
                              e5 * E5 +
                              e2e3 * E2_E3 +
                              e4 * E4 +
                              e2e2 * E2_E2 +
                              e3 * E3 +
                              e2 * E2 +
                              CONSTANT) /
                             DENOMINATOR;
        final double polyTerm = poly / (aM * FastMath.sqrt(aM) * fourM);

        return polyTerm + sum * 3;

    }

}
