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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
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
 * @since 2.0
 */
abstract class RealDuplication {

    /** Max number of iterations. */
    private static final int M_MAX = 16;

    /** Symmetric variables of the integral. */
    private final double[] v0;

    /** Mean point. */
    protected final double a0;

    /** Convergence criterion. */
    private final double q;

    /** Constructor.
     * @param v symmetric variables of the integral
     */
    RealDuplication(final double... v) {

        this.v0 = v;
        this.a0 = initialMeanPoint(v0);

        double max = 0;
        for (final double vi : v) {
            max = FastMath.max(max, FastMath.abs(a0 - vi));
        }
        this.q = convergenceCriterion(FastMath.ulp(1.0), max);

    }

    /** Get the i<sup>th</sup> symmetric variable.
     * @param i index of the variable
     * @return i<sup>th</sup> symmetric variable
     */
    protected double getVi(final int i) {
        return v0[i];
    }

    /** Compute initial mean point.
     * @param v0 symmetric variables of the integral
     * @return initial mean point
     */
    protected abstract double initialMeanPoint(double[] v0);

    /** Compute convergence criterion.
     * @param r relative tolerance
     * @param max max(|a0-v[i]|)
     * @return convergence criterion
     */
    protected abstract double convergenceCriterion(double r, double max);

    /** Compute λₘ.
     * @param m iteration index
     * @param vM reduced variables
     * @param sqrtM square roots of reduced variables
     * @param fourM 4<sup>m</sup>
     * @return λₘ
     */
    protected abstract double lambda(int m, double[] vM, double[] sqrtM, double fourM);

    /** Evaluate integral.
     * @param v0 symmetric variables of the integral
     * @param a0 initial mean point
     * @param aM reduced mean point
     * @param fourM 4<sup>m</sup>
     * @return convergence criterion
     */
    protected abstract double evaluate(double[] v0, double a0, double aM, double fourM);

    /** Compute Carlson elliptic integral.
     * @return Carlson elliptic integral
     */
    public double integral() {

        // duplication iterations
        final double[] vM    = v0.clone();
        final double[] sqrtM = v0.clone();
        double         aM    = a0;
        double fourM = 1.0;
        for (int m = 0; m < M_MAX; ++m) {

            if (m > 0 && q < fourM * FastMath.abs(aM)) {
                return evaluate(v0, a0, aM, fourM);
            }

            // apply duplication once more
            // (we know that {Field}Complex.sqrt() returns the root with nonnegative real part)
            for (int i = 0; i < vM.length; ++i) {
                sqrtM[i] = FastMath.sqrt(vM[i]);
            }
            final double lambdaN = lambda(m, vM, sqrtM, fourM);

            // update symmetric integral variables and their mean
            for (int i = 0; i < vM.length; ++i) {
                vM[i] = (vM[i] + lambdaN) * 0.25;
            }
            aM = (aM + lambdaN) * 0.25;

            fourM *= 4;

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

}
