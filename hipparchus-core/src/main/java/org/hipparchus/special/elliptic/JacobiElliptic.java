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

/** Algorithm computing Jacobi elliptic functions.
 * The Jacobi elliptic functions are related to elliptic integrals.
 * We use hare the notations from <a
 * href="https://en.wikipedia.org/wiki/Abramowitz_and_Stegun">Abramowitz and
 * Stegun</a> (Ch. 16) with a parameter {@code m}. The notations from
 * <a href="https://dlmf.nist.gov/22">Digital Library of Mathematical Functions (Ch. 22)</a>
 * are different as they use modulus {@code k} instead of parameter {@code m},
 * with {@code k² = m}.
 * @since 2.0
 */
public abstract class JacobiElliptic {

    /** Parameter of the function. */
    private final double m;

    /** Simple constructor.
     * @param m parameter of the function
     */
    protected JacobiElliptic(final double m) {
        this.m = m;
    }

    /** Get the parameter of the function.
     * @return parameter of the function
     */
    public double getM() {
        return m;
    }

    /** Evaluate the three principal Jacobi elliptic functions with pole at point n in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three principal Jacobi
     * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
     */
    public abstract CopolarN valuesN(double u);

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point s in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
     */
    public CopolarS valuesS(final double u) {
        return new CopolarS(valuesN(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point c in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
     */
    public CopolarC valuesC(final double u) {
        return new CopolarC(valuesN(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point d in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
     */
    public CopolarD valuesD(final double u) {
        return new CopolarD(valuesN(u));
    }

}
