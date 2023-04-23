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
package org.hipparchus.special.elliptic.jacobi;

import org.hipparchus.special.elliptic.carlson.CarlsonEllipticIntegral;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;

/** Algorithm computing Jacobi elliptic functions.
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

    /** Evaluate inverse of Jacobi elliptic function sn.
     * @param x value of Jacobi elliptic function {@code sn(u|m)}
     * @return u such that {@code x=sn(u|m)}
     * @since 2.1
     */
    public double arcsn(final double x) {
        // p = n, q = c, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcsp(x, -1, -getM());
    }

    /** Evaluate inverse of Jacobi elliptic function cn.
     * @param x value of Jacobi elliptic function {@code cn(u|m)}
     * @return u such that {@code x=cn(u|m)}
     * @since 2.1
     */
    public double arccn(final double x) {
        // p = c, q = n, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, 1, -getM());
    }

    /** Evaluate inverse of Jacobi elliptic function dn.
     * @param x value of Jacobi elliptic function {@code dn(u|m)}
     * @return u such that {@code x=dn(u|m)}
     * @since 2.1
     */
    public double arcdn(final double x) {
        // p = d, q = n, r = c, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, getM(), -1);
    }

    /** Evaluate inverse of Jacobi elliptic function cs.
     * @param x value of Jacobi elliptic function {@code cs(u|m)}
     * @return u such that {@code x=cs(u|m)}
     * @since 2.1
     */
    public double arccs(final double x) {
        // p = c, q = n, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcps(x, 1, 1 - getM());
    }

    /** Evaluate inverse of Jacobi elliptic function ds.
     * @param x value of Jacobi elliptic function {@code ds(u|m)}
     * @return u such that {@code x=ds(u|m)}
     * @since 2.1
     */
    public double arcds(final double x) {
        // p = d, q = c, r = n, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcps(x, getM() - 1, getM());
    }

    /** Evaluate inverse of Jacobi elliptic function ns.
     * @param x value of Jacobi elliptic function {@code ns(u|m)}
     * @return u such that {@code x=ns(u|m)}
     * @since 2.1
     */
    public double arcns(final double x) {
        // p = n, q = c, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcps(x, -1, -getM());
    }

    /** Evaluate inverse of Jacobi elliptic function dc.
     * @param x value of Jacobi elliptic function {@code dc(u|m)}
     * @return u such that {@code x=dc(u|m)}
     * @since 2.1
     */
    public double arcdc(final double x) {
        // p = d, q = c, r = n, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, getM() - 1, 1);
    }

    /** Evaluate inverse of Jacobi elliptic function nc.
     * @param x value of Jacobi elliptic function {@code nc(u|m)}
     * @return u such that {@code x=nc(u|m)}
     * @since 2.1
     */
    public double arcnc(final double x) {
        // p = n, q = c, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, -1, 1 - getM());
    }

    /** Evaluate inverse of Jacobi elliptic function sc.
     * @param x value of Jacobi elliptic function {@code sc(u|m)}
     * @return u such that {@code x=sc(u|m)}
     * @since 2.1
     */
    public double arcsc(final double x) {
        // p = c, q = n, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcsp(x, 1, 1 - getM());
    }

    /** Evaluate inverse of Jacobi elliptic function nd.
     * @param x value of Jacobi elliptic function {@code nd(u|m)}
     * @return u such that {@code x=nd(u|m)}
     * @since 2.1
     */
    public double arcnd(final double x) {
        // p = n, q = d, r = c, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, -getM(), getM() - 1);
    }

    /** Evaluate inverse of Jacobi elliptic function sd.
     * @param x value of Jacobi elliptic function {@code sd(u|m)}
     * @return u such that {@code x=sd(u|m)}
     * @since 2.1
     */
    public double arcsd(final double x) {
        // p = d, q = n, r = c, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcsp(x, getM(), getM() - 1);
    }

    /** Evaluate inverse of Jacobi elliptic function cd.
     * @param x value of Jacobi elliptic function {@code cd(u|m)}
     * @return u such that {@code x=cd(u|m)}
     * @since 2.1
     */
    public double arccd(final double x) {
        // p = c, q = d, r = n, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, 1 - getM(), getM());
    }

    /** Evaluate inverse of Jacobi elliptic function ps.
     * <p>
     * Here p, q, r are any permutation of the letters c, d, n.
     * </p>
     * @param x value of Jacobi elliptic function {@code ps(u|m)}
     * @param deltaQP Δ⁡(q, p) = q⁣s²⁡(u|m) - p⁣s²(u|m) (equation 19.5.28 of DLMF)
     * @param deltaRP Δ⁡(r, p) = r⁣s²⁡(u|m) - p⁣s²⁡(u|m) (equation 19.5.28 of DLMF)
     * @return u such that {@code x=ps(u|m)}
     * @since 2.1
     */
    private double arcps(final double x, final double deltaQP, final double deltaRP) {
        // see equation 19.25.32 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E32
        final double x2 = x * x;
        return FastMath.copySign(CarlsonEllipticIntegral.rF(x2, x2 + deltaQP, x2 + deltaRP), x);
    }

    /** Evaluate inverse of Jacobi elliptic function sp.
     * <p>
     * Here p, q, r are any permutation of the letters c, d, n.
     * </p>
     * @param x value of Jacobi elliptic function {@code sp(u|m)}
     * @param deltaQP Δ⁡(q, p) = q⁣s²⁡(u|m) - p⁣s²(u|m) (equation 19.5.28 of DLMF)
     * @param deltaRP Δ⁡(r, p) = r⁣s²⁡(u|m) - p⁣s²⁡(u|m) (equation 19.5.28 of DLMF)
     * @return u such that {@code x=sp(u|m)}
     * @since 2.1
     */
    private double arcsp(final double x, final double deltaQP, final double deltaRP) {
        // see equation 19.25.33 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E33
        final double x2 = x * x;
        return x * CarlsonEllipticIntegral.rF(1, 1 + deltaQP * x2, 1 + deltaRP * x2);
    }

    /** Evaluate inverse of Jacobi elliptic function pq.
     * <p>
     * Here p, q, r are any permutation of the letters c, d, n.
     * </p>
     * @param x value of Jacobi elliptic function {@code pq(u|m)}
     * @param deltaQP Δ⁡(q, p) = q⁣s²⁡(u|m) - p⁣s²(u|m) (equation 19.5.28 of DLMF)
     * @param deltaRQ Δ⁡(r, q) = r⁣s²⁡(u|m) - q⁣s²⁡(u|m) (equation 19.5.28 of DLMF)
     * @return u such that {@code x=pq(u|m)}
     * @since 2.1
     */
    private double arcpq(final double x, final double deltaQP, final double deltaRQ) {
        // see equation 19.25.34 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E34
        final double x2       = x * x;
        final double w        = (1 - x2) / deltaQP;
        final double positive = FastMath.sqrt(w) * CarlsonEllipticIntegral.rF(x2, 1, 1 + deltaRQ * w);
        return x < 0 ? 2 * LegendreEllipticIntegral.bigK(getM()) - positive : positive;
    }

}
