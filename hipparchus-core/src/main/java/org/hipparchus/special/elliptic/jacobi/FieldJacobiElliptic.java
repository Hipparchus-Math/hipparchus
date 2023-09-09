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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.special.elliptic.carlson.CarlsonEllipticIntegral;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;

/** Computation of Jacobi elliptic functions.
 * The Jacobi elliptic functions are related to elliptic integrals.
 * @param <T> the type of the field elements
 * @since 2.0
 */
public abstract class FieldJacobiElliptic<T extends CalculusFieldElement<T>> {

    /** Parameter of the function. */
    private final T m;

    /** Simple constructor.
     * @param m parameter of the function
     */
    protected FieldJacobiElliptic(final T m) {
        this.m = m;
    }

    /** Get the parameter of the function.
     * @return parameter of the function
     */
    public T getM() {
        return m;
    }

    /** Evaluate the three principal Jacobi elliptic functions with pole at point n in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three principal Jacobi
     * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
     */
    public abstract FieldCopolarN<T> valuesN(T u);

    /** Evaluate the three principal Jacobi elliptic functions with pole at point n in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three principal Jacobi
     * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
     */
    public FieldCopolarN<T> valuesN(final double u) {
        return valuesN(m.newInstance(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point s in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
     */
    public FieldCopolarS<T> valuesS(final T u) {
        return new FieldCopolarS<>(valuesN(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point s in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
     */
    public FieldCopolarS<T> valuesS(final double u) {
        return new FieldCopolarS<>(valuesN(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point c in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
     */
    public FieldCopolarC<T> valuesC(final T u) {
        return new FieldCopolarC<>(valuesN(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point c in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
     */
    public FieldCopolarC<T> valuesC(final double u) {
        return new FieldCopolarC<>(valuesN(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point d in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
     */
    public FieldCopolarD<T> valuesD(final T u) {
        return new FieldCopolarD<>(valuesN(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point d in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
     */
    public FieldCopolarD<T> valuesD(final double u) {
        return new FieldCopolarD<>(valuesN(u));
    }

    /** Evaluate inverse of Jacobi elliptic function sn.
     * @param x value of Jacobi elliptic function {@code sn(u|m)}
     * @return u such that {@code x=sn(u|m)}
     * @since 2.1
     */
    public T arcsn(final T x) {
        // p = n, q = c, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcsp(x, x.getField().getOne().negate(), getM().negate());
    }

    /** Evaluate inverse of Jacobi elliptic function sn.
     * @param x value of Jacobi elliptic function {@code sn(u|m)}
     * @return u such that {@code x=sn(u|m)}
     * @since 2.1
     */
    public T arcsn(final double x) {
        return arcsn(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function cn.
     * @param x value of Jacobi elliptic function {@code cn(u|m)}
     * @return u such that {@code x=cn(u|m)}
     * @since 2.1
     */
    public T arccn(final T x) {
        // p = c, q = n, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpqNoDivision(x, x.getField().getOne(), getM().negate());
    }

    /** Evaluate inverse of Jacobi elliptic function cn.
     * @param x value of Jacobi elliptic function {@code cn(u|m)}
     * @return u such that {@code x=cn(u|m)}
     * @since 2.1
     */
    public T arccn(final double x) {
        return arccn(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function dn.
     * @param x value of Jacobi elliptic function {@code dn(u|m)}
     * @return u such that {@code x=dn(u|m)}
     * @since 2.1
     */
    public T arcdn(final T x) {
        // p = d, q = n, r = c, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpqNoDivision(x, getM(), x.getField().getOne().negate());
    }

    /** Evaluate inverse of Jacobi elliptic function dn.
     * @param x value of Jacobi elliptic function {@code dn(u|m)}
     * @return u such that {@code x=dn(u|m)}
     * @since 2.1
     */
    public T arcdn(final double x) {
        return arcdn(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function cs.
     * @param x value of Jacobi elliptic function {@code cs(u|m)}
     * @return u such that {@code x=cs(u|m)}
     * @since 2.1
     */
    public T arccs(final T x) {
        // p = c, q = n, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcps(x, x.getField().getOne(), getM().subtract(1).negate());
    }

    /** Evaluate inverse of Jacobi elliptic function cs.
     * @param x value of Jacobi elliptic function {@code cs(u|m)}
     * @return u such that {@code x=cs(u|m)}
     * @since 2.1
     */
    public T arccs(final double x) {
        return arccs(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function ds.
     * @param x value of Jacobi elliptic function {@code ds(u|m)}
     * @return u such that {@code x=ds(u|m)}
     * @since 2.1
     */
    public T arcds(final T x) {
        // p = d, q = c, r = n, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcps(x, getM().subtract(1), getM());
    }

    /** Evaluate inverse of Jacobi elliptic function ds.
     * @param x value of Jacobi elliptic function {@code ds(u|m)}
     * @return u such that {@code x=ds(u|m)}
     * @since 2.1
     */
    public T arcds(final double x) {
        return arcds(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function ns.
     * @param x value of Jacobi elliptic function {@code ns(u|m)}
     * @return u such that {@code x=ns(u|m)}
     * @since 2.1
     */
    public T arcns(final T x) {
        // p = n, q = c, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcps(x, x.getField().getOne().negate(), getM().negate());
    }

    /** Evaluate inverse of Jacobi elliptic function ns.
     * @param x value of Jacobi elliptic function {@code ns(u|m)}
     * @return u such that {@code x=ns(u|m)}
     * @since 2.1
     */
    public T arcns(final double x) {
        return arcns(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function dc.
     * @param x value of Jacobi elliptic function {@code dc(u|m)}
     * @return u such that {@code x=dc(u|m)}
     * @since 2.1
     */
    public T arcdc(final T x) {
        // p = d, q = c, r = n, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, getM().subtract(1), x.getField().getOne());
    }

    /** Evaluate inverse of Jacobi elliptic function dc.
     * @param x value of Jacobi elliptic function {@code dc(u|m)}
     * @return u such that {@code x=dc(u|m)}
     * @since 2.1
     */
    public T arcdc(final double x) {
        return arcdc(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function nc.
     * @param x value of Jacobi elliptic function {@code nc(u|m)}
     * @return u such that {@code x=nc(u|m)}
     * @since 2.1
     */
    public T arcnc(final T x) {
        // p = n, q = c, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, x.getField().getOne().negate(), getM().subtract(1).negate());
    }

    /** Evaluate inverse of Jacobi elliptic function nc.
     * @param x value of Jacobi elliptic function {@code nc(u|m)}
     * @return u such that {@code x=nc(u|m)}
     * @since 2.1
     */
    public T arcnc(final double x) {
        return arcnc(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function sc.
     * @param x value of Jacobi elliptic function {@code sc(u|m)}
     * @return u such that {@code x=sc(u|m)}
     * @since 2.1
     */
    public T arcsc(final T x) {
        // p = c, q = n, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcsp(x, x.getField().getOne(), getM().subtract(1).negate());
    }

    /** Evaluate inverse of Jacobi elliptic function sc.
     * @param x value of Jacobi elliptic function {@code sc(u|m)}
     * @return u such that {@code x=sc(u|m)}
     * @since 2.1
     */
    public T arcsc(final double x) {
        return arcsc(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function nd.
     * @param x value of Jacobi elliptic function {@code nd(u|m)}
     * @return u such that {@code x=nd(u|m)}
     * @since 2.1
     */
    public T arcnd(final T x) {
        // p = n, q = d, r = c, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, getM().negate(), getM().subtract(1));
    }

    /** Evaluate inverse of Jacobi elliptic function nd.
     * @param x value of Jacobi elliptic function {@code nd(u|m)}
     * @return u such that {@code x=nd(u|m)}
     * @since 2.1
     */
    public T arcnd(final double x) {
        return arcnd(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function sd.
     * @param x value of Jacobi elliptic function {@code sd(u|m)}
     * @return u such that {@code x=sd(u|m)}
     * @since 2.1
     */
    public T arcsd(final T x) {
        // p = d, q = n, r = c, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcsp(x, getM(), getM().subtract(1));
    }

    /** Evaluate inverse of Jacobi elliptic function sd.
     * @param x value of Jacobi elliptic function {@code sd(u|m)}
     * @return u such that {@code x=sd(u|m)}
     * @since 2.1
     */
    public T arcsd(final double x) {
        return arcsd(getM().getField().getZero().newInstance(x));
    }

    /** Evaluate inverse of Jacobi elliptic function cd.
     * @param x value of Jacobi elliptic function {@code cd(u|m)}
     * @return u such that {@code x=cd(u|m)}
     * @since 2.1
     */
    public T arccd(final T x) {
        // p = c, q = d, r = n, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, getM().subtract(1).negate(), getM());
    }

    /** Evaluate inverse of Jacobi elliptic function cd.
     * @param x value of Jacobi elliptic function {@code cd(u|m)}
     * @return u such that {@code x=cd(u|m)}
     * @since 2.1
     */
    public T arccd(final double x) {
        return arccd(getM().getField().getZero().newInstance(x));
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
    private T arcps(final T x, final T deltaQP, final T deltaRP) {
        // see equation 19.25.32 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E32
        final T x2       = x.multiply(x);
        final T rf       = CarlsonEllipticIntegral.rF(x2, x2.add(deltaQP), x2.add(deltaRP));
        return FastMath.copySign(1.0, rf.getReal()) * FastMath.copySign(1.0, x.getReal()) < 0 ?
               rf.negate() : rf;
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
    private T arcsp(final T x, final T deltaQP, final T deltaRP) {
        // see equation 19.25.33 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E33
        final T x2       = x.multiply(x);
        return x.multiply(CarlsonEllipticIntegral.rF(x.getField().getOne(),
                                                     deltaQP.multiply(x2).add(1),
                                                     deltaRP.multiply(x2).add(1)));
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
    private T arcpq(final T x, final T deltaQP, final T deltaRQ) {
        // see equation 19.25.34 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E34
        final T x2       = x.multiply(x);
        final T w        = x2.subtract(1).negate().divide(deltaQP);
        final T rf       = CarlsonEllipticIntegral.rF(x2, x.getField().getOne(), deltaRQ.multiply(w).add(1));
        final T positive = w.sqrt().multiply(rf);
        return x.getReal() < 0 ? LegendreEllipticIntegral.bigK(getM()).multiply(2).subtract(positive) : positive;
    }

    /** Evaluate inverse of Jacobi elliptic function pq.
     * <p>
     * Here p, q, r are any permutation of the letters c, d, n.
     * </p>
     * <p>
     * This computed the same thing as {@link #arcpq(CalculusFieldElement, CalculusFieldElement, CalculusFieldElement)}
     * but uses the homogeneity property Rf(x, y, z) = Rf(ax, ay, az) / √a to get rid of the division
     * by deltaRQ. This division induces problems in the complex case as it may lose the sign
     * of zero for values exactly along the real or imaginary axis, hence perturbing branch cuts.
     * </p>
     * @param x value of Jacobi elliptic function {@code pq(u|m)}
     * @param deltaQP Δ⁡(q, p) = q⁣s²⁡(u|m) - p⁣s²(u|m) (equation 19.5.28 of DLMF)
     * @param deltaRQ Δ⁡(r, q) = r⁣s²⁡(u|m) - q⁣s²⁡(u|m) (equation 19.5.28 of DLMF)
     * @return u such that {@code x=pq(u|m)}
     * @since 2.1
     */
    private T arcpqNoDivision(final T x, final T deltaQP, final T deltaRQ) {
        // see equation 19.25.34 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E34
        final T x2       = x.multiply(x);
        final T wDeltaQP = x2.subtract(1).negate();
        final T rf       = CarlsonEllipticIntegral.rF(x2.multiply(deltaQP), deltaQP, deltaRQ.multiply(wDeltaQP).add(deltaQP));
        final T positive = wDeltaQP.sqrt().multiply(rf);
        return FastMath.copySign(1.0, x.getReal()) < 0 ?
               LegendreEllipticIntegral.bigK(getM()).multiply(2).subtract(positive) :
               positive;
     }

}
