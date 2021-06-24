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
package org.hipparchus.special.elliptic.jacobi;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.MathUtils;

/** Computation of Jacobi elliptic functions.
 * The Jacobi elliptic functions are related to elliptic integrals.
 * @param <T> the type of the field elements
 * @since 2.0
 */
public abstract class FieldJacobiElliptic<T extends CalculusFieldElement<T>> {

    /** Parameter of the function. */
    private final T m;

    /** Jacobi θ functions. */
    private final FieldJacobiTheta<T> jacobiTheta;

    /** Value of Jacobi θ functions at origin. */
    private final FieldTheta<T> theta0;

    /** Scaling factor. */
    private final T scaling;

    /** Simple constructor.
     * @param m parameter of the function
     */
    FieldJacobiElliptic(final T m) {

        // store parameter
        this.m = m;

        // compute nome
        final T k   = m.multiply(m);
        final T q   = LegendreEllipticIntegral.nome(k);

        // prepare underlying Jacobi θ functions
        this.jacobiTheta = new FieldJacobiTheta<>(q);
        this.theta0      = jacobiTheta.values(m.getField().getZero());
        this.scaling     = LegendreEllipticIntegral.bigK(k).reciprocal().multiply(MathUtils.SEMI_PI);

    }

    /** Get the parameter of the function.
     * @return parameter of the function
     */
    public T getM() {
        return m;
    }

    /** Get the underlying Jacobi θ functions.
     * @return underlying Jacobi θ functions
     */
    public FieldJacobiTheta<T> getJacobiThetaFunctions() {
        return jacobiTheta;
    }

    /** Evaluate the three principal Jacobi elliptic functions with pole at point n in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three principal Jacobi
     * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
     */
    public FieldCopolarN<T> valuesN(T u) {

        // evaluate Jacobi θ functions at argument
        final FieldTheta<T> thetaZ = jacobiTheta.values(u.multiply(scaling));

        // convert to Jacobi elliptic functions
        final T t02 = theta0.theta2();
        final T t03 = theta0.theta3();
        final T t04 = theta0.theta4();
        final T tz1 = thetaZ.theta1();
        final T tz2 = thetaZ.theta2();
        final T tz3 = thetaZ.theta3();
        final T tz4 = thetaZ.theta4();

        final T sn = t03.multiply(tz1).divide(t02.multiply(tz4));
        final T cn = t04.multiply(tz2).divide(t02.multiply(tz4));
        final T dn = t04.multiply(tz3).divide(t03.multiply(tz4));

        return new FieldCopolarN<>(sn, cn, dn);

    }

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

}
