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
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;

/** Algorithm for computing the principal Jacobi functions for parameter m in [0; 1].
 * @param <T> the type of the field elements
 * @since 2.0
 */
class FieldBoundedParameter<T extends CalculusFieldElement<T>> extends FieldJacobiElliptic<T> {

    /** Jacobi θ functions. */
    private final FieldJacobiTheta<T> jacobiTheta;

    /** Value of Jacobi θ functions at origin. */
    private final FieldTheta<T> t0;

    /** Scaling factor. */
    private final T scaling;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function
     */
    FieldBoundedParameter(final T m) {

        super(m);

        // compute nome
        final T q   = LegendreEllipticIntegral.nome(m);

        // prepare underlying Jacobi θ functions
        this.jacobiTheta = new FieldJacobiTheta<>(q);
        this.t0          = jacobiTheta.values(m.getField().getZero());
        this.scaling     = LegendreEllipticIntegral.bigK(m).reciprocal().multiply(m.getPi().multiply(0.5));

    }

    /** {@inheritDoc}
     * <p>
     * The algorithm for evaluating the functions is based on {@link FieldJacobiTheta
     * Jacobi theta functions}.
     * </p>
     */
    @Override
    public FieldCopolarN<T> valuesN(T u) {

        // evaluate Jacobi θ functions at argument
        final FieldTheta<T> tZ = jacobiTheta.values(u.multiply(scaling));

        // convert to Jacobi elliptic functions
        final T sn = t0.theta3().multiply(tZ.theta1()).divide(t0.theta2().multiply(tZ.theta4()));
        final T cn = t0.theta4().multiply(tZ.theta2()).divide(t0.theta2().multiply(tZ.theta4()));
        final T dn = t0.theta4().multiply(tZ.theta3()).divide(t0.theta3().multiply(tZ.theta4()));

        return new FieldCopolarN<>(sn, cn, dn);

    }

}
