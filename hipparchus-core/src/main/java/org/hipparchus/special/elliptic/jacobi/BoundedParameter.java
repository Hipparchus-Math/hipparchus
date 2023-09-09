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

import org.hipparchus.complex.Complex;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.MathUtils;

/** Algorithm for computing the principal Jacobi functions for parameter m in [0; 1].
 * @since 2.0
 */
class BoundedParameter extends JacobiElliptic {

    /** Jacobi θ functions. */
    private final JacobiTheta jacobiTheta;

    /** Value of Jacobi θ functions at origin. */
    private final Theta t0;

    /** Scaling factor. */
    private final double scaling;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function
     */
    BoundedParameter(final double m) {

        super(m);

        // compute nome
        final double q = LegendreEllipticIntegral.nome(m);

        // prepare underlying Jacobi θ functions
        this.jacobiTheta = new JacobiTheta(q);
        this.t0          = jacobiTheta.values(Complex.ZERO);
        this.scaling     = MathUtils.SEMI_PI / LegendreEllipticIntegral.bigK(m);

    }

    /** {@inheritDoc}
     * <p>
     * The algorithm for evaluating the functions is based on {@link JacobiTheta
     * Jacobi theta functions}.
     * </p>
     */
    @Override
    public CopolarN valuesN(double u) {

        // evaluate Jacobi θ functions at argument
        final Theta tZ = jacobiTheta.values(new Complex(u * scaling));

        // convert to Jacobi elliptic functions
        final double sn = t0.theta3().multiply(tZ.theta1()).divide(t0.theta2().multiply(tZ.theta4())).getRealPart();
        final double cn = t0.theta4().multiply(tZ.theta2()).divide(t0.theta2().multiply(tZ.theta4())).getRealPart();
        final double dn = t0.theta4().multiply(tZ.theta3()).divide(t0.theta3().multiply(tZ.theta4())).getRealPart();

        return new CopolarN(sn, cn, dn);

    }

}
