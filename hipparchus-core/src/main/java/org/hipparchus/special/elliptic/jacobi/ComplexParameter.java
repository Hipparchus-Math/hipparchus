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
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/** Algorithm for computing the principal Jacobi functions for complex parameter m.
 * @since 2.0
 */
class ComplexParameter extends FieldJacobiElliptic<Complex> {

    /** Jacobi θ functions. */
    private final FieldJacobiTheta<Complex> jacobiTheta;

    /** Quarter period K. */
    private final Complex bigK;

    /** Quarter period iK'. */
    private final Complex iBigKPrime;

    /** Real periodic factor for K. */
    private final double rK;

    /** Imaginary periodic factor for K. */
    private final double iK;

    /** Real periodic factor for iK'. */
    private final double rKPrime;

    /** Imaginary periodic factor for iK'. */
    private final double iKPrime;

    /** Value of Jacobi θ functions at origin. */
    private final FieldTheta<Complex> t0;

    /** Scaling factor. */
    private final Complex scaling;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function
     */
    ComplexParameter(final Complex m) {

        super(m);

        // compute nome
        final Complex q = LegendreEllipticIntegral.nome(m);

        // compute periodic factors such that
        // z = 4 K [rK Re(z) + iK Im(z)] + 4i K' [rK' Re(z) + iK' Im(z)]
        bigK                 = LegendreEllipticIntegral.bigK(m);
        iBigKPrime           = LegendreEllipticIntegral.bigKPrime(m).multiplyPlusI();
        final double inverse = 0.25 /
                               (bigK.getRealPart()      * iBigKPrime.getImaginaryPart() -
                                bigK.getImaginaryPart() * iBigKPrime.getRealPart());
        this.rK              = iBigKPrime.getImaginaryPart() *  inverse;
        this.iK              = iBigKPrime.getRealPart()      * -inverse;
        this.rKPrime         = bigK.getImaginaryPart()       * -inverse;
        this.iKPrime         = bigK.getRealPart()            *  inverse;

        // prepare underlying Jacobi θ functions
        this.jacobiTheta = new FieldJacobiTheta<>(q);
        this.t0          = jacobiTheta.values(m.getField().getZero());
        this.scaling     = bigK.reciprocal().multiply(MathUtils.SEMI_PI);

    }

    /** {@inheritDoc}
     * <p>
     * The algorithm for evaluating the functions is based on {@link FieldJacobiTheta
     * Jacobi theta functions}.
     * </p>
     */
    @Override
    public FieldCopolarN<Complex> valuesN(Complex u) {

        // perform argument reduction
        final double cK      = rK * u.getRealPart() + iK * u.getImaginaryPart();
        final double cKPrime = rKPrime * u.getRealPart() + iKPrime * u.getImaginaryPart();
        final Complex reducedU = u.linearCombination(1.0,                        u,
                                                    -4 * FastMath.rint(cK),      bigK,
                                                    -4 * FastMath.rint(cKPrime), iBigKPrime);

        // evaluate Jacobi θ functions at argument
        final FieldTheta<Complex> tZ = jacobiTheta.values(reducedU.multiply(scaling));

        // convert to Jacobi elliptic functions
        final Complex sn = t0.theta3().multiply(tZ.theta1()).divide(t0.theta2().multiply(tZ.theta4()));
        final Complex cn = t0.theta4().multiply(tZ.theta2()).divide(t0.theta2().multiply(tZ.theta4()));
        final Complex dn = t0.theta4().multiply(tZ.theta3()).divide(t0.theta3().multiply(tZ.theta4()));

        return new FieldCopolarN<>(sn, cn, dn);

    }

}
