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
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;

/** Algorithm for computing the principal Jacobi functions for complex parameter m.
 * @param <T> the type of the field elements
 * @since 2.0
 */
class FieldComplexParameter<T extends CalculusFieldElement<T>> extends FieldJacobiElliptic<FieldComplex<T>> {

    /** Jacobi θ functions. */
    private final FieldJacobiTheta<FieldComplex<T>> jacobiTheta;

    /** Quarter period K. */
    private final FieldComplex<T> bigK;

    /** Quarter period iK'. */
    private final FieldComplex<T> iBigKPrime;

    /** Real periodic factor for K. */
    private final T rK;

    /** Imaginary periodic factor for K. */
    private final T iK;

    /** Real periodic factor for iK'. */
    private final T rKPrime;

    /** Imaginary periodic factor for iK'. */
    private final T iKPrime;

    /** Value of Jacobi θ functions at origin. */
    private final FieldTheta<FieldComplex<T>> t0;

    /** Scaling factor. */
    private final FieldComplex<T> scaling;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function
     */
    FieldComplexParameter(final FieldComplex<T> m) {

        super(m);

        // compute nome
         final FieldComplex<T> q = LegendreEllipticIntegral.nome(m);

        // compute periodic factors such that
        // z = 4K [rK Re(z) + iK Im(z)] + 4K' [rK' Re(z) + iK' Im(z)]
        bigK            = LegendreEllipticIntegral.bigK(m);
        iBigKPrime      = LegendreEllipticIntegral.bigKPrime(m).multiplyPlusI();
        final T inverse = bigK.getRealPart().multiply(iBigKPrime.getImaginaryPart()).
                          subtract(bigK.getImaginaryPart().multiply(iBigKPrime.getRealPart())).
                          multiply(4).reciprocal();
        this.rK         = iBigKPrime.getImaginaryPart().multiply(inverse);
        this.iK         = iBigKPrime.getRealPart().multiply(inverse).negate();
        this.rKPrime    = bigK.getImaginaryPart().multiply(inverse).negate();
        this.iKPrime    = bigK.getRealPart().multiply(inverse);

        // prepare underlying Jacobi θ functions
        this.jacobiTheta = new FieldJacobiTheta<>(q);
        this.t0          = jacobiTheta.values(m.getField().getZero());
        this.scaling     = bigK.reciprocal().multiply(m.getPi().multiply(0.5));

    }

    /** {@inheritDoc}
     * <p>
     * The algorithm for evaluating the functions is based on {@link FieldJacobiTheta
     * Jacobi theta functions}.
     * </p>
     */
    @Override
    public FieldCopolarN<FieldComplex<T>> valuesN(FieldComplex<T> u) {

        // perform argument reduction
        final T cK      = rK.multiply(u.getRealPart()).add(iK.multiply(u.getImaginaryPart()));
        final T cKPrime = rKPrime.multiply(u.getRealPart()).add(iKPrime.multiply(u.getImaginaryPart()));
        final FieldComplex<T> reducedU = u.linearCombination(1.0,                                  u,
                                                            -4 * FastMath.rint(cK.getReal()),      bigK,
                                                            -4 * FastMath.rint(cKPrime.getReal()), iBigKPrime);

        // evaluate Jacobi θ functions at argument
        final FieldTheta<FieldComplex<T>> tZ = jacobiTheta.values(reducedU.multiply(scaling));

        // convert to Jacobi elliptic functions
        final FieldComplex<T> sn = t0.theta3().multiply(tZ.theta1()).divide(t0.theta2().multiply(tZ.theta4()));
        final FieldComplex<T> cn = t0.theta4().multiply(tZ.theta2()).divide(t0.theta2().multiply(tZ.theta4()));
        final FieldComplex<T> dn = t0.theta4().multiply(tZ.theta3()).divide(t0.theta3().multiply(tZ.theta4()));

        return new FieldCopolarN<>(sn, cn, dn);

    }

}
