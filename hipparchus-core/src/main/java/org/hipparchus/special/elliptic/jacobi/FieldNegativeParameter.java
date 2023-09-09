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
import org.hipparchus.util.FastMath;

/** Algorithm for computing the principal Jacobi functions for negative parameter m.
 * <p>
 * The rules for negative parameter change are given in Abramowitz and Stegun, section 16.10.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
class FieldNegativeParameter<T extends CalculusFieldElement<T>> extends FieldJacobiElliptic<T> {

    /** Algorithm to use for the positive parameter. */
    private final FieldJacobiElliptic<T> algorithm;

    /** Input scaling factor. */
    private final T inputScale;

    /** output scaling factor. */
    private final T outputScale;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function (must be negative here)
     */
    FieldNegativeParameter(final T m) {
        super(m);
        final T omM = m.getField().getOne().subtract(m);
        algorithm   = JacobiEllipticBuilder.build(m.negate().divide(omM));
        inputScale  = FastMath.sqrt(omM);
        outputScale = inputScale.reciprocal();
    }

    /** {@inheritDoc} */
    @Override
    public FieldCopolarN<T> valuesN(final T u) {
        final FieldCopolarD<T> trioD = new FieldCopolarD<>(algorithm.valuesN(u.multiply(inputScale)));
        return new FieldCopolarN<>(outputScale.multiply(trioD.sd()), trioD.cd(), trioD.nd());
    }

}
