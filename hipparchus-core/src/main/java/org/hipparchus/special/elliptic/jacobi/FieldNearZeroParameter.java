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
import org.hipparchus.util.FieldSinCos;

/** Algorithm for computing the principal Jacobi functions for parameters slightly above zero.
 * <p>
 * The algorithm for evaluating the functions is based on approximation
 * in terms of circular functions. It is given in Abramowitz and Stegun,
 * sections 16.13.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
class FieldNearZeroParameter<T extends CalculusFieldElement<T>> extends FieldJacobiElliptic<T> {

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function (must be zero or slightly positive here)
     */
    FieldNearZeroParameter(final T m) {
        super(m);
    }

    /** {@inheritDoc} */
    @Override
    public FieldCopolarN<T> valuesN(final T u) {
        final FieldSinCos<T> sc     = FastMath.sinCos(u);
        final T              factor = getM().multiply(u.subtract(sc.sin().multiply(sc.cos()))).multiply(0.25);
        return new FieldCopolarN<>(sc.sin().subtract(factor.multiply(sc.cos())),             // equation 16.13.1
                        sc.cos().add(factor.multiply(sc.sin())),                             // equation 16.13.2
                        getM().multiply(sc.sin()).multiply(sc.sin()).multiply(-0.5).add(1)); // equation 16.13.3
    }

}
