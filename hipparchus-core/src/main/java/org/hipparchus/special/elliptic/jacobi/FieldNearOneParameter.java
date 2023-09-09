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
import org.hipparchus.util.FieldSinhCosh;

/** Algorithm for computing the principal Jacobi functions for parameters slightly below one.
 * <p>
 * The algorithm for evaluating the functions is based on approximation
 * in terms of hyperbolic functions. It is given in Abramowitz and Stegun,
 * sections 16.15.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
class FieldNearOneParameter<T extends CalculusFieldElement<T>> extends FieldJacobiElliptic<T> {

    /** Complementary parameter of the Jacobi elliptic function. */
    private final T m1Fourth;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function (must be one or slightly below one here)
     */
    FieldNearOneParameter(final T m) {
        super(m);
        this.m1Fourth = m.getField().getOne().subtract(m).multiply(0.25);
    }

    /** {@inheritDoc} */
    @Override
    public FieldCopolarN<T> valuesN(final T u) {
        final FieldSinhCosh<T> sch    = FastMath.sinhCosh(u);
        final T                sech   = sch.cosh().reciprocal();
        final T                t      = sch.sinh().multiply(sech);
        final T                factor = sch.sinh().multiply(sch.cosh()).subtract(u).multiply(sech).multiply(m1Fourth);
        return new FieldCopolarN<>(t.add(factor.multiply(sech)),  // equation 16.15.1
                        sech.subtract(factor.multiply(t)),        // equation 16.15.2
                        sech.add(factor.multiply(t)));            // equation 16.15.3
    }

}
