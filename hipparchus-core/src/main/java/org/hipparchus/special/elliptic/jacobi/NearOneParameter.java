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

import org.hipparchus.util.FastMath;
import org.hipparchus.util.SinhCosh;

/** Algorithm for computing the principal Jacobi functions for parameters slightly below one.
 * <p>
 * The algorithm for evaluating the functions is based on approximation
 * in terms of hyperbolic functions. It is given in Abramowitz and Stegun,
 * sections 16.15.
 * </p>
 * @since 2.0
 */
class NearOneParameter extends JacobiElliptic {

    /** Complementary parameter of the Jacobi elliptic function. */
    private final double m1;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function (must be one or slightly below one here)
     */
    NearOneParameter(final double m) {
        super(m);
        this.m1 = 1.0 - m;
    }

    /** {@inheritDoc} */
    @Override
    public CopolarN valuesN(final double u) {
        final SinhCosh sch  = FastMath.sinhCosh(u);
        final double sech   =  1.0 / sch.cosh();
        final double t      = sch.sinh() * sech;
        final double factor = 0.25 * m1 * (sch.sinh() * sch.cosh()  - u) * sech;
        return new CopolarN(t + factor * sech,  // equation 16.15.1
                            sech - factor * t,  // equation 16.15.2
                            sech + factor * t); // equation 16.15.3
    }

}
