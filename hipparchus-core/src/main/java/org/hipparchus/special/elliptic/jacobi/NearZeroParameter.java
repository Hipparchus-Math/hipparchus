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
import org.hipparchus.util.SinCos;

/** Algorithm for computing the principal Jacobi functions for parameters slightly above zero.
 * <p>
 * The algorithm for evaluating the functions is based on approximation
 * in terms of circular functions. It is given in Abramowitz and Stegun,
 * sections 16.13.
 * </p>
 * @since 2.0
 */
class NearZeroParameter extends JacobiElliptic {

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function (must be zero or slightly positive here)
     */
    NearZeroParameter(final double m) {
        super(m);
    }

    /** {@inheritDoc} */
    @Override
    public CopolarN valuesN(final double u) {
        final SinCos sc     = FastMath.sinCos(u);
        final double factor = 0.25 * getM() * (u - sc.sin() * sc.cos());
        return new CopolarN(sc.sin() - factor * sc.cos(),            // equation 16.13.1
                            sc.cos() + factor * sc.sin(),            // equation 16.13.2
                            1 - 0.5 * getM() * sc.sin() * sc.sin()); // equation 16.13.3
    }

}
