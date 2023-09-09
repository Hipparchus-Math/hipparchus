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

/** Algorithm for computing the principal Jacobi functions for parameter m greater than 1.
 * <p>
 * The rules for reciprocal parameter change are given in Abramowitz and Stegun,
 * sections 16.11 and 17.4.15.
 * </p>
 * @since 2.0
 */
class BigParameter extends JacobiElliptic {

    /** Algorithm to use for the positive parameter. */
    private final JacobiElliptic algorithm;

    /** Input scaling factor. */
    private final double inputScale;

    /** output scaling factor. */
    private final double outputScale;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function (must be greater than 1 here)
     */
    BigParameter(final double m) {
        super(m);
        algorithm   = JacobiEllipticBuilder.build(1.0 / m);
        inputScale  = FastMath.sqrt(m);
        outputScale = 1.0 / inputScale;
    }

    /** {@inheritDoc} */
    @Override
    public CopolarN valuesN(final double u) {
        final CopolarN trioN = algorithm.valuesN(u * inputScale);
        return new CopolarN(outputScale * trioN.sn(), trioN.dn(), trioN.cn());
    }

}
