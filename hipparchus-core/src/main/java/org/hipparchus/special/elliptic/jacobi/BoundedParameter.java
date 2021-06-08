/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.special.elliptic.jacobi;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.SinCos;

/** Algorithm for computing the principal Jacobi functions for parameter m in [0; 1].
 * @since 2.0
 */
class BoundedParameter extends JacobiElliptic {

    /** Max number of iterations of the AGM scale. */
    private static final int N_MAX = 16;

    /** Initial value for arithmetic-geometric mean. */
    private final double b0;

    /** Initial value for arithmetic-geometric mean. */
    private final double c0;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function
     */
    BoundedParameter(final double m) {
        super(m);
        this.b0 = FastMath.sqrt(1.0 - m);
        this.c0 = FastMath.sqrt(m);
    }

    /** {@inheritDoc}
     * <p>
     * The algorithm for evaluating the functions is based on arithmetic-geometric
     * mean. It is given in Abramowitz and Stegun, sections 16.4 and 17.6.
     * </p>
     */
    @Override
    public CopolarN valuesN(double u) {

        // initialize scale
        final double[] a = new double[N_MAX];
        final double[] c = new double[N_MAX];
        a[0]      = 1.0;
        double bi = b0;
        c[0]      = c0;

        // iterate down
        double phi = u;
        for (int i = 1; i < N_MAX; ++i) {

            // 2ⁿ u
            phi += phi;
            c[i] = 0.5 * (a[i - 1] - bi);

            // arithmetic mean
            a[i] = 0.5 * (a[i - 1] + bi);

            // geometric mean
            bi = FastMath.sqrt(a[i - 1] * bi);

            // convergence (by the inequality of arithmetic and geometric means, this is non-negative)
            if (c[i] <= FastMath.ulp(a[i])) {
                // convergence has been reached

                // iterate up
                phi *= a[i];
                for (int j = i; j > 0; --j) {
                    // equation 16.4.3 in Abramowitz and Stegun
                    phi = 0.5 * (phi + FastMath.asin(c[j] * FastMath.sin(phi) / a[j]));
                }
                // using 16.1.5 rather than 16.4.4 to avoid computing another cosine
                final SinCos scPhi0 = FastMath.sinCos(phi);
                return new CopolarN(scPhi0.sin(), scPhi0.cos(),
                                    FastMath.sqrt(1 - getM() * scPhi0.sin() * scPhi0.sin()));

            }

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

}
