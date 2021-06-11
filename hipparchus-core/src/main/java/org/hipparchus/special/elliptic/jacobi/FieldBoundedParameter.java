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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.MathArrays;

/** Algorithm for computing the principal Jacobi functions for parameter m in [0; 1].
 * @param <T> the type of the field elements
 * @since 2.0
 */
class FieldBoundedParameter<T extends CalculusFieldElement<T>> extends FieldJacobiElliptic<T> {

    /** Max number of iterations of the AGM scale.
     * <p>
     * This value seems sufficient even for Dfp with high accuracy as the number
     * of digits doubles at each iteration. An experiment with 300 significant
     * digits showed we reached convergence at iteration 10
     * </p>
     */
    private static final int N_MAX = 16;

    /** Initial value for arithmetic-geometric mean. */
    private final T b0;

    /** Initial value for arithmetic-geometric mean. */
    private final T c0;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function
     */
    FieldBoundedParameter(final T m) {
        super(m);
        this.b0 = FastMath.sqrt(m.getField().getOne().subtract(m));
        this.c0 = FastMath.sqrt(m);
    }

    /** {@inheritDoc}
     * <p>
     * The algorithm for evaluating the functions is based on arithmetic-geometric
     * mean. It is given in Abramowitz and Stegun, sections 16.4 and 17.6.
     * </p>
     */
    @Override
    public FieldCopolarN<T> valuesN(T u) {

        // initialize scale
        final T one = getM().getField().getOne();
        final T[] a = MathArrays.buildArray(u.getField(), N_MAX);
        final T[] c = MathArrays.buildArray(u.getField(), N_MAX);
        a[0]        = one;
        T bi        = b0;
        c[0]        = c0;

        // iterate down
        T phi = u;
        for (int i = 1; i < N_MAX; ++i) {

            // 2ⁿ u
            phi = phi.add(phi);
            c[i] = a[i - 1].subtract(bi).multiply(0.5);

            // arithmetic mean
            a[i] = a[i - 1].add(bi).multiply(0.5);

            // geometric mean
            bi = FastMath.sqrt(a[i - 1].multiply(bi));

            // convergence (by the inequality of arithmetic and geometric means, this is non-negative)
            if (c[i].getReal() <= FastMath.ulp(a[i]).getReal()) {
                // convergence has been reached

                // iterate up
                phi = phi.multiply(a[i]);
                for (int j = i; j > 0; --j) {
                    // equation 16.4.3 in Abramowitz and Stegun
                    phi = phi.add(FastMath.asin(c[j].multiply(FastMath.sin(phi)).divide(a[j]))).multiply(0.5);
                }
                // using 16.1.5 rather than 16.4.4 to avoid computing another cosine
                final FieldSinCos<T> scPhi0 = FastMath.sinCos(phi);
                return new FieldCopolarN<>(scPhi0.sin(),
                                           scPhi0.cos(),
                                           FastMath.sqrt(one.subtract(getM().multiply(scPhi0.sin()).multiply(scPhi0.sin()))));

            }

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

}
