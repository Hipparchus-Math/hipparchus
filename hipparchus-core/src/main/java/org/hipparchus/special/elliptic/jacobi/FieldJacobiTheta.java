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

/** Algorithm computing Jacobi theta functions.
 * @param <T> the type of the field elements
 * @since 2.0
 */
public class FieldJacobiTheta<T extends CalculusFieldElement<T>> {

    /** Maximum number of terms in the Fourier series. */
    private static final int N_MAX = 100;

    /** Nome. */
    private final T q;

    /** q². */
    private final T qSquare;

    /** ∜q. */
    private final T qFourth;

    /** Simple constructor.
     * <p>
     * The nome {@code q} can be computed using ratios of complete elliptic integrals
     * ({@link org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral#nome(CalculusFieldElement)
     * LegendreEllipticIntegral.nome(m)} which are themselves defined in term of parameter m,
     * where m=k² and k is the elliptic modulus.
     * </p>
     * @param q nome
     */
    public FieldJacobiTheta(final T q) {
        this.q       = q;
        this.qSquare = q.multiply(q);
        this.qFourth = FastMath.sqrt(FastMath.sqrt(q));
    }

    /** Get the nome.
     * @return nome
     */
    public T getQ() {
        return q;
    }

    /** Evaluate the Jacobi theta functions.
     * @param z argument of the functions
     * @return container for the four Jacobi theta functions θ₁(z|τ), θ₂(z|τ), θ₃(z|τ), and θ₄(z|τ)
     */
    public FieldTheta<T> values(final T z) {

        // the computation is based on Fourier series,
        // see Digital Library of Mathematical Functions section 20.2
        // https://dlmf.nist.gov/20.2
        final T zero = q.getField().getZero();
        final T one  = q.getField().getOne();

        // base angle for Fourier Series
        final FieldSinCos<T> sc1 = FastMath.sinCos(z);

        // recursion rules initialization
        double         sgn   = 1.0;
        T              qNN   = one;
        T              qTwoN = one;
        T              qNNp1 = one;
        FieldSinCos<T> sc2n1 = sc1;
        final double   eps   = FastMath.ulp(one).getReal();

        // Fourier series
        T sum1 = sc1.sin();
        T sum2 = sc1.cos();
        T sum3 = zero;
        T sum4 = zero;
        for (int n = 1; n < N_MAX; ++n) {

            sgn   = -sgn;                            // (-1)ⁿ⁻¹     ← (-1)ⁿ
            qNN   = qNN.multiply(qTwoN).multiply(q); // q⁽ⁿ⁻¹⁾⁽ⁿ⁻¹⁾ ← qⁿⁿ
            qTwoN = qTwoN.multiply(qSquare);         // q²⁽ⁿ⁻¹⁾     ← q²ⁿ
            qNNp1 = qNNp1.multiply(qTwoN);           // q⁽ⁿ⁻¹⁾ⁿ     ← qⁿ⁽ⁿ⁺¹⁾

            sc2n1 = FieldSinCos.sum(sc2n1, sc1); // {sin|cos}([2n-1] z) ← {sin|cos}(2n z)
            sum3  = sum3.add(sc2n1.cos().multiply(qNN));
            sum4  = sum4.add(sc2n1.cos().multiply(qNN.multiply(sgn)));

            sc2n1 = FieldSinCos.sum(sc2n1, sc1); // {sin|cos}(2n z) ← {sin|cos}([2n+1] z)
            sum1  = sum1.add(sc2n1.sin().multiply(qNNp1.multiply(sgn)));
            sum2  = sum2.add(sc2n1.cos().multiply(qNNp1));

            if (qNNp1.norm() <= eps) {
                // we have reach convergence
                break;
            }

        }

        return new FieldTheta<>(sum1.multiply(qFourth.multiply(2)),
                                sum2.multiply(qFourth.multiply(2)),
                                sum3.multiply(2).add(1),
                                sum4.multiply(2).add(1));

    }

}
