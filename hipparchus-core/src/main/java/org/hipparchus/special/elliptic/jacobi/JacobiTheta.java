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

import org.hipparchus.complex.Complex;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.Precision;

/** Algorithm computing Jacobi theta functions.
 * @since 2.0
 */
public class JacobiTheta {

    /** Maximum number of terms in the Fourier series. */
    private static final int N_MAX = 100;

    /** Nome. */
    private final double q;

    /** q². */
    private final double qSquare;

    /** ∜q. */
    private final double qFourth;

    /** Simple constructor.
     * <p>
     * Beware that {@link org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral
     * Legendre elliptic integrals} are defined in terms of elliptic modulus {@code k} whereas
     * {@link JacobiEllipticBuilder#build(double) Jacobi elliptic functions} (which
     * are their inverses) are defined in terms of parameter {@code m} and {@link
     * JacobiTheta#JacobiTheta(double) Jacobi theta functions} are defined
     * in terms of the {@link
     * org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral#nome(double)
     * nome q}. All are related as {@code k² = m} and the nome can be computed from ratios of complete
     * elliptic integrals.
     * </p>
     * @param q nome
     */
    public JacobiTheta(final double q) {
        this.q       = q;
        this.qSquare = q * q;
        this.qFourth = FastMath.sqrt(FastMath.sqrt(q));
    }

    /** Get the nome.
     * @return nome
     */
    public double getQ() {
        return q;
    }

    /** Evaluate the Jacobi theta functions.
     * @param z argument of the functions
     * @return container for the four Jacobi theta functions θ₁(z|τ), θ₂(z|τ), θ₃(z|τ), and θ₄(z|τ)
     */
    public Theta values(final Complex z) {

        // the computation is based on Fourier series,
        // see Digital Library of Mathematical Functions section 20.2
        // https://dlmf.nist.gov/20.2

        // base angle for Fourier Series
        final FieldSinCos<Complex> sc1 = FastMath.sinCos(z);

        // recursion rules initialization
        double               sgn   = 1.0;
        double               qNN   = 1.0;
        double               qTwoN = 1.0;
        double               qNNp1 = 1.0;
        FieldSinCos<Complex> sc2n1 = sc1;

        // Fourier series
        Complex sum1 = sc1.sin();
        Complex sum2 = sc1.cos();
        Complex sum3 = Complex.ZERO;
        Complex sum4 = Complex.ZERO;
        for (int n = 1; n < N_MAX; ++n) {

            sgn   = -sgn;              // (-1)ⁿ⁻¹     ← (-1)ⁿ
            qNN   = qNN   * qTwoN * q; // q⁽ⁿ⁻¹⁾⁽ⁿ⁻¹⁾ ← qⁿⁿ
            qTwoN = qTwoN * qSquare;   // q²⁽ⁿ⁻¹⁾     ← q²ⁿ
            qNNp1 = qNNp1 * qTwoN;     // q⁽ⁿ⁻¹⁾ⁿ     ← qⁿ⁽ⁿ⁺¹⁾

            sc2n1 = FieldSinCos.sum(sc2n1, sc1); // {sin|cos}([2n-1] z) ← {sin|cos}(2n z)
            sum3  = sum3.add(sc2n1.cos().multiply(qNN));
            sum4  = sum4.add(sc2n1.cos().multiply(sgn * qNN));

            sc2n1 = FieldSinCos.sum(sc2n1, sc1); // {sin|cos}(2n z) ← {sin|cos}([2n+1] z)
            sum1  = sum1.add(sc2n1.sin().multiply(sgn * qNNp1));
            sum2  = sum2.add(sc2n1.cos().multiply(qNNp1));

            if (qNNp1 <= Precision.EPSILON) {
                // we have reach convergence
                return new Theta(sum1.multiply(2 * qFourth), sum2.multiply(2 * qFourth),
                                 sum3.multiply(2).add(1),    sum4.multiply(2).add(1));
            }

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

}
