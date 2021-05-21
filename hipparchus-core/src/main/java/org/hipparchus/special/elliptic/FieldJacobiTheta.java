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
package org.hipparchus.special.elliptic;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;

/** Algorithm computing Jacobi theta functions.
 * @param <T> the type of the field elements
 * @since 2.0
 */
public abstract class FieldJacobiTheta<T extends CalculusFieldElement<T>> {

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
     * Beware that {@link FieldEllipticIntegral#FieldEllipticIntegral elliptic
     * integrals} are defined in terms of elliptic modulus {@code k} whereas
     * {@link JacobiEllipticBuilder#build(CalculusFieldElement) Jacobi elliptic
     * functions} (which are their inverses) are defined in terms of parameter
     * {@code m} and {@link FieldJacobiTheta#FieldJacobiTheta Jacobi theta
     * functions} are defined in terms of the {@link FieldEllipticIntegral#getNome()
     * nome q}. All are related as {@code k² = m} and the nome can be computed
     * from ratios of complete elliptic integrals.
     * </p>
     * @param q nome
     */
    protected FieldJacobiTheta(final T q) {
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
    public FieldTheta<T> values(final FieldComplex<T> z) {

        // the computation is based on Fourier series,
        // see Digital Library of Mathematical Functions section 20.2
        // https://dlmf.nist.gov/20.2

        // base angle for Fourier Series
        final FieldSinCos<FieldComplex<T>> sc1 = FastMath.sinCos(z);

        // recursion rules initialization
        double                       sgn   = 1.0;
        T                            qTwoN = q.getField().getOne();
        T                            qNN   = q.getField().getOne();
        T                            qNNp1 = q.getField().getOne();
        FieldSinCos<FieldComplex<T>> sc2n1 = sc1;
        final double                 eps   = FastMath.ulp(q.getField().getOne()).getReal();

        // Fourier series
        FieldComplex<T> sum1 = sc1.sin();
        FieldComplex<T> sum2 = sc1.cos();
        FieldComplex<T> sum3 = FieldComplex.getZero(q.getField());
        FieldComplex<T> sum4 = FieldComplex.getZero(q.getField());
        for (int n = 1; n < N_MAX; ++n) {

            sgn   = -sgn;                            // -1ⁿ⁻¹       ← -1ⁿ
            qTwoN = qTwoN.multiply(qSquare);         // q²⁽ⁿ⁻¹⁾     ← q²ⁿ
            qNN   = qNN.multiply(qTwoN).multiply(q); // q⁽ⁿ⁻¹⁾⁽ⁿ⁻¹⁾ ← qⁿⁿ
            qNNp1 = qNNp1.multiply(qTwoN);           // q⁽ⁿ⁻¹⁾ⁿ     ← qⁿ⁽ⁿ⁺¹⁾

            sc2n1 = FieldSinCos.sum(sc2n1, sc1); // {sin|cos}([2n-1] z) ← {sin|cos}(2n z)
            sum3  = sum3.add(sc2n1.cos().multiply(qNN));
            sum4  = sum4.add(sc2n1.cos().multiply(qNN.multiply(sgn)));

            sc2n1 = FieldSinCos.sum(sc2n1, sc1); // {sin|cos}(2n z) ← {sin|cos}([2n+1] z)
            sum1  = sum1.add(sc2n1.sin().multiply(qNNp1.multiply(sgn)));
            sum2  = sum2.add(sc2n1.cos().multiply(qNNp1));

            if (qNNp1.getReal() <= eps) {
                // we have reach convergence
                return new FieldTheta<>(sum1.multiply(qFourth.multiply(2)), sum2.multiply(qFourth.multiply(2)),
                                        sum3.multiply(2).add(1),            sum4.multiply(2).add(1));
            }

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

    }

}
