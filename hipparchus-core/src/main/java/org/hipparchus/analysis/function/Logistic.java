/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.analysis.function;

import org.hipparchus.analysis.ParametricUnivariateFunction;
import org.hipparchus.analysis.differentiation.Derivative;
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * <a href="http://en.wikipedia.org/wiki/Generalised_logistic_function">
 *  Generalised logistic</a> function.
 *
 */
public class Logistic implements UnivariateDifferentiableFunction {
    /** Lower asymptote. */
    private final double a;
    /** Upper asymptote. */
    private final double k;
    /** Growth rate. */
    private final double b;
    /** Parameter that affects near which asymptote maximum growth occurs. */
    private final double oneOverN;
    /** Parameter that affects the position of the curve along the ordinate axis. */
    private final double q;
    /** Abscissa of maximum growth. */
    private final double m;

    /** Simple constructor.
     * @param k If {@code b > 0}, value of the function for x going towards +&infin;.
     * If {@code b < 0}, value of the function for x going towards -&infin;.
     * @param m Abscissa of maximum growth.
     * @param b Growth rate.
     * @param q Parameter that affects the position of the curve along the
     * ordinate axis.
     * @param a If {@code b > 0}, value of the function for x going towards -&infin;.
     * If {@code b < 0}, value of the function for x going towards +&infin;.
     * @param n Parameter that affects near which asymptote the maximum
     * growth occurs.
     * @throws MathIllegalArgumentException if {@code n <= 0}.
     */
    public Logistic(double k,
                    double m,
                    double b,
                    double q,
                    double a,
                    double n)
        throws MathIllegalArgumentException {
        if (n <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL_BOUND_EXCLUDED,
                                                   n, 0);
        }

        this.k = k;
        this.m = m;
        this.b = b;
        this.q = q;
        this.a = a;
        oneOverN = 1 / n;
    }

    /** {@inheritDoc} */
    @Override
    public double value(double x) {
        return value(m - x, k, b, q, a, oneOverN);
    }

    /**
     * Parametric function where the input array contains the parameters of
     * the {@link Logistic#Logistic(double,double,double,double,double,double)
     * logistic function}, ordered as follows:
     * <ul>
     *  <li>k</li>
     *  <li>m</li>
     *  <li>b</li>
     *  <li>q</li>
     *  <li>a</li>
     *  <li>n</li>
     * </ul>
     */
    public static class Parametric implements ParametricUnivariateFunction {

        /** Empty constructor.
         * <p>
         * This constructor is not strictly necessary, but it prevents spurious
         * javadoc warnings with JDK 18 and later.
         * </p>
         * @since 3.0
         */
        public Parametric() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
            // nothing to do
        }

        /**
         * Computes the value of the sigmoid at {@code x}.
         *
         * @param x Value for which the function must be computed.
         * @param param Values for {@code k}, {@code m}, {@code b}, {@code q},
         * {@code a} and  {@code n}.
         * @return the value of the function.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws MathIllegalArgumentException if the size of {@code param} is
         * not 6.
         * @throws MathIllegalArgumentException if {@code param[5] <= 0}.
         */
        @Override
        public double value(double x, double ... param)
            throws MathIllegalArgumentException, NullArgumentException {
            validateParameters(param);
            return Logistic.value(param[1] - x, param[0],
                                  param[2], param[3],
                                  param[4], 1 / param[5]);
        }

        /**
         * Computes the value of the gradient at {@code x}.
         * The components of the gradient vector are the partial
         * derivatives of the function with respect to each of the
         * <em>parameters</em>.
         *
         * @param x Value at which the gradient must be computed.
         * @param param Values for {@code k}, {@code m}, {@code b}, {@code q},
         * {@code a} and  {@code n}.
         * @return the gradient vector at {@code x}.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws MathIllegalArgumentException if the size of {@code param} is
         * not 6.
         * @throws MathIllegalArgumentException if {@code param[5] <= 0}.
         */
        @Override
        public double[] gradient(double x, double ... param)
            throws MathIllegalArgumentException, NullArgumentException {
            validateParameters(param);

            final double b = param[2];
            final double q = param[3];

            final double mMinusX = param[1] - x;
            final double oneOverN = 1 / param[5];
            final double exp = FastMath.exp(b * mMinusX);
            final double qExp = q * exp;
            final double qExp1 = qExp + 1;
            final double factor1 = (param[0] - param[4]) * oneOverN / FastMath.pow(qExp1, oneOverN);
            final double factor2 = -factor1 / qExp1;

            // Components of the gradient.
            final double gk = Logistic.value(mMinusX, 1, b, q, 0, oneOverN);
            final double gm = factor2 * b * qExp;
            final double gb = factor2 * mMinusX * qExp;
            final double gq = factor2 * exp;
            final double ga = Logistic.value(mMinusX, 0, b, q, 1, oneOverN);
            final double gn = factor1 * FastMath.log(qExp1) * oneOverN;

            return new double[] { gk, gm, gb, gq, ga, gn };
        }

        /**
         * Validates parameters to ensure they are appropriate for the evaluation of
         * the {@link #value(double,double[])} and {@link #gradient(double,double[])}
         * methods.
         *
         * @param param Values for {@code k}, {@code m}, {@code b}, {@code q},
         * {@code a} and {@code n}.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws MathIllegalArgumentException if the size of {@code param} is
         * not 6.
         * @throws MathIllegalArgumentException if {@code param[5] <= 0}.
         */
        private void validateParameters(double[] param)
            throws MathIllegalArgumentException, NullArgumentException {
            MathUtils.checkNotNull(param);
            MathUtils.checkDimension(param.length, 6);
            if (param[5] <= 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL_BOUND_EXCLUDED,
                                                       param[5], 0);
            }
        }
    }

    /**
     * @param mMinusX {@code m - x}.
     * @param k {@code k}.
     * @param b {@code b}.
     * @param q {@code q}.
     * @param a {@code a}.
     * @param oneOverN {@code 1 / n}.
     * @return the value of the function.
     */
    private static double value(double mMinusX,
                                double k,
                                double b,
                                double q,
                                double a,
                                double oneOverN) {
        return a + (k - a) / FastMath.pow(1 + q * FastMath.exp(b * mMinusX), oneOverN);
    }

    /** {@inheritDoc}
     */
    @Override
    public <T extends Derivative<T>> T value(T t) {
        return t.negate().add(m).multiply(b).exp().multiply(q).add(1).pow(oneOverN).reciprocal().multiply(k - a).add(a);
    }

}
