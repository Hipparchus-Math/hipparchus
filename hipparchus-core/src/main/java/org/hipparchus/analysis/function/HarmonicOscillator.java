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
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.SinCos;

/**
 * <a href="http://en.wikipedia.org/wiki/Harmonic_oscillator">
 *  simple harmonic oscillator</a> function.
 *
 */
public class HarmonicOscillator implements UnivariateDifferentiableFunction {
    /** Amplitude. */
    private final double amplitude;
    /** Angular frequency. */
    private final double omega;
    /** Phase. */
    private final double phase;

    /**
     * Harmonic oscillator function.
     *
     * @param amplitude Amplitude.
     * @param omega Angular frequency.
     * @param phase Phase.
     */
    public HarmonicOscillator(double amplitude,
                              double omega,
                              double phase) {
        this.amplitude = amplitude;
        this.omega = omega;
        this.phase = phase;
    }

    /** {@inheritDoc} */
    @Override
    public double value(double x) {
        return value(omega * x + phase, amplitude);
    }

    /**
     * Parametric function where the input array contains the parameters of
     * the harmonic oscillator function, ordered as follows:
     * <ul>
     *  <li>Amplitude</li>
     *  <li>Angular frequency</li>
     *  <li>Phase</li>
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
         * Computes the value of the harmonic oscillator at {@code x}.
         *
         * @param x Value for which the function must be computed.
         * @param param Values of norm, mean and standard deviation.
         * @return the value of the function.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws MathIllegalArgumentException if the size of {@code param} is
         * not 3.
         */
        @Override
        public double value(double x, double ... param)
            throws MathIllegalArgumentException, NullArgumentException {
            validateParameters(param);
            return HarmonicOscillator.value(x * param[1] + param[2], param[0]);
        }

        /**
         * Computes the value of the gradient at {@code x}.
         * The components of the gradient vector are the partial
         * derivatives of the function with respect to each of the
         * <em>parameters</em> (amplitude, angular frequency and phase).
         *
         * @param x Value at which the gradient must be computed.
         * @param param Values of amplitude, angular frequency and phase.
         * @return the gradient vector at {@code x}.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws MathIllegalArgumentException if the size of {@code param} is
         * not 3.
         */
        @Override
        public double[] gradient(double x, double ... param)
            throws MathIllegalArgumentException, NullArgumentException {
            validateParameters(param);

            final double amplitude = param[0];
            final double omega = param[1];
            final double phase = param[2];

            final double xTimesOmegaPlusPhase = omega * x + phase;
            final double a = HarmonicOscillator.value(xTimesOmegaPlusPhase, 1);
            final double p = -amplitude * FastMath.sin(xTimesOmegaPlusPhase);
            final double w = p * x;

            return new double[] { a, w, p };
        }

        /**
         * Validates parameters to ensure they are appropriate for the evaluation of
         * the {@link #value(double,double[])} and {@link #gradient(double,double[])}
         * methods.
         *
         * @param param Values of norm, mean and standard deviation.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws MathIllegalArgumentException if the size of {@code param} is
         * not 3.
         */
        private void validateParameters(double[] param)
            throws MathIllegalArgumentException, NullArgumentException {
            MathUtils.checkNotNull(param);
            MathUtils.checkDimension(param.length, 3);
        }
    }

    /**
     * @param xTimesOmegaPlusPhase {@code omega * x + phase}.
     * @param amplitude Amplitude.
     * @return the value of the harmonic oscillator function at {@code x}.
     */
    private static double value(double xTimesOmegaPlusPhase,
                                double amplitude) {
        return amplitude * FastMath.cos(xTimesOmegaPlusPhase);
    }

    /** {@inheritDoc}
     */
    @Override
    public <T extends Derivative<T>> T value(T t)
        throws MathIllegalArgumentException {
        final double x = t.getValue();
        double[] f = new double[t.getOrder() + 1];

        final double alpha   = omega * x + phase;
        final SinCos scAlpha = FastMath.sinCos(alpha);
        f[0] = amplitude * scAlpha.cos();
        if (f.length > 1) {
            f[1] = -amplitude * omega * scAlpha.sin();
            final double mo2 = - omega * omega;
            for (int i = 2; i < f.length; ++i) {
                f[i] = mo2 * f[i - 2];
            }
        }

        return t.compose(f);

    }

}
