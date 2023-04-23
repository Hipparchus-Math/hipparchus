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
package org.hipparchus.analysis.integration.gauss;

import org.hipparchus.dfp.DfpField;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Pair;

/**
 * Class that provides different ways to compute the nodes and weights to be
 * used by the {@link GaussIntegrator Gaussian integration rule}.
 */
public class GaussIntegratorFactory {

    /** Number of digits for Legendre high precision. */
    public static final int DEFAULT_DECIMAL_DIGITS = 40;

    /** Generator of Gauss-Legendre integrators. */
    private final RuleFactory legendre;
    /** Generator of Gauss-Legendre integrators. */
    private final RuleFactory legendreHighPrecision;
    /** Generator of Gauss-Hermite integrators. */
    private final RuleFactory hermite;
    /** Generator of Gauss-Laguerre integrators. */
    private final RuleFactory laguerre;

    /** Simple constructor.
     */
    public GaussIntegratorFactory() {
        this(DEFAULT_DECIMAL_DIGITS);
    }

    /** Simple constructor.
     * @param decimalDigits minimum number of decimal digits for {@link #legendreHighPrecision(int)}
     */
    public GaussIntegratorFactory(final int decimalDigits) {
        legendre              = new LegendreRuleFactory();
        legendreHighPrecision = new ConvertingRuleFactory<>(new FieldLegendreRuleFactory<>(new DfpField(decimalDigits)));
        hermite               = new HermiteRuleFactory();
        laguerre              = new LaguerreRuleFactory();
    }

    /**
     * Creates a Gauss-Laguerre integrator of the given order.
     * The call to the
     * {@link GaussIntegrator#integrate(org.hipparchus.analysis.UnivariateFunction)
     * integrate} method will perform an integration on the interval
     * \([0, +\infty)\): the computed value is the improper integral of
     * \(e^{-x} f(x)\)
     * where \(f(x)\) is the function passed to the
     * {@link SymmetricGaussIntegrator#integrate(org.hipparchus.analysis.UnivariateFunction)
     * integrate} method.
     *
     * @param numberOfPoints Order of the integration rule.
     * @return a Gauss-Legendre integrator.
     */
    public GaussIntegrator laguerre(int numberOfPoints) {
        return new GaussIntegrator(laguerre.getRule(numberOfPoints));
    }

    /**
     * Creates a Gauss-Legendre integrator of the given order.
     * The call to the
     * {@link GaussIntegrator#integrate(org.hipparchus.analysis.UnivariateFunction)
     * integrate} method will perform an integration on the natural interval
     * {@code [-1 , 1]}.
     *
     * @param numberOfPoints Order of the integration rule.
     * @return a Gauss-Legendre integrator.
     */
    public GaussIntegrator legendre(int numberOfPoints) {
        return new GaussIntegrator(legendre.getRule(numberOfPoints));
    }

    /**
     * Creates a Gauss-Legendre integrator of the given order.
     * The call to the
     * {@link GaussIntegrator#integrate(org.hipparchus.analysis.UnivariateFunction)
     * integrate} method will perform an integration on the given interval.
     *
     * @param numberOfPoints Order of the integration rule.
     * @param lowerBound Lower bound of the integration interval.
     * @param upperBound Upper bound of the integration interval.
     * @return a Gauss-Legendre integrator.
     * @throws MathIllegalArgumentException if number of points is not positive
     */
    public GaussIntegrator legendre(int numberOfPoints,
                                    double lowerBound,
                                    double upperBound)
        throws MathIllegalArgumentException {
        return new GaussIntegrator(transform(legendre.getRule(numberOfPoints),
                                             lowerBound, upperBound));
    }

    /**
     * Creates a Gauss-Legendre integrator of the given order.
     * The call to the
     * {@link GaussIntegrator#integrate(org.hipparchus.analysis.UnivariateFunction)
     * integrate} method will perform an integration on the natural interval
     * {@code [-1 , 1]}.
     *
     * @param numberOfPoints Order of the integration rule.
     * @return a Gauss-Legendre integrator.
     * @throws MathIllegalArgumentException if number of points is not positive
     */
    public GaussIntegrator legendreHighPrecision(int numberOfPoints)
        throws MathIllegalArgumentException {
        return new GaussIntegrator(legendreHighPrecision.getRule(numberOfPoints));
    }

    /**
     * Creates an integrator of the given order, and whose call to the
     * {@link GaussIntegrator#integrate(org.hipparchus.analysis.UnivariateFunction)
     * integrate} method will perform an integration on the given interval.
     *
     * @param numberOfPoints Order of the integration rule.
     * @param lowerBound Lower bound of the integration interval.
     * @param upperBound Upper bound of the integration interval.
     * @return a Gauss-Legendre integrator.
     * @throws MathIllegalArgumentException if number of points is not positive
     */
    public GaussIntegrator legendreHighPrecision(int numberOfPoints,
                                                 double lowerBound,
                                                 double upperBound)
        throws MathIllegalArgumentException {
        return new GaussIntegrator(transform(legendreHighPrecision.getRule(numberOfPoints),
                                             lowerBound, upperBound));
    }

    /**
     * Creates a Gauss-Hermite integrator of the given order.
     * The call to the
     * {@link SymmetricGaussIntegrator#integrate(org.hipparchus.analysis.UnivariateFunction)
     * integrate} method will perform a weighted integration on the interval
     * \([-\infty, +\infty]\): the computed value is the improper integral of
     * \(e^{-x^2}f(x)\)
     * where \(f(x)\) is the function passed to the
     * {@link SymmetricGaussIntegrator#integrate(org.hipparchus.analysis.UnivariateFunction)
     * integrate} method.
     *
     * @param numberOfPoints Order of the integration rule.
     * @return a Gauss-Hermite integrator.
     */
    public SymmetricGaussIntegrator hermite(int numberOfPoints) {
        return new SymmetricGaussIntegrator(hermite.getRule(numberOfPoints));
    }

    /**
     * Performs a change of variable so that the integration can be performed
     * on an arbitrary interval {@code [a, b]}.
     * It is assumed that the natural interval is {@code [-1, 1]}.
     *
     * @param rule Original points and weights.
     * @param a Lower bound of the integration interval.
     * @param b Lower bound of the integration interval.
     * @return the points and weights adapted to the new interval.
     */
    private Pair<double[], double[]> transform(Pair<double[], double[]> rule, double a, double b) {
        final double[] points = rule.getFirst();
        final double[] weights = rule.getSecond();

        // Scaling
        final double scale = (b - a) / 2;
        final double shift = a + scale;

        for (int i = 0; i < points.length; i++) {
            points[i] = points[i] * scale + shift;
            weights[i] *= scale;
        }

        return new Pair<>(points, weights);
    }

}
