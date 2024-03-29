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
package org.hipparchus.analysis.integration;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.UnivariateSolverUtils;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.Incrementor;
import org.hipparchus.util.MathUtils;

/**
 * Provide a default implementation for several generic functions.
 *
 */
public abstract class BaseAbstractUnivariateIntegrator implements UnivariateIntegrator {

    /** Default absolute accuracy. */
    public static final double DEFAULT_ABSOLUTE_ACCURACY = 1.0e-15;

    /** Default relative accuracy. */
    public static final double DEFAULT_RELATIVE_ACCURACY = 1.0e-6;

    /** Default minimal iteration count. */
    public static final int DEFAULT_MIN_ITERATIONS_COUNT = 3;

    /** Default maximal iteration count. */
    public static final int DEFAULT_MAX_ITERATIONS_COUNT = Integer.MAX_VALUE;

    /** The iteration count. */
    protected final Incrementor iterations;

    /** Maximum absolute error. */
    private final double absoluteAccuracy;

    /** Maximum relative error. */
    private final double relativeAccuracy;

    /** minimum number of iterations */
    private final int minimalIterationCount;

    /** The functions evaluation count. */
    private Incrementor evaluations;

    /** Function to integrate. */
    private UnivariateFunction function;

    /** Lower bound for the interval. */
    private double min;

    /** Upper bound for the interval. */
    private double max;

    /**
     * Construct an integrator with given accuracies and iteration counts.
     * <p>
     * The meanings of the various parameters are:
     * <ul>
     *   <li>relative accuracy:
     *       this is used to stop iterations if the absolute accuracy can't be
     *       achieved due to large values or short mantissa length. If this
     *       should be the primary criterion for convergence rather then a
     *       safety measure, set the absolute accuracy to a ridiculously small value,
     *       like {@link org.hipparchus.util.Precision#SAFE_MIN Precision.SAFE_MIN}.</li>
     *   <li>absolute accuracy:
     *       The default is usually chosen so that results in the interval
     *       -10..-0.1 and +0.1..+10 can be found with a reasonable accuracy. If the
     *       expected absolute value of your results is of much smaller magnitude, set
     *       this to a smaller value.</li>
     *   <li>minimum number of iterations:
     *       minimal iteration is needed to avoid false early convergence, e.g.
     *       the sample points happen to be zeroes of the function. Users can
     *       use the default value or choose one that they see as appropriate.</li>
     *   <li>maximum number of iterations:
     *       usually a high iteration count indicates convergence problems. However,
     *       the "reasonable value" varies widely for different algorithms. Users are
     *       advised to use the default value supplied by the algorithm.</li>
     * </ul>
     *
     * @param relativeAccuracy relative accuracy of the result
     * @param absoluteAccuracy absolute accuracy of the result
     * @param minimalIterationCount minimum number of iterations
     * @param maximalIterationCount maximum number of iterations
     * @exception MathIllegalArgumentException if minimal number of iterations
     * is not strictly positive
     * @exception MathIllegalArgumentException if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     */
    protected BaseAbstractUnivariateIntegrator(final double relativeAccuracy,
                                               final double absoluteAccuracy,
                                               final int minimalIterationCount,
                                               final int maximalIterationCount)
        throws MathIllegalArgumentException {

        // accuracy settings
        this.relativeAccuracy      = relativeAccuracy;
        this.absoluteAccuracy      = absoluteAccuracy;

        // iterations count settings
        if (minimalIterationCount <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL_BOUND_EXCLUDED,
                                                   minimalIterationCount, 0);
        }
        if (maximalIterationCount <= minimalIterationCount) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL_BOUND_EXCLUDED,
                                                   maximalIterationCount, minimalIterationCount);
        }
        this.minimalIterationCount = minimalIterationCount;
        this.iterations            = new Incrementor(maximalIterationCount);

        // prepare evaluations counter, but do not set it yet
        evaluations = new Incrementor();

    }

    /**
     * Construct an integrator with given accuracies.
     * @param relativeAccuracy relative accuracy of the result
     * @param absoluteAccuracy absolute accuracy of the result
     */
    protected BaseAbstractUnivariateIntegrator(final double relativeAccuracy,
                                           final double absoluteAccuracy) {
        this(relativeAccuracy, absoluteAccuracy,
             DEFAULT_MIN_ITERATIONS_COUNT, DEFAULT_MAX_ITERATIONS_COUNT);
    }

    /**
     * Construct an integrator with given iteration counts.
     * @param minimalIterationCount minimum number of iterations
     * @param maximalIterationCount maximum number of iterations
     * @exception MathIllegalArgumentException if minimal number of iterations
     * is not strictly positive
     * @exception MathIllegalArgumentException if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     */
    protected BaseAbstractUnivariateIntegrator(final int minimalIterationCount,
                                           final int maximalIterationCount)
        throws MathIllegalArgumentException {
        this(DEFAULT_RELATIVE_ACCURACY, DEFAULT_ABSOLUTE_ACCURACY,
             minimalIterationCount, maximalIterationCount);
    }

    /** {@inheritDoc} */
    @Override
    public double getRelativeAccuracy() {
        return relativeAccuracy;
    }

    /** {@inheritDoc} */
    @Override
    public double getAbsoluteAccuracy() {
        return absoluteAccuracy;
    }

    /** {@inheritDoc} */
    @Override
    public int getMinimalIterationCount() {
        return minimalIterationCount;
    }

    /** {@inheritDoc} */
    @Override
    public int getMaximalIterationCount() {
        return iterations.getMaximalCount();
    }

    /** {@inheritDoc} */
    @Override
    public int getEvaluations() {
        return evaluations.getCount();
    }

    /** {@inheritDoc} */
    @Override
    public int getIterations() {
        return iterations.getCount();
    }

    /** Get the lower bound.
     * @return the lower bound.
     */
    protected double getMin() {
        return min;
    }

    /** Get the upper bound.
     * @return the upper bound.
     */
    protected double getMax() {
        return max;
    }

    /**
     * Compute the objective function value.
     *
     * @param point Point at which the objective function must be evaluated.
     * @return the objective function value at specified point.
     * @throws MathIllegalStateException if the maximal number of function
     * evaluations is exceeded.
     */
    protected double computeObjectiveValue(final double point)
        throws MathIllegalStateException {
        evaluations.increment();
        return function.value(point);
    }

    /**
     * Prepare for computation.
     * Subclasses must call this method if they override any of the
     * {@code solve} methods.
     *
     * @param maxEval Maximum number of evaluations.
     * @param f the integrand function
     * @param lower the min bound for the interval
     * @param upper the upper bound for the interval
     * @throws NullArgumentException if {@code f} is {@code null}.
     * @throws MathIllegalArgumentException if {@code min >= max}.
     */
    protected void setup(final int maxEval,
                         final UnivariateFunction f,
                         final double lower, final double upper)
        throws MathIllegalArgumentException, NullArgumentException {

        // Checks.
        MathUtils.checkNotNull(f);
        UnivariateSolverUtils.verifyInterval(lower, upper);

        // Reset.
        min = lower;
        max = upper;
        function = f;
        evaluations = evaluations.withMaximalCount(maxEval);
        iterations.reset();
    }

    /** {@inheritDoc} */
    @Override
    public double integrate(final int maxEval, final UnivariateFunction f,
                            final double lower, final double upper)
        throws MathIllegalArgumentException, MathIllegalStateException, NullArgumentException {

        // Initialization.
        setup(maxEval, f, lower, upper);

        // Perform computation.
        return doIntegrate();

    }

    /**
     * Method for implementing actual integration algorithms in derived
     * classes.
     *
     * @return the root.
     * @throws MathIllegalStateException if the maximal number of evaluations
     * is exceeded.
     * @throws MathIllegalStateException if the maximum iteration count is exceeded
     * or the integrator detects convergence problems otherwise
     */
    protected abstract double doIntegrate()
        throws MathIllegalStateException;

}
