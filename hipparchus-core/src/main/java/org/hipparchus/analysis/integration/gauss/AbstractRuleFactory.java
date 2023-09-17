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
package org.hipparchus.analysis.integration.gauss;

import java.util.Arrays;
import java.util.SortedMap;
import java.util.TreeMap;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Incrementor;
import org.hipparchus.util.Pair;

/**
 * Base class for rules that determines the integration nodes and their
 * weights.
 * Subclasses must implement the {@link #computeRule(int) computeRule} method.
 *
 * @since 2.0
 */
public abstract class AbstractRuleFactory implements RuleFactory {

    /** List of points and weights, indexed by the order of the rule. */
    private final SortedMap<Integer, Pair<double[], double[]>> pointsAndWeights = new TreeMap<>();

    /** Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 3.0
     */
    public AbstractRuleFactory() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** {@inheritDoc} */
    @Override
    public Pair<double[], double[]> getRule(int numberOfPoints)
        throws MathIllegalArgumentException {

        if (numberOfPoints <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_POINTS,
                                                   numberOfPoints);
        }
        if (numberOfPoints > 1000) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE,
                                                   numberOfPoints, 1000);
        }

        Pair<double[], double[]> rule;
        synchronized (pointsAndWeights) {
            // Try to obtain the rule from the cache.
            rule = pointsAndWeights.get(numberOfPoints);

            if (rule == null) {
                // Rule not computed yet.

                // Compute the rule.
                rule = computeRule(numberOfPoints);

                // Cache it.
                pointsAndWeights.put(numberOfPoints, rule);
            }
        }

        // Return a copy.
        return new Pair<>(rule.getFirst().clone(), rule.getSecond().clone());

    }

    /**
     * Computes the rule for the given order.
     *
     * @param numberOfPoints Order of the rule to be computed.
     * @return the computed rule.
     * @throws MathIllegalArgumentException if the elements of the pair do not
     * have the same length.
     */
    protected abstract Pair<double[], double[]> computeRule(int numberOfPoints)
        throws MathIllegalArgumentException;

    /** Computes roots of the associated orthogonal polynomials.
     * <p>
     * The roots are found using the <a href="https://en.wikipedia.org/wiki/Aberth_method">Aberth method</a>.
     * The guess points for initializing search for degree n are fixed for degrees 1 and 2 and are
     * selected from n-1 roots of rule n-1 (the two extreme roots are used, plus the n-1 intermediate
     * points between all roots).
     * </p>
     * @param n number of roots to search for
     * @param ratioEvaluator function evaluating the ratio Pₙ(x)/Pₙ'(x)
     * @return sorted array of roots
     */
    protected double[] findRoots(final int n, final UnivariateFunction ratioEvaluator) {

        final double[] roots  = new double[n];

        // set up initial guess
        if (n == 1) {
            // arbitrary guess
            roots[0] = 0;
        } else if (n == 2) {
            // arbitrary guess
            roots[0] = -1;
            roots[1] = +1;
        } else {

            // get roots from previous rule.
            // If it has not been computed yet it will trigger a recursive call
            final double[] previousPoints = getRule(n - 1).getFirst();

            // first guess at previous first root
            roots[0] = previousPoints[0];

            // intermediate guesses between previous roots
            for (int i = 1; i < n - 1; ++i) {
                roots[i] = (previousPoints[i - 1] + previousPoints[i]) * 0.5;
            }

            // last guess at previous last root
            roots[n - 1] = previousPoints[n - 2];

        }

        // use Aberth method to find all roots simultaneously
        final double[]    ratio       = new double[n];
        final Incrementor incrementor = new Incrementor(1000);
        double            tol;
        double            maxOffset;
        do {

            // safety check that triggers an exception if too much iterations are made
            incrementor.increment();

            // find the ratio P(xᵢ)/P'(xᵢ) for all current roots approximations
            for (int i = 0; i < n; ++i) {
                ratio[i] = ratioEvaluator.value(roots[i]);
            }

            // move roots approximations all at once, using Aberth method
            maxOffset = 0;
            for (int i = 0; i < n; ++i) {
                double sum = 0;
                for (int j = 0; j < n; ++j) {
                    if (j != i) {
                        sum += 1 / (roots[i] - roots[j]);
                    }
                }
                final double offset = ratio[i] / (1 - ratio[i] * sum);
                maxOffset = FastMath.max(maxOffset, FastMath.abs(offset));
                roots[i] -= offset;
            }

            // we set tolerance to 1 ulp of the largest root
            tol = 0;
            for (final double r : roots) {
                tol = FastMath.max(tol, FastMath.ulp(r));
            }

        } while (maxOffset > tol);

        // sort the roots
        Arrays.sort(roots);

        return roots;

    }

    /** Enforce symmetry of roots.
     * @param roots roots to process in place
     */
    protected void enforceSymmetry(final double[] roots) {

        final int n = roots.length;

        // enforce symmetry
        for (int i = 0; i < n / 2; ++i) {
            final int idx = n - i - 1;
            final double c = (roots[i] - roots[idx]) * 0.5;
            roots[i]   = +c;
            roots[idx] = -c;
        }

        // If n is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (n % 2 != 0) {
            roots[n / 2] = 0;
        }

    }

}
