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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Incrementor;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.Pair;

/**
 * Base class for rules that determines the integration nodes and their
 * weights.
 * Subclasses must implement the {@link #computeRule(int) computeRule} method.
 *
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 * @since 2.0
 */
public abstract class FieldAbstractRuleFactory<T extends CalculusFieldElement<T>> implements FieldRuleFactory<T> {

    /** Field to which rule coefficients belong. */
    private final Field<T> field;

    /** List of points and weights, indexed by the order of the rule. */
    private final SortedMap<Integer, Pair<T[], T[]>> pointsAndWeights;

    /** Simple constructor
     * @param field field to which rule coefficients belong
     */
    public FieldAbstractRuleFactory(final Field<T> field) {
        this.field            = field;
        this.pointsAndWeights = new TreeMap<>();
    }

    /** Get the field to which rule coefficients belong.
     * @return field to which rule coefficients belong
     */
    public Field<T> getField() {
        return field;
    }

    /** {@inheritDoc} */
    @Override
    public Pair<T[], T[]> getRule(int numberOfPoints)
        throws MathIllegalArgumentException {

        if (numberOfPoints <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_POINTS,
                                                   numberOfPoints);
        }
        if (numberOfPoints > 1000) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE,
                                                   numberOfPoints, 1000);
        }

        Pair<T[], T[]> rule;
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
    protected abstract Pair<T[], T[]> computeRule(int numberOfPoints)
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
    protected T[] findRoots(final int n, final CalculusFieldUnivariateFunction<T> ratioEvaluator) {

        final T[] roots  = MathArrays.buildArray(field, n);

        // set up initial guess
        if (n == 1) {
            // arbitrary guess
            roots[0] = field.getZero();
        } else if (n == 2) {
            // arbitrary guess
            roots[0] = field.getOne().negate();
            roots[1] = field.getOne();
        } else {

            // get roots from previous rule.
            // If it has not been computed yet it will trigger a recursive call
            final T[] previousPoints = getRule(n - 1).getFirst();

            // first guess at previous first root
            roots[0] = previousPoints[0];

            // intermediate guesses between previous roots
            for (int i = 1; i < n - 1; ++i) {
                roots[i] = previousPoints[i - 1].add(previousPoints[i]).multiply(0.5);
            }

            // last guess at previous last root
            roots[n - 1] = previousPoints[n - 2];

        }

        // use Aberth method to find all roots simultaneously
        final T[]         ratio       = MathArrays.buildArray(field, n);
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
                T sum = field.getZero();
                for (int j = 0; j < n; ++j) {
                    if (j != i) {
                        sum = sum.add(roots[i].subtract(roots[j]).reciprocal());
                    }
                }
                final T offset = ratio[i].divide(sum.multiply(ratio[i]).negate().add(1));
                maxOffset = FastMath.max(maxOffset, FastMath.abs(offset).getReal());
                roots[i] = roots[i].subtract(offset);
            }

            // we set tolerance to 1 ulp of the largest root
            tol = 0;
            for (final T r : roots) {
                tol = FastMath.max(tol, FastMath.ulp(r.getReal()));
            }

        } while (maxOffset > tol);

        // sort the roots
        Arrays.sort(roots, (r1, r2) -> Double.compare(r1.getReal(), r2.getReal()));

        return roots;

    }

    /** Enforce symmetry of roots.
     * @param roots roots to process in place
     */
    protected void enforceSymmetry(final T[] roots) {

        final int n = roots.length;

        // enforce symmetry
        for (int i = 0; i < n / 2; ++i) {
            final int idx = n - i - 1;
            final T c = roots[i].subtract(roots[idx]).multiply(0.5);
            roots[i]   = c;
            roots[idx] = c.negate();
        }

        // If n is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (n % 2 != 0) {
            roots[n / 2] = field.getZero();
        }

    }

}
