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
package org.hipparchus.analysis.integration.gauss;

import java.util.SortedMap;
import java.util.TreeMap;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
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

    /** List of points and weights, indexed by the order of the rule. */
    private final SortedMap<Integer, Pair<T[], T[]>> pointsAndWeights = new TreeMap<>();

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

}
