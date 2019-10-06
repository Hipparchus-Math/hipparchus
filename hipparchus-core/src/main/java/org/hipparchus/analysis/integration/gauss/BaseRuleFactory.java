/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.analysis.integration.gauss;

import java.util.SortedMap;
import java.util.TreeMap;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Pair;

/**
 * Base class for rules that determines the integration nodes and their
 * weights.
 * Subclasses must implement the {@link #computeRule(int) computeRule} method.
 *
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 *
 */
public abstract class BaseRuleFactory<T extends Number> {
    /** List of points and weights, indexed by the order of the rule. */
    private final SortedMap<Integer, Pair<T[], T[]>> pointsAndWeights = new TreeMap<>();
    /** Cache for double-precision rules. */
    private final SortedMap<Integer, Pair<double[], double[]>> pointsAndWeightsDouble = new TreeMap<>();

    /**
     * Gets a copy of the quadrature rule with the given number of integration
     * points.
     * The number of points is arbitrarily limited to 1000. It prevents resources
     * exhaustion. In practice the number of points is often much lower.
     *
     * @param numberOfPoints Number of integration points.
     * @return a copy of the integration rule.
     * @throws MathIllegalArgumentException if {@code numberOfPoints < 1}.
     * @throws MathIllegalArgumentException if {@code numberOfPoints > 1000}.
     * @throws MathIllegalArgumentException if the elements of the rule pair do not
     * have the same length.
     */
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

        // Try to obtain the rule from the cache.
        Pair<double[], double[]> cached = pointsAndWeightsDouble.get(numberOfPoints);

        if (cached == null) {
            // Rule not computed yet.

            // Compute the rule.
            final Pair<T[], T[]> rule = getRuleInternal(numberOfPoints);
            cached = convertToDouble(rule);

            // Cache it.
            pointsAndWeightsDouble.put(numberOfPoints, cached);
        }

        // Return a copy.
        return new Pair<double[], double[]>(cached.getFirst().clone(),
                                            cached.getSecond().clone());
    }

    /**
     * Gets a rule.
     * Synchronization ensures that rules will be computed and added to the
     * cache at most once.
     * The returned rule is a reference into the cache.
     *
     * @param numberOfPoints Order of the rule to be retrieved.
     * @return the points and weights corresponding to the given order.
     * @throws MathIllegalArgumentException if the elements of the rule pair do not
     * have the same length.
     */
    protected Pair<T[], T[]> getRuleInternal(int numberOfPoints)
        throws MathIllegalArgumentException {
        final Pair<T[], T[]> rule;
        synchronized (pointsAndWeights) {
            rule = pointsAndWeights.get(numberOfPoints);
            if (rule == null) {
                addRule(computeRule(numberOfPoints));
                // The rule should be available now.
                return getRuleInternal(numberOfPoints);
            }
            return rule;
        }
    }

    /**
     * Stores a rule.
     *
     * @param rule Rule to be stored.
     * @throws MathIllegalArgumentException if the elements of the pair do not
     * have the same length.
     */
    protected void addRule(Pair<T[], T[]> rule) throws MathIllegalArgumentException {
        MathUtils.checkDimension(rule.getFirst().length, rule.getSecond().length);
        synchronized (pointsAndWeights) {
            pointsAndWeights.put(rule.getFirst().length, rule);
        }
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

    /**
     * Converts the from the actual {@code Number} type to {@code double}
     *
     * @param <T> Type of the number used to represent the points and
     * weights of the quadrature rules.
     * @param rule Points and weights.
     * @return points and weights as {@code double}s.
     */
    private static <T extends Number> Pair<double[], double[]> convertToDouble(Pair<T[], T[]> rule) {
        final T[] pT = rule.getFirst();
        final T[] wT = rule.getSecond();

        final int len = pT.length;
        final double[] pD = new double[len];
        final double[] wD = new double[len];

        for (int i = 0; i < len; i++) {
            pD[i] = pT[i].doubleValue();
            wD[i] = wT[i].doubleValue();
        }

        return new Pair<double[], double[]>(pD, wD);
    }
}
