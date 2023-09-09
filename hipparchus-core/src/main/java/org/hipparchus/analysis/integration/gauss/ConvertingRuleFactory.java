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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.FieldElement;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Pair;

/**
 * Factory converting {@link CalculusFieldElement field-based} {@link FieldRuleFactory} into {@link RuleFactory}.
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 * @since 2.0
 */
public class ConvertingRuleFactory<T extends FieldElement<T>> extends AbstractRuleFactory {

    /** Underlying field-based factory. */
    private final FieldRuleFactory<T> fieldFactory;

    /** Simple constructor.
     * @param fieldFactory field-based factory to convert
     */
    public ConvertingRuleFactory(final FieldRuleFactory<T> fieldFactory) {
        this.fieldFactory = fieldFactory;
    }

    /** {@inheritDoc} */
    @Override
    protected Pair<double[], double[]> computeRule(final int numberOfPoints)
        throws MathIllegalArgumentException {

        // get the field-based rule
        Pair<T[], T[]> rule = fieldFactory.getRule(numberOfPoints);

        // convert the nodes and weights
        final T[] pT = rule.getFirst();
        final T[] wT = rule.getSecond();

        final int len = pT.length;
        final double[] pD = new double[len];
        final double[] wD = new double[len];

        for (int i = 0; i < len; i++) {
            pD[i] = pT[i].getReal();
            wD[i] = wT[i].getReal();
        }

        return new Pair<>(pD, wD);

    }

}
