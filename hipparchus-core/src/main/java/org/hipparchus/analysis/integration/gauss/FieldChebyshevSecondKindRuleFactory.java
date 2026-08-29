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
import org.hipparchus.Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.Pair;

/**
 * Factory that creates field-based Gauss-Chebyshev quadrature rules of the
 * second kind.
 *
 * <p>The generated rules approximate integrals of the form:</p>
 *
 * <pre>
 *     integral from -1 to 1 of f(x) sqrt(1 - x^2) dx
 * </pre>
 *
 * <p>using:</p>
 *
 * <pre>
 *     sum from i = 1 to n of w_i f(x_i)
 * </pre>
 *
 * <p>where:</p>
 *
 * <pre>
 *     x_i = cos(i pi / (n + 1)), i = 1, ..., n
 *     w_i = pi / (n + 1) sin^2(i pi / (n + 1))
 * </pre>
 *
 * @param <T> type of the field elements
 */
public class FieldChebyshevSecondKindRuleFactory<T extends CalculusFieldElement<T>>
        extends FieldAbstractRuleFactory<T> {

    /**
     * Constructor.
     *
     * @param field field to which rule coefficients belong
     */
    public FieldChebyshevSecondKindRuleFactory(final Field<T> field) {
        super(field);
    }

    /**
     * Computes the Gauss-Chebyshev quadrature rule of the second kind.
     *
     * @param numberOfPoints order of the rule to be computed
     * @return nodes and weights of the quadrature rule
     */
    @Override
    protected Pair<T[], T[]> computeRule(final int numberOfPoints) {
        final T[] points = MathArrays.buildArray(getField(), numberOfPoints);
        final T[] weights = MathArrays.buildArray(getField(), numberOfPoints);

        final T scale = getField().getZero().newInstance(FastMath.PI / (numberOfPoints + 1.0));

        for (int i = 0; i < numberOfPoints; i++) {
            final T angle = scale.multiply(i + 1.0);
            final T sin = FastMath.sin(angle);
            final int index = numberOfPoints - 1 - i;

            points[index] = FastMath.cos(angle);
            weights[index] = scale.multiply(sin.square());
        }

        return new Pair<>(points, weights);
    }
}
