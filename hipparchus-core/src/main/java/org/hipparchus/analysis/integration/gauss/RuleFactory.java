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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Pair;

/** Interface for rules that determines the integration nodes and their weights.
 * @since 2.0
 */
public interface RuleFactory {

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
    Pair<double[], double[]> getRule(int numberOfPoints) throws MathIllegalArgumentException;

}
