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
package org.hipparchus.analysis.differentiation;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;

/** Interface for univariate functions derivatives.
 * <p>This interface represents a simple function which computes
 * both the value and the first derivative of a mathematical function.
 * The derivative is computed with respect to the input variable.</p>
 * @see UnivariateDifferentiableFunction
 * @see UnivariateFunctionDifferentiator
 */
public interface UnivariateDifferentiableFunction extends UnivariateFunction {

    /**
     * Compute the value for the function.
     * @param x the point for which the function value should be computed
     * @param <T> the type of the field elements
     * @return the value
     * @exception MathIllegalArgumentException if {@code x} does not
     * satisfy the function's constraints (argument out of bound, or unsupported
     * derivative order for example)
     */
    <T extends Derivative<T>> T value(T x) throws MathIllegalArgumentException;

}
