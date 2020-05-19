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
package org.hipparchus.analysis.differentiation;

import java.io.Serializable;

import org.hipparchus.RealFieldElement;

/** Abstract class representing both the value and the differentials of a function.
 * @since 1.7
 */
public abstract class UnivariateDerivative<T extends UnivariateDerivative<T>>
    implements RealFieldElement<T>, Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20200519L;

    /** {@inheritDoc} */
    @Override
    public double getReal() {
        return getValue();
    }

    /** Get the value part of the univariate derivative.
     * @return value part of the univariate derivative
     */
    public abstract double getValue();

    /** Get a derivative from the univariate derivative.
     * @param n derivation order (must be between 0 and {@link #getOrder()}, both inclusive)
     * @return n<sup>th</sup> derivative, or {@code Double.NaN} if n is
     * either negative or strictly larger than {@link #getOrder()}
     */
    public abstract double getDerivative(int n);

    /** Get the derivation order.
     * @return derivation order
     */
    public abstract int getOrder();

}
