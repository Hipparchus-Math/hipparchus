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
package org.hipparchus.util;

import org.hipparchus.CalculusFieldElement;

/** Holder for both sine and cosine values.
 * <p>
 * This class is a simple container, it does not provide any computational method.
 * </p>
 * @see FastMath#sinCos(org.hipparchus.CalculusFieldElement)
 * @param <T> the type of the field elements
 * @since 1.4
 */
public class FieldSinCos<T> {

    /** Value of the sine. */
    private final T sin;

    /** Value of the cosine. */
    private final T cos;

    /** Simple constructor.
     * @param sin value of the sine
     * @param cos value of the cosine
     */
    public FieldSinCos(final T sin, final T cos) {
        this.sin = sin;
        this.cos = cos;
    }

    /** Get the value of the sine.
     * @return value of the sine
     */
    public T sin() {
        return sin;
    }

    /** Get the value of the cosine.
     * @return value of the cosine
     */
    public T cos() {
        return cos;
    }

    /** Compute sine and cosine of angles sum.
     * @param scAlpha \((\sin \alpha, \cos \alpha)\)
     * @param scBeta \((\sin \beta, \cos \beta)\)
     * @param <S> the type of the field elements
     * @return \((\sin \alpha+\beta, \cos \alpha+\beta)\)
     * @since 1.8
     */
    public static <S extends CalculusFieldElement<S>> FieldSinCos<S> sum(final FieldSinCos<S> scAlpha, final FieldSinCos<S> scBeta) {
        return new FieldSinCos<>(scAlpha.sin.linearCombination(scAlpha.sin, scBeta.cos, scAlpha.cos,          scBeta.sin),
                                 scAlpha.sin.linearCombination(scAlpha.cos, scBeta.cos, scAlpha.sin.negate(), scBeta.sin));
    }

    /** Compute sine and cosine of angles difference.
     * @param scAlpha \((\sin \alpha, \cos \alpha)\)
     * @param scBeta \((\sin \beta, \cos \beta)\)
     * @param <S> the type of the field elements
     * @return \((\sin \alpha+\beta, \cos \alpha-\beta)\)
     * @since 1.8
     */
    public static <S extends CalculusFieldElement<S>> FieldSinCos<S> difference(final FieldSinCos<S> scAlpha, final FieldSinCos<S> scBeta) {
        return new FieldSinCos<>(scAlpha.sin.linearCombination(scAlpha.sin, scBeta.cos, scAlpha.cos.negate(), scBeta.sin),
                                 scAlpha.sin.linearCombination(scAlpha.cos, scBeta.cos, scAlpha.sin,          scBeta.sin));
    }

}
