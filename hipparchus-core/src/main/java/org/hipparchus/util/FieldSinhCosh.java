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

/** Holder for both hyperbolic sine and hyperbolic cosine values.
 * <p>
 * This class is a simple container, it does not provide any computational method.
 * </p>
 * @see FastMath#sinhCosh(double)
 * @param <T> the type of the field elements
 * @since 2.0
 */
public class FieldSinhCosh<T> {

    /** Value of the hyperbolic sine. */
    private final T sinh;

    /** Value of the hyperbolic cosine. */
    private final T cosh;

    /** Simple constructor.
     * @param sinh value of the hyperbolic sine
     * @param cosh value of the hyperbolic cosine
     */
    public FieldSinhCosh(final T sinh, final T cosh) {
        this.sinh = sinh;
        this.cosh = cosh;
    }

    /** Get the value of the hyperbolic sine.
     * @return value of the hyperbolic sine
     */
    public T sinh() {
        return sinh;
    }

    /** Get the value of the hyperbolic cosine.
     * @return value of the hyperbolic cosine
     */
    public T cosh() {
        return cosh;
    }

    /** Compute hyperbolic sine and hyperbolic cosine of angles sum.
     * @param schAlpha \((\sinh \alpha, \cosh \alpha)\)
     * @param schBeta \((\sinh \beta, \cosh \beta)\)
     * @param <S> the type of the field elements
     * @return \((\sinh \alpha+\beta, \cosh \alpha+\beta)\)
     */
    public static <S extends CalculusFieldElement<S>> FieldSinhCosh<S> sum(final FieldSinhCosh<S> schAlpha, final FieldSinhCosh<S> schBeta) {
        return new FieldSinhCosh<>(schAlpha.sinh.linearCombination(schAlpha.sinh, schBeta.cosh, schAlpha.cosh, schBeta.sinh),
                                   schAlpha.sinh.linearCombination(schAlpha.cosh, schBeta.cosh, schAlpha.sinh, schBeta.sinh));
    }

    /** Compute hyperbolic sine and hyperbolic cosine of angles difference.
     * @param schAlpha \((\sinh \alpha, \cosh \alpha)\)
     * @param schBeta \((\sinh \beta, \cosh \beta)\)
     * @param <S> the type of the field elements
     * @return \((\sinh \alpha+\beta, \cosh \alpha-\beta)\)
     */
    public static <S extends CalculusFieldElement<S>> FieldSinhCosh<S> difference(final FieldSinhCosh<S> schAlpha, final FieldSinhCosh<S> schBeta) {
        final S mShB = schBeta.sinh.negate();
        return new FieldSinhCosh<>(schAlpha.sinh.linearCombination(schAlpha.sinh, schBeta.cosh, schAlpha.cosh, mShB),
                                   schAlpha.sinh.linearCombination(schAlpha.cosh, schBeta.cosh, schAlpha.sinh, mShB));
    }

}
