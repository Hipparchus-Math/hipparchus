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

import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Interface representing classes that can blend with other instances of themselves using a given blending value.
 * <p> The blending value is commonly given from a
 * {@link org.hipparchus.analysis.polynomials.SmoothStepFactory.SmoothStepFunction smoothstep function}.
 *
 * @param <B> blendable class
 */
public interface Blendable<B> {

    /**
     * Blend arithmetically this instance with another one.
     *
     * @param other other instance to blend arithmetically with
     * @param blendingValue value from smoothstep function B(x). It is expected to be between [0:1] and will
     *         throw an exception otherwise.
     * @return this * (1 - B(x)) + other * B(x)
     * @throws MathIllegalArgumentException if blending value is not within [0:1]
     */
    B blendArithmeticallyWith(B other, double blendingValue) throws MathIllegalArgumentException;
}
