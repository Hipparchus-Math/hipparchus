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
package org.hipparchus.analysis.differentiation;

import java.util.HashMap;
import java.util.Map;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;

/** Field for {@link FieldUnivariateDerivative1} instances.
 * @param <T> the type of the function parameters and value
 * @since 1.7
 */
public class FieldUnivariateDerivative1Field<T extends CalculusFieldElement<T>> implements Field<FieldUnivariateDerivative1<T>> {

    /** Cached fields. */
    private static final Map<Field<?>, FieldUnivariateDerivative1Field<?>> CACHE = new HashMap<>();

    /** Zero constant. */
    private final FieldUnivariateDerivative1<T> zero;

    /** One constant. */
    private final FieldUnivariateDerivative1<T> one;

    /** Associated factory for conversions to {@link FieldDerivativeStructure}. */
    private final FDSFactory<T> factory;

    /** Private constructor for populating the cache.
     * @param valueField field for the function parameters and value
     */
    private FieldUnivariateDerivative1Field(final Field<T> valueField) {
        zero    = new FieldUnivariateDerivative1<>(valueField.getZero(), valueField.getZero());
        one     = new FieldUnivariateDerivative1<>(valueField.getOne(), valueField.getZero());
        factory = new FDSFactory<>(valueField, 1, 1);
    }

    /** Get the univariate derivative field corresponding to a value field.
     * @param valueField field for the function parameters and value
     * @param <T> the type of the function parameters and value
     * @return univariate derivative field
     */
    public static <T extends CalculusFieldElement<T>> FieldUnivariateDerivative1Field<T> getUnivariateDerivative1Field(final Field<T> valueField) {
        synchronized (CACHE) {
            FieldUnivariateDerivative1Field<?> cached = CACHE.get(valueField);
            if (cached == null) {
                cached = new FieldUnivariateDerivative1Field<>(valueField);
                CACHE.put(valueField, cached);
            }
            @SuppressWarnings("unchecked")
            final FieldUnivariateDerivative1Field<T> tCached = (FieldUnivariateDerivative1Field<T>) cached;
            return tCached;
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> getOne() {
        return one;
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> getZero() {
        return zero;
    }

    /** Get the factory for converting to {@link DerivativeStructure}.
     * <p>
     * This factory is used only for conversions. {@code UnivariateDerivative1} by
     * itself does not rely at all on {@link DSFactory}, {@link DSCompiler}
     * or {@link DerivativeStructure} for its computation. For this reason,
     * the factory here is hidden and this method is package private, so
     * only {@link UnivariateDerivative1#toDerivativeStructure()} can call it on an
     * existing {@link UnivariateDerivative1} instance
     * </p>
     * @return factory for conversion
     */
    FDSFactory<T> getConversionFactory() {
        return factory;
    }

    /** {@inheritDoc} */
    @SuppressWarnings("unchecked")
    @Override
    public Class<FieldUnivariateDerivative1<T>> getRuntimeClass() {
        return (Class<FieldUnivariateDerivative1<T>>) zero.getClass();
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object other) {
        return this == other;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return 0x712c0fc7;
    }

}
