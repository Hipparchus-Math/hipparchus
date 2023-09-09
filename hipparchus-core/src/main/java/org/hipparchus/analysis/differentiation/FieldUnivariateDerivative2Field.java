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

/** Field for {@link FieldUnivariateDerivative2} instances.
 * @param <T> the type of the function parameters and value
 * @since 1.7
 */
public class FieldUnivariateDerivative2Field<T extends CalculusFieldElement<T>> implements Field<FieldUnivariateDerivative2<T>> {

    /** Cached fields. */
    private static final Map<Field<?>, FieldUnivariateDerivative2Field<?>> CACHE = new HashMap<>();

    /** Zero constant. */
    private final FieldUnivariateDerivative2<T> zero;

    /** One constant. */
    private final FieldUnivariateDerivative2<T> one;

    /** Associated factory for conversions to {@link FieldDerivativeStructure}. */
    private final FDSFactory<T> factory;

    /** Private constructor for populating the cache.
     * @param valueField field for the function parameters and value
     */
    private FieldUnivariateDerivative2Field(final Field<T> valueField) {
        zero    = new FieldUnivariateDerivative2<>(valueField.getZero(), valueField.getZero(), valueField.getZero());
        one     = new FieldUnivariateDerivative2<>(valueField.getOne(), valueField.getZero(), valueField.getZero());
        factory = new FDSFactory<>(valueField, 1, 2);
    }

    /** Get the univariate derivative field corresponding to a value field.
     * @param valueField field for the function parameters and value
     * @param <T> the type of the function parameters and value
     * @return univariate derivative field
     */
    public static <T extends CalculusFieldElement<T>> FieldUnivariateDerivative2Field<T> getUnivariateDerivative2Field(final Field<T> valueField) {
        synchronized (CACHE) {
            FieldUnivariateDerivative2Field<?> cached = CACHE.get(valueField);
            if (cached == null) {
                cached = new FieldUnivariateDerivative2Field<>(valueField);
                CACHE.put(valueField, cached);
            }
            @SuppressWarnings("unchecked")
            final FieldUnivariateDerivative2Field<T> tCached = (FieldUnivariateDerivative2Field<T>) cached;
            return tCached;
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> getOne() {
        return one;
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> getZero() {
        return zero;
    }

    /** Get the factory for converting to {@link DerivativeStructure}.
     * <p>
     * This factory is used only for conversions. {@code UnivariateDerivative2} by
     * itself does not rely at all on {@link DSFactory}, {@link DSCompiler}
     * or {@link DerivativeStructure} for its computation. For this reason,
     * the factory here is hidden and this method is package private, so
     * only {@link UnivariateDerivative2#toDerivativeStructure()} can call it on an
     * existing {@link UnivariateDerivative2} instance
     * </p>
     * @return factory for conversions
     */
    FDSFactory<T> getConversionFactory() {
        return factory;
    }

    /** {@inheritDoc} */
    @SuppressWarnings("unchecked")
    @Override
    public Class<FieldUnivariateDerivative2<T>> getRuntimeClass() {
        return (Class<FieldUnivariateDerivative2<T>>) zero.getClass();
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object other) {
        return this == other;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return 0x3f4b793e;
    }

}
