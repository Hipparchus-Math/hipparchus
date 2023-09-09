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

import java.lang.reflect.Array;
import java.util.HashMap;
import java.util.Map;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.util.MathArrays;

/** Field for {@link Gradient} instances.
 * @param <T> the type of the function parameters and value
 * @since 1.7
 */
public class FieldGradientField<T extends CalculusFieldElement<T>> implements Field<FieldGradient<T>> {

    /** Cached fields. */
    private static final Map<Field<?>, FieldGradientField<?>[]> CACHE = new HashMap<>();

    /** Zero constant. */
    private final FieldGradient<T> zero;

    /** One constant. */
    private final FieldGradient<T> one;

    /** Associated factory for conversions to {@link DerivativeStructure}. */
    private final FDSFactory<T> factory;

    /** Private constructor.
     * @param valueField field for the function parameters and value
     * @param parameters number of free parameters
     */
    private FieldGradientField(final Field<T> valueField, final int parameters) {
        zero    = new FieldGradient<>(valueField.getZero(), MathArrays.buildArray(valueField, parameters));
        one     = new FieldGradient<>(valueField.getOne(), MathArrays.buildArray(valueField, parameters));
        factory = new FDSFactory<>(valueField, parameters, 1);
    }

    /** Get the field for number of free parameters.
     * @param valueField field for the function parameters and value
     * @param parameters number of free parameters
     * @param <T> the type of the function parameters and value
     * @return cached field
     */
    public static <T extends CalculusFieldElement<T>> FieldGradientField<T> getField(final Field<T> valueField, final int parameters) {

        FieldGradientField<?>[] cachedFields;
        synchronized (CACHE) {
            cachedFields = CACHE.get(valueField);
            if (cachedFields == null || cachedFields.length <= parameters) {
                FieldGradientField<?>[] newCachedFields =
                                (FieldGradientField<?>[]) Array.newInstance(FieldGradientField.class, parameters + 1);
                if (cachedFields != null) {
                    // preserve the already created fields
                    System.arraycopy(cachedFields, 0, newCachedFields, 0, cachedFields.length);
                }
                cachedFields = newCachedFields;
                CACHE.put(valueField, cachedFields);
            }
        }

        if (cachedFields[parameters] == null) {
            // we need to create a new field
            cachedFields[parameters] = new FieldGradientField<>(valueField, parameters);
        }

        @SuppressWarnings("unchecked")
        final FieldGradientField<T> tCached = (FieldGradientField<T>) cachedFields[parameters];
        return tCached;

    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T>  getOne() {
        return one;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T>  getZero() {
        return zero;
    }

    /** {@inheritDoc} */
    @SuppressWarnings("unchecked")
    @Override
    public Class<FieldGradient<T>> getRuntimeClass() {
        return (Class<FieldGradient<T>>) getZero().getClass();
    }

    /** Get the factory for converting to {@link DerivativeStructure}.
     * <p>
     * This factory is used only for conversions. {@code Gradient} by
     * itself does not rely at all on {@link DSFactory}, {@link DSCompiler}
     * or {@link DerivativeStructure} for its computation. For this reason,
     * the factory here is hidden and this method is package private, so
     * only {@link Gradient#toDerivativeStructure()} can call it on an
     * existing {@link Gradient} instance
     * </p>
     * @return factory for conversions
     */
    FDSFactory<T> getConversionFactory() {
        return factory;
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object other) {
        return this == other;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return 0xcd3e92ee;
    }

}
