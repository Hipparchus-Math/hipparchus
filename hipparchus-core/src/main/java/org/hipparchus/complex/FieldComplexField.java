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
package org.hipparchus.complex;

import java.util.HashMap;
import java.util.Map;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;

/**
 * Representation of the complex numbers field.
 *
 * @param <T> the type of the field elements
 * @see FieldComplex
 * @since 2.0
 */
public class FieldComplexField<T extends CalculusFieldElement<T>> implements Field<FieldComplex<T>> {

    /** Cached fields. */
    private static final Map<Field<?>, FieldComplexField<?>> CACHE = new HashMap<>();

    /** Constant 0. */
    private final FieldComplex<T> zero;

    /** Constant 1. */
    private final FieldComplex<T> one;

    /** Simple constructor.
     * @param field type of the field element
     */
    private FieldComplexField(final Field<T> field) {
        zero = FieldComplex.getZero(field);
        one  = FieldComplex.getOne(field);
    }

    /** Get the field for complex numbers.
     * @param partsField field for the real and imaginary parts
     * @param <T> the type of the field elements
     * @return cached field
     */
    public static <T extends CalculusFieldElement<T>> FieldComplexField<T> getField(final Field<T> partsField) {
        FieldComplexField<?> cachedField;
        synchronized (CACHE) {
            cachedField = CACHE.get(partsField);
            if (cachedField == null) {
                cachedField = new FieldComplexField<>(partsField);
                CACHE.put(partsField, cachedField);
            }
        }

        @SuppressWarnings("unchecked")
        final FieldComplexField<T> tCached = (FieldComplexField<T>) cachedField;
        return tCached;

    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> getOne() {
        return one;
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> getZero() {
        return zero;
    }

    /** {@inheritDoc} */
    @SuppressWarnings("unchecked")
    @Override
    public Class<FieldComplex<T>> getRuntimeClass() {
        return (Class<FieldComplex<T>>) getZero().getClass();
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object other) {
        return this == other;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return 0xd368f208;
    }

}
