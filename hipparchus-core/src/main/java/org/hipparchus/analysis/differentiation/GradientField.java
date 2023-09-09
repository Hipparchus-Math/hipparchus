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

import java.util.concurrent.atomic.AtomicReference;

import org.hipparchus.Field;
import org.hipparchus.util.FastMath;

/** Field for {@link Gradient} instances.
 * @since 1.7
 */
public class GradientField implements Field<Gradient> {

    /** Array of all fields created so far. */
    private static AtomicReference<GradientField[]> fields = new AtomicReference<>(null);

    /** Zero constant. */
    private final Gradient zero;

    /** One constant. */
    private final Gradient one;

    /** Associated factory for conversions to {@link DerivativeStructure}. */
    private final DSFactory factory;

    /** Private constructor.
     * @param parameters number of free parameters
     */
    private GradientField(final int parameters) {
        zero    = new Gradient(0.0, new double[parameters]);
        one     = new Gradient(1.0, new double[parameters]);
        factory = new DSFactory(parameters, 1);
    }

    /** Get the field for number of free parameters.
     * @param parameters number of free parameters
     * @return cached field
     */
    public static GradientField getField(int parameters) {

        // get the cached fields
        final GradientField[] cache = fields.get();
        if (cache != null && cache.length > parameters && cache[parameters] != null) {
            // the field has already been created
            return cache[parameters];
        }

        // we need to create a new field
        final int maxParameters = FastMath.max(parameters, cache == null ? 0 : cache.length);
        final GradientField[] newCache = new GradientField[maxParameters + 1];

        if (cache != null) {
            // preserve the already created fields
            System.arraycopy(cache, 0, newCache, 0, cache.length);
        }

        // create the new field
        newCache[parameters] = new GradientField(parameters);

        // atomically reset the cached fileds array
        fields.compareAndSet(cache, newCache);

        return newCache[parameters];

    }

    /** {@inheritDoc} */
    @Override
    public Gradient getOne() {
        return one;
    }

    /** {@inheritDoc} */
    @Override
    public Gradient getZero() {
        return zero;
    }

    /** {@inheritDoc} */
    @Override
    public Class<Gradient> getRuntimeClass() {
        return Gradient.class;
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
    DSFactory getConversionFactory() {
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
        return 0x26ca1af0;
    }

}
