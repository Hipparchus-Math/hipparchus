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

import java.io.Serializable;

import org.hipparchus.Field;

/** Field for {@link UnivariateDerivative1} instances.
 * <p>
 * This class is a singleton.
 * </p>
 * @since 1.7
 */
public class UnivariateDerivative1Field implements Field<UnivariateDerivative1>, Serializable {

    /** Serializable version identifier. */
    private static final long serialVersionUID = 20200519L;

    /** Zero constant. */
    private final UnivariateDerivative1 zero;

    /** One constant. */
    private final UnivariateDerivative1 one;

    /** Associated factory for conversions to {@link DerivativeStructure}. */
    private final DSFactory factory;

    /** Private constructor for the singleton.
     */
    private UnivariateDerivative1Field() {
        zero    = new UnivariateDerivative1(0.0, 0.0);
        one     = new UnivariateDerivative1(1.0, 0.0);
        factory = new DSFactory(1, 1);
    }

    /** Get the unique instance.
     * @return the unique instance
     */
    public static UnivariateDerivative1Field getInstance() {
        return LazyHolder.INSTANCE;
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 getOne() {
        return one;
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 getZero() {
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
     * @return factory for conversions
     */
    DSFactory getConversionFactory() {
        return factory;
    }

    /** {@inheritDoc} */
    @Override
    public Class<UnivariateDerivative1> getRuntimeClass() {
        return UnivariateDerivative1.class;
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object other) {
        return this == other;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return 0xd5e174f7;
    }

    // CHECKSTYLE: stop HideUtilityClassConstructor
    /** Holder for the instance.
     * <p>We use here the Initialization On Demand Holder Idiom.</p>
     */
    private static class LazyHolder {
        /** Cached field instance. */
        private static final UnivariateDerivative1Field INSTANCE = new UnivariateDerivative1Field();
    }
    // CHECKSTYLE: resume HideUtilityClassConstructor

    /** Handle deserialization of the singleton.
     * @return the singleton instance
     */
    private Object readResolve() {
        // return the singleton instance
        return LazyHolder.INSTANCE;
    }

}
