/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.util;

import java.io.Serializable;

import org.hipparchus.Field;

/**
 * Representation of real numbers with arbitrary precision field.
 * <p>
 * This class is a singleton.
 *
 * @see BigReal
 */
public class BigRealField implements Field<BigReal>, Serializable  {

    /** Serializable version identifier */
    private static final long serialVersionUID = 20160327L;

    /** Private constructor for the singleton.
     */
    private BigRealField() {
    }

    /** Get the unique instance.
     * @return the unique instance
     */
    public static BigRealField getInstance() {
        return LazyHolder.INSTANCE;
    }

    /** {@inheritDoc} */
    @Override
    public BigReal getOne() {
        return BigReal.ONE;
    }

    /** {@inheritDoc} */
    @Override
    public BigReal getZero() {
        return BigReal.ZERO;
    }

    /** {@inheritDoc} */
    @Override
    public Class<BigReal> getRuntimeClass() {
        return BigReal.class;
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object other) {
        return this == other;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return 0xf38b3541;
    }

    // CHECKSTYLE: stop HideUtilityClassConstructor
    /** Holder for the instance.
     * <p>We use here the Initialization On Demand Holder Idiom.</p>
     */
    private static class LazyHolder {
        /** Cached field instance. */
        private static final BigRealField INSTANCE = new BigRealField();
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
