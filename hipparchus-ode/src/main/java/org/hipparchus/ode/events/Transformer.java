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

package org.hipparchus.ode.events;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;


/** Transformer for {@link ODEEventHandler#g(org.hipparchus.ode.ODEStateAndDerivative) g functions}.
 * @see EventSlopeFilter
 * @see FilterType
 */
enum Transformer {

    /** Transformer computing transformed = 0.
     * <p>
     * This transformer is used when we initialize the filter, until we get at
     * least one non-zero value to select the proper transformer.
     * </p>
     */
    UNINITIALIZED {

        /**  {@inheritDoc} */
        @Override
        protected double transformed(final double g) {
            return 0;
        }

        /**  {@inheritDoc} */
        @Override
        protected <T extends CalculusFieldElement<T>> T transformed(final T g) {
            return g.getField().getZero();
        }

    },

    /** Transformer computing transformed = g.
     * <p>
     * When this transformer is applied, the roots of the original function
     * are preserved, with the same {@code increasing/decreasing} status.
     * </p>
     */
    PLUS {

        /**  {@inheritDoc} */
        @Override
        protected double transformed(final double g) {
            return g;
        }

        /**  {@inheritDoc} */
        @Override
        protected <T extends CalculusFieldElement<T>> T transformed(final T g) {
            return g;
        }

    },

    /** Transformer computing transformed = -g.
     * <p>
     * When this transformer is applied, the roots of the original function
     * are preserved, with reversed {@code increasing/decreasing} status.
     * </p>
     */
    MINUS {

        /**  {@inheritDoc} */
        @Override
        protected double transformed(final double g) {
            return -g;
        }

        /**  {@inheritDoc} */
        @Override
        protected <T extends CalculusFieldElement<T>> T transformed(final T g) {
            return g.negate();
        }

    },

    /** Transformer computing transformed = min(-{@link Precision#SAFE_MIN}, -g, +g).
     * <p>
     * When this transformer is applied, the transformed function is
     * guaranteed to be always strictly negative (i.e. there are no roots).
     * </p>
     */
    MIN {

        /**  {@inheritDoc} */
        @Override
        protected double transformed(final double g) {
            return FastMath.min(FastMath.min(-g, +g), -Precision.SAFE_MIN);
        }

        /**  {@inheritDoc} */
        @Override
        protected <T extends CalculusFieldElement<T>> T transformed(final T g) {
            return FastMath.min(FastMath.min(g.negate(), g), -Precision.SAFE_MIN);
        }

    },

    /** Transformer computing transformed = max(+{@link Precision#SAFE_MIN}, -g, +g).
     * <p>
     * When this transformer is applied, the transformed function is
     * guaranteed to be always strictly positive (i.e. there are no roots).
     * </p>
     */
    MAX {

        /**  {@inheritDoc} */
        @Override
        protected double transformed(final double g) {
            return FastMath.max(FastMath.max(-g, +g), Precision.SAFE_MIN);
        }

        /**  {@inheritDoc} */
        @Override
        protected <T extends CalculusFieldElement<T>> T transformed(final T g) {
            return FastMath.max(FastMath.max(g.negate(), g), Precision.SAFE_MIN);
        }

    };

    /** Transform value of function g.
     * @param g raw value of function g
     * @return transformed value of function g
     */
    protected abstract double transformed(double g);

    /** Transform value of function g.
     * @param g raw value of function g
     * @return transformed value of function g
     * @param <T> the type of the field elements
     * @since 2.0
     */
    protected abstract <T extends CalculusFieldElement<T>> T transformed(T g);

}
