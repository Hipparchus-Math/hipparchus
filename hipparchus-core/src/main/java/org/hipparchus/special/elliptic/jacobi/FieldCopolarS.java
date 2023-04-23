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
package org.hipparchus.special.elliptic.jacobi;

import org.hipparchus.CalculusFieldElement;

/** Copolar trio with pole at point s in Glaisher’s Notation.
 * <p>
 * This is a container for the three subsidiary Jacobi elliptic functions
 * {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
public class FieldCopolarS<T extends CalculusFieldElement<T>> {

    /** Value of the cs function. */
    private final T cs;

    /** Value of the dn function. */
    private final T ds;

    /** Value of the ns function. */
    private final T ns;

    /** Simple constructor.
     * @param trioN copolar trio with pole at point n in Glaisher’s Notation
     */
    FieldCopolarS(final FieldCopolarN<T> trioN) {
        this.ns = trioN.sn().reciprocal();
        this.cs = ns.multiply(trioN.cn());
        this.ds = ns.multiply(trioN.dn());
    }

    /** Get the value of the cs function.
     * @return cs(u|m)
     */
    public T cs() {
        return cs;
    }

    /** Get the value of the ds function.
     * @return ds(u|m)
     */
    public T ds() {
        return ds;
    }

    /** Get the value of the ns function.
     * @return ns(u|m)
     */
    public T ns() {
        return ns;
    }

}
