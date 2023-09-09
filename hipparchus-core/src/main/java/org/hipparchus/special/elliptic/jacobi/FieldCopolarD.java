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

/** Copolar trio with pole at point d in Glaisher’s Notation.
 * <p>
 * This is a container for the three subsidiary Jacobi elliptic functions
 * {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
public class FieldCopolarD<T extends CalculusFieldElement<T>> {

    /** Value of the nd function. */
    private final T nd;

    /** Value of the sd function. */
    private final T sd;

    /** Value of the cd function. */
    private final T cd;

    /** Simple constructor.
     * @param trioN copolar trio with pole at point n in Glaisher’s Notation
     */
    FieldCopolarD(final FieldCopolarN<T> trioN) {
        this.nd = trioN.dn().reciprocal();
        this.sd = nd.multiply(trioN.sn());
        this.cd = nd.multiply(trioN.cn());
    }

    /** Get the value of the nd function.
     * @return nd(u|m)
     */
    public T nd() {
        return nd;
    }

    /** Get the value of the sd function.
     * @return sd(u|m)
     */
    public T sd() {
        return sd;
    }

    /** Get the value of the cd function.
     * @return cd(u|m)
     */
    public T cd() {
        return cd;
    }

}
