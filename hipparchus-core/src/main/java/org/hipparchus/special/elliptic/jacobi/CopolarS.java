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

/** Copolar trio with pole at point s in Glaisher’s Notation.
 * <p>
 * This is a container for the three subsidiary Jacobi elliptic functions
 * {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
 * </p>
 * @since 2.0
 */
public class CopolarS {

    /** Value of the cs function. */
    private final double cs;

    /** Value of the dn function. */
    private final double ds;

    /** Value of the ns function. */
    private final double ns;

    /** Simple constructor.
     * @param trioN copolar trio with pole at point n in Glaisher’s Notation
     */
    CopolarS(final CopolarN trioN) {
        this.ns = 1.0 / trioN.sn();
        this.cs = ns  * trioN.cn();
        this.ds = ns  * trioN.dn();
    }

    /** Get the value of the cs function.
     * @return cs(u|m)
     */
    public double cs() {
        return cs;
    }

    /** Get the value of the ds function.
     * @return ds(u|m)
     */
    public double ds() {
        return ds;
    }

    /** Get the value of the ns function.
     * @return ns(u|m)
     */
    public double ns() {
        return ns;
    }

}
