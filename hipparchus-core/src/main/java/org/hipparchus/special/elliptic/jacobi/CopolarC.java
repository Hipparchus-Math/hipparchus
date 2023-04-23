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

/** Copolar trio with pole at point c in Glaisher’s Notation.
 * <p>
 * This is a container for the three subsidiary Jacobi elliptic functions
 * {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
 * </p>
 * @since 2.0
 */
public class CopolarC {

    /** Value of the dc function. */
    private final double dc;

    /** Value of the nc function. */
    private final double nc;

    /** Value of the sc function. */
    private final double sc;

    /** Simple constructor.
     * @param trioN copolar trio with pole at point n in Glaisher’s Notation
     */
    CopolarC(final CopolarN trioN) {
        this.nc = 1.0 / trioN.cn();
        this.sc = nc  * trioN.sn();
        this.dc = nc  * trioN.dn();
    }

    /** Get the value of the dc function.
     * @return dc(u|m)
     */
    public double dc() {
        return dc;
    }

    /** Get the value of the nc function.
     * @return nc(u|m)
     */
    public double nc() {
        return nc;
    }

    /** Get the value of the sc function.
     * @return sc(u|m)
     */
    public double sc() {
        return sc;
    }

}
