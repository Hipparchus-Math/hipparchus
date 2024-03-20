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

package org.hipparchus.analysis.polynomials;


/** Class for handling Jacobi polynomials keys.
 * @since 3.1
 */
public class JacobiKey {

    /** First exponent. */
    private final int v;

    /** Second exponent. */
    private final int w;

    /** Simple constructor.
     * @param v first exponent
     * @param w second exponent
     */
    public JacobiKey(final int v, final int w) {
        this.v = v;
        this.w = w;
    }

    /** Get hash code.
     * @return hash code
     */
    @Override
    public int hashCode() {
        return (v << 16) ^ w;
    }

    /** Check if the instance represent the same key as another instance.
     * @param key other key
     * @return true if the instance and the other key refer to the same polynomial
     */
    @Override
    public boolean equals(final Object key) {

        if (!(key instanceof JacobiKey)) {
            return false;
        }

        final JacobiKey otherK = (JacobiKey) key;
        return (v == otherK.v) && (w == otherK.w);
    }
}
