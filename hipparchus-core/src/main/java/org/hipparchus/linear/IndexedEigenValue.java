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
package org.hipparchus.linear;

import org.hipparchus.complex.Complex;

/** Container for index and eigenvalue pair.
 * @since 3.0
 */
class IndexedEigenvalue {

    /** Index in the diagonal matrix. */
    private int index;

    /** Eigenvalue. */
    private final Complex eigenValue;

    /** Build the container from its fields.
     * @param index index in the diagonal matrix
     * @param eigenvalue eigenvalue
     */
    IndexedEigenvalue(final int index, final Complex eigenvalue) {
        this.index      = index;
        this.eigenValue = eigenvalue;
    }

    /** Get the index in the diagonal matrix.
     * @return index in the diagonal matrix
     */
    public int getIndex() {
        return index;
    }

    /** Set the index in the diagonal matrix.
     * @param index new index in the diagonal matrix
     */
    public void setIndex(final int index) {
        this.index = index;
    }

    /** Get the eigenvalue.
     * @return eigenvalue
     */
    public Complex getEigenvalue() {
        return eigenValue;
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object other) {

        if (this == other) {
            return true;
        }

        if (other instanceof IndexedEigenvalue) {
            final IndexedEigenvalue rhs = (IndexedEigenvalue) other;
            return eigenValue.equals(rhs.eigenValue);
        }

        return false;

    }

    /**
     * Get a hashCode for the pair.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return 4563 + index + eigenValue.hashCode();
    }

}
