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

import java.util.function.Predicate;

import org.hipparchus.FieldElement;

/** Matrix decomposer using LU-decomposition.
 * @param <T> the type of the field elements
 * @since 2.2
 */
public class FieldLUDecomposer<T extends FieldElement<T>> implements FieldMatrixDecomposer<T> {

    /** Checker for zero elements. */
    private final Predicate<T> zeroChecker;

    /**
     * Creates a LU decomposer with specific zero checker for several matrices.
     * @param zeroChecker checker for zero elements
     */
    public FieldLUDecomposer(final Predicate<T> zeroChecker) {
        this.zeroChecker = zeroChecker;
    }

    /** {@inheritDoc} */
    @Override
    public FieldDecompositionSolver<T> decompose(final FieldMatrix<T> a) {
        return new FieldLUDecomposition<>(a, zeroChecker).getSolver();
    }

}
