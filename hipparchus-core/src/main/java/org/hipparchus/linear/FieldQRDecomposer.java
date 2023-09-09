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

import org.hipparchus.CalculusFieldElement;

/** Matrix decomposer using QR-decomposition.
 * @param <T> the type of the field elements
 * @since 2.2
 */
public class FieldQRDecomposer<T extends CalculusFieldElement<T>> implements FieldMatrixDecomposer<T> {

    /** Threshold under which a matrix is considered singular. */
    private final T singularityThreshold;

    /**
     * Creates a QR decomposer with specify threshold for several matrices.
     * @param singularityThreshold threshold (based on partial row norm)
     * under which a matrix is considered singular
     */
    public FieldQRDecomposer(final T singularityThreshold) {
        this.singularityThreshold = singularityThreshold;
    }

    /** {@inheritDoc} */
    @Override
    public FieldDecompositionSolver<T> decompose(final FieldMatrix<T> a) {
        return new FieldQRDecomposition<>(a, singularityThreshold).getSolver();
    }

}
