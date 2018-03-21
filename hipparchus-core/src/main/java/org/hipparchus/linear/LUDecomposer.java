/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hipparchus.linear;

/** Matrix decomposer using LU-decomposition.
 * @since 1.3
 */
public class LUDecomposer implements MatrixDecomposer {

    /** Threshold under which a matrix is considered singular. */
    private final double singularityThreshold;

    /**
     * Creates a LU decomposer with specify threshold for several matrices.
     * @param singularityThreshold threshold (based on partial row norm)
     * under which a matrix is considered singular
     */
    public LUDecomposer(final double singularityThreshold) {
        this.singularityThreshold = singularityThreshold;
    }

    /** {@inheritDoc} */
    @Override
    public DecompositionSolver decompose(final RealMatrix a) {
        return new LUDecomposition(a, singularityThreshold).getSolver();
    }

}
