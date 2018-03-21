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

import org.hipparchus.exception.MathIllegalArgumentException;

/** Interface for all algorithms providing matrix decomposition.
 * @since 1.3
 */
public interface MatrixDecomposer {

    /**
     * Get a solver for finding the A &times; X = B solution in least square sense.
     * @param a coefficient matrix A to decompose
     * @return a solver
     * @throws MathIllegalArgumentException if decomposition fails
     */
    DecompositionSolver decompose(RealMatrix a) throws MathIllegalArgumentException;

}
