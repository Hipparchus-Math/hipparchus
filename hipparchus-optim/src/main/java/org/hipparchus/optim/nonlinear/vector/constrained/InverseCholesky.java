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
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.OptimizationData;

/** Inverse of Cholesky factorization.
 * @since 4.1
 */
public class InverseCholesky implements OptimizationData {

    /** Inverse of Cholesku factorization. */
    private final RealMatrix invL;

    /** Simple constructor.
     * @param inverseL inverse of Cholesku factorization
     */
    public InverseCholesky(final RealMatrix inverseL) {
       this.invL = inverseL;
    }

    /** Get inverse of Cholesku factorization.
     * @return inverse of Cholesku factorization
     */
    public RealMatrix getInverseL() {
        return this.invL;
    }

}
