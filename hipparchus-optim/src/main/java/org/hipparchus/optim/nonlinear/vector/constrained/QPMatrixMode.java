/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.optim.OptimizationData;

/**
 * Marker OptimizationData indicating that the quadratic matrix "G" passed to the QP solver
 * is not the Hessian H, but the lower Cholesky factor L such that H = L*L^T or the inverse of Cholesky factor L^-1
 *
 * When present and set to true, QPDualActiveSolver will interpret function.getP()
 * as L or L^-1 and will NOT factorize it.
 */

public enum QPMatrixMode implements OptimizationData {

    /** QP Problem with full Hessian Matrix. */
    FULL,

    /** QP Problem with Cholesky L factor. */
    CHOLESKY,

    /** QP Problem with Inverse o Cholesky factor. */
    INVCHOLESKY

}

