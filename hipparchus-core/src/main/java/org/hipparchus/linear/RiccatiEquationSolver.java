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

/**
 *
 * An algebraic Riccati equation is a type of nonlinear equation that arises in
 * the context of infinite-horizon optimal control problems in continuous time
 * or discrete time.
 *
 * The continuous time algebraic Riccati equation (CARE):
 * \[
 * A^{T}X+XA-XBR^{-1}B^{T}X+Q=0
 * \}
 *
 * And the respective linear controller is:
 * \[
 * K = R^{-1}B^{T}P
 * \]
 *
 * A solver receives A, B, Q and R and computes P and K.
 *
 */
public interface RiccatiEquationSolver {

    /** Get the solution.
     * @return the p
     */
    RealMatrix getP();

    /**
     * Get the linear controller k.
     * @return the linear controller k
     */
    RealMatrix getK();

}
