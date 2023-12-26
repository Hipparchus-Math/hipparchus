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

import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.OptimizationData;
/** Karush–Kuhn–Tucker Solver.
 *<p>
 * Solve Equation:
 * </p>
 * \[\begin{align}
 *  |H A^{T}| &amp; = B_1\\
 *  |A  R|    &amp; = B_2
 * \end{align}\]
 * @param <T> type of the solution
 * @since 3.1
 */
 public interface KarushKuhnTuckerSolver<T> extends OptimizationData {

     /** Solve Karush–Kuhn–Tucker equation from given right hand value.
      * @param b1 first right hand vector
      * @param b2 second right hand vector
      * @return Tuple with the solution x,Lambda,value
      */
     T solve(RealVector b1, RealVector b2);

     /** Iterate Karush–Kuhn–Tucker equation from given list of Vector
      * @param b list of vectors
      * @return Tuple with the solution x,Lambda,value
      */
     T iterate(RealVector... b);

}
