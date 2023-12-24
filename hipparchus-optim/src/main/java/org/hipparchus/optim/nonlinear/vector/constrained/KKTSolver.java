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

import org.hipparchus.optim.OptimizationData;
/**
 * Karush–Kuhn–Tucker Solver.
 *<br/>
 * Solve Equation:
 *<br/>
* |H AT|=B1
* <br/>
* |A  R|=B2
 * @since 3.1
*/
 public abstract class KKTSolver<P, M, V> implements OptimizationData {

    /**
    *Solve KKT equation from given right hand value.
    *@param B1 first right hand vector
    *@param B2 second right hand vector
    *@return Tuple with the solution x,Lambda,value
    *
    */
    public abstract P solve(V B1, V B2);

     /**
    *Return KKT Matrix.
    *@param H Hessian MAtrix
    *@param A Constraint Matrix
    *@param R Wheight Matrix
    *@return M KKT MAtrix
    */
    public abstract M getKKTMatrix(M H, M A, M R);

    /**
    *Iterate KKT equation from given list of Vector
    *@param V list of Vector
    *@return P Tuple with the solution x,Lambda,value
    */
    public abstract P iterate(V... b);

}
