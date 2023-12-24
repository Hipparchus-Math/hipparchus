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

import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;

public class TESTHOCKShittkowski77 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        HockSchittkowskiFunction77 problem = new HockSchittkowskiFunction77();
       // HockSchittkowskiCostraintInequality  ineq = new HockSchittkowskiCostraintInequality (new Array2DRowRealMatrix(5,4),new ArrayRealVector(new double[]{-2,-1.5,-1.5,-1.5,-1.5}));
         HockSchittkowskiConstraintEquality77  eq = new HockSchittkowskiConstraintEquality77();
        SQPOptimizerS optimizer3 = new SQPOptimizerS();
        SQPOption option=new SQPOption();
        LagrangeSolution sol1 = optimizer3.optimize(new InitialGuess(new double[]{2,2,2,3,3}),new ObjectiveFunction(problem),eq);
    }


}
