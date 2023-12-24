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

import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;

public class TESTNEW {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        QuadraticFunction q = new QuadraticFunction(
    new double[][] { { 4.0, -2.0 }, { -2.0, 4.0 } },
    new double[] { 6.0, 0.0 },
    0.0);



LinearInequalityConstraint ineqc1 = new LinearInequalityConstraint(
    new double[][] { { 1.0, 0.0 },{ 0.0, 1.0 },{ 1.0, 1.0 } }, // constraint x > 1,
    new double[] { 0.0,0.0,2.0 });

LinearInequalityConstraint ineqc2 = new LinearInequalityConstraint(
    new double[][] { { 1.0, 0.0 },{ 1.0, 1.0 } ,{ 1.0, 0.0 }}, // constraint x > 1,
    new double[] { 0.0,2.0 ,0.0});


//    QPActiveSet optimizer1 = new QPActiveSet();
//        optimizer1.optimize(new ObjectiveFunction(q),ineqc1);
//
//    ineqc1 = new LinearInequalityConstraint(
//    new double[][] { { -1.0, 0.0 },{ 0.0, -1.0 },{ -1.0, -1.0 } }, // constraint x > 1,
//    new double[] { 0.0,0.0,-2.0 });

//     QPActiveSetNew optimizer2 = new QPActiveSetNew();
//        optimizer2.optimize(new ObjectiveFunction(q),ineqc1);

        ConstraintOptimizer optimizer = new ADMMQPOptimizer();
        LagrangeSolution sol = optimizer.optimize(new ObjectiveFunction(q),ineqc1);
        System.out.println("Converged:"+((ADMMQPOptimizer)optimizer).isConverged());
    }
}
