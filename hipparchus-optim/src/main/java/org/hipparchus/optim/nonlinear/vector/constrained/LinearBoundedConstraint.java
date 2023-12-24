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


import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.OptimizationData;

/** A set of linear inequality constraints expressed as ub>Ax>lb.
 * @since 3.1
 */
public class LinearBoundedConstraint extends BoundedConstraint implements OptimizationData {
    /** The corresponding set of individual linear constraint functions */
//    public final LinearFunction[] lcf;
    private final RealMatrix A;
//    public final RealVector UB;
//    public final RealVector LB;
    /**
     * Construct a set of linear inequality constraints from Ax &lt; B
     * @param A A matrix linear coefficient vectors
     * @param b A vector of constants
     */
    public LinearBoundedConstraint(RealMatrix A, final RealVector LB,final RealVector UB) {

        this.A = A.copy();
        this.LB = LB;
         this.UB = UB;

    }

    /**
     * Construct a set of linear inequality constraints from Ub>=Ax>=Lb.
     * @param A A matrix linear coefficient vectors
     * @param LB A vector of constants
     * @param UB A vector of constants
     */
    public LinearBoundedConstraint(final double[][] A,final double[] LB,final double[] UB) {
        this(new Array2DRowRealMatrix(A), new ArrayRealVector(LB),new ArrayRealVector(UB));
    }



    @Override
    public double[] value(double[] x) throws IllegalArgumentException {
        return this.A.operate(x);
    }

    @Override
    public int dim() {
        return A.getColumnDimension();
    }

    @Override
    public RealVector value(RealVector x) {
       return A.operate(x);
    }

    @Override
    public RealMatrix jacobian(RealVector x) {
       return A.copy();
    }

    @Override
    public int dimY() {
       return A.getRowDimension();
    }
}
