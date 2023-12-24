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

/** Set of linear inequality constraints expressed as Ax>b.
 * @since 3.1
 */
public class LinearInequalityConstraint extends InequalityConstraint implements OptimizationData {
    /** The corresponding set of individual linear constraint functions */
//    public final LinearFunction[] lcf;
    public final RealMatrix A;

    /**
     * Construct a set of linear inequality constraints from Ax &lt; B
     * @param A A matrix linear coefficient vectors
     * @param b A vector of constants
     */
    public LinearInequalityConstraint(final RealMatrix A, final RealVector b) {

        this.A = A;
        this.LB = b;
        this.UB = new ArrayRealVector(b.getDimension(),Double.POSITIVE_INFINITY);

    }

    /**
     * Construct a set of linear inequality constraints from Ax &lt; B
     * @param A A matrix linear coefficient vectors
     * @param b A vector of constants
     */
    public LinearInequalityConstraint(final double[][] A, final double[] b) {
        this(new Array2DRowRealMatrix(A), new ArrayRealVector(b));
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
       return this.A.getRowDimension();
    }
}
