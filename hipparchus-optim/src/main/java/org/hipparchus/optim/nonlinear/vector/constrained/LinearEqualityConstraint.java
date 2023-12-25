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
/** A set of linear equality constraints given as Ax = b.
 * @since 3.1
 */
public class LinearEqualityConstraint extends EqualityConstraint implements OptimizationData {

    public final RealMatrix A;

    /**
     * Construct a set of linear equality constraints Ax = b. Represents
     * equations A[i].x = b[i], for each row of A.
     *
     * @param A the matrix of linear weights
     * @param b the vector of constants
     */
    public LinearEqualityConstraint(final RealMatrix A, final RealVector b) {
        super(b);
        this.A = A;
    }

    /**
     * Construct a set of linear equality constraints Ax = b. Represents
     * equations A[i].x = b[i], for each row of A.
     *
     * @param A the matrix of linear weights
     * @param b the vector of constants
     */
    public LinearEqualityConstraint(final double[][] A, final double[] b) {
        this(new Array2DRowRealMatrix(A), new ArrayRealVector(b));
    }

    @Override
    public int dim() {
        return A.getColumnDimension();
    }

    @Override
    public RealVector value(final RealVector x) {
        return A.operate(x);
    }

    @Override
    public RealMatrix jacobian(final RealVector x) {
        return A.copy();
    }

}
