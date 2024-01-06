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

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.util.MathUtils;

/** Set of linear inequality constraints expressed as \( A x \gt B\).
 * @since 3.1
 */
public class LinearInequalityConstraint extends InequalityConstraint implements OptimizationData {

    /** Corresponding set of individual linear constraint functions. */
    private final RealMatrix a;

    /**
     * Construct a set of linear inequality constraints from \( A x \gt B\).
     * @param a A matrix linear coefficient vectors
     * @param b A vector of constants
     */
    public LinearInequalityConstraint(final RealMatrix a, final RealVector b) {
        super(b);
        MathUtils.checkDimension(b.getDimension(), a.getRowDimension());
        this.a = a;
    }

    /**
     * Construct a set of linear inequality constraints from Ax &gt; B
     * @param a A matrix linear coefficient vectors
     * @param b A vector of constants
     */
    public LinearInequalityConstraint(final double[][] a, final double[] b) {
        this(MatrixUtils.createRealMatrix(a), MatrixUtils.createRealVector(b));
    }

     @Override
    public int dim() {
        return a.getColumnDimension();
    }

    @Override
    public RealVector value(RealVector x) {
        return a.operate(x);
    }

    @Override
    public RealMatrix jacobian(RealVector x) {
        return a.copy();
    }

}
