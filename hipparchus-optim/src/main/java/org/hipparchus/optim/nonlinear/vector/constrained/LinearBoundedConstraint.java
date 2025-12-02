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

/** A set of linear inequality constraints expressed as ub&gt;Ax&gt;lb.
 * @since 3.1
 */
public class LinearBoundedConstraint extends BoundedConstraint implements OptimizationData {

    /** The corresponding set of individual linear constraint functions. */
    private final RealMatrix a;

    /** Construct a set of linear inequality constraints from Ax &lt; B
     * @param a A matrix linear coefficient vectors
     * @param lower lower bound
     * @param upper upper bound
     */
    public LinearBoundedConstraint(final RealMatrix a, final RealVector lower, final RealVector upper) {
        super(lower, upper);
        this.a = a;
    }

    /** Construct a set of linear inequality constraints from Ax &lt; B
     * @param a A matrix linear coefficient vectors
     * @param lower lower bound
     * @param upper upper bound
     */
    public LinearBoundedConstraint(final double[][] a, final double[] lower, final double[] upper) {
        this(new Array2DRowRealMatrix(a), new ArrayRealVector(lower), new ArrayRealVector(upper));
    }

    /** {@inheritDoc} */
    @Override
    public double[] value(final double[] x) {
        return this.a.operate(x);
    }

    /** {@inheritDoc} */
    @Override
    public int dim() {
        return a.getColumnDimension();
    }

    /** {@inheritDoc} */
    @Override
    public RealVector value(final RealVector x) {
        return a.operate(x);
    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix jacobian(final RealVector x) {
        return a.copy();
    }

}
