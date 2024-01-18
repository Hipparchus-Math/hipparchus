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

/** Internal Solution for ADMM QP Optimizer.
 * @since 3.1
 */
public class ADMMQPSolution extends LagrangeSolution {

    /** V-tilde auxiliary variable. */
    private final RealVector v;

    /** Z auxiliary variable. */
    private final RealVector z;

    /** Simple constructor.
     * @param x solution
     * @param lambda Lagrange multipliers
     * @param value objective function value
     */
    public ADMMQPSolution(RealVector x, RealVector lambda, Double value) {
        super(x, lambda, value);
        this.v = null;
        this.z = null;
    }

    /** Simple constructor.
     * @param x solution
     * @param v V-tilde auxiliary variable
     */
    public ADMMQPSolution(RealVector x, RealVector v) {
        super(x, null, 0.0);
        this.v = v;
        this.z = null;
    }

    /** Simple constructor.
     * @param x solution
     * @param v V-tilde auxiliary variable
     * @param y Lagrange multipliers
     * @param z Z auxiliary variable
     */
    public ADMMQPSolution(final RealVector x, final RealVector v, final RealVector y, final RealVector z) {
        super(x, y, 0.0);
        this.v = v;
        this.z = z;
    }

    /** Simple constructor.
     * @param x solution
     * @param v V-tilde auxiliary variable
     * @param y Lagrange multipliers
     * @param z Z auxiliary variable
     * @param value objective function value
     */
    public ADMMQPSolution(final RealVector x, final RealVector v, final RealVector y, final RealVector z, double value) {
        super(x, y, value);
        this.v = v;
        this.z = z;
    }

    /**
     * Returns V tilde auxiliary Variable
     *
     * @return V tilde auxiliary Variable
     */
    public RealVector getV() {
        return this.v;
    }

    /**
     * Returns Z auxiliary Variable
     *
     * @return Z auxiliary Variable
     */
    public RealVector getZ() {
        return this.z;
    }

}
