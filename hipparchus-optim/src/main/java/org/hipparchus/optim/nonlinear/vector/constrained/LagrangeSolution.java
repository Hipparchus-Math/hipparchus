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

/** Container for Lagrange t-uple.
 * @since 3.1
 */
public class LagrangeSolution {

    /** Solution vector. */
    private final RealVector x;

    /** Lagrange multipliers. */
    private final RealVector lambda;

    /** Objective function value. */
    private final double value;

    /** Simple constructor.
     * @param x solution
     * @param lambda Lagrange multipliers
     * @param value objective function value
     */
    public LagrangeSolution(final RealVector x, final RealVector lambda, final double value){
        this.x      = x;
        this.lambda = lambda;
        this.value  = value;
    }

    /**
     * Returns X solution
     *
     * @return X solution
     */
    public RealVector getX() {
        return x;
    }

    /**
     * Returns Lambda Multiplier
     *
     * @return X Lambda Multiplier
     */
    public RealVector getLambda() {
        return lambda;
    }

    /**
     * Returns min(max) evaluated function at x
     *
     * @return min(max) evaluated function at x
     */
    public double getValue() {
        return value;
    }

}
