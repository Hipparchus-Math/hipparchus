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



import org.hipparchus.analysis.MultivariateFunction;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.ArrayRealVector;

/** A MultivariateFunction that also has a defined gradient and Hessian.
 * @since 3.1
 */
public abstract class TwiceDifferentiableFunction implements MultivariateFunction {
    /**
     * Returns the dimensionality of the function domain.
     * If dim() returns (n) then this function expects an n-vector as its input.
     * @return the expected dimension of the function's domain
     */
    public abstract int dim();

    /**
     * Returns the value of this function at (x)
     *
     * @param x a point to evaluate this function at.
     * @return the value of this function at (x)
     */
    public abstract double value(RealVector x);

    /**
     * Returns the gradient of this function at (x)
     *
     * @param x a point to evaluate this gradient at
     * @return the gradient of this function at (x)
     */
    public abstract RealVector gradient(RealVector x);

    /**
     * The Hessian of this function at (x)
     *
     * @param x a point to evaluate this Hessian at
     * @return the Hessian of this function at (x)
     */
    public abstract RealMatrix hessian(RealVector x);

    /**
     * Returns the value of this function at (x)
     *
     * @param x a point to evaluate this function at.
     * @return the value of this function at (x)
     */
    @Override
    public double value(final double[] x) {
        return value(new ArrayRealVector(x, false));
    }

    /**
     * Returns the gradient of this function at (x)
     *
     * @param x a point to evaluate this gradient at
     * @return the gradient of this function at (x)
     */
    public RealVector gradient(final double[] x) {
        return gradient(new ArrayRealVector(x, false));
    }

    /**
     * The Hessian of this function at (x)
     *
     * @param x a point to evaluate this Hessian at
     * @return the Hessian of this function at (x)
     */
    public RealMatrix hessian(final double[] x) {
        return hessian(new ArrayRealVector(x, false));
    }
}
