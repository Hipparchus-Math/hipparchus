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



import org.hipparchus.analysis.MultivariateVectorFunction;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/** A MultivariateFunction that also has a defined gradient and Hessian.
 * @since 3.1
 */
public interface VectorDifferentiableFunction extends MultivariateVectorFunction {

    /**
     * Returns the dimensionality of the function domain.
     * If dim() returns (n) then this function expects an n-vector as its input.
     * @return the expected dimension of the function's domain
     */
    int dim();

     /**
     * Returns the dimensionality of the function eval.
     *
     * @return the expected dimension of the function's eval
     */
    int dimY();

    /**
     * Returns the value of this function at (x)
     *
     * @param x a point to evaluate this function at.
     * @return the value of this function at (x)
     */
    RealVector value(RealVector x);

    /**
     * Returns the value of this function at (x)
     *
     * @param x a point to evaluate this function at.
     * @return the value of this function at (x)
     */
    default double[] value(final double[] x) {
        return value(new ArrayRealVector(x, false)).toArray();
    }

    /**
     * Returns the gradient of this function at (x)
     *
     * @param x a point to evaluate this gradient at
     * @return the gradient of this function at (x)
     */
    RealMatrix jacobian(RealVector x);

    /**
     * Returns the gradient of this function at (x)
     *
     * @param x a point to evaluate this gradient at
     * @return the gradient of this function at (x)
     */
    default RealMatrix gradient(final double[] x) {
        return jacobian(new ArrayRealVector(x, false));
    }

}
