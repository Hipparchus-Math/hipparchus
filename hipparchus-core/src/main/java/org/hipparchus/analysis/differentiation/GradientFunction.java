/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.analysis.differentiation;

import org.hipparchus.analysis.MultivariateVectorFunction;

/** Class representing the gradient of a multivariate function.
 * <p>
 * The vectorial components of the function represent the derivatives
 * with respect to each function parameters.
 * </p>
 */
public class GradientFunction implements MultivariateVectorFunction {

    /** Underlying real-valued function. */
    private final MultivariateDifferentiableFunction f;

    /** Simple constructor.
     * @param f underlying real-valued function
     */
    public GradientFunction(final MultivariateDifferentiableFunction f) {
        this.f = f;
    }

    /** {@inheritDoc} */
    @Override
    public double[] value(double[] point) {

        // set up parameters
        final DSFactory factory = new DSFactory(point.length, 1);
        final DerivativeStructure[] dsX = new DerivativeStructure[point.length];
        for (int i = 0; i < point.length; ++i) {
            dsX[i] = factory.variable(i, point[i]);
        }

        // compute the derivatives
        final DerivativeStructure dsY = f.value(dsX);

        // extract the gradient
        final double[] y = new double[point.length];
        final int[] orders = new int[point.length];
        for (int i = 0; i < point.length; ++i) {
            orders[i] = 1;
            y[i] = dsY.getPartialDerivative(orders);
            orders[i] = 0;
        }

        return y;

    }

}
