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

package org.hipparchus.ode;

import org.hipparchus.util.FastMath;

/**
 * This class is used in the junit tests for the ODE integrators.

 * <p>This specific problem is the following differential equation :
 * <pre>
 *    y' = t^3 - t y
 * </pre>
 * with the initial condition y (0) = 0. The solution of this equation
 * is the following function :
 * <pre>
 *   y (t) = t^2 + 2 (exp (- t^2 / 2) - 1)
 * </pre>
 * </p>

 */
public class TestProblem2 extends TestProblemAbstract {

    /**
     * Simple constructor.
     */
    public TestProblem2() {
        super(0.0, new double[] { 0.0 }, 1.0, new double[] { 1.0 });
    }

    @Override
    public double[] doComputeDerivatives(double t, double[] y) {

        // compute the derivatives
        final  double[] yDot = new double[getDimension()];
        for (int i = 0; i < getDimension(); ++i) {
            yDot[i] = t * (t * t - y[i]);
        }
        return yDot;

    }

    @Override
    public double[] computeTheoreticalState(double t) {
        final double[] y = new double[getDimension()];
        double t2 = t * t;
        double c = t2 + 2 * (FastMath.exp (-0.5 * t2) - 1);
        for (int i = 0; i < getDimension(); ++i) {
            y[i] = c;
        }
        return y;
    }

}
