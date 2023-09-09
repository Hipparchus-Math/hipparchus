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
 * <p>This is the same as problem 1 except integration is done
 * backward in time</p>
 */
public class TestProblem5 extends TestProblemAbstract {

    /**
     * Simple constructor.
     */
    public TestProblem5() {
        super(0.0, new double[] { 1.0, 0.1 }, -4.0, new double[] { 1.0, 1.0 });
    }

    @Override
    public double[] doComputeDerivatives(double t, double[] y) {

        // compute the derivatives
        final double[] yDot = new double[getDimension()];
        for (int i = 0; i < getDimension(); ++i) {
            yDot[i] = -y[i];
        }
        return yDot;

    }

    @Override
    public double[] computeTheoreticalState(double t) {
        final double c = FastMath.exp (getInitialTime() - t);
        final double[] y = getInitialState().getPrimaryState();
        for (int i = 0; i < getDimension(); ++i) {
            y[i] *= c;
        }
        return y;
    }

}
