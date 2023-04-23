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

package org.hipparchus.migration.ode;

import org.hipparchus.ode.SecondOrderODE;

/** This interface represents a second order differential equations set.

 * <p>This interface should be implemented by all real second order
 * differential equation problems before they can be handled by the
 * {@link org.hipparchus.ode.FirstOrderConverter converter to first order}.</p>
 *
 * <p>A second order differential equations problem, as seen by an
 * integrator is the second time derivative <code>d2Y/dt^2</code> of a
 * state vector <code>Y</code>, both being one dimensional
 * arrays. From the integrator point of view, this derivative depends
 * only on the current time <code>t</code>, on the state vector
 * <code>Y</code> and on the first time derivative of the state
 * vector.</p>
 *
 * <p>For real problems, the derivative depends also on parameters
 * that do not belong to the state vector (dynamical model constants
 * for example). These constants are completely outside of the scope
 * of this interface, the classes that implement it are allowed to
 * handle them as they want.</p>
 * @deprecated as of 1.0, replaced with {@link SecondOrderODE}
 */
@Deprecated
public interface SecondOrderDifferentialEquations extends SecondOrderODE {

    /** {@inheritDoc} */
    @Override
    default double[] computeSecondDerivatives(double t, double[] y, double[] yDot) {
        final double[] yDDot = new double[y.length];
        computeSecondDerivatives(t, y, yDot, yDDot);
        return yDDot;
    }

    /** Get the current time derivative of the state vector.
     * @param t current value of the independent <I>time</I> variable
     * @param y array containing the current value of the state vector
     * @param yDot array containing the current value of the first derivative
     * of the state vector
     * @param yDDot placeholder array where to put the second time derivative
     * of the state vector
     */
    void computeSecondDerivatives(double t, double[] y, double[] yDot, double[] yDDot);

}
