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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.OrdinaryDifferentialEquation;



/** This interface represents a first order differential equations set.
 *
 * <p>This interface should be implemented by all real first order
 * differential equation problems before they can be handled by the
 * integrators {@link org.hipparchus.ode.ODEIntegrator#integrate} method.</p>
 *
 * <p>A first order differential equations problem, as seen by an
 * integrator is the time derivative <code>dY/dt</code> of a state
 * vector <code>Y</code>, both being one dimensional arrays. From the
 * integrator point of view, this derivative depends only on the
 * current time <code>t</code> and on the state vector
 * <code>Y</code>.</p>
 *
 * <p>For real problems, the derivative depends also on parameters
 * that do not belong to the state vector (dynamical model constants
 * for example). These constants are completely outside of the scope
 * of this interface, the classes that implement it are allowed to
 * handle them as they want.</p>
 *
 * @see org.hipparchus.ode.ODEIntegrator
 * @see org.hipparchus.ode.FirstOrderConverter
 * @see SecondOrderDifferentialEquations
 * @deprecated as of 1.0, replaced with {@link OrdinaryDifferentialEquation}
 */
@Deprecated
public interface FirstOrderDifferentialEquations extends OrdinaryDifferentialEquation {

    /** {@inheritDoc}
     * <p>
     * The default implementation calls {@link #computeDerivatives(double, double[], double[])}.
     * </p>
     */
    @Override
    default double[] computeDerivatives(double t, double[] y) {
        final double[] yDot = new double[y.length];
        computeDerivatives(t, y, yDot);
        return yDot;
    }

    /** Get the current time derivative of the state vector.
     * @param t current value of the independent <I>time</I> variable
     * @param y array containing the current value of the state vector
     * @param yDot placeholder array where to put the time derivative of the state vector
     * @exception MathIllegalStateException if the number of functions evaluations is exceeded
     * @exception MathIllegalArgumentException if arrays dimensions do not match equations settings
     */
    void computeDerivatives(double t, double[] y, double[] yDot)
        throws MathIllegalArgumentException, MathIllegalStateException;

}
