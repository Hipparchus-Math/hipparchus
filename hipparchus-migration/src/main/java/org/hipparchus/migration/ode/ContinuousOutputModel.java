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

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.DenseOutputModel;

/**
 * This class stores all information provided by an ODE integrator
 * during the integration process and build a continuous model of the
 * solution from this.
 *
 * <p>This class act as a step handler from the integrator point of
 * view. It is called iteratively during the integration process and
 * stores a copy of all steps information in a sorted collection for
 * later use. Once the integration process is over, the user can use
 * the {@link #setInterpolatedTime setInterpolatedTime} and {@link
 * #getInterpolatedState getInterpolatedState} to retrieve this
 * information at any time. It is important to wait for the
 * integration to be over before attempting to call {@link
 * #setInterpolatedTime setInterpolatedTime} because some internal
 * variables are set only once the last step has been handled.</p>
 *
 * <p>This is useful for example if the main loop of the user
 * application should remain independent from the integration process
 * or if one needs to mimic the behaviour of an analytical model
 * despite a numerical model is used (i.e. one needs the ability to
 * get the model value at any time or to navigate through the
 * data).</p>
 *
 * <p>If problem modeling is done with several separate
 * integration phases for contiguous intervals, the same
 * ContinuousOutputModel can be used as step handler for all
 * integration phases as long as they are performed in order and in
 * the same direction. As an example, one can extrapolate the
 * trajectory of a satellite with one model (i.e. one set of
 * differential equations) up to the beginning of a maneuver, use
 * another more complex model including thrusters modeling and
 * accurate attitude control during the maneuver, and revert to the
 * first model after the end of the maneuver. If the same continuous
 * output model handles the steps of all integration phases, the user
 * do not need to bother when the maneuver begins or ends, he has all
 * the data available in a transparent manner.</p>
 *
 * <p>An important feature of this class is that it implements the
 * <code>Serializable</code> interface. This means that the result of
 * an integration can be serialized and reused later (if stored into a
 * persistent medium like a filesystem or a database) or elsewhere (if
 * sent to another application). Only the result of the integration is
 * stored, there is no reference to the integrated problem by
 * itself.</p>
 *
 * <p>One should be aware that the amount of data stored in a
 * ContinuousOutputModel instance can be important if the state vector
 * is large, if the integration interval is long or if the steps are
 * small (which can result from small tolerance settings in {@link
 * org.hipparchus.ode.nonstiff.AdaptiveStepsizeIntegrator adaptive
 * step size integrators}).</p>
 *
 * @see org.hipparchus.migration.ode.sampling.StepHandler
 * @see org.hipparchus.ode.sampling.ODEStateInterpolator
 * @deprecated as of 1.0, replaced with {@link DenseOutputModel}
 */
@Deprecated
public class ContinuousOutputModel extends DenseOutputModel {

    /** Serializable version identifier */
    private static final long serialVersionUID = 20160403L;

    /** Interpolation time. */
    private double interpolatedTime;

    /** Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 3.0
     */
    public ContinuousOutputModel() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** Set the time of the interpolated point.
     * <p>This method should <strong>not</strong> be called before the
     * integration is over because some internal variables are set only
     * once the last step has been handled.</p>
     * <p>Setting the time outside of the integration interval is now
     * allowed, but should be used with care since the accuracy of the
     * interpolator will probably be very poor far from this interval.
     * This allowance has been added to simplify implementation of search
     * algorithms near the interval endpoints.</p>
     * <p>Note that each time this method is called, the internal arrays
     * returned in {@link #getInterpolatedState()}, {@link
     * #getInterpolatedDerivatives()} and {@link #getInterpolatedSecondaryState(int)}
     * <em>will</em> be overwritten. So if their content must be preserved
     * across several calls, user must copy them.</p>
     * @param time time of the interpolated point
     * @see #getInterpolatedState()
     * @see #getInterpolatedDerivatives()
     * @see #getInterpolatedSecondaryState(int)
     */
    public void setInterpolatedTime(final double time) {
        this.interpolatedTime = time;
    }

    /**
     * Get the time of the interpolated point.
     * If {@link #setInterpolatedTime} has not been called, it returns
     * the final integration time.
     * @return interpolation point time
     */
    public double getInterpolatedTime() {
      return interpolatedTime;
    }

    /**
     * Get the state vector of the interpolated point.
     * <p>The returned vector is a reference to a reused array, so
     * it should not be modified and it should be copied if it needs
     * to be preserved across several calls to the associated
     * {@link #setInterpolatedTime(double)} method.</p>
     * @return state vector at time {@link #getInterpolatedTime}
     * @exception MathIllegalStateException if the number of functions evaluations is exceeded
     * @see #setInterpolatedTime(double)
     * @see #getInterpolatedDerivatives()
     * @see #getInterpolatedSecondaryState(int)
     * @see #getInterpolatedSecondaryDerivatives(int)
     */
    public double[] getInterpolatedState() throws MathIllegalStateException {
      return getInterpolatedState(getInterpolatedTime()).getPrimaryState();
    }

    /**
     * Get the derivatives of the state vector of the interpolated point.
     * <p>The returned vector is a reference to a reused array, so
     * it should not be modified and it should be copied if it needs
     * to be preserved across several calls to the associated
     * {@link #setInterpolatedTime(double)} method.</p>
     * @return derivatives of the state vector at time {@link #getInterpolatedTime}
     * @exception MathIllegalStateException if the number of functions evaluations is exceeded
     * @see #setInterpolatedTime(double)
     * @see #getInterpolatedState()
     * @see #getInterpolatedSecondaryState(int)
     * @see #getInterpolatedSecondaryDerivatives(int)
     */
    public double[] getInterpolatedDerivatives() throws MathIllegalStateException {
      return getInterpolatedState(getInterpolatedTime()).getPrimaryDerivative();
    }

    /** Get the interpolated secondary state corresponding to the secondary equations.
     * <p>The returned vector is a reference to a reused array, so
     * it should not be modified and it should be copied if it needs
     * to be preserved across several calls to the associated
     * {@link #setInterpolatedTime(double)} method.</p>
     * @param secondaryStateIndex index of the secondary set, as returned by {@link
     * org.hipparchus.ode.ExpandableODE#addSecondaryEquations(org.hipparchus.ode.SecondaryODE)
     * ExpandableODE.addSecondaryEquations(secondary)}
     * @return interpolated secondary state at the current interpolation date
     * @see #setInterpolatedTime(double)
     * @see #getInterpolatedState()
     * @see #getInterpolatedDerivatives()
     * @see #getInterpolatedSecondaryDerivatives(int)
     * @exception MathIllegalStateException if the number of functions evaluations is exceeded
     */
    public double[] getInterpolatedSecondaryState(final int secondaryStateIndex)
      throws MathIllegalStateException {
      return getInterpolatedState(getInterpolatedTime()).getSecondaryState(secondaryStateIndex);
    }

    /** Get the interpolated secondary derivatives corresponding to the secondary equations.
     * <p>The returned vector is a reference to a reused array, so
     * it should not be modified and it should be copied if it needs
     * to be preserved across several calls to the associated
     * {@link #setInterpolatedTime(double)} method.</p>
     * @param secondaryStateIndex index of the secondary set, as returned by {@link
     * org.hipparchus.ode.ExpandableODE#addSecondaryEquations(org.hipparchus.ode.SecondaryODE)
     * ExpandableODE.addSecondaryEquations(secondary)}
     * @return interpolated secondary derivatives at the current interpolation date
     * @see #setInterpolatedTime(double)
     * @see #getInterpolatedState()
     * @see #getInterpolatedDerivatives()
     * @see #getInterpolatedSecondaryState(int)
     * @exception MathIllegalStateException if the number of functions evaluations is exceeded
     */
    public double[] getInterpolatedSecondaryDerivatives(final int secondaryStateIndex)
      throws MathIllegalStateException {
      return getInterpolatedState(getInterpolatedTime()).getSecondaryDerivative(secondaryStateIndex);
    }

}
