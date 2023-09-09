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

package org.hipparchus.migration.ode.sampling;

import java.io.Serializable;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.sampling.ODEStateInterpolator;

/** This interface represents an interpolator over the last step
 * during an ODE integration.
 *
 * <p>The various ODE integrators provide objects implementing this
 * interface to the step handlers. These objects are often custom
 * objects tightly bound to the integrator internal algorithms. The
 * handlers can use these objects to retrieve the state vector at
 * intermediate times between the previous and the current grid points
 * (this feature is often called dense output).</p>
 *
 * @see org.hipparchus.ode.ODEIntegrator
 * @see StepHandler
 * @deprecated as of 1.0, this class is a temporary wrapper between
 * {@link ODEStateInterpolator} and {@link MigrationStepInterpolator}
 */
@Deprecated
class MigrationStepInterpolator implements StepInterpolator {

    /** Serializable UID. */
    private static final long serialVersionUID = 20160328L;

    /** Underlying interpolator. */
    private final ODEStateInterpolator interpolator;

    /** Interpolated state. */
    private ODEStateAndDerivative interpolated;

    /** Simple constructor.
     * @param interpolator underlying interpolator
     */
    MigrationStepInterpolator(final ODEStateInterpolator interpolator) {
        this.interpolator = interpolator;
        this.interpolated = interpolator.getCurrentState();
    }

    /** {@inheritDoc} */
    @Override
    @Deprecated
    public double getPreviousTime() {
        return getPreviousState().getTime();
    }

    /** {@inheritDoc} */
    @Override
    @Deprecated
    public double getCurrentTime() {
        return getCurrentState().getTime();
    }

    /** {@inheritDoc} */
    @Override
    @Deprecated
    public double getInterpolatedTime() {
        return interpolated.getTime();
    }

    /** {@inheritDoc} */
    @Override
    @Deprecated
    public void setInterpolatedTime(final double time) {
        interpolated = getInterpolatedState(time);
    }

    /** {@inheritDoc} */
    @Override
    @Deprecated
    public double[] getInterpolatedState() throws MathIllegalStateException {
        return interpolated.getPrimaryState();
    }

    /** {@inheritDoc} */
    @Override
    @Deprecated
    public double[] getInterpolatedDerivatives() throws MathIllegalStateException {
        return interpolated.getPrimaryDerivative();
    }

    /** {@inheritDoc} */
    @Override
    @Deprecated
    public double[] getInterpolatedSecondaryState(final int index) throws MathIllegalStateException {
        return interpolated.getSecondaryState(index);
    }

    /** {@inheritDoc} */
    @Override
    @Deprecated
    public double[] getInterpolatedSecondaryDerivatives(final int index) throws MathIllegalStateException {
        return interpolated.getSecondaryDerivative(index);
    }

    /** {@inheritDoc} */
    @Override
    public boolean isForward() {
        return interpolator.isForward();
    }

    /** {@inheritDoc} */
    @Override
    public MigrationStepInterpolator copy() throws MathIllegalStateException {
        return new MigrationStepInterpolator(interpolator);
    }

    /** {@inheritDoc} */
    @Override
    public ODEStateAndDerivative getPreviousState() {
        return interpolator.getPreviousState();
    }

    /** {@inheritDoc} */
    @Override
    public boolean isPreviousStateInterpolated() {
        return interpolator.isPreviousStateInterpolated();
    }

    /** {@inheritDoc} */
    @Override
    public ODEStateAndDerivative getCurrentState() {
        return interpolator.getCurrentState();
    }

    /** {@inheritDoc} */
    @Override
    public boolean isCurrentStateInterpolated() {
        return interpolator.isCurrentStateInterpolated();
    }

    /** {@inheritDoc} */
    @Override
    public ODEStateAndDerivative getInterpolatedState(final double time) {
        return interpolator.getInterpolatedState(time);
    }

    /**
     * Replace the instance with a data transfer object for serialization.
     * @return data transfer object that will be serialized
     */
    private Object writeReplace() {
        return new DataTransferObject(interpolator, interpolated.getTime());
    }

    /** Internal class used only for serialization. */
    private static class DataTransferObject implements Serializable {

        /** Serializable UID. */
        private static final long serialVersionUID = 20160328L;

        /** Underlying interpolator.
         * @serial
         */
        private final ODEStateInterpolator interpolator;

        /** Interpolation time.
         * @serial
         */
        private final double time;

        /** Simple constructor.
         * @param interpolator underlying interpolator
         * @param time interpolation time
         */
        DataTransferObject(final ODEStateInterpolator interpolator, final double time) {
            this.interpolator = interpolator;
            this.time         = time;
        }

        /** Replace the deserialized data transfer object with a {@link MigrationStepInterpolator}.
         * @return replacement {@link MigrationStepInterpolator}
         */
        private Object readResolve() {
            final MigrationStepInterpolator msi = new MigrationStepInterpolator(interpolator);
            msi.setInterpolatedTime(time);
            return msi;
        }

    }

}
