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

package org.hipparchus.ode.nonstiff;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.AbstractFieldIntegrator;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * This abstract class holds the common part of all adaptive
 * stepsize integrators for Ordinary Differential Equations.
 *
 * <p>These algorithms perform integration with stepsize control, which
 * means the user does not specify the integration step but rather a
 * tolerance on error. The error threshold is computed as
 * </p>
 * <pre>
 * threshold_i = absTol_i + relTol_i * max (abs (ym), abs (ym+1))
 * </pre>
 * <p>
 * where absTol_i is the absolute tolerance for component i of the
 * state vector and relTol_i is the relative tolerance for the same
 * component. The user can also use only two scalar values absTol and
 * relTol which will be used for all components.
 * </p>
 * <p>
 * Note that <em>only</em> the {@link FieldODEState#getPrimaryState() main part}
 * of the state vector is used for stepsize control. The {@link
 * FieldODEState#getSecondaryState(int) secondary parts} of the state
 * vector are explicitly ignored for stepsize control.
 * </p>
 *
 * <p>If the estimated error for ym+1 is such that</p>
 * <pre>
 * sqrt((sum (errEst_i / threshold_i)^2 ) / n) &lt; 1
 * </pre>
 *
 * <p>(where n is the main set dimension) then the step is accepted,
 * otherwise the step is rejected and a new attempt is made with a new
 * stepsize.</p>
 *
 * @param <T> the type of the field elements
 *
 */

public abstract class AdaptiveStepsizeFieldIntegrator<T extends CalculusFieldElement<T>>
    extends AbstractFieldIntegrator<T> {

    /** Helper for step size control. */
    private StepsizeHelper stepsizeHelper;

    /** Build an integrator with the given stepsize bounds.
     * The default step handler does nothing.
     * @param field field to which the time and state vector elements belong
     * @param name name of the method
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scalAbsoluteTolerance allowed absolute error
     * @param scalRelativeTolerance allowed relative error
     */
    public AdaptiveStepsizeFieldIntegrator(final Field<T> field, final String name,
                                           final double minStep, final double maxStep,
                                           final double scalAbsoluteTolerance,
                                           final double scalRelativeTolerance) {
        super(field, name);
        stepsizeHelper = new StepsizeHelper(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
        resetInternalState();
    }

    /** Build an integrator with the given stepsize bounds.
     * The default step handler does nothing.
     * @param field field to which the time and state vector elements belong
     * @param name name of the method
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vecAbsoluteTolerance allowed absolute error
     * @param vecRelativeTolerance allowed relative error
     */
    public AdaptiveStepsizeFieldIntegrator(final Field<T> field, final String name,
                                           final double minStep, final double maxStep,
                                           final double[] vecAbsoluteTolerance,
                                           final double[] vecRelativeTolerance) {
        super(field, name);
        stepsizeHelper = new StepsizeHelper(minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
        resetInternalState();
    }

    /** Set the adaptive step size control parameters.
     * <p>
     * A side effect of this method is to also reset the initial
     * step so it will be automatically computed by the integrator
     * if {@link #setInitialStepSize(double) setInitialStepSize}
     * is not called by the user.
     * </p>
     * @param minimalStep minimal step (must be positive even for backward
     * integration), the last step can be smaller than this
     * @param maximalStep maximal step (must be positive even for backward
     * integration)
     * @param absoluteTolerance allowed absolute error
     * @param relativeTolerance allowed relative error
     */
    public void setStepSizeControl(final double minimalStep, final double maximalStep,
                                   final double absoluteTolerance, final double relativeTolerance) {
        stepsizeHelper = new StepsizeHelper(minimalStep, maximalStep, absoluteTolerance, relativeTolerance);
    }

    /** Set the adaptive step size control parameters.
     * <p>
     * A side effect of this method is to also reset the initial
     * step so it will be automatically computed by the integrator
     * if {@link #setInitialStepSize(double) setInitialStepSize}
     * is not called by the user.
     * </p>
     * @param minimalStep minimal step (must be positive even for backward
     * integration), the last step can be smaller than this
     * @param maximalStep maximal step (must be positive even for backward
     * integration)
     * @param absoluteTolerance allowed absolute error
     * @param relativeTolerance allowed relative error
     */
    public void setStepSizeControl(final double minimalStep, final double maximalStep,
                                   final double[] absoluteTolerance,
                                   final double[] relativeTolerance) {
        stepsizeHelper = new StepsizeHelper(minimalStep, maximalStep, absoluteTolerance, relativeTolerance);
    }

    /** Get the stepsize helper.
     * @return stepsize helper
     * @since 2.0
     */
    protected StepsizeHelper getStepSizeHelper() {
        return stepsizeHelper;
    }

    /** Set the initial step size.
     * <p>This method allows the user to specify an initial positive
     * step size instead of letting the integrator guess it by
     * itself. If this method is not called before integration is
     * started, the initial step size will be estimated by the
     * integrator.</p>
     * @param initialStepSize initial step size to use (must be positive even
     * for backward integration ; providing a negative value or a value
     * outside of the min/max step interval will lead the integrator to
     * ignore the value and compute the initial step size by itself)
     */
    public void setInitialStepSize(final double initialStepSize) {
        stepsizeHelper.setInitialStepSize(initialStepSize);
    }

    /** {@inheritDoc} */
    @Override
    protected void sanityChecks(final FieldODEState<T> initialState, final T t)
        throws MathIllegalArgumentException {
        super.sanityChecks(initialState, t);
        stepsizeHelper.setMainSetDimension(initialState.getPrimaryStateDimension());
    }

    /** Initialize the integration step.
     * @param forward forward integration indicator
     * @param order order of the method
     * @param scale scaling vector for the state vector (can be shorter than state vector)
     * @param state0 state at integration start time
     * @param mapper mapper for all the equations
     * @return first integration step
     * @exception MathIllegalStateException if the number of functions evaluations is exceeded
     * @exception MathIllegalArgumentException if arrays dimensions do not match equations settings
     */
    public double initializeStep(final boolean forward, final int order, final T[] scale,
                                 final FieldODEStateAndDerivative<T> state0,
                                 final FieldEquationsMapper<T> mapper)
        throws MathIllegalArgumentException, MathIllegalStateException {

        if (stepsizeHelper.getInitialStep() > 0) {
            // use the user provided value
            return forward ? stepsizeHelper.getInitialStep() : -stepsizeHelper.getInitialStep();
        }

        // very rough first guess : h = 0.01 * ||y/scale|| / ||y'/scale||
        // this guess will be used to perform an Euler step
        final T[] y0    = state0.getCompleteState();
        final T[] yDot0 = state0.getCompleteDerivative();
        double yOnScale2     = 0;
        double yDotOnScale2  = 0;
        for (int j = 0; j < scale.length; ++j) {
            final double ratio    = y0[j].getReal() / scale[j].getReal();
            yOnScale2            += ratio * ratio;
            final double ratioDot = yDot0[j].getReal() / scale[j].getReal();
            yDotOnScale2         += ratioDot * ratioDot;
        }

        double h = ((yOnScale2 < 1.0e-10) || (yDotOnScale2 < 1.0e-10)) ?
                   1.0e-6 : (0.01 * FastMath.sqrt(yOnScale2 / yDotOnScale2));
        if (h > getMaxStep()) {
            h = getMaxStep();
        }
        if (! forward) {
            h = -h;
        }

        // perform an Euler step using the preceding rough guess
        final T[] y1 = MathArrays.buildArray(getField(), y0.length);
        for (int j = 0; j < y0.length; ++j) {
            y1[j] = y0[j].add(yDot0[j].multiply(h));
        }
        final T[] yDot1 = computeDerivatives(state0.getTime().add(h), y1);

        // estimate the second derivative of the solution
        double yDDotOnScale = 0;
        for (int j = 0; j < scale.length; ++j) {
            final double ratioDotDot = (yDot1[j].getReal() - yDot0[j].getReal()) / scale[j].getReal();
            yDDotOnScale += ratioDotDot * ratioDotDot;
        }
        yDDotOnScale = FastMath.sqrt(yDDotOnScale) / h;

        // step size is computed such that
        // h^order * max (||y'/tol||, ||y''/tol||) = 0.01
        final double maxInv2 = FastMath.max(FastMath.sqrt(yDotOnScale2), yDDotOnScale);
        final double h1 = (maxInv2 < 1.0e-15) ?
                          FastMath.max(1.0e-6, 0.001 * FastMath.abs(h)) :
                          FastMath.pow(0.01 / maxInv2, 1.0 / order);
        h = FastMath.min(100.0 * FastMath.abs(h), h1);
        h = FastMath.max(h, 1.0e-12 * FastMath.abs(state0.getTime().getReal()));  // avoids cancellation when computing t1 - t0
        if (h < getMinStep()) {
            h = getMinStep();
        }
        if (h > getMaxStep()) {
            h = getMaxStep();
        }

        if (! forward) {
            h = -h;
        }

        return h;

    }

    /** Reset internal state to dummy values. */
    protected void resetInternalState() {
        setStepStart(null);
        setStepSize(getField().getZero().add(stepsizeHelper.getDummyStepsize()));
    }

    /** Get the minimal step.
     * @return minimal step
     */
    public double getMinStep() {
        return stepsizeHelper.getMinStep();
    }

    /** Get the maximal step.
     * @return maximal step
     */
    public double getMaxStep() {
        return stepsizeHelper.getMaxStep();
    }

}
