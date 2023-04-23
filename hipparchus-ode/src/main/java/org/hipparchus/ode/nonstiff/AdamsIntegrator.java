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

package org.hipparchus.ode.nonstiff;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.MultistepIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;


/** Base class for {@link AdamsBashforthIntegrator Adams-Bashforth} and
 * {@link AdamsMoultonIntegrator Adams-Moulton} integrators.
 */
public abstract class AdamsIntegrator extends MultistepIntegrator {

    /** Transformer. */
    private final AdamsNordsieckTransformer transformer;

    /**
     * Build an Adams integrator with the given order and step control parameters.
     * @param name name of the method
     * @param nSteps number of steps of the method excluding the one being computed
     * @param order order of the method
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scalAbsoluteTolerance allowed absolute error
     * @param scalRelativeTolerance allowed relative error
     * @exception MathIllegalArgumentException if order is 1 or less
     */
    public AdamsIntegrator(final String name, final int nSteps, final int order,
                           final double minStep, final double maxStep,
                           final double scalAbsoluteTolerance,
                           final double scalRelativeTolerance)
        throws MathIllegalArgumentException {
        super(name, nSteps, order, minStep, maxStep,
              scalAbsoluteTolerance, scalRelativeTolerance);
        transformer = AdamsNordsieckTransformer.getInstance(nSteps);
    }

    /**
     * Build an Adams integrator with the given order and step control parameters.
     * @param name name of the method
     * @param nSteps number of steps of the method excluding the one being computed
     * @param order order of the method
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vecAbsoluteTolerance allowed absolute error
     * @param vecRelativeTolerance allowed relative error
     * @exception IllegalArgumentException if order is 1 or less
     */
    public AdamsIntegrator(final String name, final int nSteps, final int order,
                           final double minStep, final double maxStep,
                           final double[] vecAbsoluteTolerance,
                           final double[] vecRelativeTolerance)
        throws IllegalArgumentException {
        super(name, nSteps, order, minStep, maxStep,
              vecAbsoluteTolerance, vecRelativeTolerance);
        transformer = AdamsNordsieckTransformer.getInstance(nSteps);
    }

    /** {@inheritDoc} */
    @Override
    public ODEStateAndDerivative integrate(final ExpandableODE equations,
                                           final ODEState initialState,
                                           final double finalTime)
        throws MathIllegalArgumentException, MathIllegalStateException {

        sanityChecks(initialState, finalTime);
        setStepStart(initIntegration(equations, initialState, finalTime));
        final boolean forward = finalTime > initialState.getTime();

        // compute the initial Nordsieck vector using the configured starter integrator
        start(equations, getStepStart(), finalTime);

        // reuse the step that was chosen by the starter integrator
        ODEStateAndDerivative stepEnd   =
                        AdamsStateInterpolator.taylor(equations.getMapper(), getStepStart(),
                                                      getStepStart().getTime() + getStepSize(),
                                                      getStepSize(), scaled, nordsieck);

        // main integration loop
        setIsLastStep(false);
        final double[] y  = getStepStart().getCompleteState();
        do {

            double[] predictedY  = null;
            final double[] predictedScaled = new double[y.length];
            Array2DRowRealMatrix predictedNordsieck = null;
            double error = 10;
            while (error >= 1.0) {

                // predict a first estimate of the state at step end
                predictedY = stepEnd.getCompleteState();

                // evaluate the derivative
                final double[] yDot = computeDerivatives(stepEnd.getTime(), predictedY);

                // predict Nordsieck vector at step end
                for (int j = 0; j < predictedScaled.length; ++j) {
                    predictedScaled[j] = getStepSize() * yDot[j];
                }
                predictedNordsieck = updateHighOrderDerivativesPhase1(nordsieck);
                updateHighOrderDerivativesPhase2(scaled, predictedScaled, predictedNordsieck);

                // evaluate error
                error = errorEstimation(y, stepEnd.getTime(), predictedY, predictedScaled, predictedNordsieck);
                if (Double.isNaN(error)) {
                    throw new MathIllegalStateException(LocalizedODEFormats.NAN_APPEARING_DURING_INTEGRATION,
                                                        stepEnd.getTime());
                }

                if (error >= 1.0) {
                    // reject the step and attempt to reduce error by stepsize control
                    final double factor = computeStepGrowShrinkFactor(error);
                    rescale(getStepSizeHelper().filterStep(getStepSize() * factor, forward, false));
                    stepEnd = AdamsStateInterpolator.taylor(equations.getMapper(), getStepStart(),
                                                            getStepStart().getTime() + getStepSize(),
                                                            getStepSize(),
                                                            scaled,
                                                            nordsieck);

                }
            }

            final AdamsStateInterpolator interpolator =
                            finalizeStep(getStepSize(), predictedY, predictedScaled, predictedNordsieck,
                                         forward, getStepStart(), stepEnd, equations.getMapper());

            // discrete events handling
            setStepStart(acceptStep(interpolator, finalTime));
            scaled    = interpolator.getScaled();
            nordsieck = interpolator.getNordsieck();

            if (!isLastStep()) {

                if (resetOccurred()) {

                    // some events handler has triggered changes that
                    // invalidate the derivatives, we need to restart from scratch
                    start(equations, getStepStart(), finalTime);

                    final double  nextT      = getStepStart().getTime() + getStepSize();
                    final boolean nextIsLast = forward ?
                                               (nextT >= finalTime) :
                                               (nextT <= finalTime);
                    final double hNew = nextIsLast ? finalTime - getStepStart().getTime() : getStepSize();

                    rescale(hNew);
                    System.arraycopy(getStepStart().getCompleteState(), 0, y, 0, y.length);

                } else {

                    // stepsize control for next step
                    final double  factor     = computeStepGrowShrinkFactor(error);
                    final double  scaledH    = getStepSize() * factor;
                    final double  nextT      = getStepStart().getTime() + scaledH;
                    final boolean nextIsLast = forward ?
                                               (nextT >= finalTime) :
                                               (nextT <= finalTime);
                    double hNew = getStepSizeHelper().filterStep(scaledH, forward, nextIsLast);

                    final double  filteredNextT      = getStepStart().getTime() + hNew;
                    final boolean filteredNextIsLast = forward ? (filteredNextT >= finalTime) : (filteredNextT <= finalTime);
                    if (filteredNextIsLast) {
                        hNew = finalTime - getStepStart().getTime();
                    }

                    rescale(hNew);
                    System.arraycopy(predictedY, 0, y, 0, y.length);

                }

                stepEnd = AdamsStateInterpolator.taylor(equations.getMapper(), getStepStart(), getStepStart().getTime() + getStepSize(),
                                                        getStepSize(), scaled, nordsieck);

            }

        } while (!isLastStep());

        final ODEStateAndDerivative finalState = getStepStart();
        setStepStart(null);
        setStepSize(Double.NaN);
        return finalState;

    }

    /** {@inheritDoc} */
    @Override
    protected Array2DRowRealMatrix initializeHighOrderDerivatives(final double h, final double[] t,
                                                                  final double[][] y,
                                                                  final double[][] yDot) {
        return transformer.initializeHighOrderDerivatives(h, t, y, yDot);
    }

    /** Update the high order scaled derivatives for Adams integrators (phase 1).
     * <p>The complete update of high order derivatives has a form similar to:
     * \[
     * r_{n+1} = (s_1(n) - s_1(n+1)) P^{-1} u + P^{-1} A P r_n
     * \]
     * this method computes the P<sup>-1</sup> A P r<sub>n</sub> part.</p>
     * @param highOrder high order scaled derivatives
     * (h<sup>2</sup>/2 y'', ... h<sup>k</sup>/k! y(k))
     * @return updated high order derivatives
     * @see #updateHighOrderDerivativesPhase2(double[], double[], Array2DRowRealMatrix)
     */
    public Array2DRowRealMatrix updateHighOrderDerivativesPhase1(final Array2DRowRealMatrix highOrder) {
        return transformer.updateHighOrderDerivativesPhase1(highOrder);
    }

    /** Update the high order scaled derivatives Adams integrators (phase 2).
     * <p>The complete update of high order derivatives has a form similar to:
     * \[
     * r_{n+1} = (s_1(n) - s_1(n+1)) P^{-1} u + P^{-1} A P r_n
     * \]
     * this method computes the (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u part.</p>
     * <p>Phase 1 of the update must already have been performed.</p>
     * @param start first order scaled derivatives at step start
     * @param end first order scaled derivatives at step end
     * @param highOrder high order scaled derivatives, will be modified
     * (h<sup>2</sup>/2 y'', ... h<sup>k</sup>/k! y(k))
     * @see #updateHighOrderDerivativesPhase1(Array2DRowRealMatrix)
     */
    public void updateHighOrderDerivativesPhase2(final double[] start,
                                                 final double[] end,
                                                 final Array2DRowRealMatrix highOrder) {
        transformer.updateHighOrderDerivativesPhase2(start, end, highOrder);
    }

    /** Estimate error.
     * @param previousState state vector at step start
     * @param predictedTime time at step end
     * @param predictedState predicted state vector at step end
     * @param predictedScaled predicted value of the scaled derivatives at step end
     * @param predictedNordsieck predicted value of the Nordsieck vector at step end
     * @return estimated normalized local discretization error
     * @since 2.0
     */
    protected abstract double errorEstimation(double[] previousState, double predictedTime,
                                              double[] predictedState, double[] predictedScaled,
                                              RealMatrix predictedNordsieck);

    /** Finalize the step.
     * @param stepSize step size used in the scaled and Nordsieck arrays
     * @param predictedState predicted state at end of step
     * @param predictedScaled predicted first scaled derivative
     * @param predictedNordsieck predicted Nordsieck vector
     * @param isForward integration direction indicator
     * @param globalPreviousState start of the global step
     * @param globalCurrentState end of the global step
     * @param equationsMapper mapper for ODE equations primary and secondary components
     * @return step interpolator
     * @since 2.0
     */
    protected abstract AdamsStateInterpolator finalizeStep(double stepSize, double[] predictedState,
                                                           double[] predictedScaled, Array2DRowRealMatrix predictedNordsieck,
                                                           boolean isForward,
                                                           ODEStateAndDerivative globalPreviousState,
                                                           ODEStateAndDerivative globalCurrentState,
                                                           EquationsMapper equationsMapper);

}
