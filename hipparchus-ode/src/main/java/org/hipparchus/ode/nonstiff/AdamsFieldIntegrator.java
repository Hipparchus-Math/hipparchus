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

import org.hipparchus.Field;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.linear.Array2DRowFieldMatrix;
import org.hipparchus.linear.FieldMatrix;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.MultistepFieldIntegrator;
import org.hipparchus.util.MathArrays;


/** Base class for {@link AdamsBashforthFieldIntegrator Adams-Bashforth} and
 * {@link AdamsMoultonFieldIntegrator Adams-Moulton} integrators.
 * @param <T> the type of the field elements
 */
public abstract class AdamsFieldIntegrator<T extends CalculusFieldElement<T>> extends MultistepFieldIntegrator<T> {

    /** Transformer. */
    private final AdamsNordsieckFieldTransformer<T> transformer;

    /**
     * Build an Adams integrator with the given order and step control parameters.
     * @param field field to which the time and state vector elements belong
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
    public AdamsFieldIntegrator(final Field<T> field, final String name,
                                final int nSteps, final int order,
                                final double minStep, final double maxStep,
                                final double scalAbsoluteTolerance,
                                final double scalRelativeTolerance)
        throws MathIllegalArgumentException {
        super(field, name, nSteps, order, minStep, maxStep,
              scalAbsoluteTolerance, scalRelativeTolerance);
        transformer = AdamsNordsieckFieldTransformer.getInstance(field, nSteps);
    }

    /**
     * Build an Adams integrator with the given order and step control parameters.
     * @param field field to which the time and state vector elements belong
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
    public AdamsFieldIntegrator(final Field<T> field, final String name,
                                final int nSteps, final int order,
                                final double minStep, final double maxStep,
                                final double[] vecAbsoluteTolerance,
                                final double[] vecRelativeTolerance)
        throws IllegalArgumentException {
        super(field, name, nSteps, order, minStep, maxStep,
              vecAbsoluteTolerance, vecRelativeTolerance);
        transformer = AdamsNordsieckFieldTransformer.getInstance(field, nSteps);
    }

    /** {@inheritDoc} */
    @Override
    public FieldODEStateAndDerivative<T> integrate(final FieldExpandableODE<T> equations,
                                                   final FieldODEState<T> initialState,
                                                   final T finalTime)
        throws MathIllegalArgumentException, MathIllegalStateException {

        sanityChecks(initialState, finalTime);
        setStepStart(initIntegration(equations, initialState, finalTime));
        final boolean forward = finalTime.subtract(initialState.getTime()).getReal() > 0;

        // compute the initial Nordsieck vector using the configured starter integrator
        start(equations, getStepStart(), finalTime);

        // reuse the step that was chosen by the starter integrator
        FieldODEStateAndDerivative<T> stepStart = getStepStart();
        FieldODEStateAndDerivative<T> stepEnd   =
                        AdamsFieldStateInterpolator.taylor(equations.getMapper(), stepStart,
                                                           stepStart.getTime().add(getStepSize()),
                                                           getStepSize(), scaled, nordsieck);

        // main integration loop
        setIsLastStep(false);
        final T[] y = stepStart.getCompleteState();
        do {

            T[] predictedY = null;
            final T[] predictedScaled = MathArrays.buildArray(getField(), y.length);
            Array2DRowFieldMatrix<T> predictedNordsieck = null;
            double error = 10;
            while (error >= 1.0) {

                // predict a first estimate of the state at step end
                predictedY = stepEnd.getCompleteState();

                // evaluate the derivative
                final T[] yDot = computeDerivatives(stepEnd.getTime(), predictedY);

                // predict Nordsieck vector at step end
                for (int j = 0; j < predictedScaled.length; ++j) {
                    predictedScaled[j] = getStepSize().multiply(yDot[j]);
                }
                predictedNordsieck = updateHighOrderDerivativesPhase1(nordsieck);
                updateHighOrderDerivativesPhase2(scaled, predictedScaled, predictedNordsieck);

                // evaluate error
                error = errorEstimation(y, stepEnd.getTime(), predictedY, predictedScaled, predictedNordsieck);
                if (Double.isNaN(error)) {
                    throw new MathIllegalStateException(LocalizedODEFormats.NAN_APPEARING_DURING_INTEGRATION,
                                                        stepEnd.getTime().getReal());
                }

                if (error >= 1.0) {
                    // reject the step and attempt to reduce error by stepsize control
                    final double factor = computeStepGrowShrinkFactor(error);
                    rescale(getStepSizeHelper().filterStep(getStepSize().multiply(factor), forward, false));
                    stepEnd = AdamsFieldStateInterpolator.taylor(equations.getMapper(), getStepStart(),
                                                                 getStepStart().getTime().add(getStepSize()),
                                                                 getStepSize(),
                                                                 scaled,
                                                                 nordsieck);

                }
            }

            final AdamsFieldStateInterpolator<T> interpolator =
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

                    final T  nextT      = getStepStart().getTime().add(getStepSize());
                    final boolean nextIsLast = forward ?
                                               nextT.subtract(finalTime).getReal() >= 0 :
                                               nextT.subtract(finalTime).getReal() <= 0;
                    final T hNew = nextIsLast ? finalTime.subtract(getStepStart().getTime()) : getStepSize();

                    rescale(hNew);
                    System.arraycopy(getStepStart().getCompleteState(), 0, y, 0, y.length);

                } else {

                    // stepsize control for next step
                    final double  factor     = computeStepGrowShrinkFactor(error);
                    final T       scaledH    = getStepSize().multiply(factor);
                    final T       nextT      = getStepStart().getTime().add(scaledH);
                    final boolean nextIsLast = forward ?
                                               nextT.subtract(finalTime).getReal() >= 0 :
                                               nextT.subtract(finalTime).getReal() <= 0;
                    T hNew = getStepSizeHelper().filterStep(scaledH, forward, nextIsLast);

                    final T       filteredNextT      = getStepStart().getTime().add(hNew);
                    final boolean filteredNextIsLast = forward ?
                                                       filteredNextT.subtract(finalTime).getReal() >= 0 :
                                                       filteredNextT.subtract(finalTime).getReal() <= 0;
                    if (filteredNextIsLast) {
                        hNew = finalTime.subtract(getStepStart().getTime());
                    }

                    rescale(hNew);
                    System.arraycopy(predictedY, 0, y, 0, y.length);

                }

                stepEnd = AdamsFieldStateInterpolator.taylor(equations.getMapper(), getStepStart(),
                                                             getStepStart().getTime().add(getStepSize()),
                                                             getStepSize(), scaled, nordsieck);

            }

        } while (!isLastStep());

        final FieldODEStateAndDerivative<T> finalState = getStepStart();
        setStepStart(null);
        setStepSize(null);
        return finalState;

    }

    /** {@inheritDoc} */
    @Override
    protected Array2DRowFieldMatrix<T> initializeHighOrderDerivatives(final T h, final T[] t,
                                                                      final T[][] y,
                                                                      final T[][] yDot) {
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
     * @see #updateHighOrderDerivativesPhase2(CalculusFieldElement[], CalculusFieldElement[], Array2DRowFieldMatrix)
     */
    public Array2DRowFieldMatrix<T> updateHighOrderDerivativesPhase1(final Array2DRowFieldMatrix<T> highOrder) {
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
     * @see #updateHighOrderDerivativesPhase1(Array2DRowFieldMatrix)
     */
    public void updateHighOrderDerivativesPhase2(final T[] start, final T[] end,
                                                 final Array2DRowFieldMatrix<T> highOrder) {
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
    protected abstract double errorEstimation(T[] previousState, T predictedTime,
                                              T[] predictedState, T[] predictedScaled,
                                              FieldMatrix<T> predictedNordsieck);

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
    protected abstract AdamsFieldStateInterpolator<T> finalizeStep(T stepSize, T[] predictedState,
                                                                   T[] predictedScaled, Array2DRowFieldMatrix<T> predictedNordsieck,
                                                                   boolean isForward,
                                                                   FieldODEStateAndDerivative<T> globalPreviousState,
                                                                   FieldODEStateAndDerivative<T> globalCurrentState,
                                                                   FieldEquationsMapper<T> equationsMapper);

}
