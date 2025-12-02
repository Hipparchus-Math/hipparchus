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
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.nonstiff.interpolators.RungeKuttaFieldStateInterpolator;
import org.hipparchus.util.MathArrays;

/**
 * This class implements the common part of all fixed step Runge-Kutta
 * integrators for Ordinary Differential Equations.
 *
 * <p>These methods are explicit Runge-Kutta methods, their Butcher
 * arrays are as follows :</p>
 * <pre>
 *    0  |
 *   c2  | a21
 *   c3  | a31  a32
 *   ... |        ...
 *   cs  | as1  as2  ...  ass-1
 *       |--------------------------
 *       |  b1   b2  ...   bs-1  bs
 * </pre>
 *
 * @see EulerFieldIntegrator
 * @see ClassicalRungeKuttaFieldIntegrator
 * @see GillFieldIntegrator
 * @see MidpointFieldIntegrator
 * @param <T> the type of the field elements
 */

public abstract class FixedStepRungeKuttaFieldIntegrator<T extends CalculusFieldElement<T>>
    extends AbstractFieldIntegrator<T> implements FieldExplicitRungeKuttaIntegrator<T> {

    /** Time steps from Butcher array (without the first zero). */
    private final T[] c;

    /** Internal weights from Butcher array (without the first empty row). */
    private final T[][] a;

    /** External weights for the high order method from Butcher array. */
    private final T[] b;

    /** Time steps from Butcher array (without the first zero). */
    private double[] realC = new double[0];

    /** Internal weights from Butcher array (without the first empty row). Real version, optional. */
    private double[][] realA = new double[0][];

    /** External weights for the high order method from Butcher array. Real version, optional. */
    private double[] realB = new double[0];

    /** Integration step. */
    private final T step;

    /** Flag setting whether coefficients in Butcher array are interpreted as Field or real numbers. */
    private boolean usingFieldCoefficients;

    /** Simple constructor.
     * Build a Runge-Kutta integrator with the given
     * step. The default step handler does nothing.
     * @param field field to which the time and state vector elements belong
     * @param name name of the method
     * @param step integration step
     */
    protected FixedStepRungeKuttaFieldIntegrator(final Field<T> field, final String name, final T step) {
        super(field, name);
        this.c    = getC();
        this.a    = getA();
        this.b    = getB();
        this.step = step.abs();
        this.usingFieldCoefficients = false;
        setStepSize(getDefaultStep());
    }

    /** Getter for the default, positive step-size assigned at constructor level.
     * @return step
     */
    public T getDefaultStep() {
        return this.step;
    }

    /**
     * Setter for the flag between real or Field coefficients in the Butcher array.
     *
     * @param usingFieldCoefficients new value for flag
     */
    public void setUsingFieldCoefficients(boolean usingFieldCoefficients) {
        this.usingFieldCoefficients = usingFieldCoefficients;
    }

    /** {@inheritDoc} */
    @Override
    public boolean isUsingFieldCoefficients() {
        return usingFieldCoefficients;
    }

    /** {@inheritDoc} */
    @Override
    public int getNumberOfStages() {
        return b.length;
    }

    /** Create an interpolator.
     * @param forward integration direction indicator
     * @param yDotK slopes at the intermediate points
     * @param globalPreviousState start of the global step
     * @param globalCurrentState end of the global step
     * @param mapper equations mapper for the all equations
     * @return external weights for the high order method from Butcher array
     */
    protected abstract RungeKuttaFieldStateInterpolator<T> createInterpolator(boolean forward, T[][] yDotK,
                                                                              FieldODEStateAndDerivative<T> globalPreviousState,
                                                                              FieldODEStateAndDerivative<T> globalCurrentState,
                                                                              FieldEquationsMapper<T> mapper);

    /** {@inheritDoc} */
    @Override
    protected FieldODEStateAndDerivative<T> initIntegration(FieldExpandableODE<T> eqn, FieldODEState<T> s0, T t) {
        if (!isUsingFieldCoefficients()) {
            realA = getRealA();
            realB = getRealB();
            realC = getRealC();
        }
        return super.initIntegration(eqn, s0, t);
    }

    /** {@inheritDoc} */
    @Override
    public FieldODEStateAndDerivative<T> integrate(final FieldExpandableODE<T> equations,
                                                   final FieldODEState<T> initialState, final T finalTime)
        throws MathIllegalArgumentException, MathIllegalStateException {

        sanityChecks(initialState, finalTime);
        setStepStart(initIntegration(equations, initialState, finalTime));
        final boolean forward = finalTime.subtract(initialState.getTime()).getReal() > 0;

        // create some internal working arrays
        final int   stages = getNumberOfStages();
        final T[][] yDotK  = MathArrays.buildArray(getField(), stages, -1);
        MathArrays.buildArray(getField(), equations.getMapper().getTotalDimension());

        // set up integration control objects
        if (forward) {
            if (getStepStart().getTime().add(step).subtract(finalTime).getReal() >= 0) {
                setStepSize(finalTime.subtract(getStepStart().getTime()));
            } else {
                setStepSize(step);
            }
        } else {
            if (getStepStart().getTime().subtract(step).subtract(finalTime).getReal() <= 0) {
                setStepSize(finalTime.subtract(getStepStart().getTime()));
            } else {
                setStepSize(step.negate());
            }
        }

        // main integration loop
        setIsLastStep(false);
        do {

            // first stage
            final T[] y = getStepStart().getCompleteState();
            yDotK[0]    = getStepStart().getCompleteDerivative();

            // next stages
            final T[] yTmp;
            if (isUsingFieldCoefficients()) {
                FieldExplicitRungeKuttaIntegrator.applyInternalButcherWeights(getEquations(), getStepStart().getTime(),
                        y, getStepSize(), a, c, yDotK);
                yTmp = FieldExplicitRungeKuttaIntegrator.applyExternalButcherWeights(y, yDotK, getStepSize(), b);
            } else {
                FieldExplicitRungeKuttaIntegrator.applyInternalButcherWeights(getEquations(), getStepStart().getTime(),
                        y, getStepSize(), realA, realC, yDotK);
                yTmp = FieldExplicitRungeKuttaIntegrator.applyExternalButcherWeights(y, yDotK, getStepSize(), realB);
            }

            incrementEvaluations(stages - 1);

            final T stepEnd   = getStepStart().getTime().add(getStepSize());
            final T[] yDotTmp = computeDerivatives(stepEnd, yTmp);
            final FieldODEStateAndDerivative<T> stateTmp = equations.getMapper().mapStateAndDerivative(stepEnd, yTmp, yDotTmp);

            // discrete events handling
            setStepStart(acceptStep(createInterpolator(forward, yDotK, getStepStart(), stateTmp, equations.getMapper()),
                                    finalTime));

            if (!isLastStep()) {

                // stepsize control for next step
                final T  nextT      = getStepStart().getTime().add(getStepSize());
                final boolean nextIsLast = forward ?
                                           (nextT.subtract(finalTime).getReal() >= 0) :
                                           (nextT.subtract(finalTime).getReal() <= 0);
                if (nextIsLast) {
                    setStepSize(finalTime.subtract(getStepStart().getTime()));
                }
            }

        } while (!isLastStep());

        final FieldODEStateAndDerivative<T> finalState = getStepStart();
        setStepStart(null);
        setStepSize(null);
        return finalState;

    }

}
