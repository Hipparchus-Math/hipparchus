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
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.util.MathArrays;

/**
 * This class implements a second order Runge-Kutta integrator for
 * Ordinary Differential Equations.
 *
 * <p>This method is an explicit Runge-Kutta method, its Butcher-array
 * is the following one :</p>
 * <pre>
 *    0  |  0    0
 *   1/2 | 1/2   0
 *       |----------
 *       |  0    1
 * </pre>
 *
 * @see EulerFieldIntegrator
 * @see ClassicalRungeKuttaFieldIntegrator
 * @see GillFieldIntegrator
 * @see ThreeEighthesFieldIntegrator
 * @see LutherFieldIntegrator
 *
 * @param <T> the type of the field elements
 */

public class MidpointFieldIntegrator<T extends CalculusFieldElement<T>> extends RungeKuttaFieldIntegrator<T> {

    /** Simple constructor.
     * Build a midpoint integrator with the given step.
     * @param field field to which the time and state vector elements belong
     * @param step integration step
     */
    public MidpointFieldIntegrator(final Field<T> field, final T step) {
        super(field, "midpoint", step);
    }

    /** {@inheritDoc} */
    @Override
    public T[] getC() {
        final T[] c = MathArrays.buildArray(getField(), 1);
        c[0] = getField().getOne().newInstance(0.5);
        return c;
    }

    /** {@inheritDoc} */
    @Override
    public T[][] getA() {
        final T[][] a = MathArrays.buildArray(getField(), 1, 1);
        a[0][0] = FieldExplicitRungeKuttaIntegrator.fraction(getField(), 1, 2);
        return a;
    }

    /** {@inheritDoc} */
    @Override
    public T[] getB() {
        final T[] b = MathArrays.buildArray(getField(), 2);
        b[0] = getField().getZero();
        b[1] = getField().getOne();
        return b;
    }

    /** {@inheritDoc} */
    @Override
    protected MidpointFieldStateInterpolator<T>
        createInterpolator(final boolean forward, T[][] yDotK,
                           final FieldODEStateAndDerivative<T> globalPreviousState,
                           final FieldODEStateAndDerivative<T> globalCurrentState,
                           final FieldEquationsMapper<T> mapper) {
        return new MidpointFieldStateInterpolator<T>(getField(), forward, yDotK,
                                                    globalPreviousState, globalCurrentState,
                                                    globalPreviousState, globalCurrentState,
                                                    mapper);
    }

}
