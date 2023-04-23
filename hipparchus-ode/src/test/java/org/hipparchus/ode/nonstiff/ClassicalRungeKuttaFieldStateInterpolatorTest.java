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


import org.hipparchus.Field;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.util.Binary64Field;
import org.junit.Test;

public class ClassicalRungeKuttaFieldStateInterpolatorTest extends RungeKuttaFieldStateInterpolatorAbstractTest {

    protected <T extends CalculusFieldElement<T>> RungeKuttaFieldStateInterpolator<T>
    createInterpolator(Field<T> field, boolean forward, T[][] yDotK,
                       FieldODEStateAndDerivative<T> globalPreviousState,
                       FieldODEStateAndDerivative<T> globalCurrentState,
                       FieldODEStateAndDerivative<T> softPreviousState,
                       FieldODEStateAndDerivative<T> softCurrentState,
                       FieldEquationsMapper<T> mapper) {
        return new ClassicalRungeKuttaFieldStateInterpolator<T>(field, forward, yDotK,
                                                               globalPreviousState, globalCurrentState,
                                                               softPreviousState, softCurrentState,
                                                               mapper);
    }

    protected <T extends CalculusFieldElement<T>> FieldButcherArrayProvider<T>
    createButcherArrayProvider(final Field<T> field) {
        return new ClassicalRungeKuttaFieldIntegrator<T>(field, field.getOne());
    }

    @Test
    public void interpolationAtBounds() {
        doInterpolationAtBounds(Binary64Field.getInstance(), 1.0e-15);
    }

    @Test
    public void interpolationInside() {
        doInterpolationInside(Binary64Field.getInstance(), 2.6e-7, 3.6e-6);
    }

    @Test
    public void nonFieldInterpolatorConsistency() {
        doNonFieldInterpolatorConsistency(Binary64Field.getInstance(), 2.8e-17, 1.2e-16, 3.4e-16, 2.1e-17);
    }

}
