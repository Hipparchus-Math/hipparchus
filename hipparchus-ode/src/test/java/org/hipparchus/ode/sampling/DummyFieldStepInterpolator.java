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

package org.hipparchus.ode.sampling;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldODEStateAndDerivative;

public class DummyFieldStepInterpolator<T extends CalculusFieldElement<T>>
    extends AbstractFieldODEStateInterpolator<T> {

    public DummyFieldStepInterpolator(final boolean forward,
                                      final FieldODEStateAndDerivative<T> globalPreviousState,
                                      final FieldODEStateAndDerivative<T> globalCurrentState,
                                      final FieldODEStateAndDerivative<T> softPreviousState,
                                      final FieldODEStateAndDerivative<T> softCurrentState,
                                      final FieldEquationsMapper<T> mapper) {
        super(forward, globalPreviousState, globalCurrentState, softPreviousState, softCurrentState, mapper);
    }

    @Override
    protected AbstractFieldODEStateInterpolator<T> create(final boolean newForward,
                                                      final FieldODEStateAndDerivative<T> newGlobalPreviousState,
                                                      final FieldODEStateAndDerivative<T> newGlobalCurrentState,
                                                      final FieldODEStateAndDerivative<T> newSoftPreviousState,
                                                      final FieldODEStateAndDerivative<T> newSoftCurrentState,
                                                      final FieldEquationsMapper<T> newMapper) {
        return new DummyFieldStepInterpolator<T>(newForward,
                                                 newGlobalPreviousState, newGlobalCurrentState,
                                                 newSoftPreviousState, newSoftCurrentState,
                                                 newMapper);
    }

    @Override
    protected FieldODEStateAndDerivative<T> computeInterpolatedStateAndDerivatives(FieldEquationsMapper<T> equationsMapper,
                                                                                   T time, T theta, T thetaH, T oneMinusThetaH) {
        return getGlobalCurrentState();
    }

}
