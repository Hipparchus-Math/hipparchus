/*
 * Licensed to the Hipparchus project under one or more
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
package org.hipparchus.ode.events;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.nonstiff.interpolators.ClassicalRungeKuttaFieldStateInterpolator;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.mockito.Mockito;

class FieldDetectorBasedEventStateTest {

    // Unit test reproducing https://gitlab.orekit.org/orekit/orekit/-/issues/1808
    @Test
    void testEpochComparisonAtLeastSignificantBit() throws NoSuchFieldException, IllegalAccessException {
        final Binary64Field field = Binary64Field.getInstance();
        final Binary64 zero = field.getZero();
        // Epoch of event
        final Binary64 eventTime = zero.add(17016.237999999998);

        // Get the interpolated state at event time
        // It will return globalCurrent at 17016.238 sec since the difference between current and previous state times is smaller than the least bit
        final Binary64[] array = MathArrays.buildArray(field, 2);
        final FieldODEStateAndDerivative<Binary64> globalCurrent = new FieldODEStateAndDerivative<>(zero.add(17016.238), array, array);
        final FieldODEStateAndDerivative<Binary64> globalPrevious = new FieldODEStateAndDerivative<>(zero.add(17016.237999999998), array, array);
        final ClassicalRungeKuttaFieldStateInterpolator<Binary64> interpolator = new ClassicalRungeKuttaFieldStateInterpolator<>(field, true, MathArrays.buildArray(field, 2, 2),
                                                                                                                                 globalPrevious, globalCurrent,
                                                                                                                                 globalPrevious, globalCurrent, null);
        final FieldODEStateAndDerivative<Binary64> interpolatedState = interpolator.getInterpolatedState(eventTime);
        Assertions.assertEquals(interpolatedState.getTime().getReal(), globalCurrent.getTime().getReal());
        Assertions.assertNotEquals(interpolatedState.getTime().getReal(), globalPrevious.getTime().getReal());

        // Configure the event state (failing before the fix)
        // Since detecting the event causing the numerical issue is tricky; we access the private field to simplify the workflow and directly set the necessary values causing the issue
        final FieldDetectorBasedEventState<Binary64> es = new FieldDetectorBasedEventState<>(new TestFieldDetector<>(field, true));
        final java.lang.reflect.Field pendingEvent = FieldDetectorBasedEventState.class.getDeclaredField("pendingEvent");
        pendingEvent.setAccessible(true);
        pendingEvent.set(es, true);
        final java.lang.reflect.Field pendingEventTime = FieldDetectorBasedEventState.class.getDeclaredField("pendingEventTime");
        pendingEventTime.setAccessible(true);
        pendingEventTime.set(es, eventTime);
        final java.lang.reflect.Field afterG = FieldDetectorBasedEventState.class.getDeclaredField("afterG");
        afterG.setAccessible(true);
        afterG.set(es, zero); // Dummy value (this value is not interesting in that case)

        // Action & verify
        Assertions.assertNotNull(es.doEvent(interpolatedState));
    }

    @Test
    void testDoEventThrowsIfTimeMismatch() throws NoSuchFieldException, IllegalAccessException {
        final Binary64Field field = Binary64Field.getInstance();
        final Binary64 zero = field.getZero();
        // Initialization
        final FieldODEEventDetector<Binary64> detector = new DummyDetector<>(field);
        final FieldDetectorBasedEventState<Binary64> eventState = new FieldDetectorBasedEventState<>(detector);
        java.lang.reflect.Field pendingEvent = FieldDetectorBasedEventState.class.getDeclaredField("pendingEvent");
        pendingEvent.setAccessible(true);
        pendingEvent.set(eventState, true);
        final java.lang.reflect.Field pendingEventTime = FieldDetectorBasedEventState.class.getDeclaredField("pendingEventTime");
        pendingEventTime.setAccessible(true);
        pendingEventTime.set(eventState, field.getOne());
        final Binary64[] array = MathArrays.buildArray(field, 1);
        final FieldODEStateAndDerivative<Binary64> state = new FieldODEStateAndDerivative<>(zero.add(1.0001), array, array);
        // Action & verify
        Assertions.assertThrows(org.hipparchus.exception.MathRuntimeException.class, () -> eventState.doEvent(state));
    }

    @ParameterizedTest
    @ValueSource(booleans = {true, false})
    void testNextCheck(final boolean isForward) {
        // GIVEN
        final TestFieldDetector<Complex> detector = new TestFieldDetector<>(ComplexField.getInstance(), isForward);
        final FieldDetectorBasedEventState<Complex> eventState = new FieldDetectorBasedEventState<>(detector);
        final FieldODEStateInterpolator<Complex> mockedInterpolator = Mockito.mock(FieldODEStateInterpolator.class);
        final FieldODEStateAndDerivative<Complex> stateAndDerivative1 = getStateAndDerivative(1);
        final FieldODEStateAndDerivative<Complex> stateAndDerivative2 = getStateAndDerivative(-1);
        if (isForward) {
            Mockito.when(mockedInterpolator.getCurrentState()).thenReturn(stateAndDerivative1);
            Mockito.when(mockedInterpolator.getPreviousState()).thenReturn(stateAndDerivative2);
        } else {
            Mockito.when(mockedInterpolator.getCurrentState()).thenReturn(stateAndDerivative2);
            Mockito.when(mockedInterpolator.getPreviousState()).thenReturn(stateAndDerivative1);
        }
        Mockito.when(mockedInterpolator.isForward()).thenReturn(isForward);
        Mockito.when(mockedInterpolator.getInterpolatedState(new Complex(0.))).thenReturn(getStateAndDerivative(0.));
        eventState.init(mockedInterpolator.getPreviousState(), mockedInterpolator.getPreviousState().getTime());
        eventState.reinitializeBegin(mockedInterpolator);
        // WHEN & THEN
        final AssertionError error = Assertions.assertThrows(AssertionError.class, () ->
                eventState.evaluateStep(mockedInterpolator));
        Assertions.assertEquals(isForward ? "forward" : "backward", error.getMessage());
    }

    private static FieldODEStateAndDerivative<Complex> getStateAndDerivative(final double time) {
        final Complex[] state = MathArrays.buildArray(ComplexField.getInstance(), 1);
        state[0] = new Complex(time);
        final Complex[] derivative = MathArrays.buildArray(ComplexField.getInstance(), 1);
        derivative[0] = Complex.ONE;
        return new FieldODEStateAndDerivative<>(state[0], state, derivative);
    }

    private static class TestFieldDetector<T extends CalculusFieldElement<T>> implements FieldODEEventDetector<T> {

        private final Field<T> field;
        private final boolean failOnForward;

        TestFieldDetector(final Field<T> field, final boolean failOnForward) {
            this.field = field;
            this.failOnForward = failOnForward;
        }

        @Override
        public FieldAdaptableInterval<T> getMaxCheckInterval() {
            return (state, isForward) -> {
                if (isForward && failOnForward) {
                    throw new AssertionError("forward");
                } else if (!isForward && !failOnForward) {
                    throw new AssertionError("backward");
                }
                return 1.;
            };
        }

        @Override
        public int getMaxIterationCount() {
            return 10;
        }

        @Override
        public BracketedRealFieldUnivariateSolver<T> getSolver() {
            return new FieldBracketingNthOrderBrentSolver<>(field.getOne(), field.getOne(), field.getOne(), 2);
        }

        @Override
        public FieldODEEventHandler<T> getHandler() {
            return (s, e, d) -> Action.CONTINUE;
        }

        @Override
        public T g(FieldODEStateAndDerivative<T> state) {
            return state.getTime();
        }
    }

    private static class DummyDetector<T extends CalculusFieldElement<T>> implements FieldODEEventDetector<T> {

        private final Field<T> field;

        public DummyDetector(final Field<T> field) {
            this.field = field;
        }

        @Override
        public FieldAdaptableInterval<T> getMaxCheckInterval() {
            return (state, isForward) -> 1.0;
        }

        @Override
        public int getMaxIterationCount() {
            return 10;
        }

        @Override
        public FieldBracketingNthOrderBrentSolver<T> getSolver() {
            return new FieldBracketingNthOrderBrentSolver<>(field.getZero(), field.getZero(), field.getZero(), 2);
        }

        @Override
        public FieldODEEventHandler<T> getHandler() {
            return (state, detector, increasing) -> Action.CONTINUE;
        }

        @Override
        public T g(FieldODEStateAndDerivative<T> state) {
            return state.getTime().getField().getZero();
        }
    }

}
