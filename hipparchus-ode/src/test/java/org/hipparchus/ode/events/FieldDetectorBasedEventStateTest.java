package org.hipparchus.ode.events;

import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.mockito.Mockito;

class FieldDetectorBasedEventStateTest {

    @ParameterizedTest
    @ValueSource(booleans = {true, false})
    void testNextCheck(final boolean isForward) {
        // GIVEN
        final TestFieldDetector detector = new TestFieldDetector(isForward);
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

    private static class TestFieldDetector implements FieldODEEventDetector<Complex> {

        private final boolean failOnForward;

        TestFieldDetector(final boolean failOnForward) {
            this.failOnForward = failOnForward;
        }

        @Override
        public FieldAdaptableInterval<Complex> getMaxCheckInterval() {
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
        public BracketedRealFieldUnivariateSolver<Complex> getSolver() {
            return new FieldBracketingNthOrderBrentSolver<>(Complex.ONE, Complex.ONE, Complex.ONE, 2);
        }

        @Override
        public FieldODEEventHandler<Complex> getHandler() {
            return (s, e, d) -> Action.CONTINUE;
        }

        @Override
        public Complex g(FieldODEStateAndDerivative<Complex> state) {
            return state.getTime();
        }
    }

}
