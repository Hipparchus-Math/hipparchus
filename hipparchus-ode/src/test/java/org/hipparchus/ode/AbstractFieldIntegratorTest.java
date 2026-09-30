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

package org.hipparchus.ode;

import org.hipparchus.Field;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.FieldAdaptableInterval;
import org.hipparchus.ode.events.FieldODEEventDetector;
import org.hipparchus.ode.events.FieldODEEventHandler;
import org.hipparchus.ode.nonstiff.EulerFieldIntegrator;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.ode.sampling.FieldODEStepHandler;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class AbstractFieldIntegratorTest {

    @Test
    void testIntegrateWithResetDerivativesAndEventDetector() {
        // GIVEN
        final Field<Binary64> field = Binary64Field.getInstance();
        final Binary64 finalTime = new Binary64(1.);
        final EulerFieldIntegrator<Binary64> integrator = new EulerFieldIntegrator<>(field, finalTime);
        final TestDetector detector = new TestDetector(0.5, Action.RESET_DERIVATIVES);
        integrator.addEventDetector(detector);
        final TestFieldProblem1<Binary64> testProblem = new TestFieldProblem1<>(field);
        final FieldODEState<Binary64> initialState = new FieldODEState<>(Binary64.ZERO, MathArrays.buildArray(field, 2));
        // WHEN
        integrator.integrate(testProblem, initialState, finalTime);
        // THEN
        assertTrue(detector.resetted);
    }

    @Test
    void testUpdateStepIsCalledOncePerStepWhileHandleStepIsCalledAtEachEvent() {
        // GIVEN
        final Field<Binary64> field = Binary64Field.getInstance();
        final Binary64 finalTime = new Binary64(3.);
        final EulerFieldIntegrator<Binary64> integrator = new EulerFieldIntegrator<>(field, finalTime);
        integrator.addStepHandler(new UpdateStepTestHandler());
        integrator.addEventDetector(new TestDetector(0.5, Action.CONTINUE));
        integrator.addEventDetector(new TestDetector(0.6, Action.CONTINUE));
        final TestFieldProblem1<Binary64> testProblem = new TestFieldProblem1<>(field);
        final FieldODEState<Binary64> initialState = new FieldODEState<>(Binary64.ZERO, MathArrays.buildArray(field, 2));
        // WHEN
        integrator.integrate(testProblem, initialState, finalTime);
        // THEN
        assertEquals(1, ((UpdateStepTestHandler) integrator.getStepHandlers().getFirst()).getUpdateStepCounter());
        assertEquals(3, ((UpdateStepTestHandler) integrator.getStepHandlers().getFirst()).getHandlerStepCounter());
    }

    private static class TestDetector implements FieldODEEventDetector<Binary64> {

        boolean resetted = false;
        double rootTime;
        Action action;

        public TestDetector(double rootTime, Action action) {
            this.rootTime = rootTime;
            this.action = action;
        }

        @Override
        public void reset(FieldODEStateAndDerivative<Binary64> intermediateState, Binary64 finalTime) {
            FieldODEEventDetector.super.reset(intermediateState, finalTime);
            resetted = true;
        }

        @Override
        public FieldAdaptableInterval<Binary64> getMaxCheckInterval() {
            return FieldAdaptableInterval.of(1);
        }

        @Override
        public int getMaxIterationCount() {
            return 10;
        }

        @Override
        public BracketedRealFieldUnivariateSolver<Binary64> getSolver() {
            return new FieldBracketingNthOrderBrentSolver<>(new Binary64(1e-14), new Binary64(1e-6), new Binary64(1e-15), 5);
        }

        @Override
        public FieldODEEventHandler<Binary64> getHandler() {
            return (state, detector, increasing) -> action;
        }

        @Override
        public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
            return state.getTime().subtract(rootTime);
        }
    }

    private static class UpdateStepTestHandler implements FieldODEStepHandler<Binary64> {

        private int handlerStepCounter = 0;
        private int updateStepCounter = 0;

        @Override
        public void handleStep(FieldODEStateInterpolator<Binary64> interpolator) {
            handlerStepCounter++;
        }

        @Override
        public void updateOnStep(FieldODEStateInterpolator<Binary64> interpolator) {
            updateStepCounter++;
        }

        public int getHandlerStepCounter() {
            return handlerStepCounter;
        }

        public int getUpdateStepCounter() {
            return updateStepCounter;
        }
    }

}