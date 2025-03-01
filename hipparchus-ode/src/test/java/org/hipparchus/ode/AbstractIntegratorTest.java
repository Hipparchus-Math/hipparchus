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

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.AdaptableInterval;
import org.hipparchus.ode.events.ODEEventDetector;
import org.hipparchus.ode.events.ODEEventHandler;
import org.hipparchus.ode.nonstiff.EulerIntegrator;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

class AbstractIntegratorTest {

    @Test
    void testIntegrateWithResetDerivativesAndEventDetector() {
        // GIVEN
        final double finalTime = 1.;
        final EulerIntegrator integrator = new EulerIntegrator(finalTime);
        final TestDetector detector = new TestDetector();
        integrator.addEventDetector(detector);
        final TestProblem1 testProblem = new TestProblem1();
        final ODEState initialState = new ODEState(0., new double[2]);
        // WHEN
        integrator.integrate(testProblem, initialState, finalTime);
        // THEN
        assertTrue(detector.resetted);
    }

    private static class TestDetector implements ODEEventDetector {
        boolean resetted = false;

        @Override
        public void reset(ODEStateAndDerivative intermediateState, double finalTime) {
            ODEEventDetector.super.reset(intermediateState, finalTime);
            resetted = true;
        }

        @Override
        public AdaptableInterval getMaxCheckInterval() {
            return AdaptableInterval.of(1);
        }

        @Override
        public int getMaxIterationCount() {
            return 10;
        }

        @Override
        public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
            return new BracketingNthOrderBrentSolver();
        }

        @Override
        public ODEEventHandler getHandler() {
            return (state, detector, increasing) -> Action.RESET_DERIVATIVES;
        }

        @Override
        public double g(ODEStateAndDerivative state) {
            return state.getTime() - 0.5;
        }
    }

}
