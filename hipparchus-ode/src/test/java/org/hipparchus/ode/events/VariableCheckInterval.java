/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.ode.events;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/** Tests for variable check interval.
 */
public class VariableCheckInterval implements OrdinaryDifferentialEquation {

    @Test
    public void testFixedInterval() {
        double tZero = 7.0;
        double width = 0.25;
        doTest(tZero, width, s -> width / 25, 710);
    }

    @Test
    public void testWidthAwareInterval() {
        double tZero = 7.0;
        double width = 0.25;
        doTest(tZero, width,
               s -> {
                   if (s.getTime() < tZero - 0.5 * width) {
                       return tZero - 0.25 * width - s.getTime();
                   } else if (s.getTime() > tZero + 0.5 * width) {
                       return s.getTime() - (tZero + 0.25 * width);
                   } else {
                       return width / 25;
                   }
               },
               21);
    }

    private void doTest(final double tZero, final double width, final AdaptableInterval checkInterval,
                        final int expectedCalls) {
        double e = 1e-15;
        ODEIntegrator integrator = new DormandPrince853Integrator(e, 100.0, 1e-7, 1e-7);
        Event evt = new Event(checkInterval, e, 999, tZero, width);
        integrator.addEventDetector(evt);
        double t = 0.0;
        double tEnd = 9.75;
        double[] y = {0.0, 0.0};
        final ODEStateAndDerivative finalState =
                        integrator.integrate(this, new ODEState(t, y), tEnd);
        t = finalState.getTime();
        Assert.assertEquals(tZero, finalState.getTime(), e);
        Assert.assertEquals(expectedCalls, evt.count);
     }

    /** {@inheritDoc} */
    public int getDimension() {
        return 2;
    }

    /** {@inheritDoc} */
    public double[] computeDerivatives(double t, double[] y) {
        return new double[] { 1.0, 2.0 };
    }

    /** State events for this unit test. */
    private class Event implements ODEEventDetector {

        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;
        private final double                        tZero;
        private final double                        width;
        private       int                           count;

        /** Constructor for the {@link Event} class.
         * @param maxCheck maximum checking interval
         * @param threshold convergence threshold (s)
         * @param maxIter maximum number of iterations in the event time search
         * @param tZero time of zero crossing event
         * @param width width of variation interval
         */
        public Event(final AdaptableInterval maxCheck, final double threshold, final int maxIter,
                     final double tZero, final double width) {
            this.maxCheck  = maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.tZero     = tZero;
            this.width     = width;
        }

        public AdaptableInterval getMaxCheckInterval() {
            return maxCheck;
        }

        public int getMaxIterationCount() {
            return maxIter;
        }

        public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
            return solver;
        }

        public ODEEventHandler getHandler() {
            return (state, detector, increasing) -> Action.STOP;
        }

        /** {@inheritDoc} */
        public double g(final ODEStateAndDerivative s) {
            ++count;
            final double t = s.getTime();
            return FastMath.max(-0.5 * width, FastMath.min(0.5 * width, t - tZero));
        }

        /** {@inheritDoc} */
        public void init(final ODEStateAndDerivative initialState, final double finalTime) {
            count = 0;
        }

    }

}
