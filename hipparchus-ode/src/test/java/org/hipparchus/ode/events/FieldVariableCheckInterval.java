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

import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853FieldIntegrator;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/** Tests for variable check interval.
 */
public class FieldVariableCheckInterval implements FieldOrdinaryDifferentialEquation<Binary64> {

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
                   if (s.getTime().getReal() < tZero - 0.5 * width) {
                       return tZero - 0.25 * width - s.getTime().getReal();
                   } else if (s.getTime().getReal() > tZero + 0.5 * width) {
                       return s.getTime().getReal() - (tZero + 0.25 * width);
                   } else {
                       return width / 25;
                   }
               },
               21);
    }

    private void doTest(final double tZero, final double width, final FieldAdaptableInterval<Binary64> checkInterval,
                        final int expectedCalls) {
        double e = 1e-15;
        FieldODEIntegrator<Binary64> integrator = new DormandPrince853FieldIntegrator<>(Binary64Field.getInstance(),
                                                                                        e, 100.0, 1e-7, 1e-7);
        Event evt = new Event(checkInterval, e, 999, tZero, width);
        integrator.addEventDetector(evt);
        Binary64 t = new Binary64(0.0);
        Binary64 tEnd = new Binary64(9.75);
        Binary64[] y = { new Binary64(0.0), new Binary64(0.0) };
        final FieldODEStateAndDerivative<Binary64> finalState =
                        integrator.integrate(new FieldExpandableODE<>(this), new FieldODEState<>(t, y), tEnd);
        t = finalState.getTime();
        Assert.assertEquals(tZero, finalState.getTime().getReal(), e);
        Assert.assertEquals(expectedCalls, evt.count);
     }

    /** {@inheritDoc} */
    public int getDimension() {
        return 2;
    }

    /** {@inheritDoc} */
    public Binary64[] computeDerivatives(Binary64 t, Binary64[] y) {
        return new Binary64[] { new Binary64(1.0), new Binary64(2.0) };
    }

    /** State events for this unit test. */
    private class Event implements FieldODEEventDetector<Binary64> {

        private final FieldAdaptableInterval<Binary64>             maxCheck;
        private final int                                          maxIter;
        private final FieldBracketingNthOrderBrentSolver<Binary64> solver;
        private final Binary64                                     tZero;
        private final Binary64                                     width;
        private       int                                          count;

        /** Constructor for the {@link Event} class.
         * @param maxCheck maximum checking interval
         * @param threshold convergence threshold (s)
         * @param maxIter maximum number of iterations in the event time search
         * @param tZero time of zero crossing event
         * @param width width of variation interval
         */
        public Event(final FieldAdaptableInterval<Binary64> maxCheck, final double threshold, final int maxIter,
                     final double tZero, final double width) {
            this.maxCheck  = maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new FieldBracketingNthOrderBrentSolver<>(Binary64.ZERO, new Binary64(threshold),
                                                                      Binary64.ZERO, 5);
            this.tZero     = new Binary64(tZero);
            this.width     = new Binary64(width);
        }

        public FieldAdaptableInterval<Binary64> getMaxCheckInterval() {
            return maxCheck;
        }

        public int getMaxIterationCount() {
            return maxIter;
        }

        public BracketedRealFieldUnivariateSolver<Binary64> getSolver() {
            return solver;
        }

        public FieldODEEventHandler<Binary64> getHandler() {
            return (state, detector, increasing) -> Action.STOP;
        }

        /** {@inheritDoc} */
        public Binary64 g(final FieldODEStateAndDerivative<Binary64> s) {
            ++count;
            final Binary64 t = s.getTime();
            return FastMath.max(width.multiply(-0.5), FastMath.min(width.multiply(0.5), t.subtract(tZero)));
        }

        /** {@inheritDoc} */
        public void init(final FieldODEStateAndDerivative<Binary64> initialState, final Binary64 finalTime) {
            count = 0;
        }

    }

}
