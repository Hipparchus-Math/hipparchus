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
package org.hipparchus.ode.events;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853FieldIntegrator;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937a;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class FieldEventSlopeFilterTest {

    @Test
    void testHistoryIncreasingForwardField() {

        // start point: g > 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                         0.5 * FastMath.PI, 30.5 * FastMath.PI, FastMath.PI, -1);

        // start point: g = 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                         0, 30.5 * FastMath.PI, FastMath.PI, -1);

        // start point: g < 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                         1.5 * FastMath.PI, 30.5 * FastMath.PI, FastMath.PI, +1);

    }

    @Test
    void testHistoryIncreasingBackwardField() {

        // start point: g > 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                         0.5 * FastMath.PI, -30.5 * FastMath.PI, FastMath.PI, -1);

        // start point: g = 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                         0, -30.5 * FastMath.PI, FastMath.PI, +1);

        // start point: g < 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                         1.5 * FastMath.PI, -30.5 * FastMath.PI, FastMath.PI, -1);

    }

    @Test
    void testHistoryDecreasingForwardField() {

        // start point: g > 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                         0.5 * FastMath.PI, 30.5 * FastMath.PI, 0, +1);

        // start point: g = 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                         0, 30.5 * FastMath.PI, 0, +1);

        // start point: g < 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                         1.5 * FastMath.PI, 30.5 * FastMath.PI, 0, +1);

    }

    @Test
    void testHistoryDecreasingBackwardField() {

        // start point: g > 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                         0.5 * FastMath.PI, -30.5 * FastMath.PI, 0, -1);

        // start point: g = 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                         0, -30.5 * FastMath.PI, 0, -1);

        // start point: g < 0
        testHistoryField(Binary64Field.getInstance(), FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                         1.5 * FastMath.PI, -30.5 * FastMath.PI, 0, +1);

    }

    private <T extends CalculusFieldElement<T>> void testHistoryField(Field<T> field, FilterType type, double t0, double t1, double refSwitch, double signEven) {
        FieldEvent<T> onlyIncreasing = new FieldEvent<>(Double.POSITIVE_INFINITY,
                                                        field.getZero().newInstance(1.0e-10), 1000, false, true);
        FieldEventSlopeFilter<FieldEvent<T>, T> eventFilter = new FieldEventSlopeFilter<>(field, onlyIncreasing, type);
        eventFilter.init(buildStateAndDerivative(field, t0), field.getZero().add(t1));

        // first pass to set up switches history for a long period
        double h = FastMath.copySign(0.05, t1 - t0);
        double n = (int) FastMath.floor((t1 - t0) / h);
        for (int i = 0; i < n; ++i) {
            double t = t0 + i * h;
            eventFilter.g(buildStateAndDerivative(field, t));
        }

        // verify old events are preserved, even if randomly accessed
        RandomGenerator rng = new Well19937a(0xb0e7401265af8cd3l);
        for (int i = 0; i < 5000; i++) {
            double t = t0 + (t1 - t0) * rng.nextDouble();
            double g = eventFilter.g(buildStateAndDerivative(field, t)).getReal();
            int turn = (int) FastMath.floor((t - refSwitch) / (2 * FastMath.PI));
            if (turn % 2 == 0) {
                assertEquals( signEven * FastMath.sin(t), g, 1.0e-10);
            } else {
                assertEquals(-signEven * FastMath.sin(t), g, 1.0e-10);
            }
        }

    }

    private <T extends CalculusFieldElement<T>> FieldODEStateAndDerivative<T> buildStateAndDerivative(Field<T> field, double t0) {
        T t0F      = field.getZero().add(t0);
        T[] y0F    = MathArrays.buildArray(field, 2);
        y0F[0]     = FastMath.sin(t0F);
        y0F[1]     = FastMath.cos(t0F);
        T[] y0DotF = MathArrays.buildArray(field, 2);
        y0DotF[0]  = FastMath.cos(t0F);
        y0DotF[1]  = FastMath.sin(t0F).negate();
        return new FieldODEStateAndDerivative<>(t0F, y0F, y0DotF);
    }

    @Test
    void testIncreasingOnlyField() {
        doTestIncreasingOnlyField(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestIncreasingOnlyField(Field<T> field) {
        FieldODEIntegrator<T> integrator = new DormandPrince853FieldIntegrator<>(field, 1.0e-3, 100.0, 1e-7, 1e-7);
        FieldEvent<T> allEvents = new FieldEvent<>(0.1, field.getZero().newInstance(1.0e-7), 100, true, true);
        integrator.addEventDetector(allEvents);
        FieldEvent<T> onlyIncreasing = new FieldEvent<>(0.1, field.getZero().newInstance(1.0e-7), 100, false, true);
        integrator.addEventDetector(new FieldEventSlopeFilter<>(field, onlyIncreasing,
                                                                FilterType.TRIGGER_ONLY_INCREASING_EVENTS));
        T t0   = field.getZero().add(0.5 * FastMath.PI);
        T tEnd = field.getZero().add(5.5 * FastMath.PI);
        T[] y  = MathArrays.buildArray(field, 2);
        y[0]   = field.getZero();
        y[1]   = field.getOne();
        assertEquals(tEnd.getReal(),
                            integrator.integrate(new FieldExpandableODE<>(new FieldSineCosine<T>()),
                                                 new FieldODEState<>(t0, y), tEnd).getTime().getReal(),
                            1.0e-7);

        assertEquals(5, allEvents.getEventCount());
        assertEquals(2, onlyIncreasing.getEventCount());

    }

    @Test
    void testDecreasingOnlyField() {
        doTestDecreasingOnlyField(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestDecreasingOnlyField(Field<T> field) {
        FieldODEIntegrator<T> integrator = new DormandPrince853FieldIntegrator<>(field, 1.0e-3, 100.0, 1e-7, 1e-7);
        FieldEvent<T> allEvents = new FieldEvent<>(0.1, field.getZero().newInstance(1.0e-7), 100, true, true);
        integrator.addEventDetector(allEvents);
        FieldEvent<T> onlyDecreasing = new FieldEvent<>(0.1, field.getZero().newInstance(1.0e-7), 100, true, false);
        integrator.addEventDetector(new FieldEventSlopeFilter<>(field, onlyDecreasing,
                                                                FilterType.TRIGGER_ONLY_DECREASING_EVENTS));
        T t0   = field.getZero().add(0.5 * FastMath.PI);
        T tEnd = field.getZero().add(5.5 * FastMath.PI);
        T[] y  = MathArrays.buildArray(field, 2);
        y[0]   = field.getZero();
        y[1]   = field.getOne();
        assertEquals(tEnd.getReal(),
                            integrator.integrate(new FieldExpandableODE<>(new FieldSineCosine<T>()),
                                                 new FieldODEState<>(t0, y), tEnd).getTime().getReal(),
                            1.0e-7);

        assertEquals(5, allEvents.getEventCount());
        assertEquals(3, onlyDecreasing.getEventCount());

    }

    @Test
    void testTwoOppositeFiltersField() {
        doestTwoOppositeFiltersField(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doestTwoOppositeFiltersField(Field<T> field)
        throws MathIllegalArgumentException, MathIllegalStateException {
        FieldODEIntegrator<T> integrator = new DormandPrince853FieldIntegrator<>(field, 1.0e-3, 100.0, 1e-7, 1e-7);
        FieldEvent<T> allEvents = new FieldEvent<>(0.1, field.getZero().newInstance(1.0e-7), 100, true, true);
        integrator.addEventDetector(allEvents);
        FieldEvent<T> onlyIncreasing = new FieldEvent<>(0.1, field.getZero().newInstance(1.0e-7), 100, false, true);
        integrator.addEventDetector(new FieldEventSlopeFilter<>(field, onlyIncreasing,
                                                                FilterType.TRIGGER_ONLY_INCREASING_EVENTS));
        FieldEvent<T> onlyDecreasing = new FieldEvent<>(0.1, field.getZero().newInstance(1.0e-7), 100, true, false);
        integrator.addEventDetector(new FieldEventSlopeFilter<>(field, onlyDecreasing,
                                                                FilterType.TRIGGER_ONLY_DECREASING_EVENTS));
        T t0   = field.getZero().add(0.5 * FastMath.PI);
        T tEnd = field.getZero().add(5.5 * FastMath.PI);
        T[] y  = MathArrays.buildArray(field, 2);
        y[0]   = field.getZero();
        y[1]   = field.getOne();
        assertEquals(tEnd.getReal(),
                            integrator.integrate(new FieldExpandableODE<>(new FieldSineCosine<T>()),
                                                 new FieldODEState<>(t0, y), tEnd).getTime().getReal(),
                            1.0e-7);

        assertEquals(5, allEvents.getEventCount());
        assertEquals(2, onlyIncreasing.getEventCount());
        assertEquals(3, onlyDecreasing.getEventCount());

    }

    private static class FieldSineCosine<T extends CalculusFieldElement<T>> implements FieldOrdinaryDifferentialEquation<T> {
        public int getDimension() {
            return 2;
        }

        public T[] computeDerivatives(T t, T[] y) {
            final T[] yDot = MathArrays.buildArray(t.getField(), getDimension());
            yDot[0] = y[1];
            yDot[1] = y[0].negate();
            return yDot;
        }
    }

    /** State events for this unit test. */
    protected static class FieldEvent<T extends CalculusFieldElement<T>> implements FieldODEEventDetector<T> {

        private final FieldAdaptableInterval<T>             maxCheck;
        private final int                                   maxIter;
        private final BracketedRealFieldUnivariateSolver<T> solver;
        private final boolean                               expectDecreasing;
        private final boolean                               expectIncreasing;
        private int                                         eventCount;

        public FieldEvent(final double maxCheck, final T threshold, final int maxIter,
                          boolean expectDecreasing, boolean expectIncreasing) {
            this.maxCheck         = (s, isForward) -> maxCheck;
            this.maxIter          = maxIter;
            this.solver           = new FieldBracketingNthOrderBrentSolver<>(threshold.getField().getZero(),
                                                                            threshold,
                                                                            threshold.getField().getZero(),
                                                                            5);
            this.expectDecreasing = expectDecreasing;
            this.expectIncreasing = expectIncreasing;
        }

        public FieldAdaptableInterval<T> getMaxCheckInterval() {
            return maxCheck;
        }

        public int getMaxIterationCount() {
            return maxIter;
        }

        public BracketedRealFieldUnivariateSolver<T> getSolver() {
            return solver;
        }

        public int getEventCount() {
            return eventCount;
        }

        @Override
        public void init(FieldODEStateAndDerivative<T> s0, T t) {
            eventCount = 0;
        }

        public T g(FieldODEStateAndDerivative<T> s) {
            return s.getPrimaryState()[0];
        }

        public FieldODEEventHandler<T> getHandler() {
            return (FieldODEStateAndDerivative<T> s, FieldODEEventDetector<T> detector, boolean increasing) -> {
                if (increasing) {
                    assertTrue(expectIncreasing);
                } else {
                    assertTrue(expectDecreasing);
                }
                eventCount++;
                return Action.RESET_STATE;
            };
        }

    }
}
