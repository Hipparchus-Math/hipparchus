/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.hipparchus.Field;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853FieldIntegrator;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Check events are detected correctly when the event times are close.
 *
 * @author Evan Ward
 */
public class FieldCloseEventsTest {

    /** type of field. */
    private static final Field<Decimal64> field = Decimal64Field.getInstance();
    private static final Decimal64 zero = field.getZero();
    private static final Decimal64 one = field.getOne();
    private static final FieldODEState<Decimal64> initialState =
            new FieldODEState<>(zero, new Decimal64[]{zero, zero});

    @Test
    public void testCloseEventsFirstOneIsReset() {
        // setup
        // a fairly rare state to reproduce this bug. Two dates, d1 < d2, that
        // are very close. Event triggers on d1 will reset state to break out of
        // event handling loop in AbstractIntegrator.acceptStep(). At this point
        // detector2 has g0Positive == true but the event time is set to just
        // before the event so g(t0) is negative. Now on processing the
        // next step the root solver checks the sign of the start, midpoint,
        // and end of the interval so we need another event less than half a max
        // check interval after d2 so that the g function will be negative at
        // all three times. Then we get a non bracketing exception.
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(Action.RESET_DERIVATIVES, 9);
        integrator.addEventHandler(detector1, 10, 1e-9, 100);
        TimeDetector detector2 = new TimeDetector(9 + 1e-15, 9 + 4.9);
        integrator.addEventHandler(detector2, 11, 1e-9, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(20));

        // verify
        List<Event> events1 = detector1.getEvents();
        Assert.assertEquals(1, events1.size());
        Assert.assertEquals(9, events1.get(0).getT(), 0.0);
        List<Event> events2 = detector2.getEvents();
        Assert.assertEquals(0, events2.size());
    }

    @Test
    public void testCloseEvents() {
        // setup
        double e = 1e-15;
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, e, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(5);
        integrator.addEventHandler(detector1, 10, 1, 100);
        TimeDetector detector2 = new TimeDetector(5.5);
        integrator.addEventHandler(detector2, 10, 1, 100);

        // action
        integrator.integrate(new Equation(), initialState, one.add(20));

        // verify
        List<Event> events1 = detector1.getEvents();
        Assert.assertEquals(1, events1.size());
        Assert.assertEquals(5, events1.get(0).getT(), 0.0);
        List<Event> events2 = detector2.getEvents();
        Assert.assertEquals(1, events2.size());
        Assert.assertEquals(5.5, events2.get(0).getT(), 0.0);
    }

    @Test
    public void testSimultaneousEvents() {
        // setup
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(5);
        integrator.addEventHandler(detector1, 10, 1, 100);
        TimeDetector detector2 = new TimeDetector(5);
        integrator.addEventHandler(detector2, 10, 1, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(20));

        // verify
        List<Event> events1 = detector1.getEvents();
        Assert.assertEquals(1, events1.size());
        Assert.assertEquals(5, events1.get(0).getT(), 0.0);
        List<Event> events2 = detector2.getEvents();
        Assert.assertEquals(1, events2.size());
        Assert.assertEquals(5, events2.get(0).getT(), 0.0);
    }

    /**
     * test the g function switching with a period shorter than the tolerance. We don't
     * need to find any of the events, but we do need to not crash. And we need to
     * preserve the alternating increasing / decreasing sequence.
     */
    @Test
    public void testFastSwitching() {
        // setup
        // step size of 10 to land in between two events we would otherwise miss
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(9.9, 10.1, 12);
        integrator.addEventHandler(detector1, 10, 0.2, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(20));

        //verify
        // finds one or three events. Not 2.
        List<Event> events1 = detector1.getEvents();
        Assert.assertEquals(1, events1.size());
        Assert.assertEquals(9.9, events1.get(0).getT(), 0.1);
        Assert.assertEquals(true, events1.get(0).isIncreasing());
    }

    /** "A Tricky Problem" from bug #239. */
    @Test
    public void testTrickyCaseLower() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = 1.0, t2 = 15, t3 = 16, t4 = 17, t5 = 18;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t3);
        TimeDetector detectorB = new TimeDetector(events, -10, t1, t2, t5);
        TimeDetector detectorC = new TimeDetector(Action.RESET_DERIVATIVES, events, t4);

        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorC, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        // but I only know one way to do that in this case.
        Assert.assertEquals(5, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(false, events.get(0).isIncreasing());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertEquals(t3, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertEquals(t4, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(t5, events.get(4).getT(), tolerance);
        Assert.assertEquals(false, events.get(4).isIncreasing());
    }

    /**
     * Test case for two event detectors. DetectorA has event at t2, DetectorB at t3, but
     * due to the root finding tolerance DetectorB's event occurs at t1. With t1 < t2 <
     * t3.
     */
    @Test
    public void testRootFindingTolerance() {
        //setup
        double maxCheck = 10;
        double t2 = 11, t3 = t2 + 1e-5;
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t2);
        TimeDetector detectorB = new FlatDetector(events, t3);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, 1e-6, 100);
        integrator.addEventHandler(detectorB, maxCheck, 0.5, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        // if these fail the event finding did its job,
        // but this test isn't testing what it is supposed to be
        Assert.assertSame(detectorB, events.get(0).getHandler());
        Assert.assertSame(detectorA, events.get(1).getHandler());
        Assert.assertTrue(events.get(0).getT() < events.get(1).getT());

        // check event detection worked
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t3, events.get(0).getT(), 0.5);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(t2, events.get(1).getT(), 1e-6);
        Assert.assertEquals(true, events.get(1).isIncreasing());
    }

    /** check when g(t < root) < 0,  g(root + convergence) < 0. */
    @Test
    public void testRootPlusToleranceHasWrongSign() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        final double toleranceB = 0.3;
        double t1 = 11, t2 = 11.1, t3 = 11.2;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t2);
        TimeDetector detectorB = new TimeDetector(events, t1, t3);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, 1e-6, 100);
        integrator.addEventHandler(detectorB, maxCheck, toleranceB, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        // we only care that the rules are satisfied, there are other solutions
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), toleranceB);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorB, events.get(0).getHandler());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getHandler());
        Assert.assertEquals(t3, events.get(2).getT(), toleranceB);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertSame(detectorB, events.get(2).getHandler());
        // chronological
        for (int i = 1; i < events.size(); i++) {
            Assert.assertTrue(events.get(i).getT() >= events.get(i - 1).getT());
        }
    }

    /** check when g(t < root) < 0,  g(root + convergence) < 0. */
    @Test
    public void testRootPlusToleranceHasWrongSignAndLessThanTb() {
        // setup
        // test is fragile w.r.t. implementation and these parameters
        double maxCheck = 10;
        double tolerance = 0.5;
        double t1 = 11, t2 = 11.4, t3 = 12.0;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorB = new FlatDetector(events, t1, t2, t3);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        // allowed to find t1 or t3.
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorB, events.get(0).getHandler());
    }

    /**
     * Check when g(t) has a multiple root. e.g. g(t < root) < 0, g(root) = 0, g(t > root)
     * < 0.
     */
    @Test
    public void testDoubleRoot() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = 11;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t1);
        TimeDetector detectorB = new TimeDetector(events, t1, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        // detector worked correctly
        Assert.assertTrue(detectorB.g(state(t1)).getReal() == 0.0);
        Assert.assertTrue(detectorB.g(state(t1 - 1e-6)).getReal() < 0);
        Assert.assertTrue(detectorB.g(state(t1 + 1e-6)).getReal() < 0);
    }

    /**
     * Check when g(t) has a multiple root. e.g. g(t < root) > 0, g(root) = 0, g(t > root)
     * > 0.
     */
    @Test
    public void testDoubleRootOppositeSign() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = 11;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t1);
        TimeDetector detectorB = new ContinuousDetector(events, -20, t1, t1);
        detectorB.g(state(t1));
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        // detector worked correctly
        Assert.assertEquals(0.0, detectorB.g(state(t1)).getReal(), 0.0);
        Assert.assertTrue(detectorB.g(state(t1 - 1e-6)).getReal() > 0);
        Assert.assertTrue(detectorB.g(state(t1 + 1e-6)).getReal() > 0);
    }

    /** check root finding when zero at both ends. */
    @Test
    public void testZeroAtBeginningAndEndOfInterval() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = 10, t2 = 20;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ContinuousDetector(events, t1, t2);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t2, events.get(1).getT(), 0.0);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getHandler());
    }

    /** check root finding when zero at both ends. */
    @Test
    public void testZeroAtBeginningAndEndOfIntervalOppositeSign() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = 10, t2 = 20;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ContinuousDetector(events, -10, t1, t2);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(false, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t2, events.get(1).getT(), 0.0);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getHandler());
    }

    /** Test where an event detector has to back up multiple times. */
    @Test
    public void testMultipleBackups() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = 11.0, t2 = 12, t3 = 13, t4 = 14, t5 = 15, t6 = 16, t7 = 17;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ContinuousDetector(events, t6);
        TimeDetector detectorB = new ContinuousDetector(events, t1, t3, t4, t7);
        TimeDetector detectorC = new ContinuousDetector(events, t2, t5);

        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorC, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        Assert.assertEquals(5, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(detectorB, events.get(0).getHandler());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertEquals(detectorC, events.get(1).getHandler());
        // reporting t3 and t4 is optional, seeing them is not.
        // we know a root was found at t3 because events are reported at t2 and t5.
        /*
        Assert.assertEquals(t3, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertEquals(detectorB, events.get(2).getHandler());
        Assert.assertEquals(t4, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(detectorB, events.get(3).getHandler());
        */
        Assert.assertEquals(t5, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertEquals(detectorC, events.get(2).getHandler());
        Assert.assertEquals(t6, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(detectorA, events.get(3).getHandler());
        Assert.assertEquals(t7, events.get(4).getT(), tolerance);
        Assert.assertEquals(false, events.get(4).isIncreasing());
        Assert.assertEquals(detectorB, events.get(4).getHandler());
    }

    /** Test a reset event triggering another event at the same time. */
    @Test
    public void testEventCausedByStateReset() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = 15.0;
        final FieldODEState<Decimal64> newState = new FieldODEState<>(
                zero.add(t1), new Decimal64[]{zero.add(-20), zero});
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ResetDetector(events, newState, t1);
        BaseDetector detectorB = new StateDetector(events, -1);

        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(40.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t1, events.get(1).getT(), tolerance);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertEquals(detectorB, events.get(1).getHandler());
        Assert.assertEquals(t1 + 19, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertEquals(detectorB, events.get(2).getHandler());
    }

    /** check when t + tolerance == t. */
    @Test
    public void testConvergenceTooTight() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-18;
        double t1 = 15;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ContinuousDetector(events, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
    }

    /** check when root finding tolerance > event finding tolerance. */
    @Test
    public void testToleranceMismatch() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-18;
        double t1 = 15.1;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new FlatDetector(events, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100,
                new FieldBracketingNthOrderBrentSolver<>(zero, zero.add(1e-3), zero, 5));

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        // use root finder tolerance instead of event finder tolerance.
        Assert.assertEquals(t1, events.get(0).getT(), 1e-3);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
    }

    /**
     * test when one event detector changes the definition of another's g function before
     * the end of the step as a result of a continue action. Not sure if this should be
     * officially supported, but it is used in Orekit's DateDetector, it's useful, and not
     * too hard to implement.
     */
    @Test
    public void testEventChangesGFunctionDefinition() {
        // setup
        double maxCheck = 5;
        double tolerance = 1e-6;
        double t1 = 11, t2 = 19;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // mutable boolean
        boolean[] swap = new boolean[1];
        final TimeDetector detectorA = new ContinuousDetector(events, t1) {
            @Override
            public Action eventOccurred(FieldODEStateAndDerivative<Decimal64> state, boolean increasing) {
                swap[0] = true;
                return super.eventOccurred(state, increasing);
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(events, t2);
        BaseDetector detectorC = new BaseDetector(Action.CONTINUE, events) {

            @Override
            public Decimal64 g(FieldODEStateAndDerivative<Decimal64> state) {
                if (swap[0]) {
                    return detectorB.g(state);
                } else {
                    return zero.add(-1);
                }
            }

        };
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorC, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorC, events.get(1).getHandler());
    }

    /** check when root finding tolerance > event finding tolerance. */
    @Test
    public void testToleranceStop() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-18; // less than 1 ulp
        double t1 = 15.1;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new FlatDetector(Action.STOP, events, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        FieldODEStateAndDerivative<Decimal64> finalState =
                integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        // use root finder tolerance instead of event finder tolerance.
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t1, finalState.getTime().getReal(), tolerance);

        // try to resume propagation
        finalState = integrator.integrate(new Equation(), finalState, zero.add(30.0));

        // verify it got to the end
        Assert.assertEquals(30.0, finalState.getTime().getReal(), 0.0);
    }

    /**
     * Test when g function is initially zero for longer than the tolerance. Can occur
     * when restarting after a stop and cancellation occurs in the g function.
     */
    @Test
    public void testLongInitialZero() {
        // setup
        double maxCheck = 10;
        double tolerance = 1;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(Action.STOP, events, -50) {
            @Override
            public Decimal64 g(FieldODEStateAndDerivative<Decimal64> state) {
                if (state.getTime().getReal() < 2) {
                    return zero;
                } else {
                    return super.g(state);
                }
            }
        };
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(0, events.size());
    }

    /**
     * The root finder requires the start point to be in the interval (a, b) which is hard
     * when there aren't many numbers between a and b. This test uses a second event
     * detector to force a very small window for the first event detector.
     */
    @Test
    public void testShortBracketingInterval() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        final double t1 = FastMath.nextUp(10.0), t2 = 10.5;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // never zero so there is no easy way out
        TimeDetector detectorA = new TimeDetector(events) {
            @Override
            public Decimal64 g(FieldODEStateAndDerivative<Decimal64> state) {
                final Decimal64 t = state.getTime();
                if (t.getReal() < t1) {
                    return one.negate();
                } else if (t.getReal() < t2) {
                    return one;
                } else {
                    return one.negate();
                }
            }
        };
        TimeDetector detectorB = new TimeDetector(events, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t1, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorB, events.get(1).getHandler());
        Assert.assertEquals(t2, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertSame(detectorA, events.get(2).getHandler());
    }



    /* The following tests are copies of the above tests, except that they propagate in
     * the reverse direction and all the signs on the time values are negated.
     */


    @Test
    public void testCloseEventsFirstOneIsResetReverse() {
        // setup
        // a fairly rare state to reproduce this bug. Two dates, d1 < d2, that
        // are very close. Event triggers on d1 will reset state to break out of
        // event handling loop in AbstractIntegrator.acceptStep(). At this point
        // detector2 has g0Positive == true but the event time is set to just
        // before the event so g(t0) is negative. Now on processing the
        // next step the root solver checks the sign of the start, midpoint,
        // and end of the interval so we need another event less than half a max
        // check interval after d2 so that the g function will be negative at
        // all three times. Then we get a non bracketing exception.
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);

        // switched for 9 to 1 to be close to the start of the step
        double t1 = -1;
        TimeDetector detector1 = new TimeDetector(Action.RESET_DERIVATIVES, t1);
        integrator.addEventHandler(detector1, 10, 1e-9, 100);
        TimeDetector detector2 = new TimeDetector(t1 - 1e-15, t1 - 4.9);
        integrator.addEventHandler(detector2, 11, 1e-9, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-20));

        // verify
        List<Event> events1 = detector1.getEvents();
        Assert.assertEquals(1, events1.size());
        Assert.assertEquals(t1, events1.get(0).getT(), 0.0);
        List<Event> events2 = detector2.getEvents();
        Assert.assertEquals(0, events2.size());
    }

    @Test
    public void testCloseEventsReverse() {
        // setup
        double e = 1e-15;
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, e, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(-5);
        integrator.addEventHandler(detector1, 10, 1, 100);
        TimeDetector detector2 = new TimeDetector(-5.5);
        integrator.addEventHandler(detector2, 10, 1, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-20));

        // verify
        List<Event> events1 = detector1.getEvents();
        Assert.assertEquals(1, events1.size());
        Assert.assertEquals(-5, events1.get(0).getT(), 0.0);
        List<Event> events2 = detector2.getEvents();
        Assert.assertEquals(1, events2.size());
        Assert.assertEquals(-5.5, events2.get(0).getT(), 0.0);
    }

    @Test
    public void testSimultaneousEventsReverse() {
        // setup
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(-5);
        integrator.addEventHandler(detector1, 10, 1, 100);
        TimeDetector detector2 = new TimeDetector(-5);
        integrator.addEventHandler(detector2, 10, 1, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-20));

        // verify
        List<Event> events1 = detector1.getEvents();
        Assert.assertEquals(1, events1.size());
        Assert.assertEquals(-5, events1.get(0).getT(), 0.0);
        List<Event> events2 = detector2.getEvents();
        Assert.assertEquals(1, events2.size());
        Assert.assertEquals(-5, events2.get(0).getT(), 0.0);
    }

    /**
     * test the g function switching with a period shorter than the tolerance. We don't
     * need to find any of the events, but we do need to not crash. And we need to
     * preserve the alternating increasing / decreasing sequence.
     */
    @Test
    public void testFastSwitchingReverse() {
        // setup
        // step size of 10 to land in between two events we would otherwise miss
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(-9.9, -10.1, -12);
        integrator.addEventHandler(detector1, 10, 0.2, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-20));

        //verify
        // finds one or three events. Not 2.
        List<Event> events1 = detector1.getEvents();
        Assert.assertEquals(1, events1.size());
        Assert.assertEquals(-9.9, events1.get(0).getT(), 0.2);
        Assert.assertEquals(true, events1.get(0).isIncreasing());
    }

    /** "A Tricky Problem" from bug #239. */
    @Test
    public void testTrickyCaseLowerReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = -1.0, t2 = -15, t3 = -16, t4 = -17, t5 = -18;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t3);
        TimeDetector detectorB = new TimeDetector(events, -50, t1, t2, t5);
        TimeDetector detectorC = new TimeDetector(Action.RESET_DERIVATIVES, events, t4);

        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorC, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        // but I only know one way to do that in this case.
        Assert.assertEquals(5, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(false, events.get(0).isIncreasing());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertEquals(t3, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertEquals(t4, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(t5, events.get(4).getT(), tolerance);
        Assert.assertEquals(false, events.get(4).isIncreasing());
    }

    /**
     * Test case for two event detectors. DetectorA has event at t2, DetectorB at t3, but
     * due to the root finding tolerance DetectorB's event occurs at t1. With t1 < t2 <
     * t3.
     */
    @Test
    public void testRootFindingToleranceReverse() {
        //setup
        double maxCheck = 10;
        double t2 = -11, t3 = t2 - 1e-5;
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t2);
        FlatDetector detectorB = new FlatDetector(events, t3);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, 1e-6, 100);
        integrator.addEventHandler(detectorB, maxCheck, 0.5, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        // if these fail the event finding did its job,
        // but this test isn't testing what it is supposed to be
        Assert.assertSame(detectorB, events.get(0).getHandler());
        Assert.assertSame(detectorA, events.get(1).getHandler());
        Assert.assertTrue(events.get(0).getT() > events.get(1).getT());

        // check event detection worked
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t3, events.get(0).getT(), 0.5);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(t2, events.get(1).getT(), 1e-6);
        Assert.assertEquals(true, events.get(1).isIncreasing());
    }

    /** check when g(t < root) < 0,  g(root + convergence) < 0. */
    @Test
    public void testRootPlusToleranceHasWrongSignReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        final double toleranceB = 0.3;
        double t1 = -11, t2 = -11.1, t3 = -11.2;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t2);
        TimeDetector detectorB = new TimeDetector(events, -50, t1, t3);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, toleranceB, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        // we only care that the rules are satisfied. There are multiple solutions.
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), toleranceB);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorB, events.get(0).getHandler());
        Assert.assertEquals(t3, events.get(1).getT(), toleranceB);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertSame(detectorB, events.get(1).getHandler());
        Assert.assertEquals(t2, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertSame(detectorA, events.get(2).getHandler());
        // ascending order
        Assert.assertTrue(events.get(0).getT() >= events.get(1).getT());
        Assert.assertTrue(events.get(1).getT() >= events.get(2).getT());
    }

    /** check when g(t < root) < 0,  g(root + convergence) < 0. */
    @Test
    public void testRootPlusToleranceHasWrongSignAndLessThanTbReverse() {
        // setup
        // test is fragile w.r.t. implementation and these parameters
        double maxCheck = 10;
        double tolerance = 0.5;
        double t1 = -11, t2 = -11.4, t3 = -12.0;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorB = new FlatDetector(events, t1, t2, t3);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        // allowed to report t1 or t3.
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorB, events.get(0).getHandler());
    }

    /**
     * Check when g(t) has a multiple root. e.g. g(t < root) < 0, g(root) = 0, g(t > root)
     * < 0.
     */
    @Test
    public void testDoubleRootReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = -11;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t1);
        TimeDetector detectorB = new TimeDetector(events, t1, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        // detector worked correctly
        Assert.assertTrue(detectorB.g(state(t1)).getReal() == 0.0);
        Assert.assertTrue(detectorB.g(state(t1 + 1e-6)).getReal() < 0);
        Assert.assertTrue(detectorB.g(state(t1 - 1e-6)).getReal() < 0);
    }

    /**
     * Check when g(t) has a multiple root. e.g. g(t < root) > 0, g(root) = 0, g(t > root)
     * > 0.
     */
    @Test
    public void testDoubleRootOppositeSignReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = -11;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(events, t1);
        TimeDetector detectorB = new ContinuousDetector(events, -50, t1, t1);
        detectorB.g(state(t1));
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        // detector worked correctly
        Assert.assertEquals(0.0, detectorB.g(state(t1)).getReal(), 0.0);
        Assert.assertTrue(detectorB.g(state(t1 + 1e-6)).getReal() > 0);
        Assert.assertTrue(detectorB.g(state(t1 - 1e-6)).getReal() > 0);
    }

    /** check root finding when zero at both ends. */
    @Test
    public void testZeroAtBeginningAndEndOfIntervalReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = -10, t2 = -20;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ContinuousDetector(events, -50, t1, t2);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t2, events.get(1).getT(), 0.0);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getHandler());
    }

    /** check root finding when zero at both ends. */
    @Test
    public void testZeroAtBeginningAndEndOfIntervalOppositeSignReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = -10, t2 = -20;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ContinuousDetector(events, t1, t2);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(false, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t2, events.get(1).getT(), 0.0);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getHandler());
    }

    /** Test where an event detector has to back up multiple times. */
    @Test
    public void testMultipleBackupsReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = -11.0, t2 = -12, t3 = -13, t4 = -14, t5 = -15, t6 = -16, t7 = -17;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ContinuousDetector(events, t6);
        TimeDetector detectorB = new ContinuousDetector(events, -50, t1, t3, t4, t7);
        TimeDetector detectorC = new ContinuousDetector(events, -50, t2, t5);

        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorC, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        Assert.assertEquals(5, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(detectorB, events.get(0).getHandler());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertEquals(detectorC, events.get(1).getHandler());
        // reporting t3 and t4 is optional, seeing them is not.
        // we know a root was found at t3 because events are reported at t2 and t5.
        /*
        Assert.assertEquals(t3, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertEquals(detectorB, events.get(2).getHandler());
        Assert.assertEquals(t4, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(detectorB, events.get(3).getHandler());
        */
        Assert.assertEquals(t5, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertEquals(detectorC, events.get(2).getHandler());
        Assert.assertEquals(t6, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(detectorA, events.get(3).getHandler());
        Assert.assertEquals(t7, events.get(4).getT(), tolerance);
        Assert.assertEquals(false, events.get(4).isIncreasing());
        Assert.assertEquals(detectorB, events.get(4).getHandler());
    }

    /** Test a reset event triggering another event at the same time. */
    @Test
    public void testEventCausedByStateResetReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = -15.0;
        final FieldODEState<Decimal64> newState =
                new FieldODEState<>(zero.add(t1), new Decimal64[]{zero.add(20), zero});
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ResetDetector(events, newState, t1);
        BaseDetector detectorB = new StateDetector(events, 1);

        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-40.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t1, events.get(1).getT(), tolerance);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertEquals(detectorB, events.get(1).getHandler());
        Assert.assertEquals(t1 - 19, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertEquals(detectorB, events.get(2).getHandler());
    }

    /** check when t + tolerance == t. */
    @Test
    public void testConvergenceTooTightReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-18;
        double t1 = -15;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ContinuousDetector(events, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
    }

    /** check when root finding tolerance > event finding tolerance. */
    @Test
    public void testToleranceMismatchReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-18;
        double t1 = -15.1;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new FlatDetector(events, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100,
                new FieldBracketingNthOrderBrentSolver<>(zero, zero.add(1e-3), zero, 5));

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        // use root finding tolerance since it is larger
        Assert.assertEquals(t1, events.get(0).getT(), 1e-3);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
    }

    /**
     * test when one event detector changes the definition of another's g function before
     * the end of the step as a result of a continue action. Not sure if this should be
     * officially supported, but it is used in Orekit's DateDetector, it's useful, and not
     * too hard to implement.
     */
    @Test
    public void testEventChangesGFunctionDefinitionReverse() {
        // setup
        double maxCheck = 5;
        double tolerance = 1e-6;
        double t1 = -11, t2 = -19;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // mutable boolean
        boolean[] swap = new boolean[1];
        final TimeDetector detectorA = new ContinuousDetector(events, t1) {
            @Override
            public Action eventOccurred(FieldODEStateAndDerivative<Decimal64> state, boolean increasing) {
                swap[0] = true;
                return super.eventOccurred(state, increasing);
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(events, t2);
        BaseDetector detectorC = new BaseDetector(Action.CONTINUE, events) {

            @Override
            public Decimal64 g(FieldODEStateAndDerivative<Decimal64> state) {
                if (swap[0]) {
                    return detectorB.g(state);
                } else {
                    return one;
                }
            }

        };
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorC, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorC, events.get(1).getHandler());
    }

    /** check when root finding tolerance > event finding tolerance. */
    @Test
    public void testToleranceStopReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-18; // less than 1 ulp
        double t1 = -15.1;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new FlatDetector(Action.STOP, events, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        FieldODEStateAndDerivative<Decimal64> finalState =
                integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        // use root finder tolerance instead of event finder tolerance.
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t1, finalState.getTime().getReal(), tolerance);

        // try to resume propagation
        finalState = integrator.integrate(new Equation(), finalState, zero.add(-30.0));

        // verify it got to the end
        Assert.assertEquals(-30.0, finalState.getTime().getReal(), 0.0);
    }

    /**
     * Test when g function is initially zero for longer than the tolerance. Can occur
     * when restarting after a stop and cancellation occurs in the g function.
     */
    @Test
    public void testLongInitialZeroReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new TimeDetector(Action.STOP, events, 50) {
            @Override
            public Decimal64 g(FieldODEStateAndDerivative<Decimal64> state) {
                if (state.getTime().getReal() > -2) {
                    return zero;
                } else {
                    return super.g(state);
                }
            }
        };
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(0, events.size());
    }

    /**
     * The root finder requires the start point to be in the interval (a, b) which is hard
     * when there aren't many numbers between a and b. This test uses a second event
     * detector to force a very small window for the first event detector.
     */
    @Test
    public void testShortBracketingIntervalReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        final double t1 = FastMath.nextDown(-10.0), t2 = -10.5;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // never zero so there is no easy way out
        TimeDetector detectorA = new TimeDetector(events) {
            @Override
            public Decimal64 g(FieldODEStateAndDerivative<Decimal64> state) {
                final Decimal64 t = state.getTime();
                if (t.getReal() > t1) {
                    return one.negate();
                } else if (t.getReal() > t2) {
                    return one;
                } else {
                    return one.negate();
                }
            }
        };
        TimeDetector detectorB = new TimeDetector(events, t1);
        FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventHandler(detectorA, maxCheck, tolerance, 100);
        integrator.addEventHandler(detectorB, maxCheck, tolerance, 100);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(false, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getHandler());
        Assert.assertEquals(t1, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorB, events.get(1).getHandler());
        Assert.assertEquals(t2, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertSame(detectorA, events.get(2).getHandler());
    }



    /* utility classes and methods */

    /**
     * Create a state at a time.
     *
     * @param t time of state.
     * @return new state.
     */
    private FieldODEStateAndDerivative<Decimal64> state(double t) {
        return new FieldODEStateAndDerivative<>(
                zero.add(t), new Decimal64[0], new Decimal64[0]);
    }

    /** Base class to record events that occured. */
    private static abstract class BaseDetector implements FieldODEEventHandler<Decimal64> {

        protected final Action action;

        /** times the event was actually triggered. */
        private final List<Event> events;

        public BaseDetector(Action action, List<Event> events) {
            this.action = action;
            this.events = events;
        }

        /**
         * Get all events that occurred.
         *
         * @return list of events in the order they occurred.
         */
        public List<Event> getEvents() {
            return events;
        }

        @Override
        public Action eventOccurred(FieldODEStateAndDerivative<Decimal64> state,
                                    boolean increasing) {
            events.add(new Event(state, increasing, this));
            return this.action;
        }

    }

    /** Trigger an event at a particular time. */
    private static class TimeDetector extends BaseDetector {

        /** time of the event to trigger. */
        protected final double[] eventTs;

        /**
         * Create a new detector.
         *
         * @param eventTs the time to trigger an event.
         */
        public TimeDetector(double... eventTs) {
            this(Action.CONTINUE, eventTs);
        }

        public TimeDetector(List<Event> events, double... eventTs) {
            this(Action.CONTINUE, events, eventTs);
        }

        public TimeDetector(Action action, double... eventTs) {
            this(action, new ArrayList<Event>(), eventTs);
        }

        public TimeDetector(Action action, List<Event> events, double... eventTs) {
            super(action, events);
            this.eventTs = eventTs.clone();
            Arrays.sort(this.eventTs);
        }

        @Override
        public Decimal64 g(final FieldODEStateAndDerivative<Decimal64> state) {
            final Decimal64 t = state.getTime();
            int i = 0;
            while (i < eventTs.length && t.getReal() > eventTs[i]) {
                i++;
            }
            i--;
            if (i < 0) {
                return t.subtract(eventTs[0]);
            } else {
                int sign = (i % 2) * 2 - 1;
                return t.subtract(eventTs[i]).multiply(-sign);
            }
        }

    }

    private static class Event {

        private final FieldODEStateAndDerivative<Decimal64> state;
        private final boolean increasing;
        private final FieldODEEventHandler<Decimal64> handler;

        public Event(FieldODEStateAndDerivative<Decimal64> state,
                     boolean increasing,
                     FieldODEEventHandler<Decimal64> handler) {
            this.increasing = increasing;
            this.state = state;
            this.handler = handler;
        }

        public boolean isIncreasing() {
            return increasing;
        }

        public double getT() {
            return state.getTime().getReal();
        }

        public FieldODEEventHandler<Decimal64> getHandler() {
            return handler;
        }

        @Override
        public String toString() {
            return "Event{" +
                    "increasing=" + increasing +
                    ", time=" + state.getTime() +
                    ", state=" + Arrays.toString(state.getCompleteState()) +
                    '}';
        }
    }

    /**
     * Same as {@link TimeDetector} except that it has a very flat g function which makes
     * root finding hard.
     */
    private static class FlatDetector extends TimeDetector {

        public FlatDetector(List<Event> events, double... eventTs) {
            super(events, eventTs);
        }

        public FlatDetector(Action action, List<Event> events, double... eventTs) {
            super(action, events, eventTs);
        }

        @Override
        public Decimal64 g(FieldODEStateAndDerivative<Decimal64> state) {
            final Decimal64 g = super.g(state);
            return g.signum();
        }

    }

    /** Linear on both ends, parabolic in the middle. */
    private static class ContinuousDetector extends TimeDetector {

        public ContinuousDetector(List<Event> events, double... eventTs) {
            super(events, eventTs);
        }

        public ContinuousDetector(Action action, List<Event> events, double... eventTs) {
            super(action, events, eventTs);
        }

        @Override
        public Decimal64 g(FieldODEStateAndDerivative<Decimal64> state) {
            final Decimal64 t = state.getTime();
            int i = 0;
            while (i < eventTs.length && t.getReal() > eventTs[i]) {
                i++;
            }
            i--;
            if (i < 0) {
                return t.subtract(eventTs[0]);
            } else if (i < eventTs.length - 1) {
                int sign = (i % 2) * 2 - 1;
                return t.subtract(eventTs[i]).multiply(t.negate().add(eventTs[i + 1])).multiply(-sign);
            } else {
                int sign = (i % 2) * 2 - 1;
                return t.subtract(eventTs[i]).multiply(-sign);
            }
        }
    }

    /** Reset the state at a particular time. */
    private static class ResetDetector extends TimeDetector {

        private final FieldODEState<Decimal64> resetState;

        public ResetDetector(List<Event> events, FieldODEState<Decimal64> state, double eventT) {
            super(Action.RESET_STATE, events, eventT);
            this.resetState = state;
        }

        @Override
        public FieldODEState<Decimal64> resetState(FieldODEStateAndDerivative<Decimal64> state) {
            return resetState;
        }

    }

    /** Switching function based on the first primary state. */
    private static class StateDetector extends BaseDetector {

        private final double triggerState;

        public StateDetector(List<Event> events, double triggerState) {
            super(Action.CONTINUE, events);
            this.triggerState = triggerState;
        }

        @Override
        public Decimal64 g(FieldODEStateAndDerivative<Decimal64> state) {
            return state.getPrimaryState()[0].subtract(this.triggerState);
        }
    }

    /** Some basic equations to integrate. */
    public static class Equation extends FieldExpandableODE<Decimal64> {
        public Equation() {
            super(new EquationODE());
        }
    }

    /** Some basic equations to integrate. */
    public static class EquationODE implements FieldOrdinaryDifferentialEquation<Decimal64> {

        public int getDimension() {
            return 2;
        }

        @Override
        public Decimal64[] computeDerivatives(Decimal64 t, Decimal64[] y) {
            return new Decimal64[]{one, one.multiply(2)};
        }

    }

}
