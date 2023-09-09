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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.hipparchus.Field;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853FieldIntegrator;
import org.hipparchus.ode.nonstiff.LutherFieldIntegrator;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.ode.sampling.FieldODEStepHandler;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
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
    private static final Field<Binary64> field = Binary64Field.getInstance();
    private static final Binary64 zero = field.getZero();
    private static final Binary64 one = field.getOne();
    private static final FieldODEState<Binary64> initialState =
            new FieldODEState<>(zero, new Binary64[]{zero, zero});

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
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(10, 1e-9, 100, Action.RESET_DERIVATIVES, 9);
        integrator.addEventDetector(detector1);
        TimeDetector detector2 = new TimeDetector(11, 1e-9, 100, 9 + 1e-15, 9 + 4.9);
        integrator.addEventDetector(detector2);

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
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, e, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(10, 1, 100, 5);
        integrator.addEventDetector(detector1);
        TimeDetector detector2 = new TimeDetector(10, 1, 100, 5.5);
        integrator.addEventDetector(detector2);

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
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(10, 1, 100, 5);
        integrator.addEventDetector(detector1);
        TimeDetector detector2 = new TimeDetector(10, 1, 100, 5);
        integrator.addEventDetector(detector2);

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
     * Previously there were some branches when tryAdvance() returned false but did not
     * set {@code t0 = t}. This allowed the order of events to not be chronological and to
     * detect events that should not have occurred, both of which are problems.
     */
    @Test
    public void testSimultaneousEventsResetReverse() {
        // setup
        double tol = 1e-10;
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);
        boolean[] firstEventOccurred = {false};
        List<Event> events = new ArrayList<>();

        TimeDetector detector1 = new TimeDetector(10, tol, 100, events, -5) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                firstEventOccurred[0] = true;
                super.getHandler().eventOccurred(state, detector, increasing);
                return Action.RESET_STATE;
                };
            }
        };
        integrator.addEventDetector(detector1);
        // this detector changes it's g function definition when detector1 fires
        TimeDetector detector2 = new TimeDetector(1, tol, 100, events, -1, -3, -5) {
            @Override
            public Binary64 g(final FieldODEStateAndDerivative<Binary64> state) {
                if (firstEventOccurred[0]) {
                    return super.g(state);
                }
                return new TimeDetector(10, tol, 100, -5).g(state);
            }
        };
        integrator.addEventDetector(detector2);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-20));

        // verify
        // order is important to make sure the test checks what it is supposed to
        Assert.assertEquals(-5, events.get(0).getT(), 0.0);
        Assert.assertTrue(events.get(0).isIncreasing());
        Assert.assertEquals(detector1, events.get(0).getDetector());
        Assert.assertEquals(-5, events.get(1).getT(), 0.0);
        Assert.assertTrue(events.get(1).isIncreasing());
        Assert.assertEquals(detector2, events.get(1).getDetector());
        Assert.assertEquals(2, events.size());
    }

    /**
     * When two event detectors have a discontinuous event caused by a {@link
     * Action#RESET_STATE} or {@link Action#RESET_DERIVATIVES}. The two event detectors
     * would each say they had an event that had to be handled before the other one, but
     * neither would actually back up at all. For Hipparchus GitHub #91.
     */
    @Test
    public void testSimultaneousDiscontinuousEventsAfterReset() {
        // setup
        double t = FastMath.PI;
        double tol = 1e-10;
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        List<Event> events = new ArrayList<>();

        TimeDetector resetDetector =
                new ResetDetector(10, tol, 100, events, new FieldODEState<>(zero.add(t), new Binary64[]{zero.add(1e100), zero}), t);
        integrator.addEventDetector(resetDetector);
        List<BaseDetector> detectors = new ArrayList<>();
        for (int i = 0; i < 2; i++) {
            BaseDetector detector1 = new StateDetector(10, tol, 100, events, 0.0);
            integrator.addEventDetector(detector1);
            detectors.add(detector1);
        }

        // action
        integrator.integrate(new Equation(), new FieldODEState<>(zero, new Binary64[]{zero.add(-1e100), zero}), zero.add(10));

        // verify
        Assert.assertEquals(t, events.get(0).getT(), tol);
        Assert.assertTrue(events.get(0).isIncreasing());
        Assert.assertEquals(resetDetector, events.get(0).getDetector());
        // next two events can occur in either order
        Assert.assertEquals(t, events.get(1).getT(), tol);
        Assert.assertTrue(events.get(1).isIncreasing());
        Assert.assertEquals(detectors.get(0), events.get(1).getDetector());
        Assert.assertEquals(t, events.get(2).getT(), tol);
        Assert.assertTrue(events.get(2).isIncreasing());
        Assert.assertEquals(detectors.get(1), events.get(2).getDetector());
        Assert.assertEquals(events.size(), 3);
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
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(10, 0.2, 100, 9.9, 10.1, 12);
        integrator.addEventDetector(detector1);

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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, events, t3);
        TimeDetector detectorB = new TimeDetector(maxCheck, tolerance, 100, events, -10, t1, t2, t5);
        TimeDetector detectorC = new TimeDetector(maxCheck, tolerance, 100, Action.RESET_DERIVATIVES, events, t4);

        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);
        integrator.addEventDetector(detectorC);

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
        TimeDetector detectorA = new TimeDetector(maxCheck, 1e-6, 100, events, t2);
        TimeDetector detectorB = new FlatDetector(maxCheck, 0.5, 100, events, t3);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        // if these fail the event finding did its job,
        // but this test isn't testing what it is supposed to be
        Assert.assertSame(detectorB, events.get(0).getDetector());
        Assert.assertSame(detectorA, events.get(1).getDetector());
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
        TimeDetector detectorA = new TimeDetector(maxCheck, 1e-6, 100, events, t2);
        TimeDetector detectorB = new TimeDetector(maxCheck, toleranceB, 100, events, t1, t3);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        // we only care that the rules are satisfied, there are other solutions
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), toleranceB);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorB, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getDetector());
        Assert.assertEquals(t3, events.get(2).getT(), toleranceB);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertSame(detectorB, events.get(2).getDetector());
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
        TimeDetector detectorB = new FlatDetector(maxCheck, tolerance, 100, events, t1, t2, t3);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        // allowed to find t1 or t3.
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorB, events.get(0).getDetector());
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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, events, t1);
        TimeDetector detectorB = new TimeDetector(maxCheck, tolerance, 100, events, t1, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, events, t1);
        TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, -20, t1, t1);
        detectorB.g(state(t1));
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, t1, t2);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), 0.0);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getDetector());
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
        TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, -10, t1, t2);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(false, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), 0.0);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getDetector());
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
        TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, t6);
        TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, t1, t3, t4, t7);
        TimeDetector detectorC = new ContinuousDetector(maxCheck, tolerance, 100, events, t2, t5);

        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        Assert.assertEquals(5, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(detectorB, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertEquals(detectorC, events.get(1).getDetector());
        // reporting t3 and t4 is optional, seeing them is not.
        // we know a root was found at t3 because events are reported at t2 and t5.
        /*
        Assert.assertEquals(t3, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertEquals(detectorB, events.get(2).getDetector());
        Assert.assertEquals(t4, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(detectorB, events.get(3).getDetector());
        */
        Assert.assertEquals(t5, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertEquals(detectorC, events.get(2).getDetector());
        Assert.assertEquals(t6, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(detectorA, events.get(3).getDetector());
        Assert.assertEquals(t7, events.get(4).getT(), tolerance);
        Assert.assertEquals(false, events.get(4).isIncreasing());
        Assert.assertEquals(detectorB, events.get(4).getDetector());
    }

    /** Test a reset event triggering another event at the same time. */
    @Test
    public void testEventCausedByStateReset() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = 15.0;
        final FieldODEState<Binary64> newState = new FieldODEState<>(
                zero.add(t1), new Binary64[]{zero.add(-20), zero});
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ResetDetector(maxCheck, tolerance, 100, events, newState, t1);
        BaseDetector detectorB = new StateDetector(maxCheck, tolerance, 100, events, -1);

        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(40.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t1, events.get(1).getT(), tolerance);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertEquals(detectorB, events.get(1).getDetector());
        Assert.assertEquals(t1 + 19, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertEquals(detectorB, events.get(2).getDetector());
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
        TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        TimeDetector detectorA = new FlatDetector(maxCheck, tolerance, 100, events, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        // use root finder tolerance instead of event finder tolerance.
        Assert.assertEquals(t1, events.get(0).getT(), 1e-3);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        final TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, t1) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                    swap[0] = true;
                    return super.getHandler().eventOccurred(state, detector, increasing);
                };
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, t2);
        BaseDetector detectorC = new BaseDetector(maxCheck, tolerance, 100, Action.CONTINUE, events) {

            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (swap[0]) {
                    return detectorB.g(state);
                } else {
                    return zero.add(-1);
                }
            }

        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorC, events.get(1).getDetector());
    }


    /**
     * test when one event detector changes the definition of another's g function before
     * the end of the step as a result of an event occurring. In this case the change
     * cancels the occurrence of the event.
     */
    @Test
    public void testEventChangesGFunctionDefinitionCancel() {
        // setup
        double maxCheck = 5;
        double tolerance = 1e-6;
        double t1 = 11, t2 = 11.1;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // mutable boolean
        boolean[] swap = new boolean[1];
        final TimeDetector detectorA =
                        new ContinuousDetector(maxCheck, tolerance, 100, Action.RESET_EVENTS, events, t1) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                    swap[0] = true;
                    return super.getHandler().eventOccurred(state, detector, increasing);
                };
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, t2);
        BaseDetector detectorC = new BaseDetector(maxCheck, tolerance, 100, Action.CONTINUE, events) {

            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (!swap[0]) {
                    return detectorB.g(state);
                } else {
                    return zero.add(-1);
                }
            }

        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
    }

    /**
     * test when one event detector changes the definition of another's g function before
     * the end of the step as a result of an event occurring. In this case the change
     * delays the occurrence of the event.
     */
    @Test
    public void testEventChangesGFunctionDefinitionDelay() {
        // setup
        double maxCheck = 5;
        double tolerance = 1e-6;
        double t1 = 11, t2 = 11.1, t3 = 11.2;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // mutable boolean
        boolean[] swap = new boolean[1];
        final TimeDetector detectorA =
                        new ContinuousDetector(maxCheck, tolerance, 100, Action.RESET_EVENTS, events, t1) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                    swap[0] = true;
                    return super.getHandler().eventOccurred(state, detector, increasing);
                };
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, t2);
        final TimeDetector detectorD = new ContinuousDetector(maxCheck, tolerance, 100, events, t3);
        BaseDetector detectorC = new BaseDetector(maxCheck, tolerance, 100, Action.CONTINUE, events) {

            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (!swap[0]) {
                    return detectorB.g(state);
                } else {
                    return detectorD.g(state);
                }
            }

        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t3, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorC, events.get(1).getDetector());
    }

    /**
     * test when one event detector changes the definition of another's g function before
     * the end of the step as a result of an event occurring. In this case the change
     * causes the event to happen sooner than originally expected.
     */
    @Test
    public void testEventChangesGFunctionDefinitionAccelerate() {
        // setup
        double maxCheck = 5;
        double tolerance = 1e-6;
        double t1 = 11, t2 = 11.1, t3 = 11.2;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // mutable boolean
        boolean[] swap = new boolean[1];
        final TimeDetector detectorA =
                        new ContinuousDetector(maxCheck, tolerance, 100, Action.RESET_EVENTS, events, t1) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                    swap[0] = true;
                    return super.getHandler().eventOccurred(state, detector, increasing);
                };
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, t2);
        final TimeDetector detectorD = new ContinuousDetector(maxCheck, tolerance, 100, events, t3);
        BaseDetector detectorC = new BaseDetector(maxCheck, tolerance, 100, Action.CONTINUE, events) {

            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (swap[0]) {
                    return detectorB.g(state);
                } else {
                    return detectorD.g(state);
                }
            }

        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorC, events.get(1).getDetector());
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
        TimeDetector detectorA = new FlatDetector(maxCheck, tolerance, 100, Action.STOP, events, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        FieldODEStateAndDerivative<Binary64> finalState =
                integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(1, events.size());
        // use root finder tolerance instead of event finder tolerance.
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, Action.STOP, events, -50) {
            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (state.getTime().getReal() < 2) {
                    return zero;
                } else {
                    return super.g(state);
                }
            }
        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, events) {
            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                final Binary64 t = state.getTime();
                if (t.getReal() < t1) {
                    return one.negate();
                } else if (t.getReal() < t2) {
                    return one;
                } else {
                    return one.negate();
                }
            }
        };
        TimeDetector detectorB = new TimeDetector(maxCheck, tolerance, 100, events, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(30.0));

        // verify
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t1, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorB, events.get(1).getDetector());
        Assert.assertEquals(t2, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertSame(detectorA, events.get(2).getDetector());
    }

    /** Check that steps are restricted correctly with a continue event. */
    @Test
    public void testEventStepHandler() {
        // setup
        double tolerance = 1e-18;
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(new TimeDetector(100.0, tolerance, 100, 5));
        StepHandler stepHandler = new StepHandler();
        integrator.addStepHandler(stepHandler);

        // action
        FieldODEStateAndDerivative<Binary64> finalState = integrator
                .integrate(new Equation(), initialState, zero.add(10));

        // verify
        Assert.assertEquals(10.0, finalState.getTime().getReal(), tolerance);
        Assert.assertEquals(0.0,
                stepHandler.initialState.getTime().getReal(), tolerance);
        Assert.assertEquals(10.0, stepHandler.finalTime.getReal(), tolerance);
        Assert.assertEquals(10.0,
                stepHandler.finalState.getTime().getReal(), tolerance);
        FieldODEStateInterpolator<Binary64> interpolator = stepHandler.interpolators.get(0);
        Assert.assertEquals(0.0,
                interpolator.getPreviousState().getTime().getReal(), tolerance);
        Assert.assertEquals(5.0,
                interpolator.getCurrentState().getTime().getReal(), tolerance);
        interpolator = stepHandler.interpolators.get(1);
        Assert.assertEquals(5.0,
                interpolator.getPreviousState().getTime().getReal(), tolerance);
        Assert.assertEquals(10.0,
                interpolator.getCurrentState().getTime().getReal(), tolerance);
        Assert.assertEquals(2, stepHandler.interpolators.size());
    }

    /** Test resetState(...) returns {@code null}. */
    @Test
    public void testEventCausedByDerivativesReset() {
        // setup
        TimeDetector detectorA = new TimeDetector(10, 1e-6, 100, Action.RESET_STATE, 15.0) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return new FieldODEEventHandler<Binary64>() {
                    @Override
                    public Action eventOccurred(FieldODEStateAndDerivative<Binary64> state,
                                                FieldODEEventDetector<Binary64> detector,
                                                boolean increasing) {
                        return Action.RESET_STATE;
                    }
                    @Override
                    public FieldODEState<Binary64> resetState(FieldODEEventDetector<Binary64> detector,
                                                               FieldODEStateAndDerivative<Binary64> state) {
                        return null;
                    }
                };
            }
        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        try {
            // action
            integrator.integrate(new Equation(), initialState, zero.add(20.0));
            Assert.fail("Expected Exception");
        } catch (NullPointerException e) {
            // expected
        }
    }

    @Test
    public void testResetChangesSign() {
        FieldOrdinaryDifferentialEquation<Binary64> equation = new FieldOrdinaryDifferentialEquation<Binary64>() {
            public int getDimension() { return 1; }
            public Binary64[] computeDerivatives(Binary64 t, Binary64[] y) { return new Binary64[] { new Binary64(1.0) }; }
        };

        LutherFieldIntegrator<Binary64> integrator = new LutherFieldIntegrator<>(Binary64Field.getInstance(), new Binary64(20.0));
        final double small = 1.0e-10;
        ResetChangesSignGenerator eventsGenerator = new ResetChangesSignGenerator(6.0, 9.0, -0.5 * small, 8.0, small, 1000);
        integrator.addEventDetector(eventsGenerator);
        final FieldODEStateAndDerivative<Binary64> end = integrator.integrate(new FieldExpandableODE<>(equation),
                new FieldODEState<>(new Binary64(0.0),
                        new Binary64[] { new Binary64(0.0) }),
                new Binary64(100.0));
        Assert.assertEquals(2,                 eventsGenerator.getCount());
        Assert.assertEquals(9.0,               end.getCompleteState()[0].getReal(), 1.0e-12);
        Assert.assertEquals(9.0 + 0.5 * small, end.getTime().getReal(),             1.0e-12);
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
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);

        // switched for 9 to 1 to be close to the start of the step
        double t1 = -1;
        TimeDetector detector1 = new TimeDetector(10, 1e-9, 100, Action.RESET_DERIVATIVES, t1);
        integrator.addEventDetector(detector1);
        TimeDetector detector2 = new TimeDetector(11, 1e-9, 100, t1 - 1e-15, t1 - 4.9);
        integrator.addEventDetector(detector2);

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
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, e, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(10, 1, 100, -5);
        integrator.addEventDetector(detector1);
        TimeDetector detector2 = new TimeDetector(10, 1, 100, -5.5);
        integrator.addEventDetector(detector2);

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
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(10, 1, 100, -5);
        integrator.addEventDetector(detector1);
        TimeDetector detector2 = new TimeDetector(10, 1, 100, -5);
        integrator.addEventDetector(detector2);

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
     * Previously there were some branches when tryAdvance() returned false but did not
     * set {@code t0 = t}. This allowed the order of events to not be chronological and to
     * detect events that should not have occurred, both of which are problems.
     */
    @Test
    public void testSimultaneousEventsReset() {
        // setup
        double tol = 1e-10;
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 100.0, 1e-7, 1e-7);
        boolean[] firstEventOccurred = {false};
        List<Event> events = new ArrayList<>();

        TimeDetector detector1 = new TimeDetector(10, tol, 100, events, 5) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                    firstEventOccurred[0] = true;
                    super.getHandler().eventOccurred(state, detector, increasing);
                    return Action.RESET_STATE;
                };
            }
        };
        integrator.addEventDetector(detector1);
        // this detector changes it's g function definition when detector1 fires
        TimeDetector detector2 = new TimeDetector(1, tol, 100, events, 1, 3, 5) {
            @Override
            public Binary64 g(final FieldODEStateAndDerivative<Binary64> state) {
                if (firstEventOccurred[0]) {
                    return super.g(state);
                }
                return new TimeDetector(1, tol, 100, 5).g(state);
            }
        };
        integrator.addEventDetector(detector2);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(20));

        // verify
        // order is important to make sure the test checks what it is supposed to
        Assert.assertEquals(5, events.get(0).getT(), 0.0);
        Assert.assertTrue(events.get(0).isIncreasing());
        Assert.assertEquals(detector1, events.get(0).getDetector());
        Assert.assertEquals(5, events.get(1).getT(), 0.0);
        Assert.assertTrue(events.get(1).isIncreasing());
        Assert.assertEquals(detector2, events.get(1).getDetector());
        Assert.assertEquals(2, events.size());
    }

    /**
     * When two event detectors have a discontinuous event caused by a {@link
     * Action#RESET_STATE} or {@link Action#RESET_DERIVATIVES}. The two event detectors
     * would each say they had an event that had to be handled before the other one, but
     * neither would actually back up at all. For Hipparchus GitHub #91.
     */
    @Test
    public void testSimultaneousDiscontinuousEventsAfterResetReverse() {
        // setup
        double t = -FastMath.PI;
        double tol = 1e-10;
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        List<Event> events = new ArrayList<>();

        TimeDetector resetDetector =
                new ResetDetector(10, tol, 100, events, new FieldODEState<>(zero.add(t), new Binary64[]{zero.add(1e100), zero}), t);
        integrator.addEventDetector(resetDetector);
        List<BaseDetector> detectors = new ArrayList<>();
        for (int i = 0; i < 2; i++) {
            BaseDetector detector1 = new StateDetector(10, tol, 100, events, 0.0);
            integrator.addEventDetector(detector1);
            detectors.add(detector1);
        }

        // action
        integrator.integrate(new Equation(), new FieldODEState<>(zero, new Binary64[]{zero.add(-1e100), zero}), zero.add(-10));

        // verify
        Assert.assertEquals(t, events.get(0).getT(), tol);
        Assert.assertTrue(events.get(0).isIncreasing());
        Assert.assertEquals(resetDetector, events.get(0).getDetector());
        // next two events can occur in either order
        Assert.assertEquals(t, events.get(1).getT(), tol);
        Assert.assertFalse(events.get(1).isIncreasing());
        Assert.assertEquals(detectors.get(0), events.get(1).getDetector());
        Assert.assertEquals(t, events.get(2).getT(), tol);
        Assert.assertFalse(events.get(2).isIncreasing());
        Assert.assertEquals(detectors.get(1), events.get(2).getDetector());
        Assert.assertEquals(events.size(), 3);
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
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(10, 0.2, 100, -9.9, -10.1, -12);
        integrator.addEventDetector(detector1);

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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, events, t3);
        TimeDetector detectorB = new TimeDetector(maxCheck, tolerance, 100, events, -50, t1, t2, t5);
        TimeDetector detectorC = new TimeDetector(maxCheck, tolerance, 100, Action.RESET_DERIVATIVES, events, t4);

        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);
        integrator.addEventDetector(detectorC);

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
        TimeDetector detectorA = new TimeDetector(maxCheck, 1e-6, 100, events, t2);
        FlatDetector detectorB = new FlatDetector(maxCheck, 0.5, 100, events, t3);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        // if these fail the event finding did its job,
        // but this test isn't testing what it is supposed to be
        Assert.assertSame(detectorB, events.get(0).getDetector());
        Assert.assertSame(detectorA, events.get(1).getDetector());
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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, events, t2);
        TimeDetector detectorB = new TimeDetector(maxCheck, toleranceB, 100, events, -50, t1, t3);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        // we only care that the rules are satisfied. There are multiple solutions.
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), toleranceB);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorB, events.get(0).getDetector());
        Assert.assertEquals(t3, events.get(1).getT(), toleranceB);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertSame(detectorB, events.get(1).getDetector());
        Assert.assertEquals(t2, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertSame(detectorA, events.get(2).getDetector());
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
        TimeDetector detectorB = new FlatDetector(maxCheck, tolerance, 100, events, t1, t2, t3);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        // allowed to report t1 or t3.
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorB, events.get(0).getDetector());
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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, events, t1);
        TimeDetector detectorB = new TimeDetector(maxCheck, tolerance, 100, events, t1, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, events, t1);
        TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, -50, t1, t1);
        detectorB.g(state(t1));
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, -50, t1, t2);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), 0.0);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getDetector());
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
        TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, t1, t2);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(false, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), 0.0);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorA, events.get(1).getDetector());
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
        TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, t6);
        TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, -50, t1, t3, t4, t7);
        TimeDetector detectorC = new ContinuousDetector(maxCheck, tolerance, 100, events, -50, t2, t5);

        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        Assert.assertEquals(5, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(detectorB, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertEquals(detectorC, events.get(1).getDetector());
        // reporting t3 and t4 is optional, seeing them is not.
        // we know a root was found at t3 because events are reported at t2 and t5.
        /*
        Assert.assertEquals(t3, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertEquals(detectorB, events.get(2).getDetector());
        Assert.assertEquals(t4, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(detectorB, events.get(3).getDetector());
        */
        Assert.assertEquals(t5, events.get(2).getT(), tolerance);
        Assert.assertEquals(false, events.get(2).isIncreasing());
        Assert.assertEquals(detectorC, events.get(2).getDetector());
        Assert.assertEquals(t6, events.get(3).getT(), tolerance);
        Assert.assertEquals(true, events.get(3).isIncreasing());
        Assert.assertEquals(detectorA, events.get(3).getDetector());
        Assert.assertEquals(t7, events.get(4).getT(), tolerance);
        Assert.assertEquals(false, events.get(4).isIncreasing());
        Assert.assertEquals(detectorB, events.get(4).getDetector());
    }

    /** Test a reset event triggering another event at the same time. */
    @Test
    public void testEventCausedByStateResetReverse() {
        // setup
        double maxCheck = 10;
        double tolerance = 1e-6;
        double t1 = -15.0;
        final FieldODEState<Binary64> newState =
                new FieldODEState<>(zero.add(t1), new Binary64[]{zero.add(20), zero});
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        TimeDetector detectorA = new ResetDetector(maxCheck, tolerance, 100, events, newState, t1);
        BaseDetector detectorB = new StateDetector(maxCheck, tolerance, 100, events, 1);

        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-40.0));

        //verify
        // really we only care that the Rules of Event Handling are not violated,
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertEquals(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t1, events.get(1).getT(), tolerance);
        Assert.assertEquals(false, events.get(1).isIncreasing());
        Assert.assertEquals(detectorB, events.get(1).getDetector());
        Assert.assertEquals(t1 - 19, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertEquals(detectorB, events.get(2).getDetector());
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
        TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), 0.0);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        TimeDetector detectorA = new FlatDetector(maxCheck, tolerance, 100, events, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        // use root finding tolerance since it is larger
        Assert.assertEquals(t1, events.get(0).getT(), 1e-3);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        final TimeDetector detectorA = new ContinuousDetector(maxCheck, tolerance, 100, events, t1) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                    swap[0] = true;
                    return super.getHandler().eventOccurred(state, detector, increasing);
                };
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, t2);
        BaseDetector detectorC = new BaseDetector(maxCheck, tolerance, 100, Action.CONTINUE, events) {

            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (swap[0]) {
                    return detectorB.g(state);
                } else {
                    return one;
                }
            }

        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorC, events.get(1).getDetector());
    }


    /**
     * test when one event detector changes the definition of another's g function before
     * the end of the step as a result of an event occurring. In this case the change
     * cancels the occurrence of the event.
     */
    @Test
    public void testEventChangesGFunctionDefinitionCancelReverse() {
        // setup
        double maxCheck = 5;
        double tolerance = 1e-6;
        double t1 = -11, t2 = -11.1;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // mutable boolean
        boolean[] swap = new boolean[1];
        final TimeDetector detectorA =
                        new ContinuousDetector(maxCheck, tolerance, 100, Action.RESET_EVENTS, events, t1) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                    swap[0] = true;
                    return super.getHandler().eventOccurred(state, detector, increasing);
                };
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, t2);
        BaseDetector detectorC = new BaseDetector(maxCheck, tolerance, 100, Action.CONTINUE, events) {

            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (!swap[0]) {
                    return detectorB.g(state);
                } else {
                    return zero.add(1);
                }
            }

        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
    }


    /**
     * test when one event detector changes the definition of another's g function before
     * the end of the step as a result of an event occurring. In this case the change
     * delays the occurrence of the event.
     */
    @Test
    public void testEventChangesGFunctionDefinitionDelayReverse() {
        // setup
        double maxCheck = 5;
        double tolerance = 1e-6;
        double t1 = -11, t2 = -11.1, t3 = -11.2;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // mutable boolean
        boolean[] swap = new boolean[1];
        final TimeDetector detectorA =
                        new ContinuousDetector(maxCheck, tolerance, 100, Action.RESET_EVENTS, events, t1) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                    swap[0] = true;
                    return super.getHandler().eventOccurred(state, detector, increasing);
                };
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, t2);
        final TimeDetector detectorD = new ContinuousDetector(maxCheck, tolerance, 100, events, t3);
        BaseDetector detectorC = new BaseDetector(maxCheck, tolerance, 100, Action.CONTINUE, events) {

            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (!swap[0]) {
                    return detectorB.g(state);
                } else {
                    return detectorD.g(state);
                }
            }

        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t3, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorC, events.get(1).getDetector());
    }

    /**
     * test when one event detector changes the definition of another's g function before
     * the end of the step as a result of an event occurring. In this case the change
     * causes the event to happen sooner than originally expected.
     */
    @Test
    public void testEventChangesGFunctionDefinitionAccelerateReverse() {
        // setup
        double maxCheck = 5;
        double tolerance = 1e-6;
        double t1 = -11, t2 = -11.1, t3 = -11.2;
        // shared event list so we know the order in which they occurred
        List<Event> events = new ArrayList<>();
        // mutable boolean
        boolean[] swap = new boolean[1];
        final TimeDetector detectorA =
                        new ContinuousDetector(maxCheck, tolerance, 100, Action.RESET_EVENTS, events, t1) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return (state, detector, increasing) -> {
                    swap[0] = true;
                    return super.getHandler().eventOccurred(state, detector, increasing);
                };
            }
        };
        final TimeDetector detectorB = new ContinuousDetector(maxCheck, tolerance, 100, events, t2);
        final TimeDetector detectorD = new ContinuousDetector(maxCheck, tolerance, 100, events, t3);
        BaseDetector detectorC = new BaseDetector(maxCheck, tolerance, 100, Action.CONTINUE, events) {

            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (swap[0]) {
                    return detectorB.g(state);
                } else {
                    return detectorD.g(state);
                }
            }

        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorC);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(2, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t2, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorC, events.get(1).getDetector());
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
        TimeDetector detectorA = new FlatDetector(maxCheck, tolerance, 100, Action.STOP, events, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        // action
        FieldODEStateAndDerivative<Binary64> finalState =
                integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(1, events.size());
        // use root finder tolerance instead of event finder tolerance.
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(true, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, Action.STOP, events, 50) {
            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                if (state.getTime().getReal() > -2) {
                    return zero;
                } else {
                    return super.g(state);
                }
            }
        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

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
        TimeDetector detectorA = new TimeDetector(maxCheck, tolerance, 100, events) {
            @Override
            public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
                final Binary64 t = state.getTime();
                if (t.getReal() > t1) {
                    return one.negate();
                } else if (t.getReal() > t2) {
                    return one;
                } else {
                    return one.negate();
                }
            }
        };
        TimeDetector detectorB = new TimeDetector(maxCheck, tolerance, 100, events, t1);
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);
        integrator.addEventDetector(detectorB);

        // action
        integrator.integrate(new Equation(), initialState, zero.add(-30.0));

        // verify
        Assert.assertEquals(3, events.size());
        Assert.assertEquals(t1, events.get(0).getT(), tolerance);
        Assert.assertEquals(false, events.get(0).isIncreasing());
        Assert.assertSame(detectorA, events.get(0).getDetector());
        Assert.assertEquals(t1, events.get(1).getT(), tolerance);
        Assert.assertEquals(true, events.get(1).isIncreasing());
        Assert.assertSame(detectorB, events.get(1).getDetector());
        Assert.assertEquals(t2, events.get(2).getT(), tolerance);
        Assert.assertEquals(true, events.get(2).isIncreasing());
        Assert.assertSame(detectorA, events.get(2).getDetector());
    }

    /** Check that steps are restricted correctly with a continue event. */
    @Test
    public void testEventStepHandlerReverse() {
        // setup
        double tolerance = 1e-18;
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(new TimeDetector(100.0, tolerance, 100, -5));
        StepHandler stepHandler = new StepHandler();
        integrator.addStepHandler(stepHandler);

        // action
        FieldODEStateAndDerivative<Binary64> finalState = integrator
                .integrate(new Equation(), initialState, zero.add(-10));

        // verify
        Assert.assertEquals(-10.0, finalState.getTime().getReal(), tolerance);
        Assert.assertEquals(0.0,
                stepHandler.initialState.getTime().getReal(), tolerance);
        Assert.assertEquals(-10.0, stepHandler.finalTime.getReal(), tolerance);
        Assert.assertEquals(-10.0,
                stepHandler.finalState.getTime().getReal(), tolerance);
        FieldODEStateInterpolator<Binary64> interpolator = stepHandler.interpolators.get(0);
        Assert.assertEquals(0.0,
                interpolator.getPreviousState().getTime().getReal(), tolerance);
        Assert.assertEquals(-5.0,
                interpolator.getCurrentState().getTime().getReal(), tolerance);
        interpolator = stepHandler.interpolators.get(1);
        Assert.assertEquals(-5.0,
                interpolator.getPreviousState().getTime().getReal(), tolerance);
        Assert.assertEquals(-10.0,
                interpolator.getCurrentState().getTime().getReal(), tolerance);
        Assert.assertEquals(2, stepHandler.interpolators.size());
    }

    /** Test resetState(...) returns {@code null}. */
    @Test
    public void testEventCausedByDerivativesResetReverse() {
        // setup
        TimeDetector detectorA = new TimeDetector(10, 1e-6, 100, Action.RESET_STATE, -15.0) {
            @Override
            public FieldODEEventHandler<Binary64> getHandler() {
                return new FieldODEEventHandler<Binary64>() {
                    @Override
                    public Action eventOccurred(FieldODEStateAndDerivative<Binary64> state,
                                                FieldODEEventDetector<Binary64> detector,
                                                boolean increasing) {
                        return Action.RESET_STATE;
                    }
                    @Override
                    public FieldODEState<Binary64> resetState(FieldODEEventDetector<Binary64> detector,
                                                               FieldODEStateAndDerivative<Binary64> state) {
                        return null;
                    }
                };
            }
        };
        FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<>(field, 10, 10, 1e-7, 1e-7);
        integrator.addEventDetector(detectorA);

        try {
            // action
            integrator.integrate(new Equation(), initialState, zero.add(-20.0));
            Assert.fail("Expected Exception");
        } catch (NullPointerException e) {
            // expected
        }
    }

    @Test
    public void testResetChangesSignReverse() {
        FieldOrdinaryDifferentialEquation<Binary64> equation = new FieldOrdinaryDifferentialEquation<Binary64>() {
            public int getDimension() { return 1; }
            public Binary64[] computeDerivatives(Binary64 t, Binary64[] y) { return new Binary64[] { new Binary64(1.0) }; }
        };

        LutherFieldIntegrator<Binary64> integrator = new LutherFieldIntegrator<>(Binary64Field.getInstance(), new Binary64(20.0));
        final double small = 1.0e-10;
        ResetChangesSignGenerator eventsGenerator = new ResetChangesSignGenerator(-6.0, -9.0, +0.5 * small, 8.0, small, 1000);
        integrator.addEventDetector(eventsGenerator);
        final FieldODEStateAndDerivative<Binary64> end = integrator.integrate(new FieldExpandableODE<>(equation),
                                                                               new FieldODEState<>(new Binary64(0.0),
                                                                                                   new Binary64[] { new Binary64(0.0) }),
                                                                               new Binary64(-100.0));
        Assert.assertEquals(2,                  eventsGenerator.getCount());
        Assert.assertEquals(-9.0,               end.getCompleteState()[0].getReal(), 1.0e-12);
        Assert.assertEquals(-9.0 - 0.5 * small, end.getTime().getReal(),             1.0e-12);
    }

    /* utility classes and methods */

    /**
     * Create a state at a time.
     *
     * @param t time of state.
     * @return new state.
     */
    private FieldODEStateAndDerivative<Binary64> state(double t) {
        return new FieldODEStateAndDerivative<>(
                zero.add(t), new Binary64[0], new Binary64[0]);
    }

    /** Base class to record events that occurred. */
    private static abstract class BaseDetector implements FieldODEEventDetector<Binary64> {

        private final FieldAdaptableInterval<Binary64>             maxCheck;
        private final int                                          maxIter;
        private final BracketedRealFieldUnivariateSolver<Binary64> solver;
        protected final Action                                     action;

        /** times the event was actually triggered. */
        private final List<Event> events;

        public BaseDetector(final double maxCheck, final double threshold, final int maxIter,
                            Action action, List<Event> events) {
            this.maxCheck  = s -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new FieldBracketingNthOrderBrentSolver<>(new Binary64(0),
                                                                      new Binary64(threshold),
                                                                      new Binary64(0),
                                                                      5);
            this.action    = action;
            this.events    = events;
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

        /**
         * Get all events that occurred.
         *
         * @return list of events in the order they occurred.
         */
        public List<Event> getEvents() {
            return events;
        }

        @Override
        public FieldODEEventHandler<Binary64> getHandler() {
            return new FieldODEEventHandler<Binary64>() {
                @Override
                public Action eventOccurred(FieldODEStateAndDerivative<Binary64> state,
                                            FieldODEEventDetector<Binary64> detector,
                                            boolean increasing) {
                    events.add(new Event(state, detector, increasing));
                    return action;
                }
            };
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
        public TimeDetector(final double maxCheck, final double threshold, final int maxIter,
                            final double... eventTs) {
            this(maxCheck, threshold, maxIter, Action.CONTINUE, eventTs);
        }

        public TimeDetector(final double maxCheck, final double threshold, final int maxIter,
                            final List<Event> events, final double... eventTs) {
            this(maxCheck, threshold, maxIter, Action.CONTINUE, events, eventTs);
        }

        public TimeDetector(final double maxCheck, final double threshold, final int maxIter,
                            final Action action, final double... eventTs) {
            this(maxCheck, threshold, maxIter, action, new ArrayList<Event>(), eventTs);
        }

        public TimeDetector(final double maxCheck, final double threshold, final int maxIter,
                            final Action action, final List<Event> events, final double... eventTs) {
            super(maxCheck, threshold, maxIter, action, events);
            this.eventTs = eventTs.clone();
            Arrays.sort(this.eventTs);
        }

        @Override
        public Binary64 g(final FieldODEStateAndDerivative<Binary64> state) {
            final Binary64 t = state.getTime();
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

        private final FieldODEStateAndDerivative<Binary64> state;
        private final boolean increasing;
        private final FieldODEEventDetector<Binary64> detector;

        public Event(FieldODEStateAndDerivative<Binary64> state,
                     FieldODEEventDetector<Binary64> detector,
                     boolean increasing) {
            this.increasing = increasing;
            this.state = state;
            this.detector = detector;
        }

        public boolean isIncreasing() {
            return increasing;
        }

        public double getT() {
            return state.getTime().getReal();
        }

        public FieldODEEventDetector<Binary64> getDetector() {
            return detector;
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

        public FlatDetector(final double maxCheck, final double threshold, final int maxIter,
                            final List<Event> events, final double... eventTs) {
            super(maxCheck, threshold, maxIter, events, eventTs);
        }

        public FlatDetector(final double maxCheck, final double threshold, final int maxIter,
                            final Action action, final List<Event> events, final double... eventTs) {
            super(maxCheck, threshold, maxIter, action, events, eventTs);
        }

        @Override
        public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
            final Binary64 g = super.g(state);
            return g.sign();
        }

    }

    /** Linear on both ends, parabolic in the middle. */
    private static class ContinuousDetector extends TimeDetector {

        public ContinuousDetector(final double maxCheck, final double threshold, final int maxIter,
                                  final List<Event> events, final double... eventTs) {
            super(maxCheck, threshold, maxIter, events, eventTs);
        }

        public ContinuousDetector(final double maxCheck, final double threshold, final int maxIter,
                                  final Action action, final List<Event> events, final double... eventTs) {
            super(maxCheck, threshold, maxIter, action, events, eventTs);
        }

        @Override
        public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
            final Binary64 t = state.getTime();
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

        private final FieldODEState<Binary64> resetState;

        public ResetDetector(final double maxCheck, final double threshold, final int maxIter,
                             final List<Event> events, final FieldODEState<Binary64> state, final double eventT) {
            super(maxCheck, threshold, maxIter, Action.RESET_STATE, events, eventT);
            this.resetState = state;
        }

        @Override
        public FieldODEEventHandler<Binary64> getHandler() {
            return new FieldODEEventHandler<Binary64>() {
                @Override
                public Action eventOccurred(FieldODEStateAndDerivative<Binary64> state,
                                            FieldODEEventDetector<Binary64> detector, boolean increasing) {
                    return ResetDetector.super.getHandler().eventOccurred(state, detector, increasing);
                }
                @Override
                public FieldODEState<Binary64> resetState(FieldODEEventDetector<Binary64> detector,
                                                           FieldODEStateAndDerivative<Binary64> state) {
                    Assert.assertEquals(eventTs[0], state.getTime().getReal(), 0);
                    return resetState;
                }
            };
        }

    }

    /** Switching function based on the first primary state. */
    private static class StateDetector extends BaseDetector {

        private final double triggerState;

        public StateDetector(final double maxCheck, final double threshold, final int maxIter,
                             final List<Event> events, final double triggerState) {
            super(maxCheck, threshold, maxIter, Action.CONTINUE, events);
            this.triggerState = triggerState;
        }

        @Override
        public Binary64 g(FieldODEStateAndDerivative<Binary64> state) {
            return state.getPrimaryState()[0].subtract(this.triggerState);
        }
    }

    /** Some basic equations to integrate. */
    public static class Equation extends FieldExpandableODE<Binary64> {
        public Equation() {
            super(new EquationODE());
        }
    }

    /** Some basic equations to integrate. */
    public static class EquationODE implements FieldOrdinaryDifferentialEquation<Binary64> {

        public int getDimension() {
            return 2;
        }

        @Override
        public Binary64[] computeDerivatives(Binary64 t, Binary64[] y) {
            return new Binary64[]{one, one.multiply(2)};
        }

    }

    private static class StepHandler implements FieldODEStepHandler<Binary64> {

        private FieldODEStateAndDerivative<Binary64> initialState;
        private Binary64 finalTime;
        private List<FieldODEStateInterpolator<Binary64>> interpolators = new ArrayList<>();
        private FieldODEStateAndDerivative<Binary64> finalState;

        @Override
        public void init(FieldODEStateAndDerivative<Binary64> initialState,
                         Binary64 finalTime) {
            this.initialState = initialState;
            this.finalTime = finalTime;
        }

        @Override
        public void handleStep(FieldODEStateInterpolator<Binary64> interpolator) {
            this.interpolators.add(interpolator);
        }

        @Override
        public void finish(FieldODEStateAndDerivative<Binary64> finalState) {
            this.finalState = finalState;
        }
    }

    private class ResetChangesSignGenerator implements FieldODEEventDetector<Binary64> {

        private final FieldAdaptableInterval<Binary64>             maxCheck;
        private final int                                          maxIter;
        private final BracketedRealFieldUnivariateSolver<Binary64> solver;
        final double                                               y1;
        final double                                               y2;
        final double                                               change;
        int                                                        count;

        public ResetChangesSignGenerator(final double y1, final double y2, final double change,
                                         final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = s -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new FieldBracketingNthOrderBrentSolver<>(new Binary64(0),
                            new Binary64(threshold),
                            new Binary64(0),
                            5);
            this.y1        = y1;
            this.y2        = y2;
            this.change    = change;
            this.count     = 0;
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
            return new FieldODEEventHandler<Binary64>() {
                public Action eventOccurred(FieldODEStateAndDerivative<Binary64> s,
                                            FieldODEEventDetector<Binary64> detector,
                                            boolean increasing) {
                    return ++count < 2 ? Action.RESET_STATE : Action.STOP;
                }

                public FieldODEState<Binary64> resetState(FieldODEEventDetector<Binary64> detector,
                                                           FieldODEStateAndDerivative<Binary64> s) {
                    return new FieldODEState<>(s.getTime(), new Binary64[] { s.getCompleteState()[0].add(change) });
                }
            };
        }

        public Binary64 g(FieldODEStateAndDerivative<Binary64> s) {
            return s.getCompleteState()[0].subtract(y1).multiply(s.getCompleteState()[0].subtract(y2));
        }

        public int getCount() {
            return count;
        }

    }

}
