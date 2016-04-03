package org.hipparchus.ode.events;

import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.junit.Assert;
import org.junit.Test;

/**
 * Unit tests for ODE integrators with near simultaneous events.
 *
 * @author Evan Ward
 */
public class CloseEventsTest {

    /** check nearly simultaneous events. */
    @Test
    public void testCloseEvents() {
        // setup
        double e = 1e-15;
        ODEIntegrator integrator = new DormandPrince853Integrator(e, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(5);
        integrator.addEventHandler(detector1, 10, 1, 100);
        TimeDetector detector2 = new TimeDetector(5.5);
        integrator.addEventHandler(detector2, 10, 1, 100);

        // action
        integrator.integrate(new Equation(), new ODEState(0, new double[2]), 20);

        // verify
        Assert.assertEquals(5, detector1.getActualT(), 0.0);
        Assert.assertEquals(5.5, detector2.getActualT(), 0.0);
    }

    /** test simultaneous events. */
    @Test
    public void testSimultaneousEvents() {
        // setup
        double e = 1e-15;
        ODEIntegrator integrator = new DormandPrince853Integrator(e, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(5);
        integrator.addEventHandler(detector1, 10, 1, 100);
        TimeDetector detector2 = new TimeDetector(5);
        integrator.addEventHandler(detector2, 10, 1, 100);

        // action
        integrator.integrate(new Equation(), new ODEState(0, new double[2]), 20);

        // verify
        Assert.assertEquals(5, detector1.getActualT(), 0.0);
        Assert.assertEquals(5, detector2.getActualT(), 0.0);
    }

    /** Mock event handler */
    private static class TimeDetector implements ODEEventHandler {

        /** requested event time. */
        private final double eventT;

        /** actual event time. */
        private double actualT;

        /**
         * Create an event detector that triggers an event at a particular time.
         *
         * @param eventT time the event should occur.
         */
        public TimeDetector(double eventT) {
            this.eventT = eventT;
        }

        /**
         * Get the time of the actual event.
         *
         * @return the actual time the event occurred.
         */
        public double getActualT() {
            return actualT;
        }

        @Override
        public double g(ODEStateAndDerivative s) {
            return s.getTime() - eventT;
        }

        @Override
        public Action eventOccurred(ODEStateAndDerivative s, boolean increasing) {
            this.actualT = s.getTime();
            return Action.CONTINUE;
        }

    }

    /** arbitrary equations. */
    public static class Equation implements OrdinaryDifferentialEquation {

        public int getDimension() {
            return 2;
        }

        public double[] computeDerivatives(double t, double[] y) {
            return new double[] { 1.0, 2.0 };
        }

    }

}
