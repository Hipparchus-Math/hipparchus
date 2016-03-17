package org.hipparchus.ode.events;

import org.hipparchus.ode.FirstOrderDifferentialEquations;
import org.hipparchus.ode.FirstOrderIntegrator;
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
        FirstOrderIntegrator integrator =
                new DormandPrince853Integrator(e, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(5);
        integrator.addEventHandler(detector1, 10, 1, 100);
        TimeDetector detector2 = new TimeDetector(5.5);
        integrator.addEventHandler(detector2, 10, 1, 100);

        // action
        integrator.integrate(new Equation(), 0, new double[2], 20, new double[2]);

        // verify
        Assert.assertEquals(5, detector1.getActualT(), 0.0);
        Assert.assertEquals(5.5, detector2.getActualT(), 0.0);
    }

    /** test simultaneous events. */
    @Test
    public void testSimultaneousEvents() {
        // setup
        double e = 1e-15;
        FirstOrderIntegrator integrator =
                new DormandPrince853Integrator(e, 100.0, 1e-7, 1e-7);

        TimeDetector detector1 = new TimeDetector(5);
        integrator.addEventHandler(detector1, 10, 1, 100);
        TimeDetector detector2 = new TimeDetector(5);
        integrator.addEventHandler(detector2, 10, 1, 100);

        // action
        integrator.integrate(new Equation(), 0, new double[2], 20, new double[2]);

        // verify
        Assert.assertEquals(5, detector1.getActualT(), 0.0);
        Assert.assertEquals(5, detector2.getActualT(), 0.0);
    }

    /** Mock event handler */
    private static class TimeDetector implements EventHandler {

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
        public void init(double t0, double[] y0, double t) {
        }

        @Override
        public double g(double t, double[] y) {
            return t - eventT;
        }

        @Override
        public Action eventOccurred(double t, double[] y, boolean increasing) {
            this.actualT = t;
            return Action.CONTINUE;
        }

        @Override
        public void resetState(double t, double[] y) {
        }

    }

    /** arbitrary equations. */
    public static class Equation implements FirstOrderDifferentialEquations {

        public int getDimension() {
            return 2;
        }

        public void computeDerivatives(double t, double[] y, double[] yDot) {
            yDot[0] = 1.0;
            yDot[1] = 2.0;
        }

    }

}
