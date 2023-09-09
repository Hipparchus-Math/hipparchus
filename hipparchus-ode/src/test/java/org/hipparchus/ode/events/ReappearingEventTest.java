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

import java.util.Arrays;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.ode.nonstiff.GraggBulirschStoerIntegrator;
import org.junit.Assert;
import org.junit.Test;

public class ReappearingEventTest {

    @Test
    public void testDormandPrince()
        throws MathIllegalArgumentException, MathIllegalStateException {
        double tEnd = test(1);
        Assert.assertEquals(10.0, tEnd, 1e-7);
    }

    @Test
    public void testGragg()
        throws MathIllegalArgumentException, MathIllegalStateException {
        double tEnd = test(2);
        Assert.assertEquals(10.0, tEnd, 1e-7);
    }

    public double test(int integratorType)
        throws MathIllegalArgumentException, MathIllegalStateException {
        double e = 1e-15;
        ODEIntegrator integrator = (integratorType == 1) ?
                                   new DormandPrince853Integrator(e, 100.0, 1e-7, 1e-7) :
                                   new GraggBulirschStoerIntegrator(e, 100.0, 1e-7, 1e-7);
        integrator.addEventDetector(new Event(0.1, e, 1000));
        double t0 = 6.0;
        double tEnd = 10.0;
        double[] y = {2.0, 2.0, 2.0, 4.0, 2.0, 7.0, 15.0};
        return integrator.integrate(new Ode(), new ODEState(t0, y), tEnd).getTime();
    }

    private static class Ode implements OrdinaryDifferentialEquation {
        public int getDimension() {
            return 7;
        }

        public double[] computeDerivatives(double t, double[] y) {
            double[] yDot = new double[y.length];
            Arrays.fill(yDot, 1.0);
            return yDot;
        }
    }

    /** State events for this unit test. */
    protected static class Event implements ODEEventDetector {

        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;

        /** Constructor for the {@link Event} class.
         * @param maxCheck maximum checking interval, must be strictly positive (s)
         * @param threshold convergence threshold (s)
         * @param maxIter maximum number of iterations in the event time search
         */
        public Event(final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = s -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
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

        public double g(ODEStateAndDerivative s) {
            return s.getPrimaryState()[6] - 15.0;
        }

    }

}
