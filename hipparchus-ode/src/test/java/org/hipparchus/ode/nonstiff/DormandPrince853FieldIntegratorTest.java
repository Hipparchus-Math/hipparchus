/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
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

package org.hipparchus.ode.nonstiff;


import java.util.Locale;

import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.util.Decimal64Field;
import org.junit.Test;

public class DormandPrince853FieldIntegratorTest extends EmbeddedRungeKuttaFieldIntegratorAbstractTest {

    protected <T extends RealFieldElement<T>> EmbeddedRungeKuttaFieldIntegrator<T>
    createIntegrator(Field<T> field, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance) {
        return new DormandPrince853FieldIntegrator<T>(field, minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
    }

    protected <T extends RealFieldElement<T>> EmbeddedRungeKuttaFieldIntegrator<T>
    createIntegrator(Field<T> field, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance) {
        return new DormandPrince853FieldIntegrator<T>(field, minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
    }

    private static class CircleODE implements OrdinaryDifferentialEquation {
        
        private double[] c;
        private double omega;
    
        public CircleODE(double[] c, double omega) {
            this.c     = c;
            this.omega = omega;
        }
    
        public int getDimension() {
            return 2;
        }
    
        public double[] computeDerivatives(double t, double[] y) {
            return new double[] {
                omega * (c[1] - y[1]),
                omega * (y[0] - c[0])
            };
        }
    
    }
    @Test
    public void testTutorial() {
        ODEIntegrator         dp853        = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
        OrdinaryDifferentialEquation ode          = new CircleODE(new double[] { 1.0, 1.0 }, 0.1);
        ODEState                     initialState = new ODEState(0.0, new double[] { 0.0, 1.0 });
        ODEStepHandler stepHandler = new ODEStepHandler() {
            public void handleStep(ODEStateInterpolator interpolator, boolean isLast) {
                double stepStart = interpolator.getPreviousState().getTime();
                double stepEnd   = interpolator.getCurrentState().getTime();
                for (int i = 0; i < 20; ++i) {
                    // we want to print 20 points for each step
                    double t = ((20 - i) * stepStart + i * stepEnd) / 20;
                    double[] y = interpolator.getInterpolatedState(t).getPrimaryState();
                    System.out.println(t + " " + y[0] + " " + y[1]);
                }
            }
        };
        dp853.addStepHandler(stepHandler);
        ODEStateAndDerivative        finalState   = dp853.integrate(ode, initialState, 16.0);
        double                       t            = finalState.getTime();
        double[]                     y            = finalState.getPrimaryState();
        System.out.format(Locale.US, "final state at %4.1f: %6.3f %6.3f%n", t, y[0], y[1]);

    }
    @Override
    public void testNonFieldIntegratorConsistency() {
        doTestNonFieldIntegratorConsistency(Decimal64Field.getInstance());
    }

    @Override
    public void testSanityChecks() {
        doTestSanityChecks(Decimal64Field.getInstance());
    }

    @Override
    public void testBackward() {
        doTestBackward(Decimal64Field.getInstance(), 8.1e-8, 1.1e-7, 1.0e-12, "Dormand-Prince 8 (5, 3)");
    }

    @Override
    public void testKepler() {
        doTestKepler(Decimal64Field.getInstance(), 4.4e-11);
    }

    @Override
    public void testForwardBackwardExceptions() {
        doTestForwardBackwardExceptions(Decimal64Field.getInstance());
    }

    @Override
    public void testMinStep() {
        doTestMinStep(Decimal64Field.getInstance());
    }

    @Override
    public void testIncreasingTolerance() {
        // the 1.3 factor is only valid for this test
        // and has been obtained from trial and error
        // there is no general relation between local and global errors
        doTestIncreasingTolerance(Decimal64Field.getInstance(), 1.3, 1.0e-12);
    }

    @Override
    public void testEvents() {
        doTestEvents(Decimal64Field.getInstance(), 2.1e-7, "Dormand-Prince 8 (5, 3)");
    }

    @Override
    public void testEventsErrors() {
        doTestEventsErrors(Decimal64Field.getInstance());
    }

    @Override
    public void testEventsNoConvergence() {
        doTestEventsNoConvergence(Decimal64Field.getInstance());
    }

    @Override
    public void testPartialDerivatives() {
        doTestPartialDerivatives(2.6e-12, new double[] { 1.3e-11, 3.6e-12, 5.2e-13, 3.6e-12, 3.6e-12 });
    }

}
