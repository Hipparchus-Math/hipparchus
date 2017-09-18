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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.ode.sampling;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.nonstiff.DormandPrince54Integrator;
import org.hipparchus.util.FastMath;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;


public class StepNormalizerTest {

    TestProblem3 pb;
    ODEIntegrator integ;
    boolean lastSeen;

    @Test
    public void testBoundaries()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        double range = pb.getFinalTime() - pb.getInitialTime();
        setLastSeen(false);
        integ.addStepHandler(new StepNormalizer(range / 10.0,
                                                new ODEFixedStepHandler() {
            private boolean firstCall = true;
            public void handleStep(ODEStateAndDerivative s, boolean isLast) {
                if (firstCall) {
                    checkValue(s.getTime(), pb.getInitialTime());
                    firstCall = false;
                }
                if (isLast) {
                    setLastSeen(true);
                    checkValue(s.getTime(), pb.getFinalTime());
                }
            }
        }));
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());
        Assert.assertTrue(lastSeen);
    }

    @Test
    public void testBeforeEnd()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        final double range = pb.getFinalTime() - pb.getInitialTime();
        setLastSeen(false);
        integ.addStepHandler(new StepNormalizer(range / 10.5,
                                                new ODEFixedStepHandler() {
            public void handleStep(ODEStateAndDerivative s, boolean isLast) {
                if (isLast) {
                    setLastSeen(true);
                    checkValue(s.getTime(), pb.getFinalTime() - range / 21.0);
                }
            }
        }));
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());
        Assert.assertTrue(lastSeen);
    }

    public void checkValue(double value, double reference) {
        Assert.assertTrue(FastMath.abs(value - reference) < 1.0e-10);
    }

    public void setLastSeen(boolean lastSeen) {
        this.lastSeen = lastSeen;
    }

    @Before
    public void setUp() {
        pb = new TestProblem3(0.9);
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        integ = new DormandPrince54Integrator(minStep, maxStep, 10.e-8, 1.0e-8);
        lastSeen = false;
    }

    @After
    public void tearDown() {
        pb    = null;
        integ = null;
    }

}
