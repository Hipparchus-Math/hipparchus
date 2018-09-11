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
package org.hipparchus.ode.sampling;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.TestFieldProblem3;
import org.hipparchus.ode.nonstiff.DormandPrince54FieldIntegrator;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;


public class FieldStepNormalizerTest {

    TestFieldProblem3<Decimal64> pb;
    FieldODEIntegrator<Decimal64> integ;
    boolean lastSeen;

    @Test
    public void testBoundaries()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        Decimal64 range = pb.getFinalTime().subtract(pb.getInitialTime());
        setLastSeen(false);
        integ.addStepHandler(new FieldStepNormalizer<>(range.divide(10.0).doubleValue(),
                                                       new FieldODEFixedStepHandler<Decimal64>() {
            private boolean firstCall = true;
            public void handleStep(FieldODEStateAndDerivative<Decimal64> s, boolean isLast) {
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
        integ.integrate(new FieldExpandableODE<>(pb), pb.getInitialState(), pb.getFinalTime());
        Assert.assertTrue(lastSeen);
    }

    @Test
    public void testBeforeEnd()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        final Decimal64 range = pb.getFinalTime().subtract(pb.getInitialTime());
        setLastSeen(false);
        integ.addStepHandler(new FieldStepNormalizer<>(range.divide(10.5).doubleValue(),
                                                       new FieldODEFixedStepHandler<Decimal64>() {
            public void handleStep(FieldODEStateAndDerivative<Decimal64> s, boolean isLast) {
                if (isLast) {
                    setLastSeen(true);
                    checkValue(s.getTime(), pb.getFinalTime().subtract(range.divide(21.0)));
                }
            }
        }));
        integ.integrate(new FieldExpandableODE<>(pb), pb.getInitialState(), pb.getFinalTime());
        Assert.assertTrue(lastSeen);
    }

    public void checkValue(Decimal64 value, Decimal64 reference) {
        Assert.assertTrue(FastMath.abs(value.doubleValue() - reference.doubleValue()) < 1.0e-10);
    }

    public void setLastSeen(boolean lastSeen) {
        this.lastSeen = lastSeen;
    }

    @Before
    public void setUp() {
        pb = new TestFieldProblem3<Decimal64>(new Decimal64(0.9));
        double minStep = 0;
        double maxStep = pb.getFinalTime().doubleValue() - pb.getInitialTime().doubleValue();
        integ = new DormandPrince54FieldIntegrator<>(Decimal64Field.getInstance(),
                                                     minStep, maxStep, 10.e-8, 1.0e-8);
        lastSeen = false;
    }

    @After
    public void tearDown() {
        pb    = null;
        integ = null;
    }

}
