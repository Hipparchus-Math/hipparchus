/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.optim.nonlinear.vector.constrained;


import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

class LinearSearchTest {

    @Test
    public void testMonotoneSuccess() {
        final LineSearch ls = new LineSearch(1.0e-7, 5, 1.0e-4, 0.5, 50, 2);
        double d = ls.search(buildMeritFunction(0.99, 0.8));
        Assertions.assertEquals(0.125, d, 1.0e-3);
        Assertions.assertFalse(ls.isBadStepDetected());
        Assertions.assertFalse(ls.isBadStepFailed());
    }

    @Test
    public void testBadStep() {
        final LineSearch ls = new LineSearch(1.0e-7, 5, 1.0e-4, 0.5, 50, 2);
        double d = ls.search(buildMeritFunction(-2.0, -0.5));
        Assertions.assertEquals(7.071e-7, d, 1.0e-10);
        Assertions.assertTrue(ls.isBadStepDetected());
        Assertions.assertFalse(ls.isBadStepFailed());
    }

    private MeritFunctionL2 buildMeritFunction(final double startX, final double startY) {
        final TwiceDifferentiableFunction objective = new RosenbrockFunction();
        final InequalityConstraint iqConstraint =
                new LinearInequalityConstraint(new double[][] {
                        { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 }
                }, new double[] {
                        -1.5, -1.5, -1.5, -1.5
                });
        final RealVector x = new ArrayRealVector(new double[] { startX, startY });
        final MeritFunctionL2 f = new MeritFunctionL2(objective, null, iqConstraint, x);
        f.update(objective.gradient(x),
                 null,
                 iqConstraint.gradient(x.toArray()),
                 x,
                 new ArrayRealVector(iqConstraint.dimY()),
                 new ArrayRealVector(new double[] { -1.0, 1.0 }),
                 new ArrayRealVector(new double[4]));
        return f;
    }

}
