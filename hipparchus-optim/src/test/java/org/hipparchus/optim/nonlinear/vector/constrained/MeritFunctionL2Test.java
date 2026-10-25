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
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealVector;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class MeritFunctionL2Test {

    /** Initial penalty value used by MeritFunctionL2. */
    private static final double SIGMA_INIT = 1.0e-2;

    @Test
    public void testRUpIsCapped() {
        final MeritFunctionL2 merit = createMerit();
        Assertions.assertEquals(1, merit.getR().getDimension());
        Assertions.assertEquals(SIGMA_INIT, merit.getR().getEntry(0), 0.0);

        merit.rUp();
        Assertions.assertEquals(10 * SIGMA_INIT, merit.getR().getEntry(0), 1.0e-15);

        // 0.1 * 10^k reaches the cap after 10 more steps, the 11th checks it stays there
        for (int k = 0; k < 11; ++k) {
            merit.rUp();
        }
        Assertions.assertEquals(merit.getRmax(), merit.getR().getEntry(0), 0.0);

        merit.resetRj();
        Assertions.assertEquals(SIGMA_INIT, merit.getR().getEntry(0), 0.0);
    }

    @Test
    public void testUpdateRjWithoutMultipliers() {
        final MeritFunctionL2 merit = createMerit();
        final RealVector dx = new ArrayRealVector(new double[] { 1.0, 1.0 });
        final RealVector u  = new ArrayRealVector(new double[] { 5.0 });
        merit.updateRj(MatrixUtils.createRealIdentityMatrix(2), null, dx, u, 0.0, 1);
        Assertions.assertEquals(SIGMA_INIT, merit.getR().getEntry(0), 0.0);
        merit.updateRj(MatrixUtils.createRealIdentityMatrix(2), new ArrayRealVector(0), dx, u, 0.0, 1);
        Assertions.assertEquals(SIGMA_INIT, merit.getR().getEntry(0), 0.0);
    }

    /** f(x) = ½ |x|², subject to x₀ ≥ 1. */
    private static MeritFunctionL2 createMerit() {
        final QuadraticFunction objective =
                new QuadraticFunction(MatrixUtils.createRealIdentityMatrix(2), new ArrayRealVector(2), 0.0);
        final LinearInequalityConstraint iq =
                new LinearInequalityConstraint(new double[][] { { 1.0, 0.0 } }, new double[] { 1.0 });
        return new MeritFunctionL2(objective, null, iq, null, new ArrayRealVector(new double[] { 2.0, 0.0 }));
    }

}
