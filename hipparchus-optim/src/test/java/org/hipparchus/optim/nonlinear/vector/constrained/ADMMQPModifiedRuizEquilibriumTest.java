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
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class ADMMQPModifiedRuizEquilibriumTest {

    @Test
    public void testScalingLimits() {
        // column 0 has norm 1e8 (clipped to 1e6), A row 0 has norm 1e-8 (below 1e-6, replaced by 1)
        // and A row 1 is zero (replaced by 1)
        final RealMatrix h = MatrixUtils.createRealDiagonalMatrix(new double[] { 1.0e8, 1.0 });
        final RealMatrix a = MatrixUtils.createRealMatrix(new double[][] {
                { 1.0e-8, 0.0 },
                { 0.0,    0.0 }
        });
        final RealVector q = new ArrayRealVector(new double[] { 1.0, 1.0 });
        final ADMMQPModifiedRuizEquilibrium ruiz = new ADMMQPModifiedRuizEquilibrium(h, a, q);
        ruiz.normalize(1.0e-9, 1);

        // D = diag(1 / √1e6, 1), E = I
        final RealVector x = new ArrayRealVector(new double[] { 2.0, 3.0 });
        Assertions.assertArrayEquals(new double[] { 2.0e-3, 3.0 }, ruiz.unscaleX(x).toArray(), 1.0e-18);
        Assertions.assertArrayEquals(new double[] { 2.0,    3.0 }, ruiz.unscaleZ(x).toArray(), 0.0);
        Assertions.assertEquals(1.0e-11, ruiz.getScaledA().getEntry(0, 0), 1.0e-25);
        Assertions.assertEquals(0.0,     ruiz.getScaledA().getEntry(1, 1), 0.0);

        // scaling and unscaling are inverse operations
        Assertions.assertArrayEquals(x.toArray(), ruiz.unscaleX(ruiz.scaleX(x)).toArray(), 1.0e-15);
    }

}
