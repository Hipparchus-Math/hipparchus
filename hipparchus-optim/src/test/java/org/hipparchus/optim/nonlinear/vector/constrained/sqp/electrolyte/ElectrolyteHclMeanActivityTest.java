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
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 6: HCl mean activity. */
public class ElectrolyteHclMeanActivityTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = 0;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 9.9999999999999995e-07;

    private static final double BETA0 = 0.1775;
    private static final double BETA1 = 0.2945;
    private static final double C_PHI = 0.00080;
    private static final double TARGET =
            ElectrolyteTestSupport.pitzerLnGammaPm(1.0, BETA0, BETA1, C_PHI);

    private static final class Objective extends ElectrolyteTestSupport.Objective {
        Objective() { super(1); }
        @Override public double value(final RealVector x) {
            final double m = x.getEntry(0);
            final double residual = ElectrolyteTestSupport.pitzerLnGammaPm(m, BETA0, BETA1, C_PHI) +
                                    FastMath.log(m) - TARGET;
            return residual * residual;
        }
        @Override public RealVector gradient(final RealVector x) {
            final double m = x.getEntry(0);
            final double residual = ElectrolyteTestSupport.pitzerLnGammaPm(m, BETA0, BETA1, C_PHI) +
                                    FastMath.log(m) - TARGET;
            final double derivative = ElectrolyteTestSupport.dPitzerLnGammaPmDm(m, BETA0, BETA1, C_PHI) +
                                      1.0 / m;
            return new ArrayRealVector(new double[] { 2.0 * residual * derivative }, false);
        }
    }

    @Test
    public void testHclMeanActivity() {
        ElectrolyteTestSupport.solve(
                "HCl mean activity",
                new Objective(),
                new double[] { 0.5 },
                new SimpleBounds(new double[] { 0.01 }, new double[] { 5.0 }),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE);
    }
}
