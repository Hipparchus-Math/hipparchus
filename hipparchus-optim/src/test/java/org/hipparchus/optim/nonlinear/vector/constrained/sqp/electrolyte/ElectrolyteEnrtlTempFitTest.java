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

/** ripopt electrolyte problem 12: temperature-dependent eNRTL fit. */
public class ElectrolyteEnrtlTempFitTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 3.0e-2;

    private static final double[] TRUE_PARAMS = { 8.045, -3987.0, -4.549, 2216.0 };
    private static final double[] TEMPERATURES = { 288.15, 298.15, 308.15, 318.15 };
    private static final double[] MOLALITIES = { 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0 };

    private static double model(final RealVector x, final double m, final double temperature) {
        final double tauCa = x.getEntry(0) + x.getEntry(1) / temperature;
        final double tauWc = x.getEntry(2) + x.getEntry(3) / temperature;
        final double sqrtM = FastMath.sqrt(m);
        return -ElectrolyteTestSupport.A_PHI * sqrtM / (1.0 + sqrtM) +
               m * tauCa * FastMath.exp(-0.2 * tauCa) +
               m * m * tauWc * FastMath.exp(-0.2 * tauWc);
    }

    private static RealVector truth() {
        return new ArrayRealVector(TRUE_PARAMS, false);
    }

    private static final class Objective extends ElectrolyteTestSupport.Objective {
        Objective() { super(4); }
        @Override public double value(final RealVector x) {
            final RealVector truth = truth();
            double value = 0.0;
            for (final double temperature : TEMPERATURES) {
                for (final double m : MOLALITIES) {
                    final double residual = model(x, m, temperature) - model(truth, m, temperature);
                    value += residual * residual;
                }
            }
            return value;
        }
        @Override public RealVector gradient(final RealVector x) {
            final RealVector truth = truth();
            final RealVector gradient = new ArrayRealVector(4);
            for (final double temperature : TEMPERATURES) {
                for (final double m : MOLALITIES) {
                    final double tauCa = x.getEntry(0) + x.getEntry(1) / temperature;
                    final double tauWc = x.getEntry(2) + x.getEntry(3) / temperature;
                    final double eCa = FastMath.exp(-0.2 * tauCa);
                    final double eWc = FastMath.exp(-0.2 * tauWc);
                    final double residual = model(x, m, temperature) - model(truth, m, temperature);
                    final double dCa = m * eCa * (1.0 - 0.2 * tauCa);
                    final double dWc = m * m * eWc * (1.0 - 0.2 * tauWc);
                    gradient.addToEntry(0, 2.0 * residual * dCa);
                    gradient.addToEntry(1, 2.0 * residual * dCa / temperature);
                    gradient.addToEntry(2, 2.0 * residual * dWc);
                    gradient.addToEntry(3, 2.0 * residual * dWc / temperature);
                }
            }
            return gradient;
        }
    }

    @Test
    public void testEnrtlTempFit() {
        ElectrolyteTestSupport.solve(
                "eNRTL T-dependent fit",
                new Objective(),
                new double[] { 7.0, -3500.0, -3.5, 1800.0 },
                new SimpleBounds(new double[] { -20.0, -10000.0, -20.0, -10000.0 },
                                 new double[] { 20.0, 10000.0, 20.0, 10000.0 }),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE);
    }
}
