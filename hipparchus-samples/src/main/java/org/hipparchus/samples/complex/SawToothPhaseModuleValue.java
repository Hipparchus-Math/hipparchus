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

package org.hipparchus.samples.complex;

import org.hipparchus.complex.Complex;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/** Domain coloring enhancing both phase and module changes.
 * <p>
 * Value encodes both phase and module using two sawtooth functions.
 * </p>
 * <p>
 * The sawtooth functions are computed from a natural logarithm and fractional parts.
 * They enhance both phase and modules changes with discontinuities and dark cells.
 * </p>
 */
public class SawToothPhaseModuleValue extends DomainColoring {

    /** Minimum brightness. */
    private final double minBrightness;

    /** Maximum brightness. */
    private final double maxBrightness;

    /** Number of lines per cycle. */
    private final int nbLines;

    /** Simple constructor.
     * @param saturation constant saturation
     * @param minBrightness minimum brightness
     * @param maxBrightness maximum brightness
     * @param nbLines number of lines per cycle
     */
    protected SawToothPhaseModuleValue(final double saturation,
                                       final double minBrightness, final double maxBrightness,
                                       final int nbLines) {
        super(saturation);
        this.minBrightness = minBrightness;
        this.maxBrightness = maxBrightness;
        this.nbLines       = nbLines;
    }

    /** Compute periodic fractional brightness.
     * @param x continuous value
     * @param s scaling factor
     * @return fractional brightness
     */
    private double fractionalBrightness(final double x, final double s) {
        double f = x * nbLines / s;
        return minBrightness + (maxBrightness - minBrightness) * (f - FastMath.floor(f));
    }

    /** {@inheritDoc} */
    @Override
    public double value(final Complex z) {
        final double module = z.norm();
        final double bM = fractionalBrightness(FastMath.log(module), MathUtils.TWO_PI);
        final double bP = fractionalBrightness(hue(z), 1.0);
        return bM * bP;
    }

}
