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

/** Base class for domain coloring.
 * <p>All methods have the following features:</p>
 * <ul>
 *   <li>hue represents phase (red for 0, then orange, yellow, green,
 *       blue at π, indigo, violet and back to red)</li>
 *   <li>saturation is constant</li>
 * </ul>
 * <p>Their differences lie on what value represents.</p>
 */
public abstract class DomainColoring {

    /** Constant saturation. */
    private final double saturation;

    /** Simple constructor.
     * @param saturation constant saturation
     */
    protected DomainColoring(final double saturation) {
        this.saturation = saturation;
    }

    /** Continuous hue.
     * @param z complex value
     * @return continuous hue
     */
    public double hue(final Complex z) {
        final double phase =  FastMath.PI + FastMath.atan2(-z.getImaginaryPart(), -z.getRealPart());
        return phase / MathUtils.TWO_PI;
    }

    /** Get saturation for a complex value.
     * @param z complex value
     * @return saturation
     */
    public double saturation(Complex z) {
        return saturation;
    }

    /** Get value for a complex value.
     * @param z complex value
     * @return value
     */
    protected abstract double value(Complex z);

}
