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

/** Classical method for domain coloring.
 * <p>
 * Value represents module.
 * </p>
 */
public class ContinuousModuleValue extends DomainColoring {

    /** Simple constructor.
     * @param saturation constant saturation
     */
    protected ContinuousModuleValue(final double saturation) {
        super(saturation);
    }

    /** {@inheritDoc} */
    @Override
    public double value(final Complex z) {
        final double module = z.norm();
        return FastMath.pow(1.0 - 1.0 / (1.0 + module * module), 0.2);
    }

}
