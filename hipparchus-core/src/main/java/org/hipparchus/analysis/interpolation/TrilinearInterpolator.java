/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Interpolate grid data using tri-linear interpolation.
 *
 * @since 4.1
 */
public class TrilinearInterpolator implements TrivariateGridInterpolator {

    /**
     * Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 4.1
     */
    public TrilinearInterpolator() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** {@inheritDoc} */
    @Override
    public TrilinearInterpolatingFunction interpolate(final double[] xval, final double[] yval, final double[] zval,
                                                      final double[][][] fval)
        throws MathIllegalArgumentException {
        return new TrilinearInterpolatingFunction(xval, yval, zval, fval);
    }

}
