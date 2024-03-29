/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.analysis.interpolation;

import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.SinCos;

/**
 * Utility class for the {@link MicrosphereProjectionInterpolator} algorithm.
 * For 2D interpolation, this class constructs the microsphere as a series of
 * evenly spaced facets (rather than generating random normals as in the
 * base implementation).
 *
 */
public class InterpolatingMicrosphere2D extends InterpolatingMicrosphere {
    /** Space dimension. */
    private static final int DIMENSION = 2;

    /**
     * Create a sphere from vectors regularly sampled around a circle.
     *
     * @param size Number of surface elements of the sphere.
     * @param maxDarkFraction Maximum fraction of the facets that can be dark.
     * If the fraction of "non-illuminated" facets is larger, no estimation
     * of the value will be performed, and the {@code background} value will
     * be returned instead.
     * @param darkThreshold Value of the illumination below which a facet is
     * considered dark.
     * @param background Value returned when the {@code maxDarkFraction}
     * threshold is exceeded.
     * @throws org.hipparchus.exception.MathIllegalArgumentException
     * if {@code size <= 0}.
     * @throws org.hipparchus.exception.MathIllegalArgumentException if
     * {@code darkThreshold < 0}.
     * @throws org.hipparchus.exception.MathIllegalArgumentException if
     * {@code maxDarkFraction} does not belong to the interval {@code [0, 1]}.
     */
    public InterpolatingMicrosphere2D(int size,
                                      double maxDarkFraction,
                                      double darkThreshold,
                                      double background) {
        super(DIMENSION, size, maxDarkFraction, darkThreshold, background);

        // Generate the microsphere normals.
        for (int i = 0; i < size; i++) {
            final double angle   = i * MathUtils.TWO_PI / size;
            final SinCos scAngle = FastMath.sinCos(angle);

            add(new double[] { scAngle.cos(),
                               scAngle.sin() },
                false);
        }
    }

    /**
     * Copy constructor.
     *
     * @param other Instance to copy.
     */
    protected InterpolatingMicrosphere2D(InterpolatingMicrosphere2D other) {
        super(other);
    }

    /**
     * Perform a copy.
     *
     * @return a copy of this instance.
     */
    @Override
    public InterpolatingMicrosphere2D copy() {
        return new InterpolatingMicrosphere2D(this);
    }
}
