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
package org.hipparchus.clustering.distance;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * Calculates the Canberra distance between two points.
 *
 */
public class CanberraDistance implements DistanceMeasure {

    /** Serializable version identifier. */
    private static final long serialVersionUID = -6972277381587032228L;

    /** Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 3.0
     */
    public CanberraDistance() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** {@inheritDoc} */
    @Override
    public double compute(double[] a, double[] b)
    throws MathIllegalArgumentException {
        MathArrays.checkEqualLength(a, b);
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            final double num = FastMath.abs(a[i] - b[i]);
            final double denom = FastMath.abs(a[i]) + FastMath.abs(b[i]);
            sum += num == 0.0 && denom == 0.0 ? 0.0 : num / denom;
        }
        return sum;
    }

}
