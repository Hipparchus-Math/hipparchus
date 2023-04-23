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
package org.hipparchus.geometry.spherical.oned;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.partitioning.Region.Location;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;


/** This class represents an arc on a circle.
 * @see ArcsSet
 */
public class Arc {

    /** The lower angular bound of the arc. */
    private final double lower;

    /** The upper angular bound of the arc. */
    private final double upper;

    /** Middle point of the arc. */
    private final double middle;

    /** Tolerance below which angles are considered identical. */
    private final double tolerance;

    /** Simple constructor.
     * <p>
     * If either {@code lower} is equals to {@code upper} or
     * the interval exceeds \( 2 \pi \), the arc is considered
     * to be the full circle and its initial defining boundaries
     * will be forgotten. {@code lower} is not allowed to be
     * greater than {@code upper} (an exception is thrown in this case).
     * {@code lower} will be canonicalized between 0 and \( 2 \pi \), and
     * upper shifted accordingly, so the {@link #getInf()} and {@link #getSup()}
     * may not return the value used at instance construction.
     * </p>
     * @param lower lower angular bound of the arc
     * @param upper upper angular bound of the arc
     * @param tolerance tolerance below which angles are considered identical
     * @exception MathIllegalArgumentException if lower is greater than upper
     * or tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    public Arc(final double lower, final double upper, final double tolerance)
        throws MathIllegalArgumentException {
        Sphere1D.checkTolerance(tolerance);
        this.tolerance = tolerance;
        if (Precision.equals(lower, upper, 0) || (upper - lower) >= MathUtils.TWO_PI) {
            // the arc must cover the whole circle
            this.lower  = 0;
            this.upper  = MathUtils.TWO_PI;
            this.middle = FastMath.PI;
        } else  if (lower <= upper) {
            this.lower  = MathUtils.normalizeAngle(lower, FastMath.PI);
            this.upper  = this.lower + (upper - lower);
            this.middle = 0.5 * (this.lower + this.upper);
        } else {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ENDPOINTS_NOT_AN_INTERVAL,
                                                lower, upper, true);
        }
    }

    /** Get the lower angular bound of the arc.
     * @return lower angular bound of the arc,
     * always between 0 and \( 2 \pi \)
     */
    public double getInf() {
        return lower;
    }

    /** Get the upper angular bound of the arc.
     * @return upper angular bound of the arc,
     * always between {@link #getInf()} and {@link #getInf()} \( + 2 \pi \)
     */
    public double getSup() {
        return upper;
    }

    /** Get the angular size of the arc.
     * @return angular size of the arc
     */
    public double getSize() {
        return upper - lower;
    }

    /** Get the barycenter of the arc.
     * @return barycenter of the arc
     */
    public double getBarycenter() {
        return middle;
    }

    /** Get the tolerance below which angles are considered identical.
     * @return tolerance below which angles are considered identical
     */
    public double getTolerance() {
        return tolerance;
    }

    /** Check a point with respect to the arc.
     * @param point point to check
     * @return a code representing the point status: either {@link
     * Location#INSIDE}, {@link Location#OUTSIDE} or {@link Location#BOUNDARY}
     */
    public Location checkPoint(final double point) {
        final double normalizedPoint = MathUtils.normalizeAngle(point, middle);
        if (normalizedPoint < lower - tolerance || normalizedPoint > upper + tolerance) {
            return Location.OUTSIDE;
        } else if (normalizedPoint > lower + tolerance && normalizedPoint < upper - tolerance) {
            return Location.INSIDE;
        } else {
            return (getSize() >= MathUtils.TWO_PI - tolerance) ? Location.INSIDE : Location.BOUNDARY;
        }
    }

    /**
     * Get the distance (arc length) from a point to the edge of the arc.
     *
     * <p>This method does not use {@link #getTolerance()}.
     *
     * @param point to test.
     * @return offset, negative if the point is inside the arc, positive if it is outside
     * the arc, or zero if {@code point} is {@link #getInf()} or {@link #getSup()}.
     */
    public double getOffset(final double point) {
        final double normalizedPoint = MathUtils.normalizeAngle(point, middle);
        if (normalizedPoint < middle) {
            return lower - normalizedPoint;
        } else {
            return normalizedPoint - upper;
        }
    }

    /**
     * Get the distance (arc length) from a point to the edge of the arc.
     *
     * <p>This method does not use {@link #getTolerance()}.
     *
     * @param point to test.
     * @return offset, negative if the point is inside the arc, positive if it is outside
     * the arc, or zero if {@code point} is {@link #getInf()} or {@link #getSup()}.
     */
    public double getOffset(final S1Point point) {
        return getOffset(point.getAlpha());
    }

}
