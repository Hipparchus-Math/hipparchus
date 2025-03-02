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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.partitioning.Hyperplane;


/** This class represents a 1D oriented hyperplane on the circle.
 * <p>An hyperplane on the 1-sphere is an angle with an orientation.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 */
public class LimitAngle implements Hyperplane<Sphere1D, S1Point, LimitAngle, SubLimitAngle> {

    /** Angle location. */
    private final S1Point location;

    /** Orientation. */
    private final boolean direct;

    /** Tolerance below which angles are considered identical. */
    private final double tolerance;

    /** Simple constructor.
     * @param location location of the hyperplane
     * @param direct if true, the plus side of the hyperplane is towards
     * angles greater than {@code location}
     * @param tolerance tolerance below which angles are considered identical
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    public LimitAngle(final S1Point location, final boolean direct, final double tolerance)
        throws MathIllegalArgumentException {
        Sphere1D.checkTolerance(tolerance);
        this.location  = location;
        this.direct    = direct;
        this.tolerance = tolerance;
    }

    /** Copy the instance.
     * <p>Since instances are immutable, this method directly returns
     * the instance.</p>
     * @return the instance itself
     */
    @Override
    public LimitAngle copySelf() {
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public double getOffset(final S1Point point) {
        final double delta = point.getAlpha() - location.getAlpha();
        return direct ? delta : -delta;
    }

    /** {@inheritDoc} */
    @Override
    public S1Point moveToOffset(final S1Point point, final double offset) {
        return new S1Point(location.getAlpha() + (direct ? offset : -offset));
    }

    /** {@inheritDoc} */
    @Override
    public S1Point arbitraryPoint() {
        return location;
    }

    /** Check if the hyperplane orientation is direct.
     * @return true if the plus side of the hyperplane is towards
     * angles greater than hyperplane location
     */
    public boolean isDirect() {
        return direct;
    }

    /** Get the reverse of the instance.
     * <p>Get a limit angle with reversed orientation with respect to the
     * instance. A new object is built, the instance is untouched.</p>
     * @return a new limit angle, with orientation opposite to the instance orientation
     */
    public LimitAngle getReverse() {
        return new LimitAngle(location, !direct, tolerance);
    }

    /** Build a region covering the whole hyperplane.
     * <p>Since this class represent zero dimension spaces which does
     * not have lower dimension sub-spaces, this method returns a dummy
     * implementation of a {@link
     * org.hipparchus.geometry.partitioning.SubHyperplane SubHyperplane}.
     * This implementation is only used to allow the {@link
     * org.hipparchus.geometry.partitioning.SubHyperplane
     * SubHyperplane} class implementation to work properly, it should
     * <em>not</em> be used otherwise.</p>
     * @return a dummy sub hyperplane
     */
    @Override
    public SubLimitAngle wholeHyperplane() {
        return new SubLimitAngle(this, null);
    }

    /** {@inheritDoc}
     * <p>Since this class represent zero dimension spaces which does
     * not have lower dimension sub-spaces, this method returns a dummy
     * implementation of a {@link
     * org.hipparchus.geometry.partitioning.SubHyperplane SubHyperplane}.
     * This implementation is only used to allow the {@link
     * org.hipparchus.geometry.partitioning.SubHyperplane
     * SubHyperplane} class implementation to work properly, it should
     * <em>not</em> be used otherwise.</p>
     */
    @Override
    public SubLimitAngle emptyHyperplane() {
        return new SubLimitAngle(this, null);
    }

    /** Build a region covering the whole space.
     * @return a region containing the instance (really an {@link
     * ArcsSet IntervalsSet} instance)
     */
    @Override
    public ArcsSet wholeSpace() {
        return new ArcsSet(tolerance);
    }

    /** {@inheritDoc} */
    @Override
    public boolean sameOrientationAs(final LimitAngle other) {
        return direct == other.direct;
    }

    /** Get the hyperplane location on the circle.
     * @return the hyperplane location
     */
    public S1Point getLocation() {
        return location;
    }

    /** {@inheritDoc} */
    @Override
    public S1Point project(S1Point point) {
        return location;
    }

    /** {@inheritDoc} */
    @Override
    public double getTolerance() {
        return tolerance;
    }

}
