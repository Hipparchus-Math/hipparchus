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
package org.hipparchus.geometry.euclidean.oned;

import org.hipparchus.geometry.partitioning.Hyperplane;

/** This class represents a 1D oriented hyperplane.
 * <p>An hyperplane in 1D is a simple point, its orientation being a
 * boolean.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 */
public class OrientedPoint implements Hyperplane<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> {

    /** Vector location. */
    private final Vector1D location;

    /** Orientation. */
    private boolean direct;

    /** Tolerance below which points are considered to belong to the hyperplane. */
    private final double tolerance;

    /** Simple constructor.
     * @param location location of the hyperplane
     * @param direct if true, the plus side of the hyperplane is towards
     * abscissas greater than {@code location}
     * @param tolerance tolerance below which points are considered to belong to the hyperplane
     */
    public OrientedPoint(final Vector1D location, final boolean direct, final double tolerance) {
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
    public OrientedPoint copySelf() {
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public double getOffset(final Vector1D point) {
        final double delta = point.getX() - location.getX();
        return direct ? delta : -delta;
    }

    /** {@inheritDoc} */
    @Override
    public Vector1D moveToOffset(final Vector1D point, final double offset) {
        return new Vector1D(location.getX() + (direct ? offset : -offset));
    }

    /** {@inheritDoc} */
    @Override
    public Vector1D arbitraryPoint() {
        return location;
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
    public SubOrientedPoint wholeHyperplane() {
        return new SubOrientedPoint(this, null);
    }

    /** {@inheritDoc}.
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
    public SubOrientedPoint emptyHyperplane() {
        return new SubOrientedPoint(this, null);
    }

    /** Build a region covering the whole space.
     * @return a region containing the instance (really an {@link
     * IntervalsSet IntervalsSet} instance)
     */
    @Override
    public IntervalsSet wholeSpace() {
        return new IntervalsSet(tolerance);
    }

    /** {@inheritDoc} */
    @Override
    public boolean sameOrientationAs(final OrientedPoint other) {
        return direct == other.direct;
    }

    /** {@inheritDoc}
     */
    @Override
    public Vector1D project(Vector1D point) {
        return location;
    }

    /** {@inheritDoc}
     */
    @Override
    public double getTolerance() {
        return tolerance;
    }

    /** Get the hyperplane location on the real line.
     * @return the hyperplane location
     */
    public Vector1D getLocation() {
        return location;
    }

    /** Check if the hyperplane orientation is direct.
     * @return true if the plus side of the hyperplane is towards
     * abscissae greater than hyperplane location
     */
    public boolean isDirect() {
        return direct;
    }

    /** Revert the instance.
     */
    public void revertSelf() {
        direct = !direct;
    }

}
