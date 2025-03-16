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
package org.hipparchus.geometry.spherical.twod;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.geometry.partitioning.Embedding;
import org.hipparchus.geometry.partitioning.Hyperplane;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.geometry.partitioning.Transform;
import org.hipparchus.geometry.spherical.oned.Arc;
import org.hipparchus.geometry.spherical.oned.ArcsSet;
import org.hipparchus.geometry.spherical.oned.LimitAngle;
import org.hipparchus.geometry.spherical.oned.S1Point;
import org.hipparchus.geometry.spherical.oned.Sphere1D;
import org.hipparchus.geometry.spherical.oned.SubLimitAngle;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.SinCos;

/** This class represents an oriented great circle on the 2-sphere.

 * <p>An oriented circle can be defined by a center point. The circle
 * is the set of points that are in the normal plan the center.</p>

 * <p>Since it is oriented the two spherical caps at its two sides are
 * unambiguously identified as a left cap and a right cap. This can be
 * used to identify the interior and the exterior in a simple way by
 * local properties only when part of a line is used to define part of
 * a spherical polygon boundary.</p>

 */
public class Circle
        implements Hyperplane<Sphere2D, S2Point, Circle, SubCircle>,
                   Embedding<Sphere2D, S2Point, Sphere1D, S1Point> {

    /** Pole or circle center. */
    private Vector3D pole;

    /** First axis in the equator plane, origin of the phase angles. */
    private Vector3D x;

    /** Second axis in the equator plane, in quadrature with respect to x. */
    private Vector3D y;

    /** Tolerance below which close sub-arcs are merged together. */
    private final double tolerance;

    /** Build a great circle from its pole.
     * <p>The circle is oriented in the trigonometric direction around pole.</p>
     * @param pole circle pole
     * @param tolerance tolerance below which close sub-arcs are merged together
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    public Circle(final Vector3D pole, final double tolerance)
        throws MathIllegalArgumentException {
        Sphere2D.checkTolerance(tolerance);
        reset(pole);
        this.tolerance = tolerance;
    }

    /** Build a great circle from two non-aligned points.
     * <p>The circle is oriented from first to second point using the path smaller than π.</p>
     * @param first first point contained in the great circle
     * @param second second point contained in the great circle
     * @param tolerance tolerance below which close sub-arcs are merged together
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    public Circle(final S2Point first, final S2Point second, final double tolerance)
        throws MathIllegalArgumentException {
        Sphere2D.checkTolerance(tolerance);
        reset(first.getVector().crossProduct(second.getVector()));
        this.tolerance = tolerance;
    }

    /** Build a circle from its internal components.
     * <p>The circle is oriented in the trigonometric direction around center.</p>
     * @param pole circle pole
     * @param x first axis in the equator plane
     * @param y second axis in the equator plane
     * @param tolerance tolerance below which close sub-arcs are merged together
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    private Circle(final Vector3D pole, final Vector3D x, final Vector3D y, final double tolerance)
        throws MathIllegalArgumentException {
        Sphere2D.checkTolerance(tolerance);
        this.pole      = pole;
        this.x         = x;
        this.y         = y;
        this.tolerance = tolerance;
    }

    /** Copy constructor.
     * <p>The created instance is completely independent from the
     * original instance, it is a deep copy.</p>
     * @param circle circle to copy
     */
    public Circle(final Circle circle) {
        this(circle.pole, circle.x, circle.y, circle.tolerance);
    }

    /** {@inheritDoc} */
    @Override
    public Circle copySelf() {
        return new Circle(this);
    }

    /** Reset the instance as if built from a pole.
     * <p>The circle is oriented in the trigonometric direction around pole.</p>
     * @param newPole circle pole
     */
    public void reset(final Vector3D newPole) {
        this.pole = newPole.normalize();
        this.x    = newPole.orthogonal();
        this.y    = Vector3D.crossProduct(newPole, x).normalize();
    }

    /** Revert the instance.
     */
    public void revertSelf() {
        // x remains the same
        y    = y.negate();
        pole = pole.negate();
    }

    /** Get the reverse of the instance.
     * <p>Get a circle with reversed orientation with respect to the
     * instance. A new object is built, the instance is untouched.</p>
     * @return a new circle, with orientation opposite to the instance orientation
     */
    public Circle getReverse() {
        return new Circle(pole.negate(), x, y.negate(), tolerance);
    }

    /** {@inheritDoc} */
    @Override
    public S2Point project(S2Point point) {
        return toSpace(toSubSpace(point));
    }

    /** {@inheritDoc} */
    @Override
    public double getTolerance() {
        return tolerance;
    }

    /** {@inheritDoc}
     * @see #getPhase(Vector3D)
     */
    @Override
    public S1Point toSubSpace(final S2Point point) {
        return new S1Point(getPhase(point.getVector()));
    }

    /** Get the phase angle of a direction.
     * <p>
     * The direction may not belong to the circle as the
     * phase is computed for the meridian plane between the circle
     * pole and the direction.
     * </p>
     * @param direction direction for which phase is requested
     * @return phase angle of the direction around the circle
     * @see #toSubSpace(S2Point)
     */
    public double getPhase(final Vector3D direction) {
        return FastMath.PI + FastMath.atan2(-direction.dotProduct(y), -direction.dotProduct(x));
    }

    /** {@inheritDoc}
     * @see #getPointAt(double)
     */
    @Override
    public S2Point toSpace(final S1Point point) {
        return new S2Point(getPointAt(point.getAlpha()));
    }

    /** Get a circle point from its phase around the circle.
     * @param alpha phase around the circle
     * @return circle point on the sphere
     * @see #toSpace(S1Point)
     * @see #getXAxis()
     * @see #getYAxis()
     */
    public Vector3D getPointAt(final double alpha) {
        final SinCos sc = FastMath.sinCos(alpha);
        return new Vector3D(sc.cos(), x, sc.sin(), y);
    }

    /** Get the X axis of the circle.
     * <p>
     * This method returns the same value as {@link #getPointAt(double)
     * getPointAt(0.0)} but it does not do any computation and always
     * return the same instance.
     * </p>
     * @return an arbitrary x axis on the circle
     * @see #getPointAt(double)
     * @see #getYAxis()
     * @see #getPole()
     */
    public Vector3D getXAxis() {
        return x;
    }

    /** Get the Y axis of the circle.
     * <p>
     * This method returns the same value as {@link #getPointAt(double)
     * getPointAt(MathUtils.SEMI_PI)} but it does not do any computation and always
     * return the same instance.
     * </p>
     * @return an arbitrary y axis point on the circle
     * @see #getPointAt(double)
     * @see #getXAxis()
     * @see #getPole()
     */
    public Vector3D getYAxis() {
        return y;
    }

    /** Get the pole of the circle.
     * <p>
     * As the circle is a great circle, the pole does <em>not</em>
     * belong to it.
     * </p>
     * @return pole of the circle
     * @see #getXAxis()
     * @see #getYAxis()
     */
    public Vector3D getPole() {
        return pole;
    }

    /** Get the arc of the instance that lies inside the other circle.
     * @param other other circle
     * @return arc of the instance that lies inside the other circle
     */
    public Arc getInsideArc(final Circle other) {
        final double alpha  = getPhase(other.pole);
        return new Arc(alpha - MathUtils.SEMI_PI, alpha + MathUtils.SEMI_PI, tolerance);
    }

    /** {@inheritDoc} */
    @Override
    public SubCircle wholeHyperplane() {
        return new SubCircle(this, new ArcsSet(tolerance));
    }

    /** {@inheritDoc} */
    @Override
    public SubCircle emptyHyperplane() {
        final RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle> factory = new RegionFactory<>();
        return new SubCircle(this, factory.getComplement(new ArcsSet(tolerance)));
    }

    /** Build a region covering the whole space.
     * @return a region containing the instance (really a {@link
     * SphericalPolygonsSet SphericalPolygonsSet} instance)
     */
    @Override
    public SphericalPolygonsSet wholeSpace() {
        return new SphericalPolygonsSet(tolerance);
    }

    /** {@inheritDoc}
     * @see #getOffset(Vector3D)
     */
    @Override
    public double getOffset(final S2Point point) {
        return getOffset(point.getVector());
    }

    /** Get the offset (oriented distance) of a direction.
     * <p>The offset is defined as the angular distance between the
     * circle center and the direction minus the circle radius. It
     * is therefore 0 on the circle, positive for directions outside of
     * the cone delimited by the circle, and negative inside the cone.</p>
     * @param direction direction to check
     * @return offset of the direction
     * @see #getOffset(S2Point)
     */
    public double getOffset(final Vector3D direction) {
        return Vector3D.angle(pole, direction) - MathUtils.SEMI_PI;
    }

    /** {@inheritDoc} */
    @Override
    public S2Point moveToOffset(final S2Point point, final double offset) {
        final SinCos scOld = FastMath.sinCos(getOffset(point));
        final SinCos scNew = FastMath.sinCos(offset);
        final double ratio = scNew.cos() / scOld.cos();
        return new S2Point(new Vector3D(ratio * scOld.sin() - scNew.sin(), pole,
                                        ratio, point.getVector()));
    }

    /** {@inheritDoc} */
    @Override
    public S2Point arbitraryPoint() {
        return new S2Point(pole.orthogonal());
    }

    /** {@inheritDoc} */
    @Override
    public boolean sameOrientationAs(final Circle other) {
        return Vector3D.dotProduct(pole, other.pole) >= 0.0;
    }

    /**
     * Get the arc on this circle between two defining points. Only the point's projection
     * on the circle matters, which is computed using {@link #getPhase(Vector3D)}.
     *
     * @param a first point.
     * @param b second point.
     * @return an arc of the circle.
     */
    public Arc getArc(final S2Point a, final S2Point b) {
        final double phaseA = getPhase(a.getVector());
        double phaseB = getPhase(b.getVector());
        if (phaseB < phaseA) {
            phaseB += 2 * FastMath.PI;
        }
        return new Arc(phaseA, phaseB, tolerance);
    }

    /** Get a {@link org.hipparchus.geometry.partitioning.Transform
     * Transform} embedding a 3D rotation.
     * @param rotation rotation to use
     * @return a new transform that can be applied to either {@link
     * org.hipparchus.geometry.Point Point}, {@link Circle Line} or {@link
     * org.hipparchus.geometry.partitioning.SubHyperplane
     * SubHyperplane} instances
     */
    public static Transform<Sphere2D, S2Point, Circle, SubCircle, Sphere1D, S1Point, LimitAngle, SubLimitAngle>
        getTransform(final Rotation rotation) {
        return new CircleTransform(rotation);
    }

    /** Class embedding a 3D rotation. */
    private static class CircleTransform
        implements Transform<Sphere2D, S2Point, Circle, SubCircle, Sphere1D, S1Point, LimitAngle, SubLimitAngle> {

        /** Underlying rotation. */
        private final Rotation rotation;

        /** Build a transform from a {@code Rotation}.
         * @param rotation rotation to use
         */
        CircleTransform(final Rotation rotation) {
            this.rotation = rotation;
        }

        /** {@inheritDoc} */
        @Override
        public S2Point apply(final S2Point point) {
            return new S2Point(rotation.applyTo(point.getVector()));
        }

        /** {@inheritDoc} */
        @Override
        public Circle apply(final Circle circle) {
            return new Circle(rotation.applyTo(circle.pole),
                              rotation.applyTo(circle.x),
                              rotation.applyTo(circle.y),
                              circle.tolerance);
        }

        /** {@inheritDoc} */
        @Override
        public SubLimitAngle apply(final SubLimitAngle sub, final Circle original, final Circle transformed) {
            // as the circle is rotated, the limit angles are rotated too
            return sub;
        }

    }

}
