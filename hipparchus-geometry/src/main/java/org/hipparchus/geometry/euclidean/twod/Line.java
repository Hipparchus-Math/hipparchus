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
package org.hipparchus.geometry.euclidean.twod;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.geometry.euclidean.oned.Euclidean1D;
import org.hipparchus.geometry.euclidean.oned.IntervalsSet;
import org.hipparchus.geometry.euclidean.oned.OrientedPoint;
import org.hipparchus.geometry.euclidean.oned.SubOrientedPoint;
import org.hipparchus.geometry.euclidean.oned.Vector1D;
import org.hipparchus.geometry.partitioning.Embedding;
import org.hipparchus.geometry.partitioning.Hyperplane;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.geometry.partitioning.Transform;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.SinCos;

/** This class represents an oriented line in the 2D plane.

 * <p>An oriented line can be defined either by prolongating a line
 * segment between two points past these points, or by one point and
 * an angular direction (in trigonometric orientation).</p>

 * <p>Since it is oriented the two half planes at its two sides are
 * unambiguously identified as a left half plane and a right half
 * plane. This can be used to identify the interior and the exterior
 * in a simple way by local properties only when part of a line is
 * used to define part of a polygon boundary.</p>

 * <p>A line can also be used to completely define a reference frame
 * in the plane. It is sufficient to select one specific point in the
 * line (the orthogonal projection of the original reference frame on
 * the line) and to use the unit vector in the line direction and the
 * orthogonal vector oriented from left half plane to right half
 * plane. We define two coordinates by the process, the
 * <em>abscissa</em> along the line, and the <em>offset</em> across
 * the line. All points of the plane are uniquely identified by these
 * two coordinates. The line is the set of points at zero offset, the
 * left half plane is the set of points with negative offsets and the
 * right half plane is the set of points with positive offsets.</p>

 */
public class Line
    implements Hyperplane<Euclidean2D, Vector2D, Line, SubLine>,
               Embedding<Euclidean2D, Vector2D, Euclidean1D, Vector1D> {

    /** Angle with respect to the abscissa axis. */
    private double angle;

    /** Cosine of the line angle. */
    private double cos;

    /** Sine of the line angle. */
    private double sin;

    /** Offset of the frame origin. */
    private double originOffset;

    /** Tolerance below which points are considered identical. */
    private final double tolerance;

    /** Reverse line. */
    private Line reverse;

    /** Build a line from two points.
     * <p>The line is oriented from p1 to p2</p>
     * @param p1 first point
     * @param p2 second point
     * @param tolerance tolerance below which points are considered identical
     */
    public Line(final Vector2D p1, final Vector2D p2, final double tolerance) {
        reset(p1, p2);
        this.tolerance = tolerance;
    }

    /** Build a line from a point and an angle.
     * @param p point belonging to the line
     * @param angle angle of the line with respect to abscissa axis
     * @param tolerance tolerance below which points are considered identical
     */
    public Line(final Vector2D p, final double angle, final double tolerance) {
        reset(p, angle);
        this.tolerance = tolerance;
    }

    /** Build a line from its internal characteristics.
     * @param angle angle of the line with respect to abscissa axis
     * @param cos cosine of the angle
     * @param sin sine of the angle
     * @param originOffset offset of the origin
     * @param tolerance tolerance below which points are considered identical
     */
    private Line(final double angle, final double cos, final double sin,
                 final double originOffset, final double tolerance) {
        this.angle        = angle;
        this.cos          = cos;
        this.sin          = sin;
        this.originOffset = originOffset;
        this.tolerance    = tolerance;
        this.reverse      = null;
    }

    /** Copy constructor.
     * <p>The created instance is completely independent of the
     * original instance, it is a deep copy.</p>
     * @param line line to copy
     */
    public Line(final Line line) {
        angle        = MathUtils.normalizeAngle(line.angle, FastMath.PI);
        cos          = line.cos;
        sin          = line.sin;
        originOffset = line.originOffset;
        tolerance    = line.tolerance;
        reverse      = null;
    }

    /** {@inheritDoc} */
    @Override
    public Line copySelf() {
        return new Line(this);
    }

    /** Reset the instance as if built from two points.
     * <p>The line is oriented from p1 to p2</p>
     * @param p1 first point
     * @param p2 second point
     */
    public void reset(final Vector2D p1, final Vector2D p2) {
        unlinkReverse();
        final double dx = p2.getX() - p1.getX();
        final double dy = p2.getY() - p1.getY();
        final double d = FastMath.hypot(dx, dy);
        if (d == 0.0) {
            angle        = 0.0;
            cos          = 1.0;
            sin          = 0.0;
            originOffset = p1.getY();
        } else {
            angle        = FastMath.PI + FastMath.atan2(-dy, -dx);
            cos          = dx / d;
            sin          = dy / d;
            originOffset = MathArrays.linearCombination(p2.getX(), p1.getY(), -p1.getX(), p2.getY()) / d;
        }
    }

    /** Reset the instance as if built from a line and an angle.
     * @param p point belonging to the line
     * @param alpha angle of the line with respect to abscissa axis
     */
    public void reset(final Vector2D p, final double alpha) {
        unlinkReverse();
        final SinCos sinCos = FastMath.sinCos(alpha);
        this.angle   = MathUtils.normalizeAngle(alpha, FastMath.PI);
        cos          = sinCos.cos();
        sin          = sinCos.sin();
        originOffset = MathArrays.linearCombination(cos, p.getY(), -sin, p.getX());
    }

    /** Revert the instance.
     */
    public void revertSelf() {
        unlinkReverse();
        if (angle < FastMath.PI) {
            angle += FastMath.PI;
        } else {
            angle -= FastMath.PI;
        }
        cos          = -cos;
        sin          = -sin;
        originOffset = -originOffset;
    }

    /** Unset the link between an instance and its reverse.
     */
    private void unlinkReverse() {
        if (reverse != null) {
            reverse.reverse = null;
        }
        reverse = null;
    }

    /** Get the reverse of the instance.
     * <p>Get a line with reversed orientation with respect to the
     * instance.</p>
     * <p>
     * As long as neither the instance nor its reverse are modified
     * (i.e. as long as none of the {@link #reset(Vector2D, Vector2D)},
     * {@link #reset(Vector2D, double)}, {@link #revertSelf()},
     * {@link #setAngle(double)} or {@link #setOriginOffset(double)}
     * methods are called), then the line and its reverse remain linked
     * together so that {@code line.getReverse().getReverse() == line}.
     * When one of the line is modified, the link is deleted as both
     * instance becomes independent.
     * </p>
     * @return a new line, with orientation opposite to the instance orientation
     */
    public Line getReverse() {
        if (reverse == null) {
            reverse = new Line((angle < FastMath.PI) ? (angle + FastMath.PI) : (angle - FastMath.PI),
                               -cos, -sin, -originOffset, tolerance);
            reverse.reverse = this;
        }
        return reverse;
    }

    /** {@inheritDoc} */
    @Override
    public Vector1D toSubSpace(final Vector2D point) {
        return new Vector1D(MathArrays.linearCombination(cos, point.getX(), sin, point.getY()));
    }

    /** {@inheritDoc} */
    @Override
    public Vector2D toSpace(final Vector1D point) {
        final double abscissa = point.getX();
        return new Vector2D(MathArrays.linearCombination(abscissa, cos, -originOffset, sin),
                            MathArrays.linearCombination(abscissa, sin,  originOffset, cos));
    }

    /** Get the intersection point of the instance and another line.
     * @param other other line
     * @return intersection point of the instance and the other line
     * or null if there are no intersection points
     */
    public Vector2D intersection(final Line other) {
        final double d = MathArrays.linearCombination(sin, other.cos, -other.sin, cos);
        if (FastMath.abs(d) < tolerance) {
            return null;
        }
        return new Vector2D(MathArrays.linearCombination(cos, other.originOffset, -other.cos, originOffset) / d,
                            MathArrays.linearCombination(sin, other.originOffset, -other.sin, originOffset) / d);
    }

    /** {@inheritDoc}
     */
    @Override
    public Vector2D project(Vector2D point) {
        return toSpace(toSubSpace(point));
    }

    /** {@inheritDoc}
     */
    @Override
    public double getTolerance() {
        return tolerance;
    }

    /** {@inheritDoc} */
    @Override
    public SubLine wholeHyperplane() {
        return new SubLine(this, new IntervalsSet(tolerance));
    }

    /** {@inheritDoc} */
    @Override
    public SubLine emptyHyperplane() {
        final RegionFactory<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> factory = new RegionFactory<>();
        return new SubLine(this, factory.getComplement(new IntervalsSet(tolerance)));
    }

    /** Build a region covering the whole space.
     * @return a region containing the instance (really a {@link
     * PolygonsSet PolygonsSet} instance)
     */
    @Override
    public PolygonsSet wholeSpace() {
        return new PolygonsSet(tolerance);
    }

    /** Get the offset (oriented distance) of a parallel line.
     * <p>This method should be called only for parallel lines otherwise
     * the result is not meaningful.</p>
     * <p>The offset is 0 if both lines are the same, it is
     * positive if the line is on the right side of the instance and
     * negative if it is on the left side, according to its natural
     * orientation.</p>
     * @param line line to check
     * @return offset of the line
     */
    public double getOffset(final Line line) {
        return originOffset +
               (MathArrays.linearCombination(cos, line.cos, sin, line.sin) > 0 ? -line.originOffset : line.originOffset);
    }

    /** {@inheritDoc} */
    @Override
    public double getOffset(final Vector2D point) {
        return MathArrays.linearCombination(sin, point.getX(), -cos, point.getY(), 1.0, originOffset);
    }

    /** {@inheritDoc} */
    @Override
    public Vector2D moveToOffset(final Vector2D point, final double offset) {
        final double delta = offset - getOffset(point);
        return new Vector2D(point.getX() + delta * sin, point.getY() - delta * cos);
    }

    /** {@inheritDoc} */
    @Override
    public Vector2D arbitraryPoint() {
        return getPointAt(Vector1D.ZERO, 0);
    }

    /** {@inheritDoc} */
    @Override
    public boolean sameOrientationAs(final Line other) {
        return MathArrays.linearCombination(sin, other.sin, cos, other.cos) >= 0.0;
    }

    /** Get one point from the plane.
     * @param abscissa desired abscissa for the point
     * @param offset desired offset for the point
     * @return one point in the plane, with given abscissa and offset
     * relative to the line
     */
    public Vector2D getPointAt(final Vector1D abscissa, final double offset) {
        final double x       = abscissa.getX();
        final double dOffset = offset - originOffset;
        return new Vector2D(MathArrays.linearCombination(x, cos,  dOffset, sin),
                            MathArrays.linearCombination(x, sin, -dOffset, cos));
    }

    /** Check if the line contains a point.
     * @param p point to check
     * @return true if p belongs to the line
     */
    public boolean contains(final Vector2D p) {
        return FastMath.abs(getOffset(p)) < tolerance;
    }

    /** Compute the distance between the instance and a point.
     * <p>This is a shortcut for invoking FastMath.abs(getOffset(p)),
     * and provides consistency with what is in the
     * org.hipparchus.geometry.euclidean.threed.Line class.</p>
     *
     * @param p to check
     * @return distance between the instance and the point
     */
    public double distance(final Vector2D p) {
        return FastMath.abs(getOffset(p));
    }

    /** Check the instance is parallel to another line.
     * @param line other line to check
     * @return true if the instance is parallel to the other line
     * (they can have either the same or opposite orientations)
     */
    public boolean isParallelTo(final Line line) {
        return FastMath.abs(MathArrays.linearCombination(sin, line.cos, -cos, line.sin)) < tolerance;
    }

    /** Translate the line to force it passing by a point.
     * @param p point by which the line should pass
     */
    public void translateToPoint(final Vector2D p) {
        originOffset = MathArrays.linearCombination(cos, p.getY(), -sin, p.getX());
    }

    /** Get the angle of the line.
     * @return the angle of the line with respect to the abscissa axis
     */
    public double getAngle() {
        return MathUtils.normalizeAngle(angle, FastMath.PI);
    }

    /** Set the angle of the line.
     * @param angle new angle of the line with respect to the abscissa axis
     */
    public void setAngle(final double angle) {
        unlinkReverse();
        this.angle = MathUtils.normalizeAngle(angle, FastMath.PI);
        final SinCos sinCos = FastMath.sinCos(this.angle);
        cos        = sinCos.cos();
        sin        = sinCos.sin();
    }

    /** Get the offset of the origin.
     * @return the offset of the origin
     */
    public double getOriginOffset() {
        return originOffset;
    }

    /** Set the offset of the origin.
     * @param offset offset of the origin
     */
    public void setOriginOffset(final double offset) {
        unlinkReverse();
        originOffset = offset;
    }

    /** Get a {@link org.hipparchus.geometry.partitioning.Transform
     * Transform} embedding an affine transform.
     * @param cXX transform factor between input abscissa and output abscissa
     * @param cYX transform factor between input abscissa and output ordinate
     * @param cXY transform factor between input ordinate and output abscissa
     * @param cYY transform factor between input ordinate and output ordinate
     * @param cX1 transform addendum for output abscissa
     * @param cY1 transform addendum for output ordinate
     * @return a new transform that can be applied to either {@link
     * Vector2D Vector2D}, {@link Line Line} or {@link
     * org.hipparchus.geometry.partitioning.SubHyperplane
     * SubHyperplane} instances
     * @exception MathIllegalArgumentException if the transform is non invertible
     */
    public static Transform<Euclidean2D, Vector2D, Line, SubLine,
                            Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> getTransform(final double cXX,
                                                                                       final double cYX,
                                                                                       final double cXY,
                                                                                       final double cYY,
                                                                                       final double cX1,
                                                                                       final double cY1)
        throws MathIllegalArgumentException {
        return new LineTransform(cXX, cYX, cXY, cYY, cX1, cY1);
    }

    /** Class embedding an affine transform.
     * <p>This class is used in order to apply an affine transform to a
     * line. Using a specific object allow to perform some computations
     * on the transform only once even if the same transform is to be
     * applied to a large number of lines (for example to a large
     * polygon)./<p>
     */
    private static class LineTransform
        implements Transform<Euclidean2D, Vector2D, Line, SubLine, Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> {

        /** Transform factor between input abscissa and output abscissa. */
        private final double cXX;

        /** Transform factor between input abscissa and output ordinate. */
        private final double cYX;

        /** Transform factor between input ordinate and output abscissa. */
        private final double cXY;

        /** Transform factor between input ordinate and output ordinate. */
        private final double cYY;

        /** Transform addendum for output abscissa. */
        private final double cX1;

        /** Transform addendum for output ordinate. */
        private final double cY1;

        /** cXY * cY1 - cYY * cX1. */
        private final double c1Y;

        /** cXX * cY1 - cYX * cX1. */
        private final double c1X;

        /** cXX * cYY - cYX * cXY. */
        private final double c11;

        /** Build an affine line transform from a n {@code AffineTransform}.
         * @param cXX transform factor between input abscissa and output abscissa
         * @param cYX transform factor between input abscissa and output ordinate
         * @param cXY transform factor between input ordinate and output abscissa
         * @param cYY transform factor between input ordinate and output ordinate
         * @param cX1 transform addendum for output abscissa
         * @param cY1 transform addendum for output ordinate
         * @exception MathIllegalArgumentException if the transform is non invertible
         */
        LineTransform(final double cXX, final double cYX, final double cXY,
                      final double cYY, final double cX1, final double cY1)
            throws MathIllegalArgumentException {

            this.cXX = cXX;
            this.cYX = cYX;
            this.cXY = cXY;
            this.cYY = cYY;
            this.cX1 = cX1;
            this.cY1 = cY1;

            c1Y = MathArrays.linearCombination(cXY, cY1, -cYY, cX1);
            c1X = MathArrays.linearCombination(cXX, cY1, -cYX, cX1);
            c11 = MathArrays.linearCombination(cXX, cYY, -cYX, cXY);

            if (FastMath.abs(c11) < 1.0e-20) {
                throw new MathIllegalArgumentException(LocalizedGeometryFormats.NON_INVERTIBLE_TRANSFORM);
            }

        }

        /** {@inheritDoc} */
        @Override
        public Vector2D apply(final Vector2D point) {
            final double  x   = point.getX();
            final double  y   = point.getY();
            return new Vector2D(MathArrays.linearCombination(cXX, x, cXY, y, cX1, 1),
                                MathArrays.linearCombination(cYX, x, cYY, y, cY1, 1));
        }

        /** {@inheritDoc} */
        @Override
        public Line apply(final Line hyperplane) {
            final double rOffset = MathArrays.linearCombination(c1X, hyperplane.cos, c1Y, hyperplane.sin, c11, hyperplane.originOffset);
            final double rCos    = MathArrays.linearCombination(cXX, hyperplane.cos, cXY, hyperplane.sin);
            final double rSin    = MathArrays.linearCombination(cYX, hyperplane.cos, cYY, hyperplane.sin);
            final double inv     = 1.0 / FastMath.sqrt(rSin * rSin + rCos * rCos);
            return new Line(FastMath.PI + FastMath.atan2(-rSin, -rCos),
                            inv * rCos, inv * rSin,
                            inv * rOffset, hyperplane.tolerance);
        }

        /** {@inheritDoc} */
        @Override
        public SubOrientedPoint apply(final SubOrientedPoint sub, final Line original, final Line transformed) {
            final OrientedPoint op     = sub.getHyperplane();
            final Vector1D      newLoc = transformed.toSubSpace(apply(original.toSpace(op.getLocation())));
            return new OrientedPoint(newLoc, op.isDirect(), original.tolerance).wholeHyperplane();
        }

    }

}
