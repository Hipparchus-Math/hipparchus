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
package org.hipparchus.geometry.euclidean.threed;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/** The class represent lines in a three dimensional space.

 * <p>Each oriented line is intrinsically associated with an abscissa
 * which is a coordinate on the line. The point at abscissa 0 is the
 * orthogonal projection of the origin on the line, another equivalent
 * way to express this is to say that it is the point of the line
 * which is closest to the origin. Abscissa increases in the line
 * direction.</p>
 * @param <T> the type of the field elements
 */
public class FieldLine<T extends CalculusFieldElement<T>> {

    /** Line direction. */
    private FieldVector3D<T> direction;

    /** Line point closest to the origin. */
    private FieldVector3D<T> zero;

    /** Tolerance below which points are considered identical. */
    private final double tolerance;

    /** Build a line from two points.
     * @param p1 first point belonging to the line (this can be any point)
     * @param p2 second point belonging to the line (this can be any point, different from p1)
     * @param tolerance tolerance below which points are considered identical
     * @exception MathIllegalArgumentException if the points are equal
     */
    public FieldLine(final FieldVector3D<T> p1, final FieldVector3D<T> p2, final double tolerance)
        throws MathIllegalArgumentException {
        reset(p1, p2);
        this.tolerance = tolerance;
    }

    /** Copy constructor.
     * <p>The created instance is completely independent from the
     * original instance, it is a deep copy.</p>
     * @param line line to copy
     */
    public FieldLine(final FieldLine<T> line) {
        this.direction = line.direction;
        this.zero      = line.zero;
        this.tolerance = line.tolerance;
    }

    /** Reset the instance as if built from two points.
     * @param p1 first point belonging to the line (this can be any point)
     * @param p2 second point belonging to the line (this can be any point, different from p1)
     * @exception MathIllegalArgumentException if the points are equal
     */
    public void reset(final FieldVector3D<T> p1, final FieldVector3D<T> p2)
        throws MathIllegalArgumentException {
        final FieldVector3D<T> delta = p2.subtract(p1);
        final T norm2 = delta.getNormSq();
        if (norm2.getReal() == 0.0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_NORM);
        }
        this.direction = new FieldVector3D<>(norm2.sqrt().reciprocal(), delta);
        zero = new FieldVector3D<>(norm2.getField().getOne(), p1,
                                   p1.dotProduct(delta).negate().divide(norm2), delta);
    }

    /** Get the tolerance below which points are considered identical.
     * @return tolerance below which points are considered identical
     */
    public double getTolerance() {
        return tolerance;
    }

    /** Get a line with reversed direction.
     * @return a new instance, with reversed direction
     */
    public FieldLine<T> revert() {
        final FieldLine<T> reverted = new FieldLine<>(this);
        reverted.direction = reverted.direction.negate();
        return reverted;
    }

    /** Get the normalized direction vector.
     * @return normalized direction vector
     */
    public FieldVector3D<T> getDirection() {
        return direction;
    }

    /** Get the line point closest to the origin.
     * @return line point closest to the origin
     */
    public FieldVector3D<T> getOrigin() {
        return zero;
    }

    /** Get the abscissa of a point with respect to the line.
     * <p>The abscissa is 0 if the projection of the point and the
     * projection of the frame origin on the line are the same
     * point.</p>
     * @param point point to check
     * @return abscissa of the point
     */
    public T getAbscissa(final FieldVector3D<T> point) {
        return point.subtract(zero).dotProduct(direction);
    }

    /** Get the abscissa of a point with respect to the line.
     * <p>The abscissa is 0 if the projection of the point and the
     * projection of the frame origin on the line are the same
     * point.</p>
     * @param point point to check
     * @return abscissa of the point
     */
    public T getAbscissa(final Vector3D point) {
        return zero.subtract(point).dotProduct(direction).negate();
    }

    /** Get one point from the line.
     * @param abscissa desired abscissa for the point
     * @return one point belonging to the line, at specified abscissa
     */
    public FieldVector3D<T> pointAt(final T abscissa) {
        return new FieldVector3D<T>(abscissa.getField().getOne(), zero,
                                    abscissa, direction);
    }

    /** Get one point from the line.
     * @param abscissa desired abscissa for the point
     * @return one point belonging to the line, at specified abscissa
     */
    public FieldVector3D<T> pointAt(final double abscissa) {
        return new FieldVector3D<T>(1, zero, abscissa, direction);
    }

    /** Check if the instance is similar to another line.
     * <p>Lines are considered similar if they contain the same
     * points. This does not mean they are equal since they can have
     * opposite directions.</p>
     * @param line line to which instance should be compared
     * @return true if the lines are similar
     */
    public boolean isSimilarTo(final FieldLine<T> line) {
        final double angle = FieldVector3D.angle(direction, line.direction).getReal();
        return ((angle < tolerance) || (angle > (FastMath.PI - tolerance))) && contains(line.zero);
    }

    /** Check if the instance contains a point.
     * @param p point to check
     * @return true if p belongs to the line
     */
    public boolean contains(final FieldVector3D<T> p) {
        return distance(p).getReal() < tolerance;
    }

    /** Check if the instance contains a point.
     * @param p point to check
     * @return true if p belongs to the line
     */
    public boolean contains(final Vector3D p) {
        return distance(p).getReal() < tolerance;
    }

    /** Compute the distance between the instance and a point.
     * @param p to check
     * @return distance between the instance and the point
     */
    public T distance(final FieldVector3D<T> p) {
        final FieldVector3D<T> d = p.subtract(zero);
        final FieldVector3D<T> n = new FieldVector3D<>(zero.getX().getField().getOne(), d,
                                                       d.dotProduct(direction).negate(), direction);
        return n.getNorm();
    }

    /** Compute the distance between the instance and a point.
     * @param p to check
     * @return distance between the instance and the point
     */
    public T distance(final Vector3D p) {
        final FieldVector3D<T> d = zero.subtract(p).negate();
        final FieldVector3D<T> n = new FieldVector3D<>(zero.getX().getField().getOne(), d,
                                                       d.dotProduct(direction).negate(), direction);
        return n.getNorm();
    }

    /** Compute the shortest distance between the instance and another line.
     * @param line line to check against the instance
     * @return shortest distance between the instance and the line
     */
    public T distance(final FieldLine<T> line) {

        final FieldVector3D<T> normal = FieldVector3D.crossProduct(direction, line.direction);
        final T n = normal.getNorm();
        if (n.getReal() < Precision.SAFE_MIN) {
            // lines are parallel
            return distance(line.zero);
        }

        // signed separation of the two parallel planes that contains the lines
        final T offset = line.zero.subtract(zero).dotProduct(normal).divide(n);

        return offset.abs();

    }

    /** Compute the point of the instance closest to another line.
     * @param line line to check against the instance
     * @return point of the instance closest to another line
     */
    public FieldVector3D<T> closestPoint(final FieldLine<T> line) {

        final T cos = direction.dotProduct(line.direction);
        final T n = cos.multiply(cos).subtract(1).negate();
        if (n.getReal() < Precision.EPSILON) {
            // the lines are parallel
            return zero;
        }

        final FieldVector3D<T> delta0 = line.zero.subtract(zero);
        final T a                     = delta0.dotProduct(direction);
        final T b                     = delta0.dotProduct(line.direction);

        return new FieldVector3D<T>(a.getField().getOne(), zero,
                                    a.subtract(b.multiply(cos)).divide(n), direction);

    }

    /** Get the intersection point of the instance and another line.
     * @param line other line
     * @return intersection point of the instance and the other line
     * or null if there are no intersection points
     */
    public FieldVector3D<T> intersection(final FieldLine<T> line) {
        final FieldVector3D<T> closest = closestPoint(line);
        return line.contains(closest) ? closest : null;
    }

}
