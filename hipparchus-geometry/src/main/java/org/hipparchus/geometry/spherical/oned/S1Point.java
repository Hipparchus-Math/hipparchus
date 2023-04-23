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

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.SinCos;

/** This class represents a point on the 1-sphere.
 * <p>Instances of this class are guaranteed to be immutable.</p>
 */
public class S1Point implements Point<Sphere1D> {

   // CHECKSTYLE: stop ConstantName
    /** A vector with all coordinates set to NaN. */
    public static final S1Point NaN = new S1Point(Double.NaN, Vector2D.NaN);
    // CHECKSTYLE: resume ConstantName

    /** Serializable UID. */
    private static final long serialVersionUID = 20131218L;

    /** Azimuthal angle \( \alpha \). */
    private final double alpha;

    /** Corresponding 2D normalized vector. */
    private final Vector2D vector;

    /** Simple constructor.
     * Build a vector from its coordinates
     * @param alpha azimuthal angle \( \alpha \)
     * @see #getAlpha()
     */
    public S1Point(final double alpha) {
        this(MathUtils.normalizeAngle(alpha, FastMath.PI), buildVector(alpha));
    }

    /** Build a point from its internal components.
     * @param alpha azimuthal angle \( \alpha \)
     * @param vector corresponding vector
     */
    private S1Point(final double alpha, final Vector2D vector) {
        this.alpha  = alpha;
        this.vector = vector;
    }

    /** Get the azimuthal angle \( \alpha \).
     * @return azimuthal angle \( \alpha \)
     * @see #S1Point(double)
     */
    public double getAlpha() {
        return alpha;
    }

    /** Get the corresponding normalized vector in the 2D euclidean space.
     * @return normalized vector
     */
    public Vector2D getVector() {
        return vector;
    }

    /** {@inheritDoc} */
    @Override
    public Space getSpace() {
        return Sphere1D.getInstance();
    }

    /** {@inheritDoc} */
    @Override
    public boolean isNaN() {
        return Double.isNaN(alpha);
    }

    /** {@inheritDoc} */
    @Override
    public double distance(final Point<Sphere1D> point) {
        return distance(this, (S1Point) point);
    }

    /** Compute the distance (angular separation) between two points.
     * @param p1 first vector
     * @param p2 second vector
     * @return the angular separation between p1 and p2
     */
    public static double distance(S1Point p1, S1Point p2) {
        return Vector2D.angle(p1.vector, p2.vector);
    }

    /**
     * Test for the equality of two points on the 1-sphere.
     * <p>
     * If all coordinates of two points are exactly the same, and none are
     * {@code Double.NaN}, the two points are considered to be equal.
     * </p>
     * <p>
     * {@code NaN} coordinates are considered to affect globally the point
     * and be equals to each other - i.e, if either (or all) coordinates of the
     * point are equal to {@code Double.NaN}, the point is equal to
     * {@link #NaN}.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two points on the 1-sphere objects are equal, false if
     *         object is null, not an instance of S1Point, or
     *         not equal to this S1Point instance
     *
     */
    @Override
    public boolean equals(Object other) {

        if (this == other) {
            return true;
        }

        if (other instanceof S1Point) {
            final S1Point rhs = (S1Point) other;
            return alpha == rhs.alpha || isNaN() && rhs.isNaN();
        }

        return false;

    }

    /**
     * Test for the equality of two points on the 1-sphere.
     * <p>
     * If all coordinates of two points are exactly the same, and none are
     * {@code Double.NaN}, the two points are considered to be equal.
     * </p>
     * <p>
     * In compliance with IEEE754 handling, if any coordinates of any of the
     * two points are {@code NaN}, then the points are considered different.
     * This implies that {@link #NaN S1Point.NaN}.equals({@link #NaN S1Point.NaN})
     * returns {@code false} despite the instance is checked against itself.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two points objects are equal, false if
     *         object is null, not an instance of S1Point, or
     *         not equal to this S1Point instance
     * @since 2.1
     */
    public boolean equalsIeee754(Object other) {

        if (this == other && !isNaN()) {
            return true;
        }

        if (other instanceof S1Point) {
            final S1Point rhs = (S1Point) other;
            return alpha == rhs.alpha;
        }

        return false;

    }

    /**
     * Get a hashCode for the point.
     * <p>
     * All NaN values have the same hash code.</p>
     *
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        if (isNaN()) {
            return 542;
        }
        return 1759 * MathUtils.hash(alpha);
    }

    /**
     * Build the 2D vector corresponding to the given angle.
     * @param alpha angle
     * @return the corresponding 2D vector
     */
    private static Vector2D buildVector(final double alpha) {
        final SinCos sc = FastMath.sinCos(alpha);
        return new Vector2D(sc.cos(), sc.sin());
    }

}
