/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.geometry.euclidean.twod;

import java.text.NumberFormat;

import org.hipparchus.Field;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * This class is a re-implementation of {@link Vector2D} using {@link CalculusFieldElement}.
 * <p>Instance of this class are guaranteed to be immutable.</p>
 * @param <T> the type of the field elements
 * @since 1.6
 */
public class FieldVector2D<T extends CalculusFieldElement<T>> {

    /** Abscissa. */
    private final T x;

    /** Ordinate. */
    private final T y;

    /** Simple constructor.
     * Build a vector from its coordinates
     * @param x abscissa
     * @param y ordinate
     * @see #getX()
     * @see #getY()
     */
    public FieldVector2D(final T x, final T y) {
        this.x = x;
        this.y = y;
    }

    /** Simple constructor.
     * Build a vector from its coordinates
     * @param v coordinates array
     * @exception MathIllegalArgumentException if array does not have 2 elements
     * @see #toArray()
     */
    public FieldVector2D(final T[] v) throws MathIllegalArgumentException {
        if (v.length != 2) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   v.length, 2);
        }
        this.x = v[0];
        this.y = v[1];
    }

    /** Multiplicative constructor
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    public FieldVector2D(final T a, final FieldVector2D<T> u) {
        this.x = a.multiply(u.x);
        this.y = a.multiply(u.y);
    }

    /** Multiplicative constructor
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    public FieldVector2D(final T a, final Vector2D u) {
        this.x = a.multiply(u.getX());
        this.y = a.multiply(u.getY());
    }

    /** Multiplicative constructor
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    public FieldVector2D(final double a, final FieldVector2D<T> u) {
        this.x = u.x.multiply(a);
        this.y = u.y.multiply(a);
    }

    /** Linear constructor
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    public FieldVector2D(final T a1, final FieldVector2D<T> u1, final T a2, final FieldVector2D<T> u2) {
        final T prototype = a1;
        this.x = prototype.linearCombination(a1, u1.getX(), a2, u2.getX());
        this.y = prototype.linearCombination(a1, u1.getY(), a2, u2.getY());
    }

    /** Linear constructor.
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    public FieldVector2D(final T a1, final Vector2D u1,
                         final T a2, final Vector2D u2) {
        final T prototype = a1;
        this.x = prototype.linearCombination(u1.getX(), a1, u2.getX(), a2);
        this.y = prototype.linearCombination(u1.getY(), a1, u2.getY(), a2);
    }

    /** Linear constructor.
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    public FieldVector2D(final double a1, final FieldVector2D<T> u1,
                         final double a2, final FieldVector2D<T> u2) {
        final T prototype = u1.getX();
        this.x = prototype.linearCombination(a1, u1.getX(), a2, u2.getX());
        this.y = prototype.linearCombination(a1, u1.getY(), a2, u2.getY());
    }

    /** Linear constructor.
     * Build a vector from three other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     */
    public FieldVector2D(final T a1, final FieldVector2D<T> u1,
                         final T a2, final FieldVector2D<T> u2,
                         final T a3, final FieldVector2D<T> u3) {
        final T prototype = a1;
        this.x = prototype.linearCombination(a1, u1.getX(), a2, u2.getX(), a3, u3.getX());
        this.y = prototype.linearCombination(a1, u1.getY(), a2, u2.getY(), a3, u3.getY());
    }

    /** Linear constructor.
     * Build a vector from three other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     */
    public FieldVector2D(final T a1, final Vector2D u1,
                         final T a2, final Vector2D u2,
                         final T a3, final Vector2D u3) {
        final T prototype = a1;
        this.x = prototype.linearCombination(u1.getX(), a1, u2.getX(), a2, u3.getX(), a3);
        this.y = prototype.linearCombination(u1.getY(), a1, u2.getY(), a2, u3.getY(), a3);
    }

    /** Linear constructor.
     * Build a vector from three other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     */
    public FieldVector2D(final double a1, final FieldVector2D<T> u1,
                         final double a2, final FieldVector2D<T> u2,
                         final double a3, final FieldVector2D<T> u3) {
        final T prototype = u1.getX();
        this.x = prototype.linearCombination(a1, u1.getX(), a2, u2.getX(), a3, u3.getX());
        this.y = prototype.linearCombination(a1, u1.getY(), a2, u2.getY(), a3, u3.getY());
    }

    /** Linear constructor.
     * Build a vector from four other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3 + a4 * u4
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     * @param a4 fourth scale factor
     * @param u4 fourth base (unscaled) vector
     */
    public FieldVector2D(final T a1, final FieldVector2D<T> u1,
                         final T a2, final FieldVector2D<T> u2,
                         final T a3, final FieldVector2D<T> u3,
                         final T a4, final FieldVector2D<T> u4) {
        final T prototype = a1;
        this.x = prototype.linearCombination(a1, u1.getX(), a2, u2.getX(), a3, u3.getX(), a4, u4.getX());
        this.y = prototype.linearCombination(a1, u1.getY(), a2, u2.getY(), a3, u3.getY(), a4, u4.getY());
    }

    /** Linear constructor.
     * Build a vector from four other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3 + a4 * u4
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     * @param a4 fourth scale factor
     * @param u4 fourth base (unscaled) vector
     */
    public FieldVector2D(final T a1, final Vector2D u1,
                         final T a2, final Vector2D u2,
                         final T a3, final Vector2D u3,
                         final T a4, final Vector2D u4) {
        final T prototype = a1;
        this.x = prototype.linearCombination(u1.getX(), a1, u2.getX(), a2, u3.getX(), a3, u4.getX(), a4);
        this.y = prototype.linearCombination(u1.getY(), a1, u2.getY(), a2, u3.getY(), a3, u4.getY(), a4);
    }

    /** Linear constructor.
     * Build a vector from four other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3 + a4 * u4
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     * @param a4 fourth scale factor
     * @param u4 fourth base (unscaled) vector
     */
    public FieldVector2D(final double a1, final FieldVector2D<T> u1,
                         final double a2, final FieldVector2D<T> u2,
                         final double a3, final FieldVector2D<T> u3,
                         final double a4, final FieldVector2D<T> u4) {
        final T prototype = u1.getX();
        this.x = prototype.linearCombination(a1, u1.getX(), a2, u2.getX(), a3, u3.getX(), a4, u4.getX());
        this.y = prototype.linearCombination(a1, u1.getY(), a2, u2.getY(), a3, u3.getY(), a4, u4.getY());
    }

    /** Build a {@link FieldVector2D} from a {@link Vector2D}.
     * @param field field for the components
     * @param v vector to convert
     */
    public FieldVector2D(final Field<T> field, final Vector2D v) {
        this.x = field.getZero().add(v.getX());
        this.y = field.getZero().add(v.getY());
    }

    /** Get null vector (coordinates: 0, 0).
     * @param field field for the components
     * @return a new vector
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldVector2D<T> getZero(final Field<T> field) {
        return new FieldVector2D<>(field, Vector2D.ZERO);
    }

    /** Get first canonical vector (coordinates: 1, 0).
     * @param field field for the components
     * @return a new vector
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldVector2D<T> getPlusI(final Field<T> field) {
        return new FieldVector2D<>(field, Vector2D.PLUS_I);
    }

    /** Get opposite of the first canonical vector (coordinates: -1).
     * @param field field for the components
     * @return a new vector
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldVector2D<T> getMinusI(final Field<T> field) {
        return new FieldVector2D<>(field, Vector2D.MINUS_I);
    }

    /** Get second canonical vector (coordinates: 0, 1).
     * @param field field for the components
     * @return a new vector
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldVector2D<T> getPlusJ(final Field<T> field) {
        return new FieldVector2D<>(field, Vector2D.PLUS_J);
    }

    /** Get opposite of the second canonical vector (coordinates: 0, -1).
     * @param field field for the components
     * @return a new vector
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldVector2D<T> getMinusJ(final Field<T> field) {
        return new FieldVector2D<>(field, Vector2D.MINUS_J);
    }

    /** Get a vector with all coordinates set to NaN.
     * @param field field for the components
     * @return a new vector
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldVector2D<T> getNaN(final Field<T> field) {
        return new FieldVector2D<>(field, Vector2D.NaN);
    }

    /** Get a vector with all coordinates set to positive infinity.
     * @param field field for the components
     * @return a new vector
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldVector2D<T> getPositiveInfinity(final Field<T> field) {
        return new FieldVector2D<>(field, Vector2D.POSITIVE_INFINITY);
    }

    /** Get a vector with all coordinates set to negative infinity.
     * @param field field for the components
     * @return a new vector
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldVector2D<T> getNegativeInfinity(final Field<T> field) {
        return new FieldVector2D<>(field, Vector2D.NEGATIVE_INFINITY);
    }

    /** Get the abscissa of the vector.
     * @return abscissa of the vector
     * @see #FieldVector2D(CalculusFieldElement, CalculusFieldElement)
     */
    public T getX() {
        return x;
    }

    /** Get the ordinate of the vector.
     * @return ordinate of the vector
    * @see #FieldVector2D(CalculusFieldElement, CalculusFieldElement)
     */
    public T getY() {
        return y;
    }

    /** Get the vector coordinates as a dimension 2 array.
     * @return vector coordinates
     * @see #FieldVector2D(CalculusFieldElement[])
     */
    public T[] toArray() {
        final T[] array = MathArrays.buildArray(x.getField(), 2);
        array[0] = x;
        array[1] = y;
        return array;
    }

    /** Convert to a constant vector without extra field parts.
     * @return a constant vector
     */
    public Vector2D toVector2D() {
        return new Vector2D(x.getReal(), y.getReal());
    }

    /** Get the L<sub>1</sub> norm for the vector.
     * @return L<sub>1</sub> norm for the vector
     */
    public T getNorm1() {
        return x.abs().add(y.abs());
    }

    /** Get the L<sub>2</sub> norm for the vector.
     * @return Euclidean norm for the vector
     */
    public T getNorm() {
        // there are no cancellation problems here, so we use the straightforward formula
        return x.multiply(x).add(y.multiply(y)).sqrt();
    }

    /** Get the square of the norm for the vector.
     * @return square of the Euclidean norm for the vector
     */
    public T getNormSq() {
        // there are no cancellation problems here, so we use the straightforward formula
        return x.multiply(x).add(y.multiply(y));
    }

    /** Get the L<sub>&infin;</sub> norm for the vector.
     * @return L<sub>&infin;</sub> norm for the vector
     */
    public T getNormInf() {
        return FastMath.max(FastMath.abs(x), FastMath.abs(y));
    }

    /** Add a vector to the instance.
     * @param v vector to add
     * @return a new vector
     */
    public FieldVector2D<T> add(final FieldVector2D<T> v) {
        return new FieldVector2D<>(x.add(v.x), y.add(v.y));
    }

    /** Add a vector to the instance.
     * @param v vector to add
     * @return a new vector
     */
    public FieldVector2D<T> add(final Vector2D v) {
        return new FieldVector2D<>(x.add(v.getX()), y.add(v.getY()));
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a new vector
     */
    public FieldVector2D<T> add(final T factor, final FieldVector2D<T> v) {
        return new FieldVector2D<>(x.getField().getOne(), this, factor, v);
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a new vector
     */
    public FieldVector2D<T> add(final T factor, final Vector2D v) {
        return new FieldVector2D<>(x.add(factor.multiply(v.getX())),
                                   y.add(factor.multiply(v.getY())));
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a new vector
     */
    public FieldVector2D<T> add(final double factor, final FieldVector2D<T> v) {
        return new FieldVector2D<>(1.0, this, factor, v);
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a new vector
     */
    public FieldVector2D<T> add(final double factor, final Vector2D v) {
        return new FieldVector2D<>(x.add(factor * v.getX()),
                                   y.add(factor * v.getY()));
    }

    /** Subtract a vector from the instance.
     * @param v vector to subtract
     * @return a new vector
     */
    public FieldVector2D<T> subtract(final FieldVector2D<T> v) {
        return new FieldVector2D<>(x.subtract(v.x), y.subtract(v.y));
    }

    /** Subtract a vector from the instance.
     * @param v vector to subtract
     * @return a new vector
     */
    public FieldVector2D<T> subtract(final Vector2D v) {
        return new FieldVector2D<>(x.subtract(v.getX()), y.subtract(v.getY()));
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a new vector
     */
    public FieldVector2D<T> subtract(final T factor, final FieldVector2D<T> v) {
        return new FieldVector2D<>(x.getField().getOne(), this, factor.negate(), v);
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a new vector
     */
    public FieldVector2D<T> subtract(final T factor, final Vector2D v) {
        return new FieldVector2D<>(x.subtract(factor.multiply(v.getX())),
                                   y.subtract(factor.multiply(v.getY())));
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a new vector
     */
    public FieldVector2D<T> subtract(final double factor, final FieldVector2D<T> v) {
        return new FieldVector2D<>(1.0, this, -factor, v);
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a new vector
     */
    public FieldVector2D<T> subtract(final double factor, final Vector2D v) {
        return new FieldVector2D<>(x.subtract(factor * v.getX()),
                                   y.subtract(factor * v.getY()));
    }

    /** Get a normalized vector aligned with the instance.
     * @return a new normalized vector
     * @exception MathRuntimeException if the norm is zero
     */
    public FieldVector2D<T> normalize() throws MathRuntimeException {
        final T s = getNorm();
        if (s.getReal() == 0) {
            throw new MathRuntimeException(LocalizedGeometryFormats.CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR);
        }
        return scalarMultiply(s.reciprocal());
    }

    /** Compute the angular separation between two vectors.
     * <p>This method computes the angular separation between two
     * vectors using the dot product for well separated vectors and the
     * cross product for almost aligned vectors. This allows to have a
     * good accuracy in all cases, even for vectors very close to each
     * other.</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return angular separation between v1 and v2
     * @exception MathRuntimeException if either vector has a null norm
     */
    public static <T extends CalculusFieldElement<T>> T angle(final FieldVector2D<T> v1, final FieldVector2D<T> v2)
        throws MathRuntimeException {

        final T normProduct = v1.getNorm().multiply(v2.getNorm());
        if (normProduct.getReal() == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_NORM);
        }

        final T dot = v1.dotProduct(v2);
        final double threshold = normProduct.getReal() * 0.9999;
        if (FastMath.abs(dot.getReal()) > threshold) {
            // the vectors are almost aligned, compute using the sine
            final T n = FastMath.abs(dot.linearCombination(v1.x, v2.y, v1.y.negate(), v2.x));
            if (dot.getReal() >= 0) {
                return FastMath.asin(n.divide(normProduct));
            }
            return FastMath.asin(n.divide(normProduct)).negate().add(dot.getPi());
        }

        // the vectors are sufficiently separated to use the cosine
        return FastMath.acos(dot.divide(normProduct));

    }

    /** Compute the angular separation between two vectors.
     * <p>This method computes the angular separation between two
     * vectors using the dot product for well separated vectors and the
     * cross product for almost aligned vectors. This allows to have a
     * good accuracy in all cases, even for vectors very close to each
     * other.</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return angular separation between v1 and v2
     * @exception MathRuntimeException if either vector has a null norm
     */
    public static <T extends CalculusFieldElement<T>> T angle(final FieldVector2D<T> v1, final Vector2D v2)
        throws MathRuntimeException {

        final T normProduct = v1.getNorm().multiply(v2.getNorm());
        if (normProduct.getReal() == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_NORM);
        }

        final T dot = v1.dotProduct(v2);
        final double threshold = normProduct.getReal() * 0.9999;
        if (FastMath.abs(dot.getReal()) > threshold) {
            // the vectors are almost aligned, compute using the sine
            final T n = FastMath.abs(dot.linearCombination(v2.getY(), v1.x, v2.getX(), v1.y.negate()));
            if (dot.getReal() >= 0) {
                return FastMath.asin(n.divide(normProduct));
            }
            return FastMath.asin(n.divide(normProduct)).negate().add(dot.getPi());
        }

        // the vectors are sufficiently separated to use the cosine
        return FastMath.acos(dot.divide(normProduct));

    }

    /** Compute the angular separation between two vectors.
     * <p>This method computes the angular separation between two
     * vectors using the dot product for well separated vectors and the
     * cross product for almost aligned vectors. This allows to have a
     * good accuracy in all cases, even for vectors very close to each
     * other.</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return angular separation between v1 and v2
     * @exception MathRuntimeException if either vector has a null norm
     */
    public static <T extends CalculusFieldElement<T>> T angle(final Vector2D v1, final FieldVector2D<T> v2)
        throws MathRuntimeException {
        return angle(v2, v1);
    }

    /** Get the opposite of the instance.
     * @return a new vector which is opposite to the instance
     */
    public FieldVector2D<T> negate() {
        return new FieldVector2D<>(x.negate(), y.negate());
    }

    /** Multiply the instance by a scalar.
     * @param a scalar
     * @return a new vector
     */
    public FieldVector2D<T> scalarMultiply(final T a) {
        return new FieldVector2D<>(x.multiply(a), y.multiply(a));
    }

    /** Multiply the instance by a scalar.
     * @param a scalar
     * @return a new vector
     */
    public FieldVector2D<T> scalarMultiply(final double a) {
        return new FieldVector2D<>(x.multiply(a), y.multiply(a));
    }

    /**
     * Returns true if any coordinate of this vector is NaN; false otherwise
     * @return  true if any coordinate of this vector is NaN; false otherwise
     */
    public boolean isNaN() {
        return Double.isNaN(x.getReal()) || Double.isNaN(y.getReal());
    }

    /**
     * Returns true if any coordinate of this vector is infinite and none are NaN;
     * false otherwise
     * @return  true if any coordinate of this vector is infinite and none are NaN;
     * false otherwise
     */
    public boolean isInfinite() {
        return !isNaN() && (Double.isInfinite(x.getReal()) || Double.isInfinite(y.getReal()));
    }

    /**
     * Test for the equality of two 2D vectors.
     * <p>
     * If all coordinates of two 2D vectors are exactly the same, and none of their
     * {@link CalculusFieldElement#getReal() real part} are <code>NaN</code>, the
     * two 2D vectors are considered to be equal.
     * </p>
     * <p>
     * <code>NaN</code> coordinates are considered to affect globally the vector
     * and be equals to each other - i.e, if either (or all) real part of the
     * coordinates of the 3D vector are <code>NaN</code>, the 2D vector is <code>NaN</code>.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two 2D vector objects are equal, false if
     *         object is null, not an instance of FieldVector2D, or
     *         not equal to this FieldVector2D instance
     *
     */
    @Override
    public boolean equals(Object other) {

        if (this == other) {
            return true;
        }

        if (other instanceof FieldVector2D) {
            @SuppressWarnings("unchecked")
            final FieldVector2D<T> rhs = (FieldVector2D<T>) other;
            if (rhs.isNaN()) {
                return this.isNaN();
            }

            return x.equals(rhs.x) && y.equals(rhs.y);

        }
        return false;
    }

    /**
     * Get a hashCode for the 3D vector.
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
        return 122 * (76 * x.hashCode() +  y.hashCode());
    }

    /** Compute the distance between the instance and another vector according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNorm1()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>1</sub> norm
     */
    public T distance1(final FieldVector2D<T> v) {
        final T dx = v.x.subtract(x).abs();
        final T dy = v.y.subtract(y).abs();
        return dx.add(dy);
    }

    /** Compute the distance between the instance and another vector according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNorm1()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>1</sub> norm
     */
    public T distance1(final Vector2D v) {
        final T dx = x.subtract(v.getX()).abs();
        final T dy = y.subtract(v.getY()).abs();
        return dx.add(dy);
    }

    /** Compute the distance between the instance and another vector according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNorm()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>2</sub> norm
     */
    public T distance(final FieldVector2D<T> v) {
        final T dx = v.x.subtract(x);
        final T dy = v.y.subtract(y);
        return dx.multiply(dx).add(dy.multiply(dy)).sqrt();
    }

    /** Compute the distance between the instance and another vector according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNorm()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>2</sub> norm
     */
    public T distance(final Vector2D v) {
        final T dx = x.subtract(v.getX());
        final T dy = y.subtract(v.getY());
        return dx.multiply(dx).add(dy.multiply(dy)).sqrt();
    }

    /** Compute the distance between the instance and another vector according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNormInf()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>&infin;</sub> norm
     */
    public T distanceInf(final FieldVector2D<T> v) {
        final T dx = FastMath.abs(x.subtract(v.x));
        final T dy = FastMath.abs(y.subtract(v.y));
        return FastMath.max(dx, dy);
    }

    /** Compute the distance between the instance and another vector according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNormInf()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>&infin;</sub> norm
     */
    public T distanceInf(final Vector2D v) {
        final T dx = FastMath.abs(x.subtract(v.getX()));
        final T dy = FastMath.abs(y.subtract(v.getY()));
        return FastMath.max(dx, dy);
    }

    /** Compute the square of the distance between the instance and another vector.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNormSq()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the square of the distance between the instance and p
     */
    public T distanceSq(final FieldVector2D<T> v) {
        final T dx = v.x.subtract(x);
        final T dy = v.y.subtract(y);
        return dx.multiply(dx).add(dy.multiply(dy));
    }

    /** Compute the square of the distance between the instance and another vector.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNormSq()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the square of the distance between the instance and p
     */
    public T distanceSq(final Vector2D v) {
        final T dx = x.subtract(v.getX());
        final T dy = y.subtract(v.getY());
        return dx.multiply(dx).add(dy.multiply(dy));
    }


    /** Compute the dot-product of the instance and another vector.
     * <p>
     * The implementation uses specific multiplication and addition
     * algorithms to preserve accuracy and reduce cancellation effects.
     * It should be very accurate even for nearly orthogonal vectors.
     * </p>
     * @see MathArrays#linearCombination(double, double, double, double, double, double)
     * @param v second vector
     * @return the dot product this.v
     */
    public T dotProduct(final FieldVector2D<T> v) {
        return x.linearCombination(x, v.getX(), y, v.getY());
    }

    /** Compute the dot-product of the instance and another vector.
     * <p>
     * The implementation uses specific multiplication and addition
     * algorithms to preserve accuracy and reduce cancellation effects.
     * It should be very accurate even for nearly orthogonal vectors.
     * </p>
     * @see MathArrays#linearCombination(double, double, double, double, double, double)
     * @param v second vector
     * @return the dot product this.v
     */
    public T dotProduct(final Vector2D v) {
        return x.linearCombination(v.getX(), x, v.getY(), y);
    }

    /**
     * Compute the cross-product of the instance and the given points.
     * <p>
     * The cross product can be used to determine the location of a point
     * with regard to the line formed by (p1, p2) and is calculated as:
     * \[
     *    P = (x_2 - x_1)(y_3 - y_1) - (y_2 - y_1)(x_3 - x_1)
     * \]
     * with \(p3 = (x_3, y_3)\) being this instance.
     * <p>
     * If the result is 0, the points are collinear, i.e. lie on a single straight line L;
     * if it is positive, this point lies to the left, otherwise to the right of the line
     * formed by (p1, p2).
     *
     * @param p1 first point of the line
     * @param p2 second point of the line
     * @return the cross-product
     *
     * @see <a href="http://en.wikipedia.org/wiki/Cross_product">Cross product (Wikipedia)</a>
     */
    public T crossProduct(final FieldVector2D<T> p1, final FieldVector2D<T> p2) {
        final T x1  = p2.getX().subtract(p1.getX());
        final T y1  = getY().subtract(p1.getY());
        final T mx2 = p1.getX().subtract(getX());
        final T y2  = p2.getY().subtract(p1.getY());
        return x1.linearCombination(x1, y1, mx2, y2);
    }

    /**
     * Compute the cross-product of the instance and the given points.
     * <p>
     * The cross product can be used to determine the location of a point
     * with regard to the line formed by (p1, p2) and is calculated as:
     * \[
     *    P = (x_2 - x_1)(y_3 - y_1) - (y_2 - y_1)(x_3 - x_1)
     * \]
     * with \(p3 = (x_3, y_3)\) being this instance.
     * <p>
     * If the result is 0, the points are collinear, i.e. lie on a single straight line L;
     * if it is positive, this point lies to the left, otherwise to the right of the line
     * formed by (p1, p2).
     *
     * @param p1 first point of the line
     * @param p2 second point of the line
     * @return the cross-product
     *
     * @see <a href="http://en.wikipedia.org/wiki/Cross_product">Cross product (Wikipedia)</a>
     */
    public T crossProduct(final Vector2D p1, final Vector2D p2) {
        final double x1  = p2.getX() - p1.getX();
        final T      y1  = getY().subtract(p1.getY());
        final T      x2 = getX().subtract(p1.getX());
        final double y2  = p2.getY() - p1.getY();
        return y1.linearCombination(x1, y1, -y2, x2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNorm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    public static <T extends CalculusFieldElement<T>> T  distance1(final FieldVector2D<T> p1, final FieldVector2D<T> p2) {
        return p1.distance1(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNorm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    public static <T extends CalculusFieldElement<T>> T  distance1(final FieldVector2D<T> p1, final Vector2D p2) {
        return p1.distance1(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNorm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    public static <T extends CalculusFieldElement<T>> T  distance1(final Vector2D p1, final FieldVector2D<T> p2) {
        return p2.distance1(p1);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNorm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    public static <T extends CalculusFieldElement<T>> T distance(final FieldVector2D<T> p1, final FieldVector2D<T> p2) {
        return p1.distance(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNorm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    public static <T extends CalculusFieldElement<T>> T distance(final FieldVector2D<T> p1, final Vector2D p2) {
        return p1.distance(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNorm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    public static <T extends CalculusFieldElement<T>> T distance( final Vector2D p1, final FieldVector2D<T> p2) {
        return p2.distance(p1);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNormInf()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>&infin;</sub> norm
     */
    public static <T extends CalculusFieldElement<T>> T distanceInf(final FieldVector2D<T> p1, final FieldVector2D<T> p2) {
        return p1.distanceInf(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNormInf()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>&infin;</sub> norm
     */
    public static <T extends CalculusFieldElement<T>> T distanceInf(final FieldVector2D<T> p1, final Vector2D p2) {
        return p1.distanceInf(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNormInf()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>&infin;</sub> norm
     */
    public static <T extends CalculusFieldElement<T>> T distanceInf(final Vector2D p1, final FieldVector2D<T> p2) {
        return p2.distanceInf(p1);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNormSq()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the square of the distance between p1 and p2
     */
    public static <T extends CalculusFieldElement<T>> T distanceSq(final FieldVector2D<T> p1, final FieldVector2D<T> p2) {
        return p1.distanceSq(p2);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNormSq()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the square of the distance between p1 and p2
     */
    public static <T extends CalculusFieldElement<T>> T distanceSq(final FieldVector2D<T> p1, final Vector2D p2) {
        return p1.distanceSq(p2);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).getNormSq()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the square of the distance between p1 and p2
     */
    public static <T extends CalculusFieldElement<T>> T distanceSq(final Vector2D p1, final FieldVector2D<T> p2) {
        return p2.distanceSq(p1);
    }

    /** Compute the orientation of a triplet of points.
     * @param p first vector of the triplet
     * @param q second vector of the triplet
     * @param r third vector of the triplet
     * @param <T> the type of the field elements
     * @return a positive value if (p, q, r) defines a counterclockwise oriented
     * triangle, a negative value if (p, q, r) defines a clockwise oriented
     * triangle, and 0 if (p, q, r) are collinear or some points are equal
     * @since 1.2
     */
    public static <T extends CalculusFieldElement<T>> T orientation(final FieldVector2D<T> p, final FieldVector2D<T> q, final FieldVector2D<T> r) {
        final T prototype = p.getX();
        final T[] a = MathArrays.buildArray(prototype.getField(), 6);
        a[0] = p.getX();
        a[1] = p.getX().negate();
        a[2] = q.getX();
        a[3] = q.getX().negate();
        a[4] = r.getX();
        a[5] = r.getX().negate();
        final T[] b = MathArrays.buildArray(prototype.getField(), 6);
        b[0] = q.getY();
        b[1] = r.getY();
        b[2] = r.getY();
        b[3] = p.getY();
        b[4] = p.getY();
        b[5] = q.getY();
        return prototype.linearCombination(a, b);
    }

    /** Get a string representation of this vector.
     * @return a string representation of this vector
     */
    @Override
    public String toString() {
        return Vector2DFormat.getVector2DFormat().format(toVector2D());
    }

    /** Get a string representation of this vector.
     * @param format the custom format for components
     * @return a string representation of this vector
     */
    public String toString(final NumberFormat format) {
        return new Vector2DFormat(format).format(toVector2D());
    }

}
