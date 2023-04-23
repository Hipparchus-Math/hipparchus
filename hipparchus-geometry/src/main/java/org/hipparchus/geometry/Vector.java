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
package org.hipparchus.geometry;

import org.hipparchus.analysis.polynomials.SmoothStepFactory;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.Blendable;

import java.text.NumberFormat;

/** This interface represents a generic vector in a vectorial space or a point in an affine space.
 * @param <S> Type of the space.
 * @param <V> Type of vector implementing this interface.
 * @see Space
 * @see Point
 */
public interface Vector<S extends Space, V extends Vector<S,V>> extends Point<S>, Blendable<Vector<S,V>> {

    /** Get the null vector of the vectorial space or origin point of the affine space.
     * @return null vector of the vectorial space or origin point of the affine space
     */
    V getZero();

    /** Get the L<sub>1</sub> norm for the vector.
     * @return L<sub>1</sub> norm for the vector
     */
    double getNorm1();

    /** Get the L<sub>2</sub> norm for the vector.
     * @return Euclidean norm for the vector
     */
    double getNorm();

    /** Get the square of the norm for the vector.
     * @return square of the Euclidean norm for the vector
     */
    double getNormSq();

    /** Get the L<sub>&infin;</sub> norm for the vector.
     * @return L<sub>&infin;</sub> norm for the vector
     */
    double getNormInf();

    /** Add a vector to the instance.
     * @param v vector to add
     * @return a new vector
     */
    V add(Vector<S,V> v);

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a new vector
     */
    V add(double factor, Vector<S,V> v);

    /** Subtract a vector from the instance.
     * @param v vector to subtract
     * @return a new vector
     */
    V subtract(Vector<S,V> v);

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a new vector
     */
    V subtract(double factor, Vector<S,V> v);

    /** Get the opposite of the instance.
     * @return a new vector which is opposite to the instance
     */
    V negate();

    /** Get a normalized vector aligned with the instance.
     * @return a new normalized vector
     * @exception MathRuntimeException if the norm is zero
     */
    default V normalize() throws MathRuntimeException{
        double s = getNorm();
        if (s == 0) {
            throw new MathRuntimeException(LocalizedGeometryFormats.CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR);
        }
        return scalarMultiply(1 / s);
    }

    /** Multiply the instance by a scalar.
     * @param a scalar
     * @return a new vector
     */
    V scalarMultiply(double a);

    /**
     * Returns true if any coordinate of this vector is infinite and none are NaN;
     * false otherwise
     * @return  true if any coordinate of this vector is infinite and none are NaN;
     * false otherwise
     */
    boolean isInfinite();

    /** Compute the distance between the instance and another vector according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNorm1()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>1</sub> norm
     */
    double distance1(Vector<S,V> v);

    /** Compute the distance between the instance and another vector according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNormInf()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>&infin;</sub> norm
     */
    double distanceInf(Vector<S,V> v);

    /** Compute the square of the distance between the instance and another vector.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).getNormSq()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the square of the distance between the instance and p
     */
    double distanceSq(Vector<S,V> v);

    /** Compute the dot-product of the instance and another vector.
     * @param v second vector
     * @return the dot product this.v
     */
    double dotProduct(Vector<S,V> v);

    /** Get a string representation of this vector.
     * @param format the custom format for components
     * @return a string representation of this vector
     */
    String toString(NumberFormat format);

    /** {@inheritDoc} */
    @Override
    default V blendArithmeticallyWith(Vector<S,V> other, double blendingValue)
            throws MathIllegalArgumentException {
        SmoothStepFactory.checkBetweenZeroAndOneIncluded(blendingValue);
        return this.scalarMultiply(1 - blendingValue).add(other.scalarMultiply(blendingValue));
    }
}
