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

import java.io.Serializable;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.geometry.Space;
import org.hipparchus.util.FastMath;

/**
 * This class implements a one-dimensional sphere (i.e. a circle).
 * <p>
 * We use here the topologists definition of the 1-sphere (see
 * <a href="http://mathworld.wolfram.com/Sphere.html">Sphere</a> on
 * MathWorld), i.e. the 1-sphere is the one-dimensional closed curve
 * defined in 2D as x<sup>2</sup>+y<sup>2</sup>=1.
 * </p>
 */
public class Sphere1D implements Serializable, Space {

    /** Smallest tolerance that can be managed.
     * <p>
     * Tolerances smaller than this value will generate exceptions.
     * </p>
     * @since 1.4
     */
    public static final double SMALLEST_TOLERANCE = FastMath.ulp(2 * FastMath.PI);

    /** Serializable version identifier. */
    private static final long serialVersionUID = 20131218L;

    /** Private constructor for the singleton.
     */
    private Sphere1D() {
    }

    /** Get the unique instance.
     * @return the unique instance
     */
    public static Sphere1D getInstance() {
        return LazyHolder.INSTANCE;
    }

    /** Check tolerance against {@link #SMALLEST_TOLERANCE}.
     * @param tolerance tolerance to check
     * @exception MathIllegalArgumentException if tolerance is smaller
     * than {@link #SMALLEST_TOLERANCE}
     */
    public static void checkTolerance(final double tolerance)
        throws MathIllegalArgumentException {
        if (tolerance < Sphere1D.SMALLEST_TOLERANCE) {
            throw new MathIllegalArgumentException(LocalizedGeometryFormats.TOO_SMALL_TOLERANCE,
                                                   tolerance, "Sphere1D.SMALLEST_TOLERANCE", SMALLEST_TOLERANCE);
        }
    }

    /** {@inheritDoc} */
    @Override
    public int getDimension() {
        return 1;
    }

    /** {@inheritDoc}
     * <p>
     * As the 1-dimension sphere does not have proper sub-spaces,
     * this method always throws a {@link NoSubSpaceException}
     * </p>
     * @return nothing
     * @throws NoSubSpaceException in all cases
     */
    @Override
    public Space getSubSpace() throws NoSubSpaceException {
        throw new NoSubSpaceException();
    }

    // CHECKSTYLE: stop HideUtilityClassConstructor
    /** Holder for the instance.
     * <p>We use here the Initialization On Demand Holder Idiom.</p>
     */
    private static class LazyHolder {
        /** Cached field instance. */
        private static final Sphere1D INSTANCE = new Sphere1D();
    }
    // CHECKSTYLE: resume HideUtilityClassConstructor

    /** Handle deserialization of the singleton.
     * @return the singleton instance
     */
    private Object readResolve() {
        // return the singleton instance
        return LazyHolder.INSTANCE;
    }

    /** Specialized exception for inexistent sub-space.
     * <p>
     * This exception is thrown when attempting to get the sub-space of a one-dimensional space
     * </p>
     */
    public static class NoSubSpaceException extends MathRuntimeException {

        /** Serializable UID. */
        private static final long serialVersionUID = 20140225L;

        /** Simple constructor.
         */
        public NoSubSpaceException() {
            super(LocalizedGeometryFormats.NOT_SUPPORTED_IN_DIMENSION_N, 1);
        }

    }

}
