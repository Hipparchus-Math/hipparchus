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
package org.hipparchus.complex;

import java.io.Serializable;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.SinCos;

/**
 * A helper class for the computation and caching of the {@code n}-th roots
 * of unity.
 */
public class RootsOfUnity implements Serializable {

    /** Serializable version id. */
    private static final long serialVersionUID = 20120201L;

    /** Number of roots of unity. */
    private int omegaCount;

    /** Real part of the roots. */
    private double[] omegaReal;

    /**
     * Imaginary part of the {@code n}-th roots of unity, for positive values
     * of {@code n}. In this array, the roots are stored in counter-clockwise
     * order.
     */
    private double[] omegaImaginaryCounterClockwise;

    /**
     * Imaginary part of the {@code n}-th roots of unity, for negative values
     * of {@code n}. In this array, the roots are stored in clockwise order.
     */
    private double[] omegaImaginaryClockwise;

    /**
     * {@code true} if {@link #computeRoots(int)} was called with a positive
     * value of its argument {@code n}. In this case, counter-clockwise ordering
     * of the roots of unity should be used.
     */
    private boolean isCounterClockWise;

    /**
     * Build an engine for computing the {@code n}-th roots of unity.
     */
    public RootsOfUnity() {
        omegaCount = 0;
        omegaReal = null;
        omegaImaginaryCounterClockwise = null;
        omegaImaginaryClockwise = null;
        isCounterClockWise = true;
    }

    /**
     * Returns {@code true} if {@link #computeRoots(int)} was called with a
     * positive value of its argument {@code n}. If {@code true}, then
     * counter-clockwise ordering of the roots of unity should be used.
     *
     * @return {@code true} if the roots of unity are stored in counter-clockwise order
     * @throws MathIllegalStateException if no roots of unity have been computed yet
     */
    public boolean isCounterClockWise()
            throws MathIllegalStateException {
        synchronized (this) {
            if (omegaCount == 0) {
                throw new MathIllegalStateException(LocalizedCoreFormats.ROOTS_OF_UNITY_NOT_COMPUTED_YET);
            }
            return isCounterClockWise;
        }
    }

    /**
     * Computes the {@code n}-th roots of unity.
     * <p>
     * The roots are stored in {@code omega[]}, such that {@code omega[k] = w ^ k},
     * where {@code k = 0, ..., n - 1}, {@code w = exp(2 * pi * i / n)} and
     * {@code i = sqrt(-1)}.
     * <p>
     * Note that {@code n} can be positive of negative
     * <ul>
     * <li>{@code abs(n)} is always the number of roots of unity.</li>
     * <li>If {@code n > 0}, then the roots are stored in counter-clockwise order.</li>
     * <li>If {@code n < 0}, then the roots are stored in clockwise order.</li>
     * </ul>
     *
     * @param n the (signed) number of roots of unity to be computed
     * @throws MathIllegalArgumentException if {@code n = 0}
     */
    public void computeRoots(int n) throws MathIllegalArgumentException {

        if (n == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.CANNOT_COMPUTE_0TH_ROOT_OF_UNITY);
        }

        synchronized (this) {
            isCounterClockWise = n > 0;

            // avoid repetitive calculations
            final int absN = FastMath.abs(n);

            if (absN == omegaCount) {
                return;
            }

            // calculate everything from scratch
            final double t  = 2.0 * FastMath.PI / absN;
            final SinCos sc = FastMath.sinCos(t);
            omegaReal = new double[absN];
            omegaImaginaryCounterClockwise = new double[absN];
            omegaImaginaryClockwise = new double[absN];
            omegaReal[0] = 1.0;
            omegaImaginaryCounterClockwise[0] = 0.0;
            omegaImaginaryClockwise[0] = 0.0;
            for (int i = 1; i < absN; i++) {
                omegaReal[i] = omegaReal[i - 1] * sc.cos() -
                                omegaImaginaryCounterClockwise[i - 1] * sc.sin();
                omegaImaginaryCounterClockwise[i] = omegaReal[i - 1] * sc.sin() +
                                omegaImaginaryCounterClockwise[i - 1] * sc.cos();
                omegaImaginaryClockwise[i] = -omegaImaginaryCounterClockwise[i];
            }
            omegaCount = absN;
        }
    }

    /**
     * Get the real part of the {@code k}-th {@code n}-th root of unity.
     *
     * @param k index of the {@code n}-th root of unity
     * @return real part of the {@code k}-th {@code n}-th root of unity
     * @throws MathIllegalStateException if no roots of unity have been computed yet
     * @throws MathIllegalArgumentException if {@code k} is out of range
     */
    public double getReal(int k)
            throws MathIllegalArgumentException, MathIllegalStateException {

        synchronized (this) {
            if (omegaCount == 0) {
                throw new MathIllegalStateException(LocalizedCoreFormats.ROOTS_OF_UNITY_NOT_COMPUTED_YET);
            }
            if ((k < 0) || (k >= omegaCount)) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE_ROOT_OF_UNITY_INDEX,
                                                       k, 0, omegaCount - 1);
            }

            return omegaReal[k];
        }
    }

    /**
     * Get the imaginary part of the {@code k}-th {@code n}-th root of unity.
     *
     * @param k index of the {@code n}-th root of unity
     * @return imaginary part of the {@code k}-th {@code n}-th root of unity
     * @throws MathIllegalStateException if no roots of unity have been computed yet
     * @throws MathIllegalArgumentException if {@code k} is out of range
     */
    public double getImaginary(int k)
            throws MathIllegalArgumentException, MathIllegalStateException {

        synchronized (this) {
            if (omegaCount == 0) {
                throw new MathIllegalStateException(LocalizedCoreFormats.ROOTS_OF_UNITY_NOT_COMPUTED_YET);
            }
            if ((k < 0) || (k >= omegaCount)) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE_ROOT_OF_UNITY_INDEX,
                                                       k, 0, omegaCount - 1);
            }

            return isCounterClockWise ?
                   omegaImaginaryCounterClockwise[k] :
                   omegaImaginaryClockwise[k];
        }
    }

    /**
     * Returns the number of roots of unity currently stored.
     * <p>
     * If {@link #computeRoots(int)} was called with {@code n}, then this method
     * returns {@code abs(n)}. If no roots of unity have been computed yet, this
     * method returns 0.
     *
     * @return the number of roots of unity currently stored
     */
    public int getNumberOfRoots() {
        synchronized (this) {
            return omegaCount;
        }
    }
}
