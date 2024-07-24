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
package org.hipparchus.transform;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.Random;

/**
 * Abstract test for classes implementing the {@link RealTransformer} interface.
 * This abstract test handles the automatic generation of random data of various
 * sizes. For each generated data array, actual values (returned by the
 * transformer to be tested) are compared to expected values, returned by the
 * {@link #transform(double[], TransformType)} (to be implemented by the user:
 * a naive method may be used). Methods are also provided to test that invalid
 * parameters throw the expected exceptions.
 *
 */
public abstract class RealTransformerAbstractTest <T> {

    /** The common seed of all random number generators used in this test. */
    private final static long SEED = 20110119L;

    /**
     * Returns a new instance of the {@link RealTransformer} to be tested.
     *
     * @return a the transformer to be tested
     */
    abstract RealTransformer createRealTransformer();

    /**
     * Returns an invalid data size. Transforms with this data size should
     * trigger a {@link MathIllegalArgumentException}.
     *
     * @param i the index of the invalid data size ({@code 0 <= i <}
     * {@link #getNumberOfInvalidDataSizes()}
     * @return an invalid data size
     */
    abstract int getInvalidDataSize(int i);

    /**
     * Returns the total number of invalid data sizes to be tested. If data
     * array of any
     * size can be handled by the {@link RealTransformer} to be tested, this
     * method should return {@code 0}.
     *
     * @return the total number of invalid data sizes
     */
    abstract int getNumberOfInvalidDataSizes();

    /**
     * Returns the total number of valid data sizes to be tested.
     *
     * @return the total number of valid data sizes
     */
    abstract int getNumberOfValidDataSizes();

    /**
     * Returns the expected relative accuracy for data arrays of size
     * {@code getValidDataSize(i)}.
     *
     * @param i the index of the valid data size
     * @return the expected relative accuracy
     */
    abstract double getRelativeTolerance(int i);

    /**
     * Returns a valid data size. This method allows for data arrays of various
     * sizes to be automatically tested (by allowing multiple values of the
     * specified index).
     *
     * @param i the index of the valid data size ({@code 0 <= i <}
     * {@link #getNumberOfValidDataSizes()}
     * @return a valid data size
     */
    abstract int getValidDataSize(int i);

    /**
     * Returns a function for the accuracy check of
     * {@link RealTransformer#transform(UnivariateFunction, double, double, int, TransformType)}.
     * This function should be valid. In other words, none of the above methods
     * should throw an exception when passed this function.
     *
     * @return a valid function
     */
    abstract UnivariateFunction getValidFunction();

    /**
     * Returns a sampling lower bound for the accuracy check of
     * {@link RealTransformer#transform(UnivariateFunction, double, double, int, TransformType)}.
     * This lower bound should be valid. In other words, none of the above
     * methods should throw an exception when passed this bound.
     *
     * @return a valid lower bound
     */
    abstract double getValidLowerBound();

    /**
     * Returns a sampling upper bound for the accuracy check of
     * {@link RealTransformer#transform(UnivariateFunction, double, double, int, TransformType)}.
     * This upper bound should be valid. In other words, none of the above
     * methods should throw an exception when passed this bound.
     *
     * @return a valid bound
     */
    abstract double getValidUpperBound();

    /**
     * Returns the expected transform of the specified real data array.
     *
     * @param x the real data array to be transformed
     * @param type the type of transform (forward, inverse) to be performed
     * @return the expected transform
     */
    abstract double[] transform(double[] x, TransformType type);

    /**
     * Initialize subclass parameters.
     * @param normalization normalization to use
     */
    abstract void init(T normalization);

    /*
     * Check of preconditions.
     */

    /**
     * {@link RealTransformer#transform(double[], TransformType)} should throw a
     * {@link MathIllegalArgumentException} if data size is invalid.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testTransformRealInvalidDataSize(T normalization) {
        init(normalization);
        final TransformType[] type = TransformType.values();
        final RealTransformer transformer = createRealTransformer();
        for (int i = 0; i < getNumberOfInvalidDataSizes(); i++) {
            final int n = getInvalidDataSize(i);
            for (int j = 0; j < type.length; j++) {
                try {
                    transformer.transform(createRealData(n), type[j]);
                    Assertions.fail(type[j] + ", " + n);
                } catch (MathIllegalArgumentException e) {
                    // Expected: do nothing
                }
            }
        }
    }

    /**
     * {@link RealTransformer#transform(UnivariateFunction, double, double, int, TransformType)}
     * should throw a {@link MathIllegalArgumentException} if number of samples
     * is invalid.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testTransformFunctionInvalidDataSize(T normalization) {
        init(normalization);
        final TransformType[] type = TransformType.values();
        final RealTransformer transformer = createRealTransformer();
        final UnivariateFunction f = getValidFunction();
        final double a = getValidLowerBound();
        final double b = getValidUpperBound();
        for (int i = 0; i < getNumberOfInvalidDataSizes(); i++) {
            final int n = getInvalidDataSize(i);
            for (int j = 0; j < type.length; j++) {
                try {
                    transformer.transform(f, a, b, n, type[j]);
                    Assertions.fail(type[j] + ", " + n);
                } catch (MathIllegalArgumentException e) {
                    // Expected: do nothing
                }
            }
        }
    }

    /**
     * {@link RealTransformer#transform(UnivariateFunction, double, double, int, TransformType)}
     * should throw a {@link MathIllegalArgumentException} if number of samples
     * is not strictly positive.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testTransformFunctionNotStrictlyPositiveNumberOfSamples(T normalization) {
        init(normalization);
        final TransformType[] type = TransformType.values();
        final RealTransformer transformer = createRealTransformer();
        final UnivariateFunction f = getValidFunction();
        final double a = getValidLowerBound();
        final double b = getValidUpperBound();
        for (int i = 0; i < getNumberOfValidDataSizes(); i++) {
            final int n = getValidDataSize(i);
            for (int j = 0; j < type.length; j++) {
                try {
                    transformer.transform(f, a, b, -n, type[j]);
                    Assertions.fail(type[j] + ", " + (-n));
                } catch (MathIllegalArgumentException e) {
                    // Expected: do nothing
                }
            }
        }
    }

    /**
     * {@link RealTransformer#transform(UnivariateFunction, double, double, int, TransformType)}
     * should throw a {@link MathIllegalArgumentException} if sampling bounds are
     * not correctly ordered.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testTransformFunctionInvalidBounds(T normalization) {
        init(normalization);
        final TransformType[] type = TransformType.values();
        final RealTransformer transformer = createRealTransformer();
        final UnivariateFunction f = getValidFunction();
        final double a = getValidLowerBound();
        final double b = getValidUpperBound();
        for (int i = 0; i < getNumberOfValidDataSizes(); i++) {
            final int n = getValidDataSize(i);
            for (int j = 0; j < type.length; j++) {
                try {
                    transformer.transform(f, b, a, n, type[j]);
                    Assertions.fail(type[j] + ", " + b + ", " + a);
                } catch (MathIllegalArgumentException e) {
                    // Expected: do nothing
                }
            }
        }
    }

    /*
     * Accuracy tests of transform of valid data.
     */

    /**
     * Accuracy check of {@link RealTransformer#transform(double[], TransformType)}.
     * For each valid data size returned by
     * {@link #getValidDataSize(int) getValidDataSize(i)},
     * a random data array is generated with
     * {@link #createRealData(int) createRealData(i)}. The actual
     * transform is computed and compared to the expected transform, return by
     * {@link #transform(double[], TransformType)}. Actual and expected values
     * should be equal to within the relative error returned by
     * {@link #getRelativeTolerance(int) getRelativeTolerance(i)}.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testTransformReal(T normalization) {
        init(normalization);
        final TransformType[] type = TransformType.values();
        for (int i = 0; i < getNumberOfValidDataSizes(); i++) {
            final int n = getValidDataSize(i);
            final double tol = getRelativeTolerance(i);
            for (int j = 0; j < type.length; j++) {
                doTestTransformReal(n, tol, type[j]);
            }
        }
    }

    /**
     * Accuracy check of
     * {@link RealTransformer#transform(UnivariateFunction, double, double, int, TransformType)}.
     * For each valid data size returned by
     * {@link #getValidDataSize(int) getValidDataSize(i)},
     * the {@link UnivariateFunction} returned by {@link #getValidFunction()} is
     * sampled. The actual transform is computed and compared to the expected
     * transform, return by {@link #transform(double[], TransformType)}. Actual
     * and expected values should be equal to within the relative error returned
     * by {@link #getRelativeTolerance(int) getRelativeTolerance(i)}.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testTransformFunction(T normalization) {
        init(normalization);
        final TransformType[] type = TransformType.values();
        for (int i = 0; i < getNumberOfValidDataSizes(); i++) {
            final int n = getValidDataSize(i);
            final double tol = getRelativeTolerance(i);
            for (int j = 0; j < type.length; j++) {
                doTestTransformFunction(n, tol, type[j]);
            }
        }
    }

    /*
     * Utility methods.
     */

    /**
     * Returns a random array of doubles. Random generator always uses the same
     * seed.
     *
     * @param n the size of the array to be returned
     * @return a random array of specified size
     */
    double[] createRealData(final int n) {
        final Random random = new Random(SEED);
        final double[] data = new double[n];
        for (int i = 0; i < n; i++) {
            data[i] = 2.0 * random.nextDouble() - 1.0;
        }
        return data;
    }

    /*
     * The tests per se.
     */

    private void doTestTransformReal(final int n, final double tol,
        final TransformType type) {
        final RealTransformer transformer = createRealTransformer();
        final double[] x = createRealData(n);
        final double[] expected = transform(x, type);
        final double[] actual = transformer.transform(x, type);
        for (int i = 0; i < n; i++) {
            final String msg = String.format("%d, %d", n, i);
            final double delta = tol * FastMath.abs(expected[i]);
            Assertions.assertEquals(expected[i], actual[i], delta, msg);
        }
    }

    private void doTestTransformFunction(final int n, final double tol,
        final TransformType type) {
        final RealTransformer transformer = createRealTransformer();
        final UnivariateFunction f = getValidFunction();
        final double a = getValidLowerBound();
        final double b = getValidUpperBound();
        final double[] x = createRealData(n);
        for (int i = 0; i < n; i++) {
            final double t = a + i * (b - a) / n;
            x[i] = f.value(t);
        }
        final double[] expected = transform(x, type);
        final double[] actual = transformer.transform(f, a, b, n, type);
        for (int i = 0; i < n; i++) {
            final String msg = String.format("%d, %d", n, i);
            final double delta = tol * FastMath.abs(expected[i]);
            Assertions.assertEquals(expected[i], actual[i], delta, msg);
        }
    }
}
