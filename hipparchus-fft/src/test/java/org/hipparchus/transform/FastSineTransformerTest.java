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
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.analysis.function.Sinc;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.Arrays;
import java.util.Collection;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test case for fast sine transformer.
 * <p>
 * FST algorithm is exact, the small tolerance number is used only
 * to account for round-off errors.
 *
 */
public final class FastSineTransformerTest extends RealTransformerAbstractTest<DstNormalization> {

    private DstNormalization normalization;

    private int[] invalidDataSize;

    private double[] relativeTolerance;

    private int[] validDataSize;

    public void initFastSineTransformerTest(final DstNormalization normalization) {
        this.normalization = normalization;
        this.validDataSize = new int[] {
            1, 2, 4, 8, 16, 32, 64, 128
        };
        this.invalidDataSize = new int[] {
            129
        };
        this.relativeTolerance = new double[] {
            1E-15, 1E-15, 1E-14, 1E-14, 1E-13, 1E-12, 1E-11, 1E-11
        };
    }

    /**
     * Returns an array containing {@code true, false} in order to check both
     * standard and orthogonal DSTs.
     *
     * @return an array of parameters for this parameterized test
     */
    public static Collection<Object[]> data() {
        final DstNormalization[] normalization = DstNormalization.values();
        final Object[][] data = new DstNormalization[normalization.length][1];
        for (int i = 0; i < normalization.length; i++) {
            data[i][0] = normalization[i];
        }
        return Arrays.asList(data);
    }

    /**
     * {@inheritDoc}
     *
     * Overriding the default implementation allows to ensure that the first
     * element of the data set is zero.
     */
    @Override
    double[] createRealData(final int n) {
        final double[] data = super.createRealData(n);
        data[0] = 0.0;
        return data;
    }

    @Override
    RealTransformer createRealTransformer() {
        return new FastSineTransformer(normalization);
    }

    @Override
    int getInvalidDataSize(final int i) {
        return invalidDataSize[i];
    }

    @Override
    int getNumberOfInvalidDataSizes() {
        return invalidDataSize.length;
    }

    @Override
    int getNumberOfValidDataSizes() {
        return validDataSize.length;
    }

    @Override
    double getRelativeTolerance(final int i) {
        return relativeTolerance[i];
    }

    @Override
    int getValidDataSize(final int i) {
        return validDataSize[i];
    }

    @Override
    UnivariateFunction getValidFunction() {
        return new Sinc();
    }

    @Override
    double getValidLowerBound() {
        return 0.0;
    }

    @Override
    double getValidUpperBound() {
        return FastMath.PI;
    }

    @Override
    double[] transform(final double[] x, final TransformType type) {
        final int n = x.length;
        final double[] y = new double[n];
        final double[] sin = new double[2 * n];
        for (int i = 0; i < sin.length; i++) {
            sin[i] = FastMath.sin(FastMath.PI * i / n);
        }
        for (int j = 0; j < n; j++) {
            double yj = 0.0;
            for (int i = 0; i < n; i++) {
                yj += x[i] * sin[(i * j) % sin.length];
            }
            y[j] = yj;
        }
        final double s;
        if (type == TransformType.FORWARD) {
            if (normalization == DstNormalization.STANDARD_DST_I) {
                s = 1.0;
            } else if (normalization == DstNormalization.ORTHOGONAL_DST_I) {
                s = FastMath.sqrt(2.0 / n);
            } else {
                throw new MathIllegalStateException(LocalizedCoreFormats.ILLEGAL_STATE);
            }
        } else if (type == TransformType.INVERSE) {
            if (normalization == DstNormalization.STANDARD_DST_I) {
                s = 2.0 / n;
            } else if (normalization == DstNormalization.ORTHOGONAL_DST_I) {
                s = FastMath.sqrt(2.0 / n);
            } else {
                throw new MathIllegalStateException(LocalizedCoreFormats.ILLEGAL_STATE);
            }
        } else {
             // Should never occur. This clause is a safeguard in case other
             // types are used to TransformType (which should not be done).
            throw MathRuntimeException.createInternalError();
        }
        TransformUtils.scaleArray(y, s);
        return y;
    }

    @Override
    void init(DstNormalization normalization) {
        initFastSineTransformerTest(normalization);
    }

    /*
     * Additional tests.
     */
    @MethodSource("data")
    @ParameterizedTest
    void testTransformRealFirstElementNotZero(final DstNormalization normalization) {
        initFastSineTransformerTest(normalization);
        final TransformType[] type = TransformType.values();
        final double[] data = new double[] {
            1.0, 1.0, 1.0, 1.0
        };
        final RealTransformer transformer = createRealTransformer();
        for (int j = 0; j < type.length; j++) {
            try {
                transformer.transform(data, type[j]);
                fail(type[j].toString());
            } catch (MathIllegalArgumentException e) {
                // Expected: do nothing
            }
        }
    }

    /*
     * Additional (legacy) tests.
     */

    /**
     * Test of transformer for the ad hoc data.
     */
    @MethodSource("data")
    @ParameterizedTest
    void testAdHocData(final DstNormalization normalization) {
        initFastSineTransformerTest(normalization);
        FastSineTransformer transformer;
        transformer = new FastSineTransformer(DstNormalization.STANDARD_DST_I);
        double[] result;
        double tolerance = 1E-12;

        double[] x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
        double[] y = { 0.0, 20.1093579685034, -9.65685424949238,
                       5.98642305066196, -4.0, 2.67271455167720,
                      -1.65685424949238, 0.795649469518633 };

        result = transformer.transform(x, TransformType.FORWARD);
        for (int i = 0; i < result.length; i++) {
            assertEquals(y[i], result[i], tolerance);
        }

        result = transformer.transform(y, TransformType.INVERSE);
        for (int i = 0; i < result.length; i++) {
            assertEquals(x[i], result[i], tolerance);
        }

        TransformUtils.scaleArray(x, FastMath.sqrt(x.length / 2.0));
        transformer = new FastSineTransformer(DstNormalization.ORTHOGONAL_DST_I);

        result = transformer.transform(y, TransformType.FORWARD);
        for (int i = 0; i < result.length; i++) {
            assertEquals(x[i], result[i], tolerance);
        }

        result = transformer.transform(x, TransformType.INVERSE);
        for (int i = 0; i < result.length; i++) {
            assertEquals(y[i], result[i], tolerance);
        }
    }

    /**
     * Test of transformer for the sine function.
     */
    @MethodSource("data")
    @ParameterizedTest
    void testSinFunction(final DstNormalization normalization) {
        initFastSineTransformerTest(normalization);
        UnivariateFunction f = new Sin();
        FastSineTransformer transformer;
        transformer = new FastSineTransformer(DstNormalization.STANDARD_DST_I);
        double min;
        double max;
        double[] result;
        double tolerance = 1E-12;
        int N = 1 << 8;

        min = 0.0; max = 2.0 * FastMath.PI;
        result = transformer.transform(f, min, max, N, TransformType.FORWARD);
        assertEquals(N >> 1, result[2], tolerance);
        for (int i = 0; i < N; i += (i == 1 ? 2 : 1)) {
            assertEquals(0.0, result[i], tolerance);
        }

        min = -FastMath.PI; max = FastMath.PI;
        result = transformer.transform(f, min, max, N, TransformType.FORWARD);
        assertEquals(-(N >> 1), result[2], tolerance);
        for (int i = 0; i < N; i += (i == 1 ? 2 : 1)) {
            assertEquals(0.0, result[i], tolerance);
        }
    }

    /**
     * Test of parameters for the transformer.
     */
    @MethodSource("data")
    @ParameterizedTest
    void testParameters(final DstNormalization normalization) throws Exception {
        initFastSineTransformerTest(normalization);
        UnivariateFunction f = new Sin();
        FastSineTransformer transformer;
        transformer = new FastSineTransformer(DstNormalization.STANDARD_DST_I);

        try {
            // bad interval
            transformer.transform(f, 1, -1, 64, TransformType.FORWARD);
            fail("Expecting MathIllegalArgumentException - bad interval");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // bad samples number
            transformer.transform(f, -1, 1, 0, TransformType.FORWARD);
            fail("Expecting MathIllegalArgumentException - bad samples number");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // bad samples number
            transformer.transform(f, -1, 1, 100, TransformType.FORWARD);
            fail("Expecting MathIllegalArgumentException - bad samples number");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }
}
