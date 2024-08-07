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
import org.hipparchus.util.FastMath;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.Arrays;
import java.util.Collection;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test case for fast cosine transformer.
 * <p>
 * FCT algorithm is exact, the small tolerance number is used only to account
 * for round-off errors.
 *
 */
public final class FastCosineTransformerTest
    extends RealTransformerAbstractTest<DctNormalization> {

    private DctNormalization normalization;

    private int[] invalidDataSize;

    private double[] relativeTolerance;

    private int[] validDataSize;

    public void initFastCosineTransformerTest(final DctNormalization normalization) {
        this.normalization = normalization;
        this.validDataSize = new int[] {
            2, 3, 5, 9, 17, 33, 65, 129
        };
        this.invalidDataSize = new int[] {
            128
        };
        this.relativeTolerance = new double[] {
            1E-15, 1E-15, 1E-14, 1E-13, 1E-13, 1E-12, 1E-11, 1E-10
        };
    }

    @Override
    void init(DctNormalization normalization) {
        initFastCosineTransformerTest(normalization);
    }

    /**
     * Returns an array containing {@code true, false} in order to check both
     * standard and orthogonal DCTs.
     *
     * @return an array of parameters for this parameterized test
     */
    public static Collection<Object[]> data() {
        final DctNormalization[] normalization = DctNormalization.values();
        final Object[][] data = new DctNormalization[normalization.length][1];
        for (int i = 0; i < normalization.length; i++){
            data[i][0] = normalization[i];
        }
        return Arrays.asList(data);
    }

    @Override
    RealTransformer createRealTransformer() {
        return new FastCosineTransformer(normalization);
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
        final double[] cos = new double[2 * (n - 1)];
        for (int i = 0; i < cos.length; i++) {
            cos[i] = FastMath.cos(FastMath.PI * i / (n - 1.0));
        }
        int sgn = 1;
        for (int j = 0; j < n; j++) {
            double yj = 0.5 * (x[0] + sgn * x[n - 1]);
            for (int i = 1; i < n - 1; i++) {
                yj += x[i] * cos[(i * j) % cos.length];
            }
            y[j] = yj;
            sgn *= -1;
        }
        final double s;
        if (type == TransformType.FORWARD) {
            if (normalization == DctNormalization.STANDARD_DCT_I) {
                s = 1.0;
            } else if (normalization == DctNormalization.ORTHOGONAL_DCT_I) {
                s = FastMath.sqrt(2.0 / (n - 1.0));
            } else {
                throw new MathIllegalStateException(LocalizedCoreFormats.ILLEGAL_STATE);
            }
        } else if (type == TransformType.INVERSE) {
            if (normalization == DctNormalization.STANDARD_DCT_I) {
                s = 2.0 / (n - 1.0);
            } else if (normalization == DctNormalization.ORTHOGONAL_DCT_I) {
                s = FastMath.sqrt(2.0 / (n - 1.0));
            } else {
                throw new MathIllegalStateException(LocalizedCoreFormats.ILLEGAL_STATE);
            }
        } else {
            /*
             * Should never occur. This clause is a safeguard in case other
             * types are used to TransformType (which should not be done).
             */
            throw new MathIllegalStateException(LocalizedCoreFormats.ILLEGAL_STATE);
        }
        TransformUtils.scaleArray(y, s);
        return y;
    }

    /*
     * Additional tests.
     */

    /** Test of parameters for the transformer. */
    @MethodSource("data")
    @ParameterizedTest
    void testParameters(final DctNormalization normalization)
        throws Exception {
        initFastCosineTransformerTest(normalization);
        UnivariateFunction f = new Sin();
        FastCosineTransformer transformer;
        transformer = new FastCosineTransformer(DctNormalization.STANDARD_DCT_I);

        try {
            // bad interval
            transformer.transform(f, 1, -1, 65, TransformType.FORWARD);
            fail("Expecting MathIllegalArgumentException - bad interval");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // bad samples number
            transformer.transform(f, -1, 1, 1, TransformType.FORWARD);
            fail("Expecting MathIllegalArgumentException - bad samples number");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // bad samples number
            transformer.transform(f, -1, 1, 64, TransformType.FORWARD);
            fail("Expecting MathIllegalArgumentException - bad samples number");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    /** Test of transformer for the sine function. */
    @MethodSource("data")
    @ParameterizedTest
    void testSinFunction(final DctNormalization normalization) {
        initFastCosineTransformerTest(normalization);
        UnivariateFunction f = new Sin();
        FastCosineTransformer transformer;
        transformer = new FastCosineTransformer(DctNormalization.STANDARD_DCT_I);
        double min;
        double max;
        double[] result;
        double tolerance = 1E-12;
        int N = 9;

        double[] expected =
            {
                0.0, 3.26197262739567, 0.0, -2.17958042710327, 0.0,
                -0.648846697642915, 0.0, -0.433545502649478, 0.0
            };
        min = 0.0;
        max = 2.0 * FastMath.PI * N / (N - 1);
        result = transformer.transform(f, min, max, N, TransformType.FORWARD);
        for (int i = 0; i < N; i++) {
            assertEquals(expected[i], result[i], tolerance);
        }

        min = -FastMath.PI;
        max = FastMath.PI * (N + 1) / (N - 1);
        result = transformer.transform(f, min, max, N, TransformType.FORWARD);
        for (int i = 0; i < N; i++) {
            assertEquals(-expected[i], result[i], tolerance);
        }
    }

    /** Test of transformer for the ad hoc data. */
    @MethodSource("data")
    @ParameterizedTest
    void testAdHocData(final DctNormalization normalization) {
        initFastCosineTransformerTest(normalization);
        FastCosineTransformer transformer;
        transformer = new FastCosineTransformer(DctNormalization.STANDARD_DCT_I);
        double[] result;
        double tolerance = 1E-12;

        double[] x = {
                        0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0
        };
        double[] y =
                        {
                                        172.0, -105.096569476353, 27.3137084989848, -12.9593152353742,
                                        8.0, -5.78585076868676, 4.68629150101524, -4.15826451958632,
                                        4.0
                        };

        result = transformer.transform(x, TransformType.FORWARD);
        for (int i = 0; i < result.length; i++) {
            assertEquals(y[i], result[i], tolerance);
        }

        result = transformer.transform(y, TransformType.INVERSE);
        for (int i = 0; i < result.length; i++) {
            assertEquals(x[i], result[i], tolerance);
        }

        TransformUtils.scaleArray(x, FastMath.sqrt(0.5 * (x.length - 1)));

        transformer = new FastCosineTransformer(DctNormalization.ORTHOGONAL_DCT_I);
        result = transformer.transform(y, TransformType.FORWARD);
        for (int i = 0; i < result.length; i++) {
            assertEquals(x[i], result[i], tolerance);
        }

        result = transformer.transform(x, TransformType.INVERSE);
        for (int i = 0; i < result.length; i++) {
            assertEquals(y[i], result[i], tolerance);
        }
    }
}
