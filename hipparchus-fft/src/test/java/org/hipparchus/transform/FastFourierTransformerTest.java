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
import org.hipparchus.complex.Complex;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test case for fast Fourier transformer.
 * <p>
 * FFT algorithm is exact, the small tolerance number is used only
 * to account for round-off errors.
 *
 */
final class FastFourierTransformerTest {
    /** The common seed of all random number generators used in this test. */
    private final static long SEED = 20110111L;

    /*
     * Precondition checks.
     */

    @Test
    void testTransformComplexSizeNotAPowerOfTwo() {
        final int n = 127;
        final Complex[] x = createComplexData(n);
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                final FastFourierTransformer fft;
                fft = new FastFourierTransformer(norm[i]);
                try {
                    fft.transform(x, type[j]);
                    fail(norm[i] + ", " + type[j] +
                        ": MathIllegalArgumentException was expected");
                } catch (MathIllegalArgumentException e) {
                    // Expected behaviour
                }
            }
        }
    }

    @Test
    void testTransformRealSizeNotAPowerOfTwo() {
        final int n = 127;
        final double[] x = createRealData(n);
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                final FastFourierTransformer fft;
                fft = new FastFourierTransformer(norm[i]);
                try {
                    fft.transform(x, type[j]);
                    fail(norm[i] + ", " + type[j] +
                        ": MathIllegalArgumentException was expected");
                } catch (MathIllegalArgumentException e) {
                    // Expected behaviour
                }
            }
        }
    }

    @Test
    void testTransformFunctionSizeNotAPowerOfTwo() {
        final int n = 127;
        final UnivariateFunction f = new Sin();
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                final FastFourierTransformer fft;
                fft = new FastFourierTransformer(norm[i]);
                try {
                    fft.transform(f, 0.0, Math.PI, n, type[j]);
                    fail(norm[i] + ", " + type[j] +
                        ": MathIllegalArgumentException was expected");
                } catch (MathIllegalArgumentException e) {
                    // Expected behaviour
                }
            }
        }
    }

    @Test
    void testTransformFunctionNotStrictlyPositiveNumberOfSamples() {
        final int n = -128;
        final UnivariateFunction f = new Sin();
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                final FastFourierTransformer fft;
                fft = new FastFourierTransformer(norm[i]);
                try {
                    fft.transform(f, 0.0, Math.PI, n, type[j]);
                    fft.transform(f, 0.0, Math.PI, n, type[j]);
                    fail(norm[i] + ", " + type[j] +
                        ": MathIllegalArgumentException was expected");
                } catch (MathIllegalArgumentException e) {
                    // Expected behaviour
                }
            }
        }
    }

    @Test
    void testTransformFunctionInvalidBounds() {
        final int n = 128;
        final UnivariateFunction f = new Sin();
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                final FastFourierTransformer fft;
                fft = new FastFourierTransformer(norm[i]);
                try {
                    fft.transform(f, Math.PI, 0.0, n, type[j]);
                    fail(norm[i] + ", " + type[j] +
                        ": MathIllegalArgumentException was expected");
                } catch (MathIllegalArgumentException e) {
                    // Expected behaviour
                }
            }
        }
    }

    @Test
    void testTransformOnePoint() {
        double[][] data = new double[][] { { 1.0 }, { 2.0} };
        FastFourierTransformer.transformInPlace(data,
                                                DftNormalization.STANDARD,
                                                TransformType.FORWARD);
        assertEquals(1.0, data[0][0], 1.0e-15);
        assertEquals(2.0, data[1][0], 1.0e-15);
    }

    /*
     * Utility methods for checking (successful) transforms.
     */

    private static Complex[] createComplexData(final int n) {
        final Random random = new Random(SEED);
        final Complex[] data = new Complex[n];
        for (int i = 0; i < n; i++) {
            final double re = 2.0 * random.nextDouble() - 1.0;
            final double im = 2.0 * random.nextDouble() - 1.0;
            data[i] = new Complex(re, im);
        }
        return data;
    }

    private static double[] createRealData(final int n) {
        final Random random = new Random(SEED);
        final double[] data = new double[n];
        for (int i = 0; i < n; i++) {
            data[i] = 2.0 * random.nextDouble() - 1.0;
        }
        return data;
    }

    /** Naive implementation of DFT, for reference. */
    private static Complex[] dft(final Complex[] x, final int sgn) {
        final int n = x.length;
        final double[] cos = new double[n];
        final double[] sin = new double[n];
        final Complex[] y = new Complex[n];
        for (int i = 0; i < n; i++) {
            final double arg = 2.0 * FastMath.PI * i / n;
            cos[i] = FastMath.cos(arg);
            sin[i] = FastMath.sin(arg);
        }
        for (int i = 0; i < n; i++) {
            double yr = 0.0;
            double yi = 0.0;
            for (int j = 0; j < n; j++) {
                final int index = (i * j) % n;
                final double c = cos[index];
                final double s = sin[index];
                final double xr = x[j].getReal();
                final double xi = x[j].getImaginary();
                yr += c * xr - sgn * s * xi;
                yi += sgn * s * xr + c * xi;
            }
            y[i] = new Complex(yr, yi);
        }
        return y;
    }

    private static void doTestTransformComplex(final int n, final double tol,
        final DftNormalization normalization,
        final TransformType type) {
        final FastFourierTransformer fft;
        fft = new FastFourierTransformer(normalization);
        final Complex[] x = createComplexData(n);
        final Complex[] expected;
        final double s;
        if (type==TransformType.FORWARD) {
            expected = dft(x, -1);
            if (normalization == DftNormalization.STANDARD){
                s = 1.0;
            } else {
                s = 1.0 / FastMath.sqrt(n);
            }
        } else {
            expected = dft(x, 1);
            if (normalization == DftNormalization.STANDARD) {
                s = 1.0 / n;
            } else {
                s = 1.0 / FastMath.sqrt(n);
            }
        }
        final Complex[] actual = fft.transform(x, type);
        for (int i = 0; i < n; i++) {
            final String msg;
            msg = String.format("%s, %s, %d, %d", normalization, type, n, i);
            final double re = s * expected[i].getReal();
            assertEquals(re, actual[i].getReal(),
                tol * FastMath.abs(re),
                msg);
            final double im = s * expected[i].getImaginary();
            assertEquals(im, actual[i].getImaginary(), tol *
                FastMath.abs(re), msg);
        }
    }

    private static void doTestTransformReal(final int n, final double tol,
        final DftNormalization normalization,
        final TransformType type) {
        final FastFourierTransformer fft;
        fft = new FastFourierTransformer(normalization);
        final double[] x = createRealData(n);
        final Complex[] xc = new Complex[n];
        for (int i = 0; i < n; i++) {
            xc[i] = new Complex(x[i], 0.0);
        }
        final Complex[] expected;
        final double s;
        if (type == TransformType.FORWARD) {
            expected = dft(xc, -1);
            if (normalization == DftNormalization.STANDARD) {
                s = 1.0;
            } else {
                s = 1.0 / FastMath.sqrt(n);
            }
        } else {
            expected = dft(xc, 1);
            if (normalization == DftNormalization.STANDARD) {
                s = 1.0 / n;
            } else {
                s = 1.0 / FastMath.sqrt(n);
            }
        }
        final Complex[] actual = fft.transform(x, type);
        for (int i = 0; i < n; i++) {
            final String msg;
            msg = String.format("%s, %s, %d, %d", normalization, type, n, i);
            final double re = s * expected[i].getReal();
            assertEquals(re, actual[i].getReal(),
                tol * FastMath.abs(re),
                msg);
            final double im = s * expected[i].getImaginary();
            assertEquals(im, actual[i].getImaginary(), tol *
                FastMath.abs(re), msg);
        }
    }

    private static void doTestTransformFunction(final UnivariateFunction f,
        final double min, final double max, int n, final double tol,
        final DftNormalization normalization,
        final TransformType type) {
        final FastFourierTransformer fft;
        fft = new FastFourierTransformer(normalization);
        final Complex[] x = new Complex[n];
        for (int i = 0; i < n; i++) {
            final double t = min + i * (max - min) / n;
            x[i] = new Complex(f.value(t));
        }
        final Complex[] expected;
        final double s;
        if (type == TransformType.FORWARD) {
            expected = dft(x, -1);
            if (normalization == DftNormalization.STANDARD) {
                s = 1.0;
            } else {
                s = 1.0 / FastMath.sqrt(n);
            }
        } else {
            expected = dft(x, 1);
            if (normalization == DftNormalization.STANDARD) {
                s = 1.0 / n;
            } else {
                s = 1.0 / FastMath.sqrt(n);
            }
        }
        final Complex[] actual = fft.transform(f, min, max, n, type);
        for (int i = 0; i < n; i++) {
            final String msg = String.format("%d, %d", n, i);
            final double re = s * expected[i].getReal();
            assertEquals(re, actual[i].getReal(),
                tol * FastMath.abs(re),
                msg);
            final double im = s * expected[i].getImaginary();
            assertEquals(im, actual[i].getImaginary(), tol *
                FastMath.abs(re), msg);
        }
    }

    /*
     * Tests of standard transform (when data is valid).
     */

    @Test
    void testTransformComplex() {
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                doTestTransformComplex(2, 1.0E-15, norm[i], type[j]);
                doTestTransformComplex(4, 1.0E-14, norm[i], type[j]);
                doTestTransformComplex(8, 1.0E-14, norm[i], type[j]);
                doTestTransformComplex(16, 1.0E-13, norm[i], type[j]);
                doTestTransformComplex(32, 1.0E-13, norm[i], type[j]);
                doTestTransformComplex(64, 1.0E-12, norm[i], type[j]);
                doTestTransformComplex(128, 1.0E-12, norm[i], type[j]);
            }
        }
    }

    @Test
    void testStandardTransformReal() {
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                doTestTransformReal(2, 1.0E-15, norm[i], type[j]);
                doTestTransformReal(4, 1.0E-14, norm[i], type[j]);
                doTestTransformReal(8, 1.0E-14, norm[i], type[j]);
                doTestTransformReal(16, 1.0E-13, norm[i], type[j]);
                doTestTransformReal(32, 1.0E-13, norm[i], type[j]);
                doTestTransformReal(64, 1.0E-13, norm[i], type[j]);
                doTestTransformReal(128, 1.0E-11, norm[i], type[j]);
            }
        }
    }

    @Test
    void testStandardTransformFunction() {
        final UnivariateFunction f = new Sinc();
        final double min = -FastMath.PI;
        final double max = FastMath.PI;
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                doTestTransformFunction(f, min, max, 2, 1.0E-15, norm[i], type[j]);
                doTestTransformFunction(f, min, max, 4, 1.0E-14, norm[i], type[j]);
                doTestTransformFunction(f, min, max, 8, 1.0E-14, norm[i], type[j]);
                doTestTransformFunction(f, min, max, 16, 1.0E-13, norm[i], type[j]);
                doTestTransformFunction(f, min, max, 32, 1.0E-13, norm[i], type[j]);
                doTestTransformFunction(f, min, max, 64, 1.0E-12, norm[i], type[j]);
                doTestTransformFunction(f, min, max, 128, 1.0E-11, norm[i], type[j]);
            }
        }
    }

    /*
     * Additional tests for 1D data.
     */

    /**
     * Test of transformer for the ad hoc data taken from Mathematica.
     */
    @Test
    void testAdHocData() {
        FastFourierTransformer transformer;
        transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] result; double tolerance = 1E-12;

        double[] x = {1.3, 2.4, 1.7, 4.1, 2.9, 1.7, 5.1, 2.7};
        Complex[] y = {
            new Complex(21.9, 0.0),
            new Complex(-2.09497474683058, 1.91507575950825),
            new Complex(-2.6, 2.7),
            new Complex(-1.10502525316942, -4.88492424049175),
            new Complex(0.1, 0.0),
            new Complex(-1.10502525316942, 4.88492424049175),
            new Complex(-2.6, -2.7),
            new Complex(-2.09497474683058, -1.91507575950825)};

        result = transformer.transform(x, TransformType.FORWARD);
        for (int i = 0; i < result.length; i++) {
            assertEquals(y[i].getReal(), result[i].getReal(), tolerance);
            assertEquals(y[i].getImaginary(), result[i].getImaginary(), tolerance);
        }

        result = transformer.transform(y, TransformType.INVERSE);
        for (int i = 0; i < result.length; i++) {
            assertEquals(x[i], result[i].getReal(), tolerance);
            assertEquals(0.0, result[i].getImaginary(), tolerance);
        }

        double[] x2 = {10.4, 21.6, 40.8, 13.6, 23.2, 32.8, 13.6, 19.2};
        TransformUtils.scaleArray(x2, 1.0 / FastMath.sqrt(x2.length));
        Complex[] y2 = y;

        transformer = new FastFourierTransformer(DftNormalization.UNITARY);
        result = transformer.transform(y2, TransformType.FORWARD);
        for (int i = 0; i < result.length; i++) {
            assertEquals(x2[i], result[i].getReal(), tolerance);
            assertEquals(0.0, result[i].getImaginary(), tolerance);
        }

        result = transformer.transform(x2, TransformType.INVERSE);
        for (int i = 0; i < result.length; i++) {
            assertEquals(y2[i].getReal(), result[i].getReal(), tolerance);
            assertEquals(y2[i].getImaginary(), result[i].getImaginary(), tolerance);
        }
    }

    /**
     * Test of transformer for the sine function.
     */
    @Test
    void testSinFunction() {
        UnivariateFunction f = new Sin();
        FastFourierTransformer transformer;
        transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] result; int N = 1 << 8;
        double min, max, tolerance = 1E-12;

        min = 0.0; max = 2.0 * FastMath.PI;
        result = transformer.transform(f, min, max, N, TransformType.FORWARD);
        assertEquals(0.0, result[1].getReal(), tolerance);
        assertEquals(-(N >> 1), result[1].getImaginary(), tolerance);
        assertEquals(0.0, result[N-1].getReal(), tolerance);
        assertEquals(N >> 1, result[N-1].getImaginary(), tolerance);
        for (int i = 0; i < N-1; i += (i == 0 ? 2 : 1)) {
            assertEquals(0.0, result[i].getReal(), tolerance);
            assertEquals(0.0, result[i].getImaginary(), tolerance);
        }

        min = -FastMath.PI; max = FastMath.PI;
        result = transformer.transform(f, min, max, N, TransformType.INVERSE);
        assertEquals(0.0, result[1].getReal(), tolerance);
        assertEquals(-0.5, result[1].getImaginary(), tolerance);
        assertEquals(0.0, result[N-1].getReal(), tolerance);
        assertEquals(0.5, result[N-1].getImaginary(), tolerance);
        for (int i = 0; i < N-1; i += (i == 0 ? 2 : 1)) {
            assertEquals(0.0, result[i].getReal(), tolerance);
            assertEquals(0.0, result[i].getImaginary(), tolerance);
        }
    }

}
