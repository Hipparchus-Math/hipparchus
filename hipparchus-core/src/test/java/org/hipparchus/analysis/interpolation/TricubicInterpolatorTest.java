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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.CalculusFieldTrivariateFunction;
import org.hipparchus.analysis.FieldTrivariateFunction;
import org.hipparchus.analysis.TrivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test case for the {@link TricubicInterpolator tricubic interpolator}.
 */
public final class TricubicInterpolatorTest {
    /**
     * Test preconditions.
     */
    @Test
    void testPreconditions() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2.5};
        double[] zval = new double[] {-12, -8, -5.5, -3, 0, 2.5};
        double[][][] fval = new double[xval.length][yval.length][zval.length];

        @SuppressWarnings("unused")
        TrivariateFunction tcf = new TricubicInterpolator().interpolate(xval, yval, zval, fval);

        double[] wxval = new double[] {3, 2, 5, 6.5};
        try {
            tcf = new TricubicInterpolator().interpolate(wxval, yval, zval, fval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[] wyval = new double[] {-4, -1, -1, 2.5};
        try {
            tcf = new TricubicInterpolator().interpolate(xval, wyval, zval, fval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[] wzval = new double[] {-12, -8, -9, -3, 0, 2.5};
        try {
            tcf = new TricubicInterpolator().interpolate(xval, yval, wzval, fval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[][][] wfval = new double[xval.length - 1][yval.length][zval.length];
        try {
            tcf = new TricubicInterpolator().interpolate(xval, yval, zval, wfval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        wfval = new double[xval.length][yval.length - 1][zval.length];
        try {
            tcf = new TricubicInterpolator().interpolate(xval, yval, zval, wfval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        wfval = new double[xval.length][yval.length][zval.length - 1];
        try {
            tcf = new TricubicInterpolator().interpolate(xval, yval, zval, wfval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
    }

    public void testIsValid() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2.5};
        double[] zval = new double[] {-12, -8, -5.5, -3, 0, 2.5};
        double[][][] fval = new double[xval.length][yval.length][zval.length];

        TricubicInterpolatingFunction tcf = new TricubicInterpolator().interpolate(xval, yval, zval, fval);

        // Valid.
        assertTrue(tcf.isValidPoint(4, -3, -8));
        assertTrue(tcf.isValidPoint(5, -3, -8));
        assertTrue(tcf.isValidPoint(4, -1, -8));
        assertTrue(tcf.isValidPoint(5, -1, -8));
        assertTrue(tcf.isValidPoint(4, -3, 0));
        assertTrue(tcf.isValidPoint(5, -3, 0));
        assertTrue(tcf.isValidPoint(4, -1, 0));
        assertTrue(tcf.isValidPoint(5, -1, 0));

        // Invalid.
        assertFalse(tcf.isValidPoint(3.5, -3, -8));
        assertFalse(tcf.isValidPoint(4.5, -3.1, -8));
        assertFalse(tcf.isValidPoint(4.5, -2, 0));
        assertFalse(tcf.isValidPoint(4.5, 0, -3.5));
        assertFalse(tcf.isValidPoint(-10, 4.1, -1));
    }

    /**
     * Test for a plane.
     * <p>
     *  f(x, y, z) = 2 x - 3 y - 4 z + 5
     * </p>
     */
    @Test
    void testPlane() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2, 2.5};
        double[] zval = new double[] {-12, -8, -5.5, -3, 0, 2.5};

        // Function values
        TrivariateFunction f = (x, y, z) -> 2 * x - 3 * y - 4 * z + 5;
        CalculusFieldTrivariateFunction<Binary64> fT = new FieldTrivariateFunction() {
            @Override
            public <T extends CalculusFieldElement<T>> T value(T x, T y, T z) {
                return x.twice().add(y.multiply(-3)).add(z.multiply(-4)).add(5);
            }
        }.toCalculusFieldTrivariateFunction(Binary64Field.getInstance());

        double[][][] fval = new double[xval.length][yval.length][zval.length];

        for (int i = 0; i < xval.length; i++) {
            for (int j = 0; j < yval.length; j++) {
                for (int k = 0; k < zval.length; k++) {
                    fval[i][j][k] = f.value(xval[i], yval[j], zval[k]);
                }
            }
        }

        TricubicInterpolatingFunction tcf = new TricubicInterpolator().interpolate(xval,
                yval,
                zval,
                fval);
        double x, y, z;
        double expected, result;
        Binary64 x64, y64, z64;
        Binary64 expectedT, resultT;

        x = 4;
        y = -3;
        z = 0;
        expected = f.value(x, y, z);
        result = tcf.value(x, y, z);
        assertEquals(expected, result, 1e-15, "On sample point");

        x64 = new Binary64(x);
        y64 = new Binary64(y);
        z64 = new Binary64(z);
        expectedT = fT.value(x64, y64, z64);
        resultT = tcf.value(x64, y64, z64);
        assertEquals(expectedT.getReal(), resultT.getReal(), 1e-15, "On sample point with CalculusFieldElement");

        x = 4.5;
        y = -1.5;
        z = -4.25;
        expected = f.value(x, y, z);
        result = tcf.value(x, y, z);
        assertEquals(expected, result, 1e-14, "Half-way between sample points (middle of the patch)");

        x64 = new Binary64(x);
        y64 = new Binary64(y);
        z64 = new Binary64(z);
        expectedT = fT.value(x64, y64, z64);
        resultT = tcf.value(x64, y64, z64);
        assertEquals(expectedT.getReal(), resultT.getReal(), 1e-14, "Half-way between sample points (middle of the patch) with CalculusFieldElement");
    }

    /**
     * Sine wave.
     * <p>
     *  f(x, y, z) = a cos [&omega; z - k<sub>y</sub> x - k<sub>y</sub> y]
     * </p>
     * with A = 0.2, &omega; = 0.5, k<sub>x</sub> = 2, k<sub>y</sub> = 1.
     */
    @Test
    void testWave() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2, 2.5};
        double[] zval = new double[] {-12, -8, -5.5, -3, 0, 4};

        final double a = 0.2;
        final double omega = 0.5;
        final double kx = 2;
        final double ky = 1;

        // Function values
        TrivariateFunction f = (x, y, z) -> a * FastMath.cos(omega * z - kx * x - ky * y);
        CalculusFieldTrivariateFunction<Binary64> fT = (x, y, z) ->
            z.multiply(omega).subtract(x.multiply(kx)).subtract(y.multiply(ky)).cos().multiply(a);

        double[][][] fval = new double[xval.length][yval.length][zval.length];
        for (int i = 0; i < xval.length; i++) {
            for (int j = 0; j < yval.length; j++) {
                for (int k = 0; k < zval.length; k++) {
                    fval[i][j][k] = f.value(xval[i], yval[j], zval[k]);
                }
            }
        }

        TricubicInterpolatingFunction tcf = new TricubicInterpolator().interpolate(xval,
                yval,
                zval,
                fval);

        double x, y, z;
        double expected, result;
        Binary64 x64, y64, z64;
        Binary64 expectedT, resultT;

        x = 4;
        y = -3;
        z = 0;
        expected = f.value(x, y, z);
        result = tcf.value(x, y, z);
        assertEquals(expected, result, 1e-14, "On sample point");

        x64 = new Binary64(x);
        y64 = new Binary64(y);
        z64 = new Binary64(z);
        expectedT = fT.value(x64, y64, z64);
        resultT = tcf.value(x64, y64, z64);
        assertEquals(expectedT.getReal(), resultT.getReal(), 1e-14, "On sample point with CalculusFieldElement");

        x = 4.5;
        y = -1.5;
        z = -4.25;
        expected = f.value(x, y, z);
        result = tcf.value(x, y, z);
        assertEquals(expected, result, 1e-1, "Half-way between sample points (middle of the patch)"); // XXX Too high tolerance!

        x64 = new Binary64(x);
        y64 = new Binary64(y);
        z64 = new Binary64(z);
        expectedT = fT.value(x64, y64, z64);
        resultT = tcf.value(x64, y64, z64);
        assertEquals(expectedT.getReal(), resultT.getReal(), 1e-1, "Half-way between sample points (middle of the patch) with CalculusFieldElement"); // XXX Too high tolerance!
    }
}
