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

package org.hipparchus.analysis;

import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.Derivative;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.MultivariateDifferentiableFunction;
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.analysis.function.Add;
import org.hipparchus.analysis.function.Constant;
import org.hipparchus.analysis.function.Cos;
import org.hipparchus.analysis.function.Cosh;
import org.hipparchus.analysis.function.Divide;
import org.hipparchus.analysis.function.Identity;
import org.hipparchus.analysis.function.Inverse;
import org.hipparchus.analysis.function.Log;
import org.hipparchus.analysis.function.Max;
import org.hipparchus.analysis.function.Min;
import org.hipparchus.analysis.function.Minus;
import org.hipparchus.analysis.function.Multiply;
import org.hipparchus.analysis.function.Pow;
import org.hipparchus.analysis.function.Power;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.analysis.function.Sinc;
import org.hipparchus.analysis.function.Subtract;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test for {@link FunctionUtils}.
 */
public class FunctionUtilsTest {
    private final double EPS = FastMath.ulp(1d);

    @Test
    public void testCompose() {
        UnivariateFunction id = new Identity();
        Assertions.assertEquals(3, FunctionUtils.compose(id, id, id).value(3), EPS);

        UnivariateFunction c = new Constant(4);
        Assertions.assertEquals(4, FunctionUtils.compose(id, c).value(3), EPS);
        Assertions.assertEquals(4, FunctionUtils.compose(c, id).value(3), EPS);

        UnivariateFunction m = new Minus();
        Assertions.assertEquals(-3, FunctionUtils.compose(m).value(3), EPS);
        Assertions.assertEquals(3, FunctionUtils.compose(m, m).value(3), EPS);

        UnivariateFunction inv = new Inverse();
        Assertions.assertEquals(-0.25, FunctionUtils.compose(inv, m, c, id).value(3), EPS);

        UnivariateFunction pow = new Power(2);
        Assertions.assertEquals(81, FunctionUtils.compose(pow, pow).value(3), EPS);
    }

    @Test
    public void testComposeDifferentiable() {
        DSFactory factory = new DSFactory(1, 1);
        UnivariateDifferentiableFunction id = new Identity();
        Assertions.assertEquals(1, FunctionUtils.compose(id, id, id).value(factory.variable(0, 3)).getPartialDerivative(1), EPS);
        Assertions.assertEquals(1.5, FunctionUtils.compose(id, id, id).value(1.5), EPS);

        UnivariateDifferentiableFunction c = new Constant(4);
        Assertions.assertEquals(0, FunctionUtils.compose(id, c).value(factory.variable(0, 3)).getPartialDerivative(1), EPS);
        Assertions.assertEquals(0, FunctionUtils.compose(c, id).value(factory.variable(0, 3)).getPartialDerivative(1), EPS);

        UnivariateDifferentiableFunction m = new Minus();
        Assertions.assertEquals(-1, FunctionUtils.compose(m).value(factory.variable(0, 3)).getPartialDerivative(1), EPS);
        Assertions.assertEquals(1, FunctionUtils.compose(m, m).value(factory.variable(0, 3)).getPartialDerivative(1), EPS);

        UnivariateDifferentiableFunction inv = new Inverse();
        Assertions.assertEquals(0.25, FunctionUtils.compose(inv, m, id).value(factory.variable(0, 2)).getPartialDerivative(1), EPS);

        UnivariateDifferentiableFunction pow = new Power(2);
        Assertions.assertEquals(108, FunctionUtils.compose(pow, pow).value(factory.variable(0, 3)).getPartialDerivative(1), EPS);

        UnivariateDifferentiableFunction log = new Log();
        double a = 9876.54321;
        Assertions.assertEquals(pow.value(factory.variable(0, a)).getPartialDerivative(1) / pow.value(a),
                            FunctionUtils.compose(log, pow).value(factory.variable(0, a)).getPartialDerivative(1), EPS);
    }

    @Test
    public void testAdd() {
        UnivariateFunction id = new Identity();
        UnivariateFunction c = new Constant(4);
        UnivariateFunction m = new Minus();
        UnivariateFunction inv = new Inverse();

        Assertions.assertEquals(4.5, FunctionUtils.add(inv, m, c, id).value(2), EPS);
        Assertions.assertEquals(4 + 2, FunctionUtils.add(c, id).value(2), EPS);
        Assertions.assertEquals(4 - 2, FunctionUtils.add(c, FunctionUtils.compose(m, id)).value(2), EPS);
    }

    @Test
    public void testAddDifferentiable() {
        UnivariateDifferentiableFunction sin = new Sin();
        UnivariateDifferentiableFunction c = new Constant(4);
        UnivariateDifferentiableFunction m = new Minus();
        UnivariateDifferentiableFunction inv = new Inverse();

        final double a = 123.456;
        DSFactory factory = new DSFactory(1, 1);
        Assertions.assertEquals(- 1 / (a * a) -1 + FastMath.cos(a),
                            FunctionUtils.add(inv, m, c, sin).value(factory.variable(0, a)).getPartialDerivative(1),
                            EPS);
        Assertions.assertEquals(4 + FastMath.sin(1.2), FunctionUtils.add(sin, c).value(1.2), EPS);
    }

    @Test
    public void testMultiply() {
        UnivariateFunction c = new Constant(4);
        Assertions.assertEquals(16, FunctionUtils.multiply(c, c).value(12345), EPS);

        UnivariateFunction inv = new Inverse();
        UnivariateFunction pow = new Power(2);
        Assertions.assertEquals(1, FunctionUtils.multiply(FunctionUtils.compose(inv, pow), pow).value(3.5), EPS);
    }

    @Test
    public void testMultiplyDifferentiable() {
        UnivariateDifferentiableFunction c = new Constant(4);
        UnivariateDifferentiableFunction id = new Identity();
        DSFactory factory = new DSFactory(1, 1);
        final double a = 1.2345678;
        Assertions.assertEquals(8 * a, FunctionUtils.multiply(c, id, id).value(factory.variable(0, a)).getPartialDerivative(1), EPS);

        UnivariateDifferentiableFunction inv = new Inverse();
        UnivariateDifferentiableFunction pow = new Power(2.5);
        UnivariateDifferentiableFunction cos = new Cos();
        Assertions.assertEquals(1.5 * FastMath.sqrt(a) * FastMath.cos(a) - FastMath.pow(a, 1.5) * FastMath.sin(a),
                            FunctionUtils.multiply(inv, pow, cos).value(factory.variable(0, a)).getPartialDerivative(1), EPS);

        UnivariateDifferentiableFunction cosh = new Cosh();
        Assertions.assertEquals(1.5 * FastMath.sqrt(a) * FastMath.cosh(a) + FastMath.pow(a, 1.5) * FastMath.sinh(a),
                            FunctionUtils.multiply(inv, pow, cosh).value(factory.variable(0, a)).getPartialDerivative(1), 8 * EPS);
        Assertions.assertEquals(16, FunctionUtils.multiply(c, c).value(FastMath.PI), EPS);
    }

    @Test
    public void testCombine() {
        BivariateFunction bi = new Subtract();
        UnivariateFunction id = new Identity();
        UnivariateFunction m = new Minus();
        UnivariateFunction c = FunctionUtils.combine(bi, id, m);
        Assertions.assertEquals(4.6912, c.value(2.3456), EPS);

        bi = new Multiply();
        UnivariateFunction inv = new Inverse();
        c = FunctionUtils.combine(bi, id, inv);
        Assertions.assertEquals(1, c.value(2.3456), EPS);
    }

    @Test
    public void testCollector() {
        BivariateFunction bi = new Add();
        MultivariateFunction coll = FunctionUtils.collector(bi, 0);
        Assertions.assertEquals(10, coll.value(new double[] {1, 2, 3, 4}), EPS);

        bi = new Multiply();
        coll = FunctionUtils.collector(bi, 1);
        Assertions.assertEquals(24, coll.value(new double[] {1, 2, 3, 4}), EPS);

        bi = new Max();
        coll = FunctionUtils.collector(bi, Double.NEGATIVE_INFINITY);
        Assertions.assertEquals(10, coll.value(new double[] {1, -2, 7.5, 10, -24, 9.99}), 0);

        bi = new Min();
        coll = FunctionUtils.collector(bi, Double.POSITIVE_INFINITY);
        Assertions.assertEquals(-24, coll.value(new double[] {1, -2, 7.5, 10, -24, 9.99}), 0);
    }

    @Test
    public void testSinc() {
        BivariateFunction div = new Divide();
        UnivariateFunction sin = new Sin();
        UnivariateFunction id = new Identity();
        UnivariateFunction sinc1 = FunctionUtils.combine(div, sin, id);
        UnivariateFunction sinc2 = new Sinc();

        for (int i = 0; i < 10; i++) {
            double x = FastMath.random();
            Assertions.assertEquals(sinc1.value(x), sinc2.value(x), EPS);
        }
    }

    @Test
    public void testFixingArguments() {
        UnivariateFunction scaler = FunctionUtils.fix1stArgument(new Multiply(), 10);
        Assertions.assertEquals(1.23456, scaler.value(0.123456), EPS);

        UnivariateFunction pow1 = new Power(2);
        UnivariateFunction pow2 = FunctionUtils.fix2ndArgument(new Pow(), 2);

        for (int i = 0; i < 10; i++) {
            double x = FastMath.random() * 10;
            Assertions.assertEquals(pow1.value(x), pow2.value(x), 0);
        }
    }

    @Test
    public void testSampleWrongBounds(){
        assertThrows(MathIllegalArgumentException.class, () -> {
            FunctionUtils.sample(new Sin(), FastMath.PI, 0.0, 10);
        });
    }

    @Test
    public void testSampleNegativeNumberOfPoints(){
        assertThrows(MathIllegalArgumentException.class, () -> {
            FunctionUtils.sample(new Sin(), 0.0, FastMath.PI, -1);
        });
    }

    @Test
    public void testSampleNullNumberOfPoints(){
        assertThrows(MathIllegalArgumentException.class, () -> {
            FunctionUtils.sample(new Sin(), 0.0, FastMath.PI, 0);
        });
    }

    @Test
    public void testSample() {
        final int n = 11;
        final double min = 0.0;
        final double max = FastMath.PI;
        final double[] actual = FunctionUtils.sample(new Sin(), min, max, n);
        for (int i = 0; i < n; i++) {
            final double x = min + (max - min) / n * i;
            Assertions.assertEquals(FastMath.sin(x), actual[i], 0.0, "x = " + x);
        }
    }

    @Test
    public void testToDifferentiableUnivariate() {

        final UnivariateFunction f0 = new UnivariateFunction() {
            @Override
            public double value(final double x) {
                return x * x;
            }
        };
        final UnivariateFunction f1 = new UnivariateFunction() {
            @Override
            public double value(final double x) {
                return 2 * x;
            }
        };
        final UnivariateFunction f2 = new UnivariateFunction() {
            @Override
            public double value(final double x) {
                return 2;
            }
        };
        final UnivariateDifferentiableFunction f = FunctionUtils.toDifferentiable(f0, f1, f2);

        DSFactory factory = new DSFactory(1, 2);
        for (double t = -1.0; t < 1; t += 0.01) {
            // x = sin(t)
            DerivativeStructure dsT = factory.variable(0, t);
            DerivativeStructure y = f.value(dsT.sin());
            Assertions.assertEquals(FastMath.sin(t) * FastMath.sin(t),               f.value(FastMath.sin(t)),  1.0e-15);
            Assertions.assertEquals(FastMath.sin(t) * FastMath.sin(t),               y.getValue(),              1.0e-15);
            Assertions.assertEquals(2 * FastMath.cos(t) * FastMath.sin(t),           y.getPartialDerivative(1), 1.0e-15);
            Assertions.assertEquals(2 * (1 - 2 * FastMath.sin(t) * FastMath.sin(t)), y.getPartialDerivative(2), 1.0e-15);
        }

        try {
            f.value(new DSFactory(1, 3).constant(0.0));
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.NUMBER_TOO_LARGE, e.getSpecifier());
            Assertions.assertEquals(2, ((Integer) e.getParts()[1]).intValue());
            Assertions.assertEquals(3, ((Integer) e.getParts()[0]).intValue());
        }
    }

    @Test
    public void testToDifferentiableMultivariate() {

        final double a = 1.5;
        final double b = 0.5;
        final MultivariateFunction f = new MultivariateFunction() {
            @Override
            public double value(final double[] point) {
                return a * point[0] + b * point[1];
            }
        };
        final MultivariateVectorFunction gradient = new MultivariateVectorFunction() {
            @Override
            public double[] value(final double[] point) {
                return new double[] { a, b };
            }
        };
        final MultivariateDifferentiableFunction mdf = FunctionUtils.toDifferentiable(f, gradient);

        DSFactory factory11 = new DSFactory(1, 1);
        for (double t = -1.0; t < 1; t += 0.01) {
            // x = sin(t), y = cos(t), hence the method really becomes univariate
            DerivativeStructure dsT = factory11.variable(0, t);
            DerivativeStructure y = mdf.value(new DerivativeStructure[] { dsT.sin(), dsT.cos() });
            Assertions.assertEquals(a * FastMath.sin(t) + b * FastMath.cos(t), y.getValue(),              1.0e-15);
            Assertions.assertEquals(a * FastMath.cos(t) - b * FastMath.sin(t), y.getPartialDerivative(1), 1.0e-15);
        }

        DSFactory factory21 = new DSFactory(2, 1);
        for (double u = -1.0; u < 1; u += 0.01) {
            DerivativeStructure dsU = factory21.variable(0, u);
            for (double v = -1.0; v < 1; v += 0.01) {
                DerivativeStructure dsV = factory21.variable(1, v);
                DerivativeStructure y = mdf.value(new DerivativeStructure[] { dsU, dsV });
                Assertions.assertEquals(a * u + b * v, mdf.value(new double[] { u, v }), 1.0e-15);
                Assertions.assertEquals(a * u + b * v, y.getValue(),                     1.0e-15);
                Assertions.assertEquals(a,             y.getPartialDerivative(1, 0),     1.0e-15);
                Assertions.assertEquals(b,             y.getPartialDerivative(0, 1),     1.0e-15);
            }
        }

        DSFactory factory13 = new DSFactory(1, 3);
        try {
            mdf.value(new DerivativeStructure[] { factory13.constant(0.0), factory13.constant(0.0) });
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.NUMBER_TOO_LARGE, e.getSpecifier());
            Assertions.assertEquals(1, ((Integer) e.getParts()[1]).intValue());
            Assertions.assertEquals(3, ((Integer) e.getParts()[0]).intValue());
        }
    }

    @Test
    public void testToDifferentiableMultivariateInconsistentGradient() {

        final double a = 1.5;
        final double b = 0.5;
        final MultivariateFunction f = new MultivariateFunction() {
            @Override
            public double value(final double[] point) {
                return a * point[0] + b * point[1];
            }
        };
        final MultivariateVectorFunction gradient = new MultivariateVectorFunction() {
            @Override
            public double[] value(final double[] point) {
                return new double[] { a, b, 0.0 };
            }
        };
        final MultivariateDifferentiableFunction mdf = FunctionUtils.toDifferentiable(f, gradient);

        DSFactory factory = new DSFactory(1, 1);
        try {
            DerivativeStructure dsT = factory.variable(0, 0.0);
            mdf.value(new DerivativeStructure[] { dsT.sin(), dsT.cos() });
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            Assertions.assertEquals(3, ((Integer) e.getParts()[0]).intValue());
            Assertions.assertEquals(2, ((Integer) e.getParts()[1]).intValue());
        }
    }

    @Test
    public void testDerivativeUnivariate() {

        final UnivariateDifferentiableFunction f = new UnivariateDifferentiableFunction() {

            @Override
            public double value(double x) {
                return x * x;
            }

            @Override
            public <T extends Derivative<T>> T value(T x) {
                return x.square();
            }

        };

        final UnivariateFunction f0 = FunctionUtils.derivative(f, 0);
        final UnivariateFunction f1 = FunctionUtils.derivative(f, 1);
        final UnivariateFunction f2 = FunctionUtils.derivative(f, 2);

        for (double t = -1.0; t < 1; t += 0.01) {
            Assertions.assertEquals(t * t, f0.value(t), 1.0e-15);
            Assertions.assertEquals(2 * t, f1.value(t), 1.0e-15);
            Assertions.assertEquals(2,     f2.value(t), 1.0e-15);
        }

    }

    @Test
    public void testDerivativeMultivariate() {

        final double a = 1.5;
        final double b = 0.5;
        final double c = 0.25;
        final MultivariateDifferentiableFunction mdf = new MultivariateDifferentiableFunction() {

            @Override
            public double value(double[] point) {
                return a * point[0] * point[0] + b * point[1] * point[1] + c * point[0] * point[1];
            }

            @Override
            public DerivativeStructure value(DerivativeStructure[] point) {
                DerivativeStructure x  = point[0];
                DerivativeStructure y  = point[1];
                DerivativeStructure x2 = x.square();
                DerivativeStructure y2 = y.square();
                DerivativeStructure xy = x.multiply(y);
                return x2.multiply(a).add(y2.multiply(b)).add(xy.multiply(c));
            }

        };

        final MultivariateFunction f       = FunctionUtils.derivative(mdf, new int[] { 0, 0 });
        final MultivariateFunction dfdx    = FunctionUtils.derivative(mdf, new int[] { 1, 0 });
        final MultivariateFunction dfdy    = FunctionUtils.derivative(mdf, new int[] { 0, 1 });
        final MultivariateFunction d2fdx2  = FunctionUtils.derivative(mdf, new int[] { 2, 0 });
        final MultivariateFunction d2fdy2  = FunctionUtils.derivative(mdf, new int[] { 0, 2 });
        final MultivariateFunction d2fdxdy = FunctionUtils.derivative(mdf, new int[] { 1, 1 });

        for (double x = -1.0; x < 1; x += 0.01) {
            for (double y = -1.0; y < 1; y += 0.01) {
                Assertions.assertEquals(a * x * x + b * y * y + c * x * y, f.value(new double[]       { x, y }), 1.0e-15);
                Assertions.assertEquals(2 * a * x + c * y,                 dfdx.value(new double[]    { x, y }), 1.0e-15);
                Assertions.assertEquals(2 * b * y + c * x,                 dfdy.value(new double[]    { x, y }), 1.0e-15);
                Assertions.assertEquals(2 * a,                             d2fdx2.value(new double[]  { x, y }), 1.0e-15);
                Assertions.assertEquals(2 * b,                             d2fdy2.value(new double[]  { x, y }), 1.0e-15);
                Assertions.assertEquals(c,                                 d2fdxdy.value(new double[] { x, y }), 1.0e-15);
            }
        }

    }

}
