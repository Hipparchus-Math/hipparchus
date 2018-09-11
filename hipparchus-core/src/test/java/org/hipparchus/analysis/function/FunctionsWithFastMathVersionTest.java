/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.analysis.function;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class FunctionsWithFastMathVersionTest {

    @Test
    public void testAbs() {
        doTestF0(new Abs(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testAcos() {
        doTestF0(new Acos(), -1.25, -0.25, 0.0, 0.25, 1.25);
        doTestFn(new Acos(), -1.25, -0.25, 0.0, 0.25, 1.25);
    }

    @Test
    public void testAcosH() {
        doTestF0(new Acosh(), -10.0, -5.0, -0.5, 0.5, 5.0, 10.0);
        doTestFn(new Acosh(), -10.0, -5.0, -0.5, 0.5, 5.0, 10.0);
    }

    @Test
    public void testAsin() {
        doTestF0(new Asin(), -1.25, -0.25, 0.0, 0.25, 1.25);
        doTestFn(new Asin(), -1.25, -0.25, 0.0, 0.25, 1.25);
    }

    @Test
    public void testAsinh() {
        doTestF0(new Asinh(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Asinh(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testAtan() {
        doTestF0(new Atan(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Atan(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testAtanh() {
        doTestF0(new Atanh(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Atanh(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testCbrt() {
        doTestF0(new Cbrt(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Cbrt(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testCeil() {
        doTestF0(new Ceil(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testCos() {
        doTestF0(new Cos(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Cos(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testCosh() {
        doTestF0(new Cosh(), -10.0, -1.25, -0.5, 0.0, 0.5, 1.25, 10.0);
        doTestFn(new Cosh(), -10.0, -1.25, -0.5, 0.0, 0.5, 1.25, 10.0);
    }

    @Test
    public void testExp() {
        doTestF0(new Exp(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Exp(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testExpm1() {
        doTestF0(new Expm1(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Expm1(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testFloor() {
        doTestF0(new Floor(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testLog() {
        doTestF0(new Log(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Log(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testLog10() {
        doTestF0(new Log10(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Log10(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testLog1p() {
        doTestF0(new Log1p(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Log1p(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testRint() {
        doTestF0(new Rint(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testSignum() {
        doTestF0(new Signum(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testSin() {
        doTestF0(new Sin(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Sin(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testSinh() {
        doTestF0(new Sinh(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Sinh(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testSqrt() {
        doTestF0(new Sqrt(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Sqrt(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testTan() {
        doTestF0(new Tan(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Tan(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testTanh() {
        doTestF0(new Tanh(), -10.0, -1.25, 0.0, 1.25, 10.0);
        doTestFn(new Tanh(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testUlp() {
        doTestF0(new Ulp(), -10.0, -1.25, 0.0, 1.25, 10.0);
    }

    @Test
    public void testAtan2() {
        final double[] d = { Double.NEGATIVE_INFINITY, -13.4, -0.4, 0.0, 0.4, 13.4, Double.POSITIVE_INFINITY, Double.NaN};
        for (double x : d) {
            for (double y : d) {
                double aRef = FastMath.atan2(y, x);
                double a    = new Atan2().value(y, x);
                if (Double.isNaN(aRef)) {
                    Assert.assertTrue(Double.isNaN(a));
                } else if (Double.isInfinite(aRef)) {
                    Assert.assertTrue(Double.isInfinite(a));
                    Assert.assertTrue(a * aRef > 0.0);
                } else {
                    Assert.assertEquals(aRef, a, FastMath.ulp(aRef));
                }
            }
        }
    }

    private void doTestF0(final UnivariateFunction f, double... x) {
        try {
            final Method fastMathVersion = FastMath.class.getMethod(f.getClass().getSimpleName().toLowerCase(), Double.TYPE);
            for (double xi : x) {
                checkF0Equality(f, fastMathVersion, xi);
            }
            checkF0Equality(f, fastMathVersion, Double.NEGATIVE_INFINITY);
            checkF0Equality(f, fastMathVersion, Double.POSITIVE_INFINITY);
            checkF0Equality(f, fastMathVersion, Double.NaN);
        } catch (NoSuchMethodException e) {
            Assert.fail(e.getLocalizedMessage());
        }
    }

    private void checkF0Equality(final UnivariateFunction f, final Method ref, final double x) {
        try {
            double yRef = ((Double) ref.invoke(null, x)).doubleValue();
            double y    = f.value(x);
            if (Double.isNaN(yRef)) {
                Assert.assertTrue(Double.isNaN(y));
            } else if (Double.isInfinite(yRef)) {
                Assert.assertTrue(Double.isInfinite(y));
                Assert.assertTrue(y * yRef > 0.0);
            } else {
                Assert.assertEquals(yRef, y, FastMath.ulp(yRef));
            }
        } catch (InvocationTargetException | IllegalAccessException | IllegalArgumentException e) {
            Assert.fail(e.getLocalizedMessage());
        }
    }

    private void doTestFn(final UnivariateDifferentiableFunction f, double... x) {
        try {
            final Method fastMathVersion = FastMath.class.getMethod(f.getClass().getSimpleName().toLowerCase(),
                                                                    RealFieldElement.class);
            for (double xi : x) {
                checkFnEqualities(f, fastMathVersion, xi);
            }
            checkFnEqualities(f, fastMathVersion, Double.NEGATIVE_INFINITY);
            checkFnEqualities(f, fastMathVersion, Double.POSITIVE_INFINITY);
            checkFnEqualities(f, fastMathVersion, Double.NaN);
        } catch (NoSuchMethodException e) {
            Assert.fail(e.getLocalizedMessage());
        }
    }

    private void checkFnEqualities(final UnivariateDifferentiableFunction f, final Method ref, final double x) {
        try {
            DSFactory           factory = new DSFactory(1, 5);
            DerivativeStructure xDS     = factory.variable(0, x);
            DerivativeStructure yRef    = (DerivativeStructure) ref.invoke(null, xDS);
            DerivativeStructure y       = f.value(xDS);
            for (int order = 0; order < factory.getCompiler().getOrder(); ++order) {
                if (Double.isNaN(yRef.getPartialDerivative(order))) {
                    Assert.assertTrue(Double.isNaN(y.getPartialDerivative(order)));
                } else if (Double.isInfinite(yRef.getPartialDerivative(order))) {
                    Assert.assertTrue(Double.isInfinite(y.getPartialDerivative(order)));
                    Assert.assertTrue(y.getPartialDerivative(order) * yRef.getPartialDerivative(order) > 0.0);
                } else {
                    Assert.assertEquals(yRef.getPartialDerivative(order), y.getPartialDerivative(order), FastMath.ulp(yRef.getPartialDerivative(order)));
                }
            }
        } catch (InvocationTargetException | IllegalAccessException | IllegalArgumentException e) {
            Assert.fail(e.getLocalizedMessage());
        }
    }

}
