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

package org.hipparchus.analysis.differentiation;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;


/**
 * Test for class {@link GradientFunction}.
 */
class GradientFunctionTest {

    @Test
    void test2DDistance() {
        EuclideanDistance f = new EuclideanDistance();
        GradientFunction g = new GradientFunction(f);
        for (double x = -10; x < 10; x += 0.5) {
            for (double y = -10; y < 10; y += 0.5) {
                double[] point = new double[] { x, y };
                UnitTestUtils.customAssertEquals(f.gradient(point), g.value(point), 1.0e-15);
            }
        }
    }

    @Test
    void test3DDistance() {
        EuclideanDistance f = new EuclideanDistance();
        GradientFunction g = new GradientFunction(f);
        for (double x = -10; x < 10; x += 0.5) {
            for (double y = -10; y < 10; y += 0.5) {
                for (double z = -10; z < 10; z += 0.5) {
                    double[] point = new double[] { x, y, z };
                    UnitTestUtils.customAssertEquals(f.gradient(point), g.value(point), 1.0e-15);
                }
            }
        }
    }

    private static class EuclideanDistance implements MultivariateDifferentiableFunction {

        @Override
        public double value(double[] point) {
            double d2 = 0;
            for (double x : point) {
                d2 += x * x;
            }
            return FastMath.sqrt(d2);
        }

        @Override
        public DerivativeStructure value(DerivativeStructure[] point)
            throws MathIllegalArgumentException {
            DerivativeStructure d2 = point[0].getField().getZero();
            for (DerivativeStructure x : point) {
                d2 = d2.add(x.square());
            }
            return d2.sqrt();
        }

        public double[] gradient(double[] point) {
            double[] gradient = new double[point.length];
            double d = value(point);
            for (int i = 0; i < point.length; ++i) {
                gradient[i] = point[i] / d;
            }
            return gradient;
        }

    }

}
