/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.analysis.polynomials;

import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.junit.Assert;
import org.junit.Test;

public class SmoothStepFactoryTest {

    final double THRESHOLD = 1e-15;

    @Test(expected = MathIllegalArgumentException.class)
    public void testExceptionBelowBoundary() {
        // Given
        final double                               x          = 2;
        final SmoothStepFactory.SmoothStepFunction smoothstep = SmoothStepFactory.getGeneralOrder(1);

        // When
        smoothstep.value(x);
    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testExceptionOverBoundary() {
        // Given
        final double                               x          = 17;
        final SmoothStepFactory.SmoothStepFunction smoothstep = SmoothStepFactory.getGeneralOrder(1);

        // When
        smoothstep.value(x);
    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testEdgesConsistency() {
        // Given
        final double                               leftEdge   = 5;
        final double                               rightEdge  = 2;
        final double                               x          = 3;
        final SmoothStepFactory.SmoothStepFunction smoothstep = SmoothStepFactory.getGeneralOrder(1);

        // When
        smoothstep.value(leftEdge, rightEdge, x);
    }

    @Test
    public void testBoundaries() {
        // Given
        final double leftEdge  = 5;
        final double rightEdge = 10;
        final double x1        = 2;
        final double x2        = 11;

        final SmoothStepFactory.SmoothStepFunction clamp = SmoothStepFactory.getClamp();

        // When
        final double computedResult1 = clamp.value(leftEdge, rightEdge, x1);
        final double computedResult2 = clamp.value(leftEdge, rightEdge, x2);

        // Then
        Assert.assertEquals(0, computedResult1, THRESHOLD);
        Assert.assertEquals(1, computedResult2, THRESHOLD);
    }

    @Test
    public void testNormalizedInput() {

        // Given
        final double                               x     = 0.4;
        final SmoothStepFactory.SmoothStepFunction cubic = SmoothStepFactory.getCubic();

        // When
        final double computedResult = cubic.value(x);

        // Then
        Assert.assertEquals(0.352, computedResult, THRESHOLD);

    }

    @Test
    public void testClampFunction() {

        // Given
        final double leftEdge  = 5;
        final double rightEdge = 10;
        final double x         = 7;

        final SmoothStepFactory.SmoothStepFunction clamp = SmoothStepFactory.getClamp();

        // When
        final double computedResult = clamp.value(leftEdge, rightEdge, x);

        // Then
        Assert.assertEquals(0.4, computedResult, THRESHOLD);

    }

    @Test
    public void testQuadraticFunction1() {

        // Given
        final double leftEdge  = 5;
        final double rightEdge = 10;
        final double x         = 7;

        final SmoothStepFactory.SmoothStepFunction quadratic = SmoothStepFactory.getQuadratic();

        // When
        final double computedResult = quadratic.value(leftEdge, rightEdge, x);

        // Then
        Assert.assertEquals(0.32, computedResult, THRESHOLD);

    }

    @Test
    public void testQuadraticFunction2() {

        // Given
        final double leftEdge  = 5;
        final double rightEdge = 10;
        final double x         = 8;

        final SmoothStepFactory.SmoothStepFunction quadratic = SmoothStepFactory.getQuadratic();

        // When
        final double computedResult = quadratic.value(leftEdge, rightEdge, x);

        // Then
        Assert.assertEquals(0.68, computedResult, THRESHOLD);

    }

    @Test
    public void testCubicFunction() {

        // Given
        final double leftEdge  = 5;
        final double rightEdge = 10;
        final double x         = 7;

        final SmoothStepFactory.SmoothStepFunction cubic = SmoothStepFactory.getCubic();

        // When
        final double computedResult = cubic.value(leftEdge, rightEdge, x);

        // Then
        Assert.assertEquals(0.352, computedResult, THRESHOLD);

    }

    @Test
    public void testQuinticFunction() {

        // Given
        final double leftEdge  = 5;
        final double rightEdge = 10;
        final double x         = 7;

        final SmoothStepFactory.SmoothStepFunction quintic = SmoothStepFactory.getQuintic();

        // When
        final double computedResult = quintic.value(leftEdge, rightEdge, x);

        // Then
        Assert.assertEquals(0.31744, computedResult, THRESHOLD);

    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testFieldEdgesConsistency() {
        // Given
        final Field<Binary64> field = Binary64Field.getInstance();

        final double   leftEdge  = 5;
        final double   rightEdge = 2;
        final Binary64 x         = new Binary64(3);
        final SmoothStepFactory.FieldSmoothStepFunction<Binary64> smoothstep =
                SmoothStepFactory.getFieldGeneralOrder(field, 1);

        // When
        smoothstep.value(leftEdge, rightEdge, x);
    }

    @Test
    public void testFieldBoundaries() {
        // Given
        final Field<Binary64> field = Binary64Field.getInstance();

        final double   leftEdge  = 5;
        final double   rightEdge = 10;
        final Binary64 x1        = new Binary64(2);
        final Binary64 x2        = new Binary64(11);

        final SmoothStepFactory.FieldSmoothStepFunction<Binary64> clamp = SmoothStepFactory.getClamp(field);

        // When
        final Binary64 computedResult1 = clamp.value(leftEdge, rightEdge, x1);
        final Binary64 computedResult2 = clamp.value(leftEdge, rightEdge, x2);

        // Then
        Assert.assertEquals(0, computedResult1.getReal(), THRESHOLD);
        Assert.assertEquals(1, computedResult2.getReal(), THRESHOLD);
    }

    @Test
    public void testFieldNormalizedInput() {

        // Given
        final Field<Binary64> field = Binary64Field.getInstance();

        final double                                              x     = 0.4;
        final SmoothStepFactory.FieldSmoothStepFunction<Binary64> cubic = SmoothStepFactory.getCubic(field);

        // When
        final Binary64 computedResult = cubic.value(x);

        // Then
        Assert.assertEquals(0.352, computedResult.getReal(), THRESHOLD);

    }

    @Test
    public void testFieldClampFunction() {

        // Given
        final Field<Binary64> field     = Binary64Field.getInstance();
        final double          leftEdge  = 5;
        final double          rightEdge = 10;
        final Binary64        x         = new Binary64(7);

        final SmoothStepFactory.FieldSmoothStepFunction<Binary64> clamp = SmoothStepFactory.getClamp(field);

        // When
        final Binary64 computedResult = clamp.value(leftEdge, rightEdge, x);

        // Then
        Assert.assertEquals(0.4, computedResult.getReal(), THRESHOLD);

    }

    @Test
    public void testFieldQuadraticFunction1() {

        // Given
        final Field<Binary64> field     = Binary64Field.getInstance();
        final double          leftEdge  = 5;
        final double          rightEdge = 10;
        final double          x         = 7;
        final Binary64        xField    = new Binary64(x);

        final SmoothStepFactory.FieldSmoothStepFunction<Binary64> quadratic = SmoothStepFactory.getQuadratic(field);

        // When
        final Binary64 computedResult  = quadratic.value(leftEdge, rightEdge, xField);
        final Binary64 computedResult2 = quadratic.value((x - leftEdge) / (rightEdge - leftEdge));

        // Then
        Assert.assertEquals(0.32, computedResult.getReal(), THRESHOLD);
        Assert.assertEquals(computedResult.getReal(), computedResult2.getReal(), THRESHOLD);

    }

    @Test
    public void testFieldQuadraticFunction2() {

        // Given
        final Field<Binary64> field     = Binary64Field.getInstance();
        final double          leftEdge  = 5;
        final double          rightEdge = 10;
        final double          x         = 8;
        final Binary64        xField    = new Binary64(x);

        final SmoothStepFactory.FieldSmoothStepFunction<Binary64> quadratic = SmoothStepFactory.getQuadratic(field);

        final Binary64 computedResult  = quadratic.value(leftEdge, rightEdge, xField);
        final Binary64 computedResult2 = quadratic.value((x - leftEdge) / (rightEdge - leftEdge));

        // Then
        Assert.assertEquals(0.68, computedResult.getReal(), THRESHOLD);
        Assert.assertEquals(computedResult.getReal(), computedResult2.getReal(), THRESHOLD);

    }

    @Test
    public void testFieldCubicFunction() {

        // Given
        final Field<Binary64> field     = Binary64Field.getInstance();
        final double          leftEdge  = 5;
        final double          rightEdge = 10;
        final Binary64        x         = new Binary64(7);

        final SmoothStepFactory.FieldSmoothStepFunction<Binary64> cubic = SmoothStepFactory.getCubic(field);

        // When
        final Binary64 computedResult = cubic.value(leftEdge, rightEdge, x);

        // Then
        Assert.assertEquals(0.352, computedResult.getReal(), THRESHOLD);

    }

    @Test
    public void testFieldQuinticFunction() {

        // Given
        final Field<Binary64> field       = Binary64Field.getInstance();
        final double          leftEdge    = 5;
        final double          rightEdge   = 10;
        final double          x           = 7;
        final double          xNormalized = (x - leftEdge) / (rightEdge - leftEdge);
        final Binary64        xField      = new Binary64(x);

        final SmoothStepFactory.FieldSmoothStepFunction<Binary64> quintic = SmoothStepFactory.getQuintic(field);

        // When
        final Binary64 computedResult  = quintic.value(leftEdge, rightEdge, xField);
        final Binary64 computedResult2 = quintic.value(xNormalized);
        final Binary64 computedResult3 = quintic.value(new Binary64(xNormalized));

        // Then
        Assert.assertEquals(0.31744, computedResult.getReal(), THRESHOLD);
        Assert.assertEquals(computedResult.getReal(), computedResult2.getReal(), THRESHOLD);
        Assert.assertEquals(computedResult2.getReal(), computedResult3.getReal(), THRESHOLD);
    }

}
