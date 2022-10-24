package org.hipparchus.analysis.polynomials;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.Assert;
import org.junit.Test;

public class SmoothStepFactoryTest {

    final double THRESHOLD = 1e-15;

    @Test(expected = MathIllegalArgumentException.class)
    public void testExceptionLeftAndRightEdge() {
        // Given
        final double leftEdge  = 5;
        final double rightEdge = leftEdge - 1;

        final SmoothStepFactory.SmoothStepFunction smoothstep = SmoothStepFactory.getGeneralOrder(1);

        // When
        smoothstep.value(leftEdge, rightEdge, 0);
    }

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

    @Test
    public void testBoundaries() {
        // Given
        final double leftEdge     = 5;
        final double rightEdge    = 10;
        final double x1           = 2;
        final double x2           = 11;

        final SmoothStepFactory.SmoothStepFunction clamp = SmoothStepFactory.getClamp();

        // When
        final double computedResult1           = clamp.value(leftEdge, rightEdge, x1);
        final double computedResult2           = clamp.value(leftEdge, rightEdge, x2);

        // Then
        Assert.assertEquals(0, computedResult1, THRESHOLD);
        Assert.assertEquals(1, computedResult2, THRESHOLD);
    }

    @Test
    public void testNormalizedInput() {

        // Given
        final double                               x     = 0.4;
        final SmoothStepFactory.SmoothStepFunction clamp = SmoothStepFactory.getCubic();

        // When
        final double computedResult = clamp.value(x);

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

}
