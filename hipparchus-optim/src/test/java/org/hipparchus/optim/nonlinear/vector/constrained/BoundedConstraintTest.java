package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.junit.Assert;
import org.junit.Test;

public class BoundedConstraintTest {

    @Test
    public void testConstructor() {
        // GIVEN
        final RealVector lowerBound = MatrixUtils.createRealVector(new double[] { -1. });
        final RealVector upperBound = MatrixUtils.createRealVector(new double[] { 1. });
        // WHEN
        final BoundedConstraint boundedConstraint = new TestBoundedConstraint(lowerBound, upperBound);
        // THEN
        Assert.assertEquals(lowerBound.getEntry(0), boundedConstraint.getLowerBound().getEntry(0), 0.);
        Assert.assertEquals(upperBound.getEntry(0), boundedConstraint.getUpperBound().getEntry(0), 0.);
    }

    @Test
    public void testConstructorNullLowerBound() {
        // GIVEN
        final RealVector upperBound = MatrixUtils.createRealVector(new double[1]);
        // WHEN
        final BoundedConstraint boundedConstraint = new TestBoundedConstraint(null, upperBound);
        // THEN
        Assert.assertEquals(boundedConstraint.getLowerBound().getDimension(),
                boundedConstraint.getUpperBound().getDimension());
    }

    @Test
    public void testConstructorNullUpperBound() {
        // GIVEN
        final RealVector lowerBound = MatrixUtils.createRealVector(new double[1]);
        // WHEN
        final BoundedConstraint boundedConstraint = new TestBoundedConstraint(lowerBound, null);
        // THEN
        Assert.assertEquals(boundedConstraint.getLowerBound().getDimension(),
                boundedConstraint.getUpperBound().getDimension());
    }

    private static class TestBoundedConstraint extends BoundedConstraint {

        TestBoundedConstraint(RealVector lower, RealVector upper) {
            super(lower, upper);
        }

        @Override
        public int dim() {
            return 0;
        }

        @Override
        public RealVector value(RealVector x) {
            return null;
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            return null;
        }
    }

}
