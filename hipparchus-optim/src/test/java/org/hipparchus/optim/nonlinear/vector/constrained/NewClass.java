//package org.hipparchus.optim.nonlinear.vector.constrained;
//
//
//
//import org.hipparchus.linear.*;
//import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
//import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
//
//import org.junit.jupiter.api.Test;
//
//import static org.junit.jupiter.api.Assertions.assertEquals;
//
//class HS310Test {
//
//    /** Objective function: (x1*x2)^2 * (10 - x1)^2 * (10 - x1 - x2*(10 - x1)^5)^2 */
//    static final class HS310Objective extends TwiceDifferentiableFunction {
//        @Override public int dim() { return 2; }
//
//        @Override
//        public double value(RealVector x) {
//            double x1 = x.getEntry(0);
//            double x2 = x.getEntry(1);
//
//            double A = x1 * x2;
//            double B = 10.0 - x1;
//            double C = B - x2 * Math.pow(B, 5);
//
//            return (A * A) * (B * B) * (C * C);
//        }
//
//        @Override
//        public RealVector gradient(RealVector x) {
//            double x1 = x.getEntry(0);
//            double x2 = x.getEntry(1);
//
//            double A = x1 * x2;
//            double B = 10.0 - x1;
//            double B4 = Math.pow(B, 4);
//            double B5 = B * B4;
//            double C = B - x2 * B5;
//
//            double common = 2.0 * A * B * C;
//            double g1 = common * (x2 - 1.0 - 5.0 * x2 * B4);
//            double g2 = common * (x1 - B5);
//
//            return new ArrayRealVector(new double[]{ g1, g2 }, false);
//        }
//    }
//
//    @Test
//    void testHS310() {
//        // --- Setup optimizer ---
//        SQPOptimizerS2 opt = new SQPOptimizerS2();
//        opt.setDebugPrinter(System.out::println);
//
//        HS310Objective f = new HS310Objective();
//
//        RealVector x0 = new ArrayRealVector(new double[]{ -12.0, 10.0 });
//        RealVector expected = new ArrayRealVector(new double[]{ 0.0, 0.0 });
//        double fEx = 0.0;
//
//        // --- Optimize ---
//        var result = opt.optimize(f, x0, null, null, null);
//
//        // --- Check results ---
//        double f = f.value(result.getPoint());
//        double tol = 1.0e-6 * (Math.abs(fEx) + 1.0);
//
//        assertEquals(fEx, f, tol, "objective mismatch at optimum");
//        assertEquals(expected.getEntry(0), result.getPoint().getEntry(0), tol, "x1 mismatch");
//        assertEquals(expected.getEntry(1), result.getPoint().getEntry(1), tol, "x2 mismatch");
//    }
//}
