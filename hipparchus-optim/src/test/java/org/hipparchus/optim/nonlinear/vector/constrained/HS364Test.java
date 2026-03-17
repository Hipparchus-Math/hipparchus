/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class HS364Test {

    private static final class HS364Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 6; }

        @Override
        public double value(final RealVector x) {
            return tp364a(x);
        }

        @Override public RealVector gradient(final RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(final RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static final class HS364Ineq extends InequalityConstraint {
        HS364Ineq() { super(new ArrayRealVector(new double[] {0.0, 0.0, 0.0, 0.0})); }

        @Override
        public RealVector value(final RealVector x) {
            final double xmu1 = 0.7853981633;
            final double xmu2 = 2.356194491;
            return new ArrayRealVector(new double[] {
                -x.getEntry(0) + x.getEntry(1) + x.getEntry(2) - x.getEntry(3),
                -x.getEntry(0) - x.getEntry(1) + x.getEntry(2) + x.getEntry(3),
                -x.getEntry(1) * x.getEntry(1) - x.getEntry(2) * x.getEntry(2) +
                (x.getEntry(3) - x.getEntry(0)) * (x.getEntry(3) - x.getEntry(0)) +
                2.0 * x.getEntry(1) * x.getEntry(2) * FastMath.cos(xmu1),
                x.getEntry(1) * x.getEntry(1) + x.getEntry(2) * x.getEntry(2) -
                (x.getEntry(3) + x.getEntry(0)) * (x.getEntry(3) + x.getEntry(0)) -
                2.0 * x.getEntry(1) * x.getEntry(2) * FastMath.cos(xmu2)
            }, false);
        }

        @Override public RealMatrix jacobian(final RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 6; }
    }

//    @Test
//    void testHS364ReferencePoint() {
//        final RealVector x = new ArrayRealVector(new double[] {
//            0.99616882, 4.1960616, 2.9771652, 3.9631949, 1.6536702, 1.2543998
//        }, false);
//
//        final HS364Obj obj = new HS364Obj();
//        final HS364Ineq ineq = new HS364Ineq();
//
//        assertTrue(FastMath.abs(obj.value(x) - 0.0606002) < 5.0e-3);
//        final RealVector g = ineq.value(x);
//        assertTrue(g.getEntry(0) >= -1.0e-6 && g.getEntry(1) >= -1.0e-6 &&
//                   g.getEntry(2) >= -1.0e-5 && g.getEntry(3) >= -1.0e-5);
//    }

    private static double tp364a(final RealVector x) {
        final double[] x1 = new double[31];
        final double[] y1 = new double[31];
        final double[] phi = new double[31];
        final double[] x1a = new double[31];
        final double[] y1a = new double[31];

        final double xInc = 2.0 * FastMath.PI / 30.0;
        for (int i = 0; i < 31; i++) {
            phi[i] = xInc * i;
        }
        tp364b(phi, x1, y1);

        double sum = 0.0;
        for (int i = 0; i < 31; i++) {
            final double coss = tp364c(x, phi[i]);
            final double sins = FastMath.sqrt(FastMath.max(0.0, 1.0 - coss * coss));
            final double cosy = (x.getEntry(3) + x.getEntry(2) * coss - x.getEntry(0) * FastMath.cos(phi[i])) / x.getEntry(1);
            final double siny = (x.getEntry(2) * sins - x.getEntry(0) * FastMath.sin(phi[i])) / x.getEntry(1);
            x1a[i] = x.getEntry(0) * FastMath.cos(phi[i]) + x.getEntry(4) * cosy - x.getEntry(5) * siny;
            y1a[i] = x.getEntry(0) * FastMath.sin(phi[i]) + x.getEntry(4) * siny + x.getEntry(5) * cosy;
            sum += (x1a[i] - x1[i]) * (x1a[i] - x1[i]) + (y1a[i] - y1[i]) * (y1a[i] - y1[i]);
        }

        return FastMath.sqrt(FastMath.max(0.0, sum / 31.0));
    }

    private static void tp364b(final double[] phi, final double[] x1, final double[] y1) {
        for (int i = 0; i < 31; i++) {
            x1[i] = 0.4 + FastMath.sin((2.0 * FastMath.PI) * ((FastMath.PI - phi[i]) / (2.0 * FastMath.PI) - 0.16));
            y1[i] = 2.0 + 0.9 * FastMath.sin(FastMath.PI - phi[i]);
        }
    }

    private static double tp364c(final RealVector x, final double phi) {
        final double m = 2.0 * x.getEntry(0) * x.getEntry(2) * FastMath.sin(phi);
        final double l = 2.0 * x.getEntry(2) * x.getEntry(3) - 2.0 * x.getEntry(0) * x.getEntry(2) * FastMath.cos(phi);
        final double k = x.getEntry(0) * x.getEntry(0) - x.getEntry(1) * x.getEntry(1) +
                         x.getEntry(2) * x.getEntry(2) + x.getEntry(3) * x.getEntry(3) -
                         2.0 * x.getEntry(3) * x.getEntry(0) * FastMath.cos(phi);
        final double a = l * l + m * m;
        final double b = 2.0 * k * l;
        final double c = k * k - m * m;

        double term = FastMath.sqrt(FastMath.abs(b * b - 4.0 * a * c));
        if ((FastMath.PI - phi) < 0.0) {
            term = -term;
        }
        return (-b + term) / (2.0 * a);
    }
}
