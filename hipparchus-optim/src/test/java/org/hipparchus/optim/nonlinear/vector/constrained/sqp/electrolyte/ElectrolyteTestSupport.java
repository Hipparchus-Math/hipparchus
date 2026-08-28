/*
 * Licensed to the Hipparchus project under one or more contributor license
 * agreements. See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership. The Hipparchus project licenses
 * this file to You under the Apache License, Version 2.0.
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import java.util.Locale;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;

import static org.junit.jupiter.api.Assertions.assertTrue;

/** Shared formulas and execution helpers for the ripopt electrolyte suite. */
final class ElectrolyteTestSupport {

    static final double A_DH = 0.5091;
    static final double B_DH = 0.3283;
    static final double A_PHI = 0.3915;
    static final double LN_KW = -32.2387;
    static final double LN10 = 2.302585093;
    static final double M_W = 0.018015;

    private ElectrolyteTestSupport() {
    }

    /** Base class supplying a numerical Hessian from the declared gradient. */
    abstract static class Objective extends TwiceDifferentiableFunction {

        private final int dimension;

        Objective(final int dimension) {
            this.dimension = dimension;
        }

        @Override
        public int dim() {
            return dimension;
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            final RealMatrix hessian = MatrixUtils.createRealMatrix(dimension, dimension);
            for (int j = 0; j < dimension; ++j) {
                final double xj = x.getEntry(j);
                final double h = 1.0e-6 * FastMath.max(FastMath.abs(xj), 1.0e-4);
                final RealVector xp = x.copy();
                final RealVector xm = x.copy();
                xp.setEntry(j, xj + h);
                xm.setEntry(j, xj - h);
                final RealVector gp = gradient(xp);
                final RealVector gm = gradient(xm);
                for (int i = 0; i < dimension; ++i) {
                    hessian.setEntry(i, j, (gp.getEntry(i) - gm.getEntry(i)) / (2.0 * h));
                }
            }
            // Numerical differentiation may introduce a small asymmetry.
            return hessian.add(hessian.transpose()).scalarMultiply(0.5);
        }
    }

    /** Gibbs-energy objective used by the log-molality speciation problems. */
    static final class GibbsObjective extends Objective {

        private final double[] mu0;
        private final double[] charges;
        private final double[] dhA;
        private final double[] dhB;

        GibbsObjective(final double[] mu0,
                       final double[] charges,
                       final double[] dhA,
                       final double[] dhB) {
            super(mu0.length);
            this.mu0 = mu0.clone();
            this.charges = charges.clone();
            this.dhA = dhA.clone();
            this.dhB = dhB.clone();
        }

        @Override
        public double value(final RealVector x) {
            final double ionicStrength = ionicStrength(x, charges);
            double value = 0.0;
            for (int i = 0; i < mu0.length; ++i) {
                final double molality = FastMath.exp(x.getEntry(i));
                final double lnGamma = lnGammaDh(charges[i], dhA[i], dhB[i], ionicStrength);
                value += molality * (mu0[i] + lnGamma + x.getEntry(i));
            }
            return value;
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final int n = mu0.length;
            final double ionicStrength = ionicStrength(x, charges);
            double s1 = 0.0;
            for (int i = 0; i < n; ++i) {
                final double molality = FastMath.exp(x.getEntry(i));
                s1 += molality * dLnGammaDhDi(charges[i], dhA[i], dhB[i], ionicStrength);
            }

            final RealVector gradient = new ArrayRealVector(n);
            for (int j = 0; j < n; ++j) {
                final double molality = FastMath.exp(x.getEntry(j));
                final double lnGamma = lnGammaDh(charges[j], dhA[j], dhB[j], ionicStrength);
                final double a = mu0[j] + lnGamma + x.getEntry(j) + 1.0;
                final double dIdx = 0.5 * charges[j] * charges[j] * molality;
                gradient.setEntry(j, molality * a + dIdx * s1);
            }
            return gradient;
        }
    }

    static double lnGammaDh(final double z,
                            final double a,
                            final double bDot,
                            final double ionicStrength) {
        if (ionicStrength <= 0.0) {
            return 0.0;
        }
        final double sqrtI = FastMath.sqrt(ionicStrength);
        final double denominator = 1.0 + a * B_DH * sqrtI;
        return -A_DH * z * z * sqrtI / denominator + bDot * ionicStrength;
    }

    static double dLnGammaDhDi(final double z,
                               final double a,
                               final double bDot,
                               final double ionicStrength) {
        if (ionicStrength <= 1.0e-30) {
            return -A_DH * z * z * 1.0e15 * 0.5;
        }
        final double sqrtI = FastMath.sqrt(ionicStrength);
        final double denominator = 1.0 + a * B_DH * sqrtI;
        return -A_DH * z * z / (2.0 * sqrtI * denominator * denominator) + bDot;
    }

    static double d2LnGammaDhDi2(final double z,
                                 final double a,
                                 final double ionicStrength) {
        if (ionicStrength <= 1.0e-30) {
            return 0.0;
        }
        final double sqrtI = FastMath.sqrt(ionicStrength);
        final double ab = a * B_DH;
        final double denominator = 1.0 + ab * sqrtI;
        final double term1 = 1.0 /
                (4.0 * ionicStrength * sqrtI * denominator * denominator);
        final double term2 = ab /
                (2.0 * ionicStrength * denominator * denominator * denominator);
        return A_DH * z * z * (term1 + term2);
    }

    static double ionicStrength(final RealVector x, final double[] charges) {
        double ionicStrength = 0.0;
        for (int i = 0; i < charges.length; ++i) {
            if (charges[i] != 0.0) {
                ionicStrength += charges[i] * charges[i] * FastMath.exp(x.getEntry(i));
            }
        }
        return 0.5 * ionicStrength;
    }

    static double pitzerLnGammaPm(final double m,
                                  final double beta0,
                                  final double beta1,
                                  final double cPhi) {
        if (m <= 0.0) {
            return 0.0;
        }
        final double sqrtM = FastMath.sqrt(m);
        final double fGamma = -A_PHI *
                (sqrtM / (1.0 + 1.2 * sqrtM) +
                 (2.0 / 1.2) * FastMath.log(1.0 + 1.2 * sqrtM));
        final double x = 2.0 * sqrtM;
        final double exp = FastMath.exp(-x);
        final double bGamma = 2.0 * beta0 +
                2.0 * beta1 / (4.0 * m) * (1.0 - (1.0 + x) * exp);
        return fGamma + m * bGamma + 1.5 * cPhi * m * m;
    }

    static double dPitzerLnGammaPmDm(final double m,
                                     final double beta0,
                                     final double beta1,
                                     final double cPhi) {
        if (m <= 1.0e-30) {
            return 0.0;
        }
        final double sqrtM = FastMath.sqrt(m);
        final double denominator = 1.0 + 1.2 * sqrtM;
        final double dfGamma = -A_PHI *
                (0.5 / (sqrtM * denominator * denominator) +
                 1.0 / (sqrtM * denominator));
        final double dMb = 2.0 * beta0 + beta1 * FastMath.exp(-2.0 * sqrtM);
        return dfGamma + dMb + 3.0 * cPhi * m;
    }

    static double pitzerOsmotic(final double m,
                                final double beta0,
                                final double beta1,
                                final double cPhi) {
        if (m <= 0.0) {
            return 1.0;
        }
        final double sqrtM = FastMath.sqrt(m);
        return 1.0 - A_PHI * sqrtM / (1.0 + 1.2 * sqrtM) +
               m * (beta0 + beta1 * FastMath.exp(-2.0 * sqrtM)) +
               m * m * cPhi;
    }

    static double dPitzerOsmoticDm(final double m,
                                   final double beta0,
                                   final double beta1,
                                   final double cPhi) {
        if (m <= 1.0e-30) {
            return 0.0;
        }
        final double sqrtM = FastMath.sqrt(m);
        final double denominator = 1.0 + 1.2 * sqrtM;
        final double dDh = -A_PHI * 0.5 /
                (sqrtM * denominator * denominator);
        final double exp = FastMath.exp(-2.0 * sqrtM);
        final double dB = beta0 + beta1 * exp * (1.0 - sqrtM);
        return dDh + dB + 2.0 * m * cPhi;
    }

    static double[] pitzerOsmoticPartials(final double m) {
        if (m <= 0.0) {
            return new double[] { 0.0, 0.0, 0.0 };
        }
        final double sqrtM = FastMath.sqrt(m);
        return new double[] { m, m * FastMath.exp(-2.0 * sqrtM), m * m };
    }

    static double[] nrtlBinary(final double x1,
                               final double tau12,
                               final double tau21,
                               final double alpha) {
        final double x2 = 1.0 - x1;
        final double g12 = FastMath.exp(-alpha * tau12);
        final double g21 = FastMath.exp(-alpha * tau21);
        final double lnGamma1 = x2 * x2 *
                (tau21 * FastMath.pow(g21 / (x1 + x2 * g21), 2) +
                 tau12 * g12 / FastMath.pow(x2 + x1 * g12, 2));
        final double lnGamma2 = x1 * x1 *
                (tau12 * FastMath.pow(g12 / (x2 + x1 * g12), 2) +
                 tau21 * g21 / FastMath.pow(x1 + x2 * g21, 2));
        return new double[] { lnGamma1, lnGamma2 };
    }

    static double dNrtlLnGamma1Dx1(final double x1,
                                   final double tau12,
                                   final double tau21,
                                   final double alpha) {
        final double eps = 1.0e-8;
        return (nrtlBinary(x1 + eps, tau12, tau21, alpha)[0] -
                nrtlBinary(x1 - eps, tau12, tau21, alpha)[0]) / (2.0 * eps);
    }

    static double dNrtlLnGamma2Dx1(final double x1,
                                   final double tau12,
                                   final double tau21,
                                   final double alpha) {
        final double eps = 1.0e-8;
        return (nrtlBinary(x1 + eps, tau12, tau21, alpha)[1] -
                nrtlBinary(x1 - eps, tau12, tau21, alpha)[1]) / (2.0 * eps);
    }

    static SimpleBounds bounds(final double[] lower, final double[] upper) {
        return new SimpleBounds(lower, upper);
    }

    static double[] filled(final int dimension, final double value) {
        final double[] array = new double[dimension];
        for (int i = 0; i < dimension; ++i) {
            array[i] = value;
        }
        return array;
    }

    static LagrangeSolution solve(final String name,
                                  final TwiceDifferentiableFunction objective,
                                  final double[] initial,
                                  final SimpleBounds bounds,
                                  final double expectedObjective,
                                  final double objectiveTolerance) {
        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        final SQPOptimizerS2 optimizer = ElectrolyteProblemTestUtils.newOptimizer();
        final long start = System.nanoTime();
        final LagrangeSolution solution = optimizer.optimize(
                new MaxIter(3000),
                new InitialGuess(initial),
                new ObjectiveFunction(objective),
                bounds,
                option);
        printResult(name, solution, expectedObjective, 0.0, start);
        assertExpectedObjective(name, expectedObjective, objectiveTolerance, solution);
        return solution;
    }

    static LagrangeSolution solve(final String name,
                                  final TwiceDifferentiableFunction objective,
                                  final EqualityConstraint equality,
                                  final double[] initial,
                                  final SimpleBounds bounds,
                                  final double expectedObjective,
                                  final double objectiveTolerance,
                                  final double feasibilityTolerance) {
        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        final SQPOptimizerS2 optimizer = ElectrolyteProblemTestUtils.newOptimizer();
        final long start = System.nanoTime();
        final LagrangeSolution solution = optimizer.optimize(
                new MaxIter(3000),
                new InitialGuess(initial),
                new ObjectiveFunction(objective),
                equality,
                bounds,
                option);
        final double residual = equality.value(solution.getX()).getLInfNorm();
        printResult(name, solution, expectedObjective, residual, start);
        assertExpectedObjective(name, expectedObjective, objectiveTolerance, solution);
        assertTrue(residual <= feasibilityTolerance,
                   name + " equality residual = " + residual);
        return solution;
    }

    private static void assertExpectedObjective(final String name,
                                                final double expectedObjective,
                                                final double objectiveTolerance,
                                                final LagrangeSolution solution) {
        final double actualObjective = solution.getValue();
        final double absoluteError = FastMath.abs(actualObjective - expectedObjective);

        assertTrue(Double.isFinite(actualObjective),
                   name + " returned a non-finite objective.");

        assertTrue(absoluteError <= objectiveTolerance,
                   name + " expected objective = " + expectedObjective +
                   ", actual objective = " + actualObjective +
                   ", absolute error = " + absoluteError +
                   ", tolerance = " + objectiveTolerance);
    }

    private static void printResult(final String name,
                                    final LagrangeSolution solution,
                                    final double expectedObjective,
                                    final double residual,
                                    final long start) {
        final double absoluteError =
                FastMath.abs(solution.getValue() - expectedObjective);

        System.out.printf(Locale.US,
                          "%n%s%n" +
                          "objective          : %.16e%n" +
                          "expected objective : %.16e%n" +
                          "absolute error     : %.16e%n" +
                          "constraint LInf    : %.16e%n" +
                          "time               : %.6f s%n",
                          name,
                          solution.getValue(),
                          expectedObjective,
                          absoluteError,
                          residual,
                          (System.nanoTime() - start) * 1.0e-9);
    }

}
