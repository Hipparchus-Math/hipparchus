package org.hipparchus.optim.nonlinear.vector.constrained.sqp.CUTEst;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

/**
 * CUTEst CHEBYQAD with the original dimension N = 10.
 */
public class CHEBYQADTest {

    private static final int N = 10;
    private static final double EXPECTED_OBJECTIVE = 4.772713e-3;

    private static final class Objective extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            final double[] residual = residuals(point);
            double value = 0.0;
            for (final double r : residual) {
                value += r * r;
            }
            return value;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            final double[] x = point.toArray();
            final double[] residual = residuals(point);
            final double[] gradient = new double[N];

            for (int variable = 0; variable < N; ++variable) {
                final double shifted = 2.0 * x[variable] - 1.0;

                double tPrevious = 1.0;
                double tCurrent = shifted;
                double dPrevious = 0.0;
                double dCurrent = 2.0;

                gradient[variable] +=
                        2.0 * residual[0] * dCurrent / N;

                for (int degree = 2; degree <= N; ++degree) {
                    final double tNext =
                            2.0 * shifted * tCurrent - tPrevious;
                    final double dNext =
                            4.0 * tCurrent +
                            2.0 * shifted * dCurrent -
                            dPrevious;

                    gradient[variable] +=
                            2.0 * residual[degree - 1] * dNext / N;

                    tPrevious = tCurrent;
                    tCurrent = tNext;
                    dPrevious = dCurrent;
                    dCurrent = dNext;
                }
            }

            return new ArrayRealVector(gradient, false);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            return new Array2DRowRealMatrix(N, N);
        }

        private static double[] residuals(final RealVector point) {
            final double[] residual = new double[N];

            for (int variable = 0; variable < N; ++variable) {
                final double shifted =
                        2.0 * point.getEntry(variable) - 1.0;

                double tPrevious = 1.0;
                double tCurrent = shifted;

                residual[0] += tCurrent / N;

                for (int degree = 2; degree <= N; ++degree) {
                    final double tNext =
                            2.0 * shifted * tCurrent - tPrevious;

                    residual[degree - 1] += tNext / N;
                    tPrevious = tCurrent;
                    tCurrent = tNext;
                }
            }

            for (int degree = 2; degree <= N; degree += 2) {
                residual[degree - 1] +=
                        1.0 / (degree * degree - 1.0);
            }

            return residual;
        }
    }

    private static double[] initialPoint() {
        final double[] start = new double[N];
        for (int i = 0; i < N; ++i) {
            start[i] = (i + 1.0) / (N + 1.0);
        }
        return start;
    }

    @Test
    public void testCHEBYQAD() {
        final SQPOptimizerS2 optimizer = CUTEstProblemUtils.newOptimizer();
        final SQPOption option = CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution = optimizer.optimize(
                new MaxIter(5000),
                new InitialGuess(initialPoint()),
                new ObjectiveFunction(new Objective()),
                new SimpleBounds(new double[N], ones(N)),
                option);

        CUTEstProblemUtils.assertExpectedObjective(EXPECTED_OBJECTIVE, solution);
    }

    private static double[] ones(final int dimension) {
        final double[] values = new double[dimension];
        for (int i = 0; i < dimension; ++i) {
            values[i] = 1.0;
        }
        return values;
    }
}