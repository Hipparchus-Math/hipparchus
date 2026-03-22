package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.*;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS099Test {

    /** Shared state (P,Q,R,S and their partials). */
    private static final class State {
        final double[] P = new double[8]; // 1..8 used (index 0..7)
        final double[] Q = new double[8];
        final double[] R = new double[8];
        final double[] S = new double[8];
        final double[][] DP = new double[8][7]; // DP(i,j) Fortran -> DP[i-1][j-1]
        final double[][] DQ = new double[8][7];
        final double[][] DR = new double[8][7];
        final double[][] DS = new double[8][7];
    }

    /** Problem data A(1..8), T(1..8). */
    private static final double[] A = new double[] {
            0.0, 50.0, 50.0, 75.0, 75.0, 75.0, 100.0, 100.0
    };
    private static final double[] T = new double[] {
            0.0, 25.0, 50.0, 100.0, 150.0, 200.0, 290.0, 380.0
    };

    /** Compute P,Q,R,S and DP,DQ,DR,DS exactly as in TP99. */
    private static State computeState(final RealVector x) {
        // x is size 7, Fortran X(1..7). Java index iJava = i-1
        final State st = new State();

        // Init as in Fortran: P(1)=Q(1)=R(1)=S(1)=0
        st.P[0] = 0.0; st.Q[0] = 0.0; st.R[0] = 0.0; st.S[0] = 0.0;

        // Init derivative arrays to 0 (Java default already 0).

        // Forward sweep for i = 2..8
        for (int i = 2; i <= 8; i++) {
            final int i1 = i - 1;
            final double xi1 = x.getEntry(i1 - 1); // X(i-1)

            final double V1 = A[i - 1] * Math.sin(xi1) - 32.0;
            final double V2 = A[i - 1] * Math.cos(xi1);
            final double V3 = T[i - 1] - T[i1 - 1];
            final double V4 = 0.5 * V3 * V3;

            st.P[i - 1] = V2 * V4 + V3 * st.R[i1 - 1] + st.P[i1 - 1];
            st.Q[i - 1] = V1 * V4 + V3 * st.S[i1 - 1] + st.Q[i1 - 1];
            st.R[i - 1] = V2 * V3 + st.R[i1 - 1];
            st.S[i - 1] = V1 * V3 + st.S[i1 - 1];
        }

        // Derivatives DP,DQ,DR,DS (i=2..8, j=1..7)
        for (int i = 2; i <= 8; i++) {
            for (int j = 1; j <= 7; j++) {

                final int i1 = i - 1;

                if (j == i - 1) {
                    // J == I-1 branch
                    final double xi1 = x.getEntry(i1 - 1);
                    final double V1 = A[i - 1] * Math.sin(xi1);
                    final double V2 = A[i - 1] * Math.cos(xi1);
                    final double V3 = T[i - 1] - T[i1 - 1];
                    final double V4 = 0.5 * V3 * V3;

                    st.DP[i - 1][i1 - 1] = -V1 * V4 + V3 * st.DR[i1 - 1][i1 - 1] + st.DP[i1 - 1][i1 - 1];
                    st.DQ[i - 1][i1 - 1] =  V2 * V4 + V3 * st.DS[i1 - 1][i1 - 1] + st.DQ[i1 - 1][i1 - 1];
                    st.DR[i - 1][i1 - 1] = -V1 * V3 + st.DR[i1 - 1][i1 - 1];
                    st.DS[i - 1][i1 - 1] =  V2 * V3 + st.DS[i1 - 1][i1 - 1];

                } else if (j < i - 1) {
                    // J < I-1 branch
                    final double V3 = T[i - 1] - T[i1 - 1];
                    st.DP[i - 1][j - 1] = V3 * st.DR[i1 - 1][j - 1] + st.DP[i1 - 1][j - 1];
                    st.DQ[i - 1][j - 1] = V3 * st.DS[i1 - 1][j - 1] + st.DQ[i1 - 1][j - 1];
                    st.DR[i - 1][j - 1] = st.DR[i1 - 1][j - 1];
                    st.DS[i - 1][j - 1] = st.DS[i1 - 1][j - 1];

                } // else j > i-1: remains zero (as in Fortran)
            }
        }

        return st;
    }

    /** Objective: f(x) = -R(8)^2, gradient: -2*R(8)*DR(8,·). */
    private static final class HS099Objective extends TwiceDifferentiableFunction {
        @Override public int dim() { return 7; }

        @Override public double value(RealVector x) {
            final State st = computeState(x);
            final double r8 = st.R[7];
            return -r8 * r8;
        }

        @Override public RealVector gradient(RealVector x) {
            final State st = computeState(x);
            final double r8 = st.R[7];
            final double coeff = -2.0 * r8;
            final double[] g = new double[7];
            for (int j = 0; j < 7; j++) {
                g[j] = coeff * st.DR[7][j];
            }
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            // Not required by the tests/optimizer (let it approximate if needed).
            throw new UnsupportedOperationException();
        }
    }

    /** Inequalities: G1=Q(8)-1e5 >= 0, G2=S(8)-1e3 >= 0. */
    private static final class HS099Ineq extends InequalityConstraint {
        HS099Ineq() {
            super(new ArrayRealVector(new double[] {0.0, 0.0}));
        }

        @Override public RealVector value(RealVector x) {
            final State st = computeState(x);
            return new ArrayRealVector(new double[] { st.Q[7] - 1.0e5, st.S[7] - 1.0e3 }, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            final State st = computeState(x);
            final double[][] J = new double[2][7];
            // Row 0 = ∂Q(8)/∂x_j
            System.arraycopy(st.DQ[7], 0, J[0], 0, 7);
            // Row 1 = ∂S(8)/∂x_j
            System.arraycopy(st.DS[7], 0, J[1], 0, 7);
            return new Array2DRowRealMatrix(J, false);
        }

        @Override public int dim() { return 7; }
    }

    private static InitialGuess guess() {
        // X(i)=0.5 for i=1..7 in Fortran
        return new InitialGuess(new double[] {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5});
    }

    private static SimpleBounds bounds() {
        final double[] lo = {0,0,0,0,0,0,0};
        final double[] up = {1.58,1.58,1.58,1.58,1.58,1.58,1.58};
        return new SimpleBounds(lo, up);
    }

    @Test
    public void testHS099() {
        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);

        final LagrangeSolution sol = optimizer.optimize(
                guess(),
                new ObjectiveFunction(new HS099Objective()),
                new HS099Ineq(),
                bounds()
        );

        // Fortran FEX: -0.831079891516D+09
        assertEquals(-0.831079891516e9, sol.getValue(), 1e5); // tolerant: large magnitude
    }
}
