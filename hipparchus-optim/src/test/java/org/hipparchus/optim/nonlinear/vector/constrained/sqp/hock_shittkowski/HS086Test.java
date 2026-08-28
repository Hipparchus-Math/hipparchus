package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS086Test {

    /** Funzione obiettivo: T = Σ(Ei Xi + Di Xi³ + Σ Cji Xi Xj). */
    private static class HS086Obj extends TwiceDifferentiableFunction {
        private final double[] E = {-15, -27, -36, -18, -12};
        private final double[] D = {4, 8, 10, 6, 2};
        private final double[][] C = buildC();

        @Override public int dim() { return 5; }

        @Override
        public double value(RealVector x) {
            double T = 0.0;
            for (int i = 0; i < 5; i++) {
                double T1 = 0.0;
                for (int j = 0; j < 5; j++) {
                    T1 += C[j][i] * x.getEntry(i) * x.getEntry(j);
                }
                T += E[i] * x.getEntry(i) + D[i] * FastMath.pow(x.getEntry(i), 3) + T1;
            }
            return T;
        }

        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }

        private static double[][] buildC() {
            double[][] C = new double[5][5];
            C[0][0] = 30;  C[0][1] = -20;  C[0][2] = -10;  C[0][3] = 32;  C[0][4] = -10;
            C[1][1] = 39;  C[1][2] = -6;   C[1][3] = -31;  C[1][4] = 32;
            C[2][2] = 10;  C[2][3] = -6;   C[2][4] = -10;
            C[3][3] = 39;  C[3][4] = -20;
            C[4][4] = 30;
            // simmetrizza (Cji = Cij)
            for (int i = 0; i < 5; i++) {
                for (int j = i + 1; j < 5; j++) {
                    C[j][i] = C[i][j];
                }
            }
            return C;
        }
    }

   /** 10 vincoli lineari di uguaglianza: G = A·x − B = 0. */
private static class HS086Ineq extends InequalityConstraint {
    private final RealMatrix A;
    private final RealVector B;

    HS086Ineq() {
        super(new ArrayRealVector(new double[]{ // 10 moltiplicatori iniziali = 0
            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
        }));
        this.A = new Array2DRowRealMatrix(buildA());
        this.B = new ArrayRealVector(buildB());
    }

    @Override public RealVector value(RealVector x) {
        // G(x) = A*x - B  (Fortran: G(I)=A(I,1)*X1+...+A(I,5)*X5 - B(I))
        return A.operate(x).subtract(B);
    }

    // Jacobiano: costante = A (10x5)
    @Override public RealMatrix jacobian(RealVector x) { return A; }

    // Dim = numero di variabili (coerente con gli altri test)
    @Override public int dim() { return 5; }

    private static double[][] buildA() {
        return new double[][]{
            {-16,  2,  0,   1,   0},
            {  0, -2,  0, 0.4,   2},
            { -3.5,0,  2,   0,   0},
            {  0, -2,  0,  -4,  -1},
            {  0, -9, -2,   1, -2.8},
            {  2,  0, -4,   0,   0},
            { -1, -1, -1,  -1,  -1},
            { -1, -2, -3,  -2,  -1},
            {  1,  2,  3,   4,   5},
            {  1,  1,  1,   1,   1}
        };
    }

    private static double[] buildB() {
        return new double[]{ -40, -2, -0.25, -4, -4, -1, -40, -60, 5, 1 };
    }
}


    @Test
    public void testHS086() {
        // Guess (MODE=1)
        InitialGuess guess = new InitialGuess(new double[]{ 0.0, 0.0, 0.0, 0.0, 1 });

        // Bounds: x ≥ 0
        SimpleBounds bounds = new SimpleBounds(
                new double[]{0.0, 0.0, 0.0, 0.0, 0.0},
                new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
                             Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}
        );

        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
      

        double expected = -32.3486789716;

        LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new HS086Obj()),
                new HS086Ineq(),
                bounds
        );

        HSProblemTestUtils.assertExpectedObjective(expected, sol);
    }
}
