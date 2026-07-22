package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

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
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS117Test {

    // e(1..5)
    private static final double[] E = { -15.0, -27.0, -36.0, -18.0, -12.0 };

    // d(1..5)
    private static final double[] D = { 4.0, 8.0, 10.0, 6.0, 2.0 };

    // c(1..5,1..5) simmetrica
    private static final double[][] C = new double[5][5];
    static {
        C[0][0]=30;   C[0][1]=-20; C[0][2]=-10; C[0][3]=32;  C[0][4]=-10;
        C[1][0]=-20;  C[1][1]=39;  C[1][2]=-6;  C[1][3]=-31; C[1][4]=32;
        C[2][0]=-10;  C[2][1]=-6;  C[2][2]=10;  C[2][3]=-6;  C[2][4]=-10;
        C[3][0]=32;   C[3][1]=-31; C[3][2]=-6;  C[3][3]=39;  C[3][4]=-20;
        C[4][0]=-10;  C[4][1]=32;  C[4][2]=-10; C[4][3]=-20; C[4][4]=30;
    }

    // a(1..10,1..5)
    private static final double[][] A = new double[10][5];
    static {
        // row 1
        A[0][0]=-16;  A[0][1]= 2;   A[0][2]= 0;   A[0][3]= 1;   A[0][4]= 0;
        // row 2
        A[1][0]= 0;   A[1][1]=-2;   A[1][2]= 0;   A[1][3]= 0.4; A[1][4]= 2;  // <-- 0.4
        // row 3
        A[2][0]=-3.5; A[2][1]= 0;   A[2][2]= 2;   A[2][3]= 0;   A[2][4]= 0;
        // row 4
        A[3][0]= 0;   A[3][1]=-2;   A[3][2]= 0;   A[3][3]=-4;   A[3][4]=-1;
        // row 5
        A[4][0]= 0;   A[4][1]=-9;   A[4][2]=-2;   A[4][3]= 1;   A[4][4]=-2.8;
        // row 6
        A[5][0]= 2;   A[5][1]= 0;   A[5][2]=-4;   A[5][3]= 0;   A[5][4]= 0;
        // row 7: -1
        for (int j=0;j<5;j++) A[6][j] = -1;
        // row 8
        A[7][0]=-1;   A[7][1]=-2;   A[7][2]=-3;   A[7][3]=-2;   A[7][4]=-1;
        // row 9: 1..5
        for (int j=0;j<5;j++) A[8][j] = j+1;
        // row 10: 1
        for (int j=0;j<5;j++) A[9][j] = 1;
    }

    // b(1..10)
    private static final double[] B = { -40, -2, -0.25, -4, -4, -1, -40, -60, 5, 1 };

    // Bounds: x_i >= 0 (nessun upper bound)
    private static final double[] LB = new double[15];
    private static final double[] UB = new double[15];
    static {
        for (int i=0;i<15;i++){ LB[i]=0.0; UB[i]=Double.POSITIVE_INFINITY; }
    }

    // Obiettivo: -sum b_i x_i + y^T C y + 2 sum d_j y_j^3, con y_j=x_{10+j}
    private static class Obj extends TwiceDifferentiableFunction {
        @Override public int dim(){ return 15; }
        @Override public double value(RealVector X){
            double t3 = 0.0;
            for (int i=0;i<10;i++) t3 += B[i]*X.getEntry(i);
            double t1 = 0.0;
            for (int i=0;i<5;i++){
                final double yi = X.getEntry(10+i);
                for (int j=0;j<5;j++){
                    t1 += C[i][j]*yi*X.getEntry(10+j);
                }
            }
            double t2 = 0.0;
            for (int j=0;j<5;j++){
                final double y = X.getEntry(10+j);
                t2 += D[j]*y*y*y;
            }
            return -t3 + t1 + 2.0*t2;
        }
        @Override public RealVector gradient(RealVector x){ throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x){ throw new UnsupportedOperationException(); }
    }

    // 5 disuguaglianze: 2*C*y + 3*d .* y.^2 + e - A^T x >= 0
    private static class Ineq extends InequalityConstraint {
        Ineq(){ super(new ArrayRealVector(new double[5])); }
        @Override public int dim(){ return 15; }
        @Override public RealVector value(RealVector X){
            final double[] y = new double[5];
            for (int j=0;j<5;j++) y[j] = X.getEntry(10+j);

            // T4 = C*y
            final double[] T4 = new double[5];
            for (int j=0;j<5;j++){
                double s=0;
                for (int i=0;i<5;i++) s += C[i][j]*y[i];
                T4[j]=s;
            }

            // T5_j = sum_i A(i,j) x_i
            final double[] T5 = new double[5];
            for (int j=0;j<5;j++){
                double s=0;
                for (int i=0;i<10;i++) s += A[i][j]*X.getEntry(i);
                T5[j]=s;
            }

            final double[] g = new double[5];
            for (int j=0;j<5;j++){
                g[j] = 2.0*T4[j] + 3.0*D[j]*y[j]*y[j] + E[j] - T5[j];
            }
            return new ArrayRealVector(g);
        }
        @Override public RealMatrix jacobian(RealVector x){ throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS117_fixed(){
        // Guess come nel model
        final double[] x0 = new double[15];
        for (int i=0;i<15;i++) x0[i]=0.001;
        x0[6] = 60.0; // x7

        final InitialGuess guess = new InitialGuess(x0);
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = opt.optimize(
                guess,
                new ObjectiveFunction(new Obj()),
                new Ineq(),
                bounds
        );

        // Best known objective ≈ 32.34867897
        HSProblemTestUtils.assertExpectedObjective(32.34867897, sol);
    }
}
