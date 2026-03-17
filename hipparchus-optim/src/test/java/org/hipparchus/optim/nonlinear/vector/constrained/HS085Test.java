package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS085Test {

    // A(2..17) e B(2..17) come nel Fortran (1-based; indice 0 inutilizzato)
    private static final double[] A = new double[18];
    private static final double[] B = new double[18];
    static {
        A[2]=17.505;    A[3]=11.275;     A[4]=214.228;   A[5]=7.458;
        A[6]=0.961;     A[7]=1.612;      A[8]=0.146;     A[9]=107.99;
        A[10]=922.693;  A[11]=926.832;   A[12]=18.766;   A[13]=1.072163e3;
        A[14]=8.961448e3; A[15]=0.063;   A[16]=7.108433e4; A[17]=2.802713e6;

        B[2]=1.0536667e3; B[3]=35.03;    B[4]=665.585;   B[5]=584.463;
        B[6]=265.916;     B[7]=7.046;    B[8]=0.222;     B[9]=273.366;
        B[10]=1.286105e3; B[11]=1.444046e3; B[12]=537.141; B[13]=3.247039e3;
        B[14]=2.6844086e4; B[15]=0.386;  B[16]=1.4e5;    B[17]=1.2146108e7;
    }

    // Bounds TP85
    private static final double[] LB = { 704.4148, 68.6, 0.0, 193.0, 25.0 };
    private static final double[] UB = { 906.3855, 288.88, 134.75, 287.0966, 84.1988 };

    // --------- Obiettivo (MODE=2): FX = -( ... ) --------------
    private static class Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 5; }
        @Override public double value(RealVector X) {
            final T t = new T(X);
            // FX = -( 5.843e-7*Y17 - 1.17e-4*Y14 - 0.1365 - 2.358e-5*Y13 - 1.502e-6*Y16
            //         - 0.0321*Y12 - 4.324e-3*Y5 - 1e-4*C15/C16 - 37.48*Y2/C12 )
            return -( 5.843e-7 * t.y17
                    - 1.17e-4  * t.y14
                    - 0.1365
                    - 2.358e-5 * t.y13
                    - 1.502e-6 * t.y16
                    - 0.0321   * t.y12
                    - 4.324e-3 * t.y5
                    - 1e-4     * (t.c15 / t.c16)
                    - 37.48    * (t.y2 / t.c12) );
        }
        @Override public RealVector gradient(RealVector x){ throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x){ throw new UnsupportedOperationException(); }
    }

    // --------- Disuguaglianze (MODE=4): 38 vincoli g(x) >= 0 --------------
    private static class Ineq extends InequalityConstraint {
        Ineq(){ super(new ArrayRealVector(new double[38])); }
        @Override public int dim(){ return 5; }
        @Override public RealVector value(RealVector X) {
            final T t = new T(X);
            final double[] g = new double[38];
            int k=0;

            g[k++] = 1.5*t.x2 - t.x3;         // 1
            g[k++] = t.y1 - 213.1;           // 2
            g[k++] = 405.23 - t.y1;          // 3

            // 4..19: Y(i+1) - A(i+1), i=1..16  (Y2..Y17)
            g[k++] = t.y2  - A[2];
            g[k++] = t.y3  - A[3];
            g[k++] = t.y4  - A[4];
            g[k++] = t.y5  - A[5];
            g[k++] = t.y6  - A[6];
            g[k++] = t.y7  - A[7];
            g[k++] = t.y8  - A[8];
            g[k++] = t.y9  - A[9];
            g[k++] = t.y10 - A[10];
            g[k++] = t.y11 - A[11];
            g[k++] = t.y12 - A[12];
            g[k++] = t.y13 - A[13];
            g[k++] = t.y14 - A[14];
            g[k++] = t.y15 - A[15];
            g[k++] = t.y16 - A[16];
            g[k++] = t.y17 - A[17];

            // 20..35: B(i+1) - Y(i+1), i=1..16  (Y2..Y17)
            g[k++] = B[2]  - t.y2;
            g[k++] = B[3]  - t.y3;
            g[k++] = B[4]  - t.y4;
            g[k++] = B[5]  - t.y5;
            g[k++] = B[6]  - t.y6;
            g[k++] = B[7]  - t.y7;
            g[k++] = B[8]  - t.y8;
            g[k++] = B[9]  - t.y9;
            g[k++] = B[10] - t.y10;
            g[k++] = B[11] - t.y11;
            g[k++] = B[12] - t.y12;
            g[k++] = B[13] - t.y13;
            g[k++] = B[14] - t.y14;
            g[k++] = B[15] - t.y15;
            g[k++] = B[16] - t.y16;
            g[k++] = B[17] - t.y17;

            // 36..38: vincoli finali
            g[k++] = t.y4 - 0.28 * t.y5 / 0.72;                 // 36
            g[k++] = 21.0 - 3496.0 * t.y2 / t.c12;              // 37
            g[k++] = 62212.0 / t.c17 - 110.6 - t.y1;            // 38
            return new ArrayRealVector(g);
        }
        @Override public RealMatrix jacobian(RealVector x){ throw new UnsupportedOperationException(); }
    }

    // ---------- Intermedi TP85 (MODE=17): calcolo Y(1..17), C(1..17) ----------
    private static final class T {
        final double x1,x2,x3,x4,x5;
        final double y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17;
        final double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17;

        T(RealVector X){
            x1=X.getEntry(0); x2=X.getEntry(1); x3=X.getEntry(2); x4=X.getEntry(3); x5=X.getEntry(4);

            y1 = x2 + x3 + 41.6;
            c1 = 0.024*x4 - 4.62;
            y2 = 12.5/c1 + 12.0;

            final double V3 = y2 * x1;
            c2 = (3.535e-4*x1 + 0.5311)*x1 + 0.08705*V3;
            c3 = 0.052*x1 + 78.0 + 2.377e-3*V3;

            y3 = c2 / c3;
            y4 = 19.0 * y3;

            final double V1 = x1 - y3;
            c4 = (0.1956*V1/x2 + 0.04782)*V1 + 0.6376*y4 + 1.594*y3;
            c5 = 100.0 * x2;
            c6 = V1 - y4;
            c7 = 0.95 - c4 / c5;

            y5 = c6 * c7;
            final double V2 = y5 + y4;
            y6 = V1 - V2;

            c8 = 0.995 * V2;
            y7 = c8 / y1;
            y8 = c8 / 3798.0;

            c9 = y7 - 0.0663*y7/y8 - 0.3153;
            y9 = 96.82/c9 + 0.321*y1;

            y10 = 1.29*y5 + 1.258*y4 + 2.29*y3 + 1.71*y6;
            y11 = 1.71*x1 - 0.452*y4 + 0.58*y3;

            c10 = 12.3/752.3;
            c11 = 1.74125 * V3;                     // 1.75 * 0.995 == 1.74125
            c12 = 0.995*y10 + 1998.0;

            y12 = c10*x1 + c11/c12;
            y13 = c12 - 1.75*y2;

            final double V4 = y9 + x5;
            y14 = 3623.0 + 64.4*x2 + 58.4*x3 + 1.46312e5 / V4;

            c13 = 0.995*y10 + 60.8*x2 + 48.0*x4 - 0.1121*y14 - 5095.0;
            y15 = y13 / c13;

            y16 = 1.48e5 - 3.31e5*y15 + 40.0*y13 - 61.0*y15*y13;

            c14 = 2.324e3*y10 - 2.874e7*y2;
            y17 = 1.413e7 - 1.328e3*y10 - 531.0*y11 + c14/c12;

            c15 = y13/y15 - y13/0.52;
            c16 = 1.104 - 0.72*y15;
            c17 = V4;
        }
    }

    @Test
    public void testHS085() {
        final InitialGuess guess = new InitialGuess(new double[]{900.0, 80.0, 115.0, 267.0, 27.0});
        final SimpleBounds bounds = new SimpleBounds(LB, UB);
        SQPOption sqpOption=new SQPOption();
        sqpOption.setMaxLineSearchIteration(20);
       
        sqpOption.setGradientMode(GradientMode.FORWARD);
        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = opt.optimize(
                guess,
                new ObjectiveFunction(new Obj()),
                new Ineq(),
                bounds,
                sqpOption
        );

        // FEX dal TP85: -0.19051553D+01
        HSProblemTestUtils.assertExpectedObjective(-1.9051553, sol);
    }
}
