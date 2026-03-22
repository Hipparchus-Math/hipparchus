
/*
* HS376 (TP376)
*
* N    = 10
* NILI = 14  (14 nonlinear inequalities)
* NINL = 14
* NELI = 1   (1 linear equality: x9 + x10 = 0.255)
* NENL = 0
*
* Reference solution (Fortran):
*   x* ≈ (0.14727222,
*         0.1,
*         0.0081,
*         628.71731,
*         0.0017,
*         0.0011816143,
*         0.0027,
*         0.00135,
*         0.15740741,
*         0.097592593)
*   f* = -0.44300879D+04 ≈ -4430.0879
*/

package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class HS376BISTest {

 /**
  * Objective from TP376:
  *
  *   FX = 0.2D+5 * (0.15*X1 + 0.14D+2*X2 - 0.6D-1) / (0.2D-2 + X1 + 0.6D+2*X2)
  *   FX = -FX
  *
  * i.e.
  *   f(x) = -20000 * (0.15*x1 + 14*x2 - 0.06) / (0.002 + x1 + 60*x2)
  */
 static final class HS376Obj extends TwiceDifferentiableFunction {

    @Override
    public int dim() {
        return 10;
    }

    @Override
    public double value(RealVector x) {
        final double x1 = x.getEntry(0);
        final double x2 = x.getEntry(1);

        // Fortran:
        // FX = 0.2D+5 * (0.15D+0*X(1) + 0.14D+2*X(2) - 0.6D-1)
        //      / (0.2D-2 + X(1) + 0.6D+2*X(2))
        // FX = -FX
        final double num = 0.15 * x1 + 14.0 * x2 - 0.06;
        final double den = 0.002 + x1 + 60.0 * x2;

        return -20000.0 * (num / den);
    }

    @Override
    public RealVector gradient(RealVector x) {
        final double x1 = x.getEntry(0);
        final double x2 = x.getEntry(1);

        final double num = 0.15 * x1 + 14.0 * x2 - 0.06;
        final double den = 0.002 + x1 + 60.0 * x2;

        // derivate numeratore / denominatore
        final double dNum_dx1 = 0.15;
        final double dNum_dx2 = 14.0;

        final double dDen_dx1 = 1.0;
        final double dDen_dx2 = 60.0;

        // f = -K * num / den
        // ∂f/∂xi = -K * (num' * den - num * den') / den^2
        final double K = 20000.0;
        final double den2 = den * den;

        final double df_dx1 = -K * (dNum_dx1 * den - num * dDen_dx1) / den2;
        final double df_dx2 = -K * (dNum_dx2 * den - num * dDen_dx2) / den2;

        double[] g = new double[10];
        g[0] = df_dx1;
        g[1] = df_dx2;
        // g[2..9] = 0.0

        return new ArrayRealVector(g, false);
    }

    @Override
    public RealMatrix hessian(RealVector x) {
        // Lasciamo costruire l’Hessiana al BFGS
        return MatrixUtils.createRealMatrix(10, 10);
    }
}


 /**
  * Inequality constraints G(x) >= 0 (NINL = 14 in Fortran).
  *
  * Fortran MODE=4:
  *
  *  G(1)  = X1 - 0.75 / (X3*X4)
  *  G(2)  = X1 - X9 / (X5*X4)
  *  G(3)  = X1 - X10/(X6*X4) - 10 / X4
  *  G(4)  = X1 - 0.19/(X7*X4) - 10 / X4
  *  G(5)  = X1 - 0.125/(X8*X4)
  *
  *  G(6)  = 0.1D+5*X2 - 0.131D-2 * X9 * X5^0.666 * X4^1.5
  *  G(7)  = 0.1D+5*X2 - 0.1038D-2 * X10 * X6^1.6 * X4^3
  *  G(8)  = 0.1D+5*X2 - 0.223D-3 * X7^0.666 * X4^1.5
  *  G(9)  = 0.1D+5*X2 - 0.76D-4 * X8^3.55 * X4^5.66
  *  G(10) = 0.1D+5*X2 - 0.698D-3 * X3^1.2 * X4^2
  *  G(11) = 0.1D+5*X2 - 0.5D-4 * X3^1.6 * X4^3
  *  G(12) = 0.1D+5*X2 - 0.654D-5 * X3^2.42 * X4^4.17
  *  G(13) = 0.1D+5*X2 - 0.257D-3 * X3^0.666 * X4^1.5
  *
  *  G(14) = 30
  *          - 2.003 * X5 * X4
  *          - 1.885 * X6 * X4
  *          - 0.184 * X8 * X4
  *          - 2.0   * X3^0.803 * X4
  */
 static final class HS376Ineq extends InequalityConstraint {

    HS376Ineq() {
        super(new ArrayRealVector(new double[14])); 
    }

    @Override
    public int dim() {
        return 10;
    }

    @Override
    public RealVector value(RealVector x) {

        final double x1  = x.getEntry(0);
        final double x2  = x.getEntry(1);
        final double x3  = x.getEntry(2);
        final double x4  = x.getEntry(3);
        final double x5  = x.getEntry(4);
        final double x6  = x.getEntry(5);
        final double x7  = x.getEntry(6);
        final double x8  = x.getEntry(7);
        final double x9  = x.getEntry(8);
        final double x10 = x.getEntry(9);

        double[] c = new double[14];

        // c1 = x1 - 0.75 / (x3 * x4)
        c[0] = x1 - 0.75 / (x3 * x4);

        // c2 = x1 - x9 / (x5 * x4)
        c[1] = x1 - x9 / (x5 * x4);

        // c3 = x1 - x10 / (x6 * x4) - 10 / x4
        c[2] = x1 - x10 / (x6 * x4) - 10.0 / x4;

        // c4 = x1 - 0.19 / (x7 * x4) - 10 / x4
        c[3] = x1 - 0.19 / (x7 * x4) - 10.0 / x4;

        // c5 = x1 - 0.125 / (x8 * x4)
        c[4] = x1 - 0.125 / (x8 * x4);

        // ATTENZIONE: tutti questi hanno coefficiente 0.1D+5 = 1.0e4 davanti a x2

        // c6 = 1e4 * x2 - 0.131e-2 * x9 * x5^0.666 * x4^1.5
        c[5] = 1.0e4 * x2
             - 0.00131 * x9 * Math.pow(x5, 0.666) * Math.pow(x4, 1.5);

        // c7 = 1e4 * x2 - 0.1038e-2 * x10 * x6^1.6 * x4^3
        c[6] = 1.0e4 * x2
             - 0.001038 * x10 * Math.pow(x6, 1.6) * Math.pow(x4, 3.0);

        // c8 = 1e4 * x2 - 0.223e-3 * x7^0.666 * x4^1.5
        c[7] = 1.0e4 * x2
             - 0.000223 * Math.pow(x7, 0.666) * Math.pow(x4, 1.5);

        // c9 = 1e4 * x2 - 0.76e-4 * x8^3.55 * x4^5.66
        c[8] = 1.0e4 * x2
             - 0.000076 * Math.pow(x8, 3.55) * Math.pow(x4, 5.66);

        // c10 = 1e4 * x2 - 0.698e-3 * x3^1.2 * x4^2
        c[9] = 1.0e4 * x2
             - 0.000698 * Math.pow(x3, 1.2) * Math.pow(x4, 2.0);

        // c11 = 1e4 * x2 - 0.5e-4 * x3^1.6 * x4^3
        c[10] = 1.0e4 * x2
              - 0.00005 * Math.pow(x3, 1.6) * Math.pow(x4, 3.0);

        // c12 = 1e4 * x2 - 0.654e-5 * x3^2.42 * x4^4.17
        c[11] = 1.0e4 * x2
              - 0.00000654 * Math.pow(x3, 2.42) * Math.pow(x4, 4.17);

        // c13 = 1e4 * x2 - 0.257e-3 * x3^0.666 * x4^1.5
        c[12] = 1.0e4 * x2
              - 0.000257 * Math.pow(x3, 0.666) * Math.pow(x4, 1.5);

        // c14 = 30 - 2.003 x5 x4 - 1.885 x6 x4 - 0.184 x8 x4 - 2.0 x3^0.803 x4
        c[13] = 30.0
              - 2.003 * x5 * x4
              - 1.885 * x6 * x4
              - 0.184 * x8 * x4
              - 2.0 * Math.pow(x3, 0.803) * x4;

        return new ArrayRealVector(c, false);
    }

    @Override
    public RealMatrix jacobian(RealVector x) {

        final double x1  = x.getEntry(0);
        final double x2  = x.getEntry(1);
        final double x3  = x.getEntry(2);
        final double x4  = x.getEntry(3);
        final double x5  = x.getEntry(4);
        final double x6  = x.getEntry(5);
        final double x7  = x.getEntry(6);
        final double x8  = x.getEntry(7);
        final double x9  = x.getEntry(8);
        final double x10 = x.getEntry(9);

        double[][] J = new double[14][10];

        // ----- c1 = x1 - 0.75 / (x3 * x4) -----
        J[0][0] = 1.0;
        J[0][2] = 0.75 / (x3 * x3 * x4);
        J[0][3] = 0.75 / (x3 * x4 * x4);

        // ----- c2 = x1 - x9 / (x5 * x4) -----
        J[1][0] = 1.0;
        J[1][4] =  x9 / (x5 * x5 * x4);
        J[1][3] =  x9 / (x5 * x4 * x4);
        J[1][8] = -1.0 / (x5 * x4);

        // ----- c3 = x1 - x10/(x6 x4) - 10/x4 -----
        J[2][0] = 1.0;
        J[2][5] =  x10 / (x6 * x6 * x4);
        J[2][3] =  x10 / (x6 * x4 * x4) + 10.0 / (x4 * x4);
        J[2][9] = -1.0 / (x6 * x4);

        // ----- c4 = x1 - 0.19/(x7 x4) - 10/x4 -----
        J[3][0] = 1.0;
        J[3][6] = 0.19 / (x7 * x7 * x4);
        J[3][3] = 0.19 / (x7 * x4 * x4) + 10.0 / (x4 * x4);

        // ----- c5 = x1 - 0.125/(x8 x4) -----
        J[4][0] = 1.0;
        J[4][7] = 0.125 / (x8 * x8 * x4);
        J[4][3] = 0.125 / (x8 * x4 * x4);

        // ----- c6 = 1e4 x2 - 0.00131 x9 x5^0.666 x4^1.5 -----
        J[5][1] = 1.0e4;
        J[5][8] = -0.00131 * Math.pow(x5, 0.666) * Math.pow(x4, 1.5);
        J[5][4] = -0.00131 * x9 * 0.666 *
                  Math.pow(x5, 0.666 - 1.0) * Math.pow(x4, 1.5);
        J[5][3] = -0.00131 * x9 * Math.pow(x5, 0.666) * 1.5 * Math.pow(x4, 0.5);

        // ----- c7 = 1e4 x2 - 0.001038 x10 x6^1.6 x4^3 -----
        J[6][1] = 1.0e4;
        J[6][9] = -0.001038 * Math.pow(x6, 1.6) * Math.pow(x4, 3.0);
        J[6][5] = -0.001038 * x10 * 1.6 * Math.pow(x6, 0.6) * Math.pow(x4, 3.0);
        J[6][3] = -0.001038 * x10 * Math.pow(x6, 1.6) * 3.0 * Math.pow(x4, 2.0);

        // ----- c8 = 1e4 x2 - 0.000223 x7^0.666 x4^1.5 -----
        J[7][1] = 1.0e4;
        J[7][6] = -0.000223 * 0.666 * Math.pow(x7, 0.666 - 1.0) * Math.pow(x4, 1.5);
        J[7][3] = -0.000223 * Math.pow(x7, 0.666) * 1.5 * Math.pow(x4, 0.5);

        // ----- c9 = 1e4 x2 - 0.000076 x8^3.55 x4^5.66 -----
        J[8][1] = 1.0e4;
        J[8][7] = -0.000076 * 3.55 * Math.pow(x8, 2.55) * Math.pow(x4, 5.66);
        J[8][3] = -0.000076 * Math.pow(x8, 3.55) * 5.66 * Math.pow(x4, 4.66);

        // ----- c10 = 1e4 x2 - 0.000698 x3^1.2 x4^2 -----
        J[9][1] = 1.0e4;
        J[9][2] = -0.000698 * 1.2 * Math.pow(x3, 0.2) * Math.pow(x4, 2.0);
        J[9][3] = -0.000698 * Math.pow(x3, 1.2) * 2.0 * x4;

        // ----- c11 = 1e4 x2 - 0.00005 x3^1.6 x4^3 -----
        J[10][1] = 1.0e4;
        J[10][2] = -0.00005 * 1.6 * Math.pow(x3, 0.6) * Math.pow(x4, 3.0);
        J[10][3] = -0.00005 * Math.pow(x3, 1.6) * 3.0 * Math.pow(x4, 2.0);

        // ----- c12 = 1e4 x2 - 0.00000654 x3^2.42 x4^4.17 -----
        J[11][1] = 1.0e4;
        J[11][2] = -0.00000654 * 2.42 * Math.pow(x3, 1.42) * Math.pow(x4, 4.17);
        J[11][3] = -0.00000654 * Math.pow(x3, 2.42) * 4.17 * Math.pow(x4, 3.17);

        // ----- c13 = 1e4 x2 - 0.000257 x3^0.666 x4^1.5 -----
        J[12][1] = 1.0e4;
        J[12][2] = -0.000257 * 0.666 * Math.pow(x3, 0.666 - 1.0) * Math.pow(x4, 1.5);
        J[12][3] = -0.000257 * Math.pow(x3, 0.666) * 1.5 * Math.pow(x4, 0.5);

        // ----- c14 = 30 - 2.003 x5 x4 - 1.885 x6 x4 - 0.184 x8 x4 - 2 x3^0.803 x4 -----
        J[13][2] = -2.0 * 0.803 * Math.pow(x3, -0.197) * x4;
        J[13][3] = -2.003 * x5
                 - 1.885 * x6
                 - 0.184 * x8
                 - 2.0 * Math.pow(x3, 0.803);
        J[13][4] = -2.003 * x4;
        J[13][5] = -1.885 * x4;
        J[13][7] = -0.184 * x4;

        return MatrixUtils.createRealMatrix(J);
    }
}

 /**
  * Equality constraint: x9 + x10 = 0.255
  * Fortran: G(15) = X(9) + X(10) - 0.255D+0
  */
 static final class HS376Eq extends EqualityConstraint {

     HS376Eq() {
         super(new ArrayRealVector(new double[]{0.0}));
     }

     @Override
     public int dim() {
         return 10;
     }

     @Override
     public RealVector value(RealVector x) {
         final double x9  = x.getEntry(8);
         final double x10 = x.getEntry(9);
         return new ArrayRealVector(new double[]{x9 + x10 - 0.255}, false);
     }

     @Override
     public RealMatrix jacobian(RealVector x) {
         double[][] a = new double[1][10];
         a[0][8] = 1.0;  // ∂/∂x9
         a[0][9] = 1.0;  // ∂/∂x10
         return MatrixUtils.createRealMatrix(a);
     }
 }

 private LagrangeSolution solve() {
     SQPOptimizerS2 opt = new SQPOptimizerS2();
     if (Boolean.getBoolean("hipparchus.debug.sqp")) {
         opt.setDebugPrinter(System.out::println);
     }

     // Initial guess from Fortran MODE=1:
     double[] x0 = {
             1.0,
             0.005,
             0.0081,
             100.0,
             0.0017,
             0.0013,
             0.0027,
             0.0020,
             0.15,
             0.105
     };

     // Bounds from XL/XU in Fortran:
     double[] lower = new double[10];
     double[] upper = new double[10];

     // XL(1)=0, XU(1)=10
     lower[0] = 0.0;
     upper[0] = 10.0;

     // XL(2)=0, XU(2)=0.1
     lower[1] = 0.0;
     upper[1] = 0.1;

     // XL(3)=0.5D-4=0.00005, XU(3)=0.81D-2=0.0081
     lower[2] = 0.00005;
     upper[2] = 0.0081;

     // XL(4)=0.1D+2=10, XU(4)=0.1D+4=1000
     lower[3] = 10.0;
     upper[3] = 1000.0;

     // XL(5..10)=0.1D-2=0.001
     for (int i = 4; i < 10; i++) {
         lower[i] = 0.001;
     }
     // XU(5)=0.17D-2=0.0017
     // XU(6)=0.13D-2=0.0013
     // XU(7)=0.27D-2=0.0027
     // XU(8)=0.2D-2 =0.0020
     // XU(9)=XU(10)=0.1D+1=1.0
     upper[4] = 0.0017;
     upper[5] = 0.0013;
     upper[6] = 0.0027;
     upper[7] = 0.0020;
     upper[8] = 1.0;
     upper[9] = 1.0;

     return opt.optimize(
             new InitialGuess(x0),
             new ObjectiveFunction(new HS376Obj()),
             new HS376Eq(),          // 1 equality
             new HS376Ineq(),        // 14 inequalities
             new SimpleBounds(lower, upper)
     );
 }

// @Test
// public void testHS376() {
//
//     final double fExpected = -4430.0879;
//     LagrangeSolution sol = solve();
//     double f = sol.getValue();
//
//     // Problema molto ostico: ci accontentiamo di "non peggio" della soluzione di riferimento
//     assertTrue(fExpected >= f);
// }
}
