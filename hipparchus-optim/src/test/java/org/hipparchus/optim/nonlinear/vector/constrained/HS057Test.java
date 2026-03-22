/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.hipparchus.linear.MatrixUtils;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS057Test {

    /** Obiettivo: somma dei quadrati F(i)^2 come in TP57 MODE=2, con F(i) definita in MODE=17. */
    private static class HS057Obj extends TwiceDifferentiableFunction {
        private final double[] A; // length 44
        private final double[] B; // length 44

        HS057Obj() {
            double[][] ab = buildAB();
            this.A = ab[0];
            this.B = ab[1];
        }

        @Override
        public int dim() {
            return 2;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            double sum = 0.0;
            for (int i = 0; i < 44; i++) {
                double Fi = B[i] - x1
                            - (0.49 - x1) * FastMath.exp(-x2 * (A[i] - 8.0));
                sum += Fi * Fi;
            }
            return sum;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // Traduzione diretta del blocco MODE=3 in TP57
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            double g1 = 0.0;
            double g2 = 0.0;

            for (int i = 0; i < 44; i++) {
                final double ai = A[i];
                final double bi = B[i];

                final double v1 = FastMath.exp(-x2 * (ai - 8.0));                // V1
                final double Fi =
                        bi - x1 - (0.49 - x1) * v1;                             // F(I)

                final double dFdx1 = -1.0 + v1;                                  // DF(I,1)
                final double dFdx2 = (ai - 8.0) * (0.49 - x1) * v1;              // DF(I,2)

                // S(J) = sum 2*F(I)*DF(I,J)  -> GF(J) = S(J)
                g1 += 2.0 * Fi * dFdx1;
                g2 += 2.0 * Fi * dFdx2;
            }

            return new ArrayRealVector(new double[] { g1, g2 });
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Puoi implementare l’Hessiana analitica, ma per HS057 non è strettamente necessario
            throw new UnsupportedOperationException("Hessian not implemented for HS057.");
        }
    }


    /** Un solo vincolo di disuguaglianza NONLINEARE: g(x) = -x1*x2 + 0.49*x2 - 0.09. */
    private static class HS057Ineq extends InequalityConstraint {
        HS057Ineq() {
            super(new ArrayRealVector(new double[] { 0.0 }));
        }

        @Override
        public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            double g = -x1 * x2 + 0.49 * x2 - 0.09;
            // Vincolo nella forma g(x) >= 0 come in Fortran (A*X + B >= 0)
            return new ArrayRealVector(new double[] { g });
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            // Traduzione diretta di GG in MODE=5:
            // GG(1,1) = -X(2)
            // GG(1,2) = -X(1) + 0.49D0
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            double[][] data = new double[1][2];
            data[0][0] = -x2;
            data[0][1] = -x1 + 0.49;

            return MatrixUtils.createRealMatrix(data);
        }

        @Override
        public int dim() {
            return 2;
        }
    }

    @Test
    public void testHS057() {
        // Guess da MODE=1
        InitialGuess guess = new InitialGuess(new double[] { 0.42, 5.0 });

        // Bounds separati (MODE=1): LXL(1)=TRUE, XL(1)=0.4 ; LXL(2)=TRUE, XL(2)=-4 ; senza upper
        SimpleBounds bounds = new SimpleBounds(
                new double[] { 0.4, -4.0 },
                new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY }
        );

        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);
        SQPOption sqpOpt=new SQPOption();
        sqpOpt.setGradientMode(GradientMode.EXTERNAL);
        sqpOpt.setEps(1.0e-11);
        double expected = 0.0284596697213; // FEX in Fortran

        LagrangeSolution sol = optimizer.optimize(
                sqpOpt,
                guess,
                new ObjectiveFunction(new HS057Obj()),
                new HS057Ineq(),
                bounds
        );

        assertEquals(expected, sol.getValue(), 1e-8);
    }

    /**
     * Ricostruisce gli array A(1..44) e B(1..44) come in TP57 MODE=18 (Fortran),
     * restituendo versioni 0-based: A[0..43], B[0..43].
     */
    private static double[][] buildAB() {
        double[] A = new double[44];
        double[] B = new double[44];

        // Helper per tradurre gli indici 1-based Fortran -> 0-based Java
        java.util.function.BiConsumer<Integer, Double> setA = (i, v) -> A[i - 1] = v;
        java.util.function.BiConsumer<Integer, Double> setB = (i, v) -> B[i - 1] = v;

        // DO 20 I=1,2
        for (int I = 1; I <= 2; I++) {
            setA.accept(I,       8.0);
            setA.accept(16 + I, 18.0);
            setA.accept(30 + I, 28.0);
            setA.accept(35 + I, 32.0);
            setA.accept(38 + I, 36.0);
            setA.accept(40 + I, 38.0);

            setB.accept(I,        0.49);
            setB.accept(6 + I,    0.46);
            setB.accept(11 + I,   0.43);
            setB.accept(14 + I,   0.43);
            setB.accept(18 + I,   0.42);
            setB.accept(21 + I,   0.41);
            setB.accept(25 + I,   0.40);
            setB.accept(29 + I,   0.41);
            setB.accept(36 + I,   0.40);
            setB.accept(40 + I,   0.40);
            setB.accept(42 + I,   0.39);
        }

        // DO 21 I=1,3
        for (int I = 1; I <= 3; I++) {
            setA.accept(10 + I, 14.0);  // 11..13
            setA.accept(13 + I, 16.0);  // 14..16
            setA.accept(18 + I, 20.0);  // 19..21
            setA.accept(21 + I, 22.0);  // 22..24
            setA.accept(24 + I, 24.0);  // 25..27
            setA.accept(27 + I, 26.0);  // 28..30
            setA.accept(32 + I, 30.0);  // 33..35

            setB.accept(31 + I, 0.40);  // 32..34
        }

        // DO 22 I=1,4
        for (int I = 1; I <= 4; I++) {
            setA.accept(2 + I, 10.0);   // 3..6
            setA.accept(6 + I, 12.0);   // 7..10
        }

        setA.accept(38, 34.0);
        setA.accept(43, 40.0);
        setA.accept(44, 42.0);

        setB.accept(3,  0.48);
        setB.accept(4,  0.47);
        setB.accept(5,  0.48);
        setB.accept(6,  0.47);
        setB.accept(9,  0.45);
        setB.accept(10, 0.43);
        setB.accept(11, 0.45);
        setB.accept(14, 0.44);
        setB.accept(17, 0.46);
        setB.accept(18, 0.45);
        setB.accept(21, 0.43);
        setB.accept(24, 0.40);
        setB.accept(25, 0.42);
        setB.accept(28, 0.41);
        setB.accept(29, 0.40);
        setB.accept(35, 0.38);
        setB.accept(36, 0.41);
        setB.accept(39, 0.41);
        setB.accept(40, 0.38);

        return new double[][] { A, B };
    }
}
