/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Disabled;

/**
 * HS362 / TP362 – time-to-speed car model with 5 gear ratios (N = 5).
 * <p>
 * Obiettivo: minimizzare il tempo T/100 necessario a passare da 5 mph a 100 mph,
 * simulato dalla funzione esterna tp362A(x).
 * <p>
 * Vincoli:
 * <ul>
 *   <li>Box bounds: 15 ≤ x₁ ≤ 20; 3 ≤ x₂,x₃,x₄ ≤ 20; 2 ≤ x₅ ≤ 20</li>
 *   <li>Linear inequalities: x₁ ≥ x₂ ≥ x₃ ≥ x₄ ≥ x₅</li>
 * </ul>
 */
public class HS362Test {

    /**
     * Obiettivo: TP362A(X) – tempo di accelerazione T/100.
     * Traduzione diretta della routine di simulazione (nessun gradiente analitico).
     */
    private static class TP362Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return 5;
        }

        @Override
        public double value(RealVector x) {
            return tp362A(x.toArray());
        }

        @Override
        public RealVector gradient(RealVector x) {
            throw new UnsupportedOperationException("TP362: no analytic gradient");
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("TP362: no analytic Hessian");
        }

        /**
         * Calcola la coppia motore in funzione dei giri/min (rpm) tramite
         * una curva a tratti cubica.
         */
        private double calculateTorque(double rpm) {
            if (rpm >= 600.0 && rpm < 1900.0) {
                return  0.3846154e-7  * Math.pow(rpm, 3)
                      - 0.2108974359e-3 * Math.pow(rpm, 2)
                      + 0.42455128205133 * rpm
                      - 187.11538461540295;
            } else if (rpm >= 1900.0 && rpm < 3000.0) {
                return -0.492424e-8    * Math.pow(rpm, 3)
                      + 0.1867424242e-4 * Math.pow(rpm, 2)
                      + 0.1229545454547e-1 * rpm
                      + 64.999999999986;
            } else if (rpm >= 3000.0 && rpm < 4500.0) {
                return -0.26667e-9  * Math.pow(rpm, 3)
                      + 0.3e-5      * Math.pow(rpm, 2)
                      - 0.1263333333336e-1 * rpm
                      + 155.10000000002947;
            } else if (rpm >= 4500.0 && rpm < 5600.0) {
                return -0.664141e-8    * Math.pow(rpm, 3)
                      + 0.8337626263e-4 * Math.pow(rpm, 2)
                      - 0.34351868688129 * rpm
                      + 597.363636363847145;
            } else if (rpm >= 5600.0 && rpm < 6000.0) {
                return -0.2539683e-7    * Math.pow(rpm, 3)
                      + 0.38158730157e-3 * Math.pow(rpm, 2)
                      - 1.9223492062348  * rpm
                      + 3380.66666645715304;
            }
            return 0.0;
        }

        /**
         * Simula l'accelerazione del veicolo da 5 mph a 100 mph.
         * Ritorna T/100 (come nella routine originale).
         */
        private double tp362A(double[] X) {

            // Costanti del modello
            final double RAD    = 1.085;
            final double CON1   = 1.466667;
            final double CON2   = 12.90842;
            final double RPMIN  = 600.0;
            final double RPMAX  = 5700.0;
            final double EI     = 0.6;
            final double VI     = 98.0;
            final double DT     = 0.01;
            final double VMAX   = 100.0;
            final double V0     = 5.0;
            final double TSHIFT = 0.25;
            final double TMAX   = 100.0;
            final int    DIM    = 5;

            int it = 0;
            double acc  = 0.0;
            double acc0;
            double v    = V0;
            double t    = 0.0;
            int gearIdx = 0;

            while (true) {
                // 302: forza resistente (aerodinamica + rotolamento)
                double force = 0.0239 * v * v + 31.2;

                // 301: calcolo RPM per il rapporto corrente
                if (gearIdx >= DIM) {
                    // Nessun rapporto disponibile → termina
                    return t / 100.0;
                }

                double xi  = X[gearIdx];
                double rpm = v * CON2 * xi;

                // 300: regime troppo basso → penalizza con TMAX
                if (rpm < RPMIN) {
                    return TMAX;  // nessuna divisione per 100 qui, come nell'originale
                }

                // 305: superato regime massimo → cambio marcia
                if (rpm >= RPMAX) {
                    gearIdx++;
                    if (gearIdx >= DIM) {
                        return t / 100.0;
                    }

                    if (t == 0.0) {
                        // all'inizio, passa direttamente al nuovo rapporto
                        continue;
                    }

                    double tt = t + TSHIFT;

                    // 306–307: fase di cambio marcia (decelerazione)
                    while (true) {
                        // 306: decelerazione per sola resistenza
                        double accShift = -force * RAD * RAD / VI;
                        it++;
                        t = DT * it;
                        v = v + accShift * DT / CON1;

                        if (t < tt) {
                            // 307: durante il cambio, cambia modello di forza
                            force = 0.0293 * v * v + 31.2;
                            continue;
                        } else {
                            // Esci dalla fase di cambio e torna al ciclo principale (302)
                            break;
                        }
                    }

                    if (t > TMAX || v >= VMAX) {
                        return t / 100.0;
                    }

                    // Ricomincia con il nuovo rapporto
                    continue;
                }

                // Fase di accelerazione normale
                double torque = calculateTorque(rpm);

                acc0 = acc;
                acc  = RAD * (xi * torque - force * RAD) / (EI * xi * xi + VI);

                it++;
                t = DT * it;

                // Integrazione trapezoidale su v
                v = v + (acc0 + acc) * 0.5 * DT / CON1;

                if (t > TMAX || v >= VMAX) {
                    return t / 100.0;
                }
            }
        }
    }

    /**
     * Vincoli lineari: x₁ ≥ x₂ ≥ x₃ ≥ x₄ ≥ x₅.
     */
    private static class TP362Ineq extends InequalityConstraint {

        TP362Ineq() {
            super(new ArrayRealVector(new double[] {0.0, 0.0, 0.0, 0.0}));
        }

        @Override
        public int dim() {
            return 5;
        }

        @Override
        public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[] {
                x.getEntry(0) - x.getEntry(1),
                x.getEntry(1) - x.getEntry(2),
                x.getEntry(2) - x.getEntry(3),
                x.getEntry(3) - x.getEntry(4)
            }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix j = new Array2DRowRealMatrix(4, 5);
            j.setEntry(0, 0,  1.0); j.setEntry(0, 1, -1.0);
            j.setEntry(1, 1,  1.0); j.setEntry(1, 2, -1.0);
            j.setEntry(2, 2,  1.0); j.setEntry(2, 3, -1.0);
            j.setEntry(3, 3,  1.0); j.setEntry(3, 4, -1.0);
            return j;
        }
    }

    /**
     * Test di ottimizzazione completo per HS362/TP362.
     * LEX = .FALSE. → si verifica solo che il valore trovato sia ≤ FEX.
     */
//    @Disabled
    @Test
    public void testHS362_optimization() {

        // Punto iniziale (X0)
        double[] x0 = {15.1, 9.05, 6.14, 4.55, 3.61};

        // Bound: XL(i)=3, XU(i)=20; con override XL(1)=15, XL(5)=2
        double[] lower = {15.0, 3.0, 3.0, 3.0, 2.0};
        double[] upper = {20.0, 20.0, 20.0, 20.0, 20.0};

        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new TP362Obj()),
            new TP362Ineq(),
            bounds
        );

        double f = sol.getValue();

        // FEX = 0.418D-01; LEX = .FALSE. → si richiede solo FEX >= f
        final double fExpected = 0.0418;
        assertTrue(fExpected >= f,
                   "HS362: expected F <= " + fExpected + " but got F = " + f);
    }
}
