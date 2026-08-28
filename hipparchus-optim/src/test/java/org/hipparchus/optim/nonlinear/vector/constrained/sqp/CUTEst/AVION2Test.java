/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.CUTEst;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code AVION2}.
 *
 * <p>Dassault France airplane-design problem supplied by A. R. Conn.
 * The problem has 49 bounded variables and 15 linear equality constraints.
 * Its nonlinear objective is the sum of 17 squared engineering residuals.</p>
 *
 * <p>The formulas, variable order, bounds, starting point and reference
 * objective are translated directly from {@code AVION2.SIF}. Objective,
 * constraint and Jacobian derivatives are deliberately evaluated by forward
 * finite differences, consistently with the other SQPOptimizerS2 CUTEst
 * tests.</p>
 */
public class AVION2Test {

    /** Number of variables. */
    private static final int N = 49;

    /** Number of equality constraints. */
    private static final int M = 15;

    /** Best objective value stated in the SIF file. */
    private static final double EXPECTED_OBJECTIVE =
            9.46801297093018e+07;

    // Original SIF variable order.
    private static final int SR       = 0;
    private static final int LR       = 1;
    private static final int PK       = 2;
    private static final int EF       = 3;
    private static final int SX       = 4;
    private static final int LX       = 5;
    private static final int SD       = 6;
    private static final int SK       = 7;
    private static final int ST       = 8;
    private static final int SF       = 9;
    private static final int LF       = 10;
    private static final int AM       = 11;
    private static final int CA       = 12;
    private static final int CB       = 13;
    private static final int SO       = 14;
    private static final int SS       = 15;
    private static final int IMPDER   = 16;
    private static final int IMPK     = 17;
    private static final int IMPFUS   = 18;
    private static final int QI       = 19;
    private static final int PT       = 20;
    private static final int MV       = 21;
    private static final int MC       = 22;
    private static final int MD       = 23;
    private static final int PD       = 24;
    private static final int NS       = 25;
    private static final int VS       = 26;
    private static final int CR       = 27;
    private static final int PM       = 28;
    private static final int DV       = 29;
    private static final int MZ       = 30;
    private static final int VN       = 31;
    private static final int QV       = 32;
    private static final int QF       = 33;
    private static final int IMPTRAIN = 34;
    private static final int IMPMOT   = 35;
    private static final int IMPNMOT  = 36;
    private static final int IMPPET   = 37;
    private static final int IMPPIL   = 38;
    private static final int IMPCAN   = 39;
    private static final int IMPSNA   = 40;
    private static final int MS       = 41;
    private static final int EL       = 42;
    private static final int DE       = 43;
    private static final int DS       = 44;
    private static final int IMPVOIL  = 45;
    private static final int NM       = 46;
    private static final int NP       = 47;
    private static final int NG       = 48;

    /** Official SIF starting point. */
    private static final double[] START = {
        2.7452e+01, 1.5000e+00, 1.0000e+01, 0.0000e+00,
        1.9217e+01, 1.5000e+00, 3.5688e+00, 4.0696e+00,
        3.4315e+01, 8.8025e+01, 5.1306e+00, 0.0000e+00,
       -1.4809e-01, 7.5980e-01, 0.0000e+00, 0.0000e+00,
        1.1470e+02, 5.0000e+02, 1.7605e+03, 2.3256e+03,
        5.6788e+00, 1.4197e+04, 1.2589e+04, 2.8394e+04,
        2.0000e-01, 1.0000e+00, 0.0000e+00, 1.0000e+02,
        1.5000e+01, 0.0000e+00, 5.0000e+02, 1.0000e+01,
        8.1490e+02, 3.1405e+03, 1.9450e+03, 1.9085e+02,
        3.5000e+01, 1.0000e+02, 2.0000e+02, 1.2000e+02,
        7.0000e+02, 1.0000e+03, 4.9367e+00, 0.0000e+00,
        0.0000e+00, 5.0000e+03, 1.0000e+00, 1.0000e+00,
        1.0000e+00
    };

    /** Variable lower bounds in the original SIF order. */
    private static final double[] LOWER = {
        10.0, 0.0, 0.0, 0.0, 7.0, 1.5, 2.0, 2.0, 30.0, 20.0,
        0.001, 0.0, -0.2, 0.1, 0.0, 0.0, 100.0, 500.0, 500.0,
        1000.0, 2.0, 2000.0, 3000.0, 5000.0, 0.2, 1.0, 0.0,
        100.0, 4.0, 0.0, 500.0, 10.0, 250.0, 750.0, 250.0,
        10.0, 35.0, 100.0, 200.0, 120.0, 700.0, 100.0, 2.0,
        0.0, 0.0, 500.0, 1.0, 1.0, 1.0
    };

    /** Variable upper bounds in the original SIF order. */
    private static final double[] UPPER = {
        150.0, 10.0, 10.0, 5.0, 120.0, 8.0, 20.0, 30.0, 500.0,
        200.0, 20.0, 10.0, -0.001, 2.0, 1.0, 2.0, 1000.0,
        5000.0, 5000.0, 20000.0, 30.0, 20000.0, 30000.0,
        50000.0, 0.8, 5.0, 20.0, 400.0, 15.0, 10.0, 10000.0,
        50.0, 5000.0, 15000.0, 3000.0, 5000.0, 70.0, 3000.0,
        400.0, 240.0, 1900.0, 1000.0, 20.0, 1.0, 2.0, 5000.0,
        2.0, 2.0, 2.0
    };

    /**
     * Sum of the 17 squared nonlinear engineering residuals.
     */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {

            final double sr       = point.getEntry(SR);
            final double pk       = point.getEntry(PK);
            final double ef       = point.getEntry(EF);
            final double sx       = point.getEntry(SX);
            final double lx       = point.getEntry(LX);
            final double sd       = point.getEntry(SD);
            final double sk       = point.getEntry(SK);
            final double st       = point.getEntry(ST);
            final double lf       = point.getEntry(LF);
            final double am       = point.getEntry(AM);
            final double ca       = point.getEntry(CA);
            final double cb       = point.getEntry(CB);
            final double so       = point.getEntry(SO);
            final double ss       = point.getEntry(SS);
            final double impder   = point.getEntry(IMPDER);
            final double impk     = point.getEntry(IMPK);
            final double qi       = point.getEntry(QI);
            final double pt       = point.getEntry(PT);
            final double mv       = point.getEntry(MV);
            final double mc       = point.getEntry(MC);
            final double md       = point.getEntry(MD);
            final double pd       = point.getEntry(PD);
            final double vs       = point.getEntry(VS);
            final double cr       = point.getEntry(CR);
            final double pm       = point.getEntry(PM);
            final double dv       = point.getEntry(DV);
            final double mz       = point.getEntry(MZ);
            final double vn       = point.getEntry(VN);
            final double qv       = point.getEntry(QV);
            final double qf       = point.getEntry(QF);
            final double impmot   = point.getEntry(IMPMOT);
            final double ms       = point.getEntry(MS);
            final double el       = point.getEntry(EL);
            final double de       = point.getEntry(DE);
            final double ds       = point.getEntry(DS);
            final double impvoil  = point.getEntry(IMPVOIL);
            final double nm       = point.getEntry(NM);

            // E4: SK - 0.01 * PK * SR.
            final double e4 =
                    sk -
                    0.01 * pk * sr;

            // E6: CA - (SS - SO - CB * LF) / LF^2.
            final double lfSquared =
                    lf * lf;

            final double e6 =
                    ca -
                    (ss - so - cb * lf) /
                    lfSquared;

            // E7: -2 AM + SO + SS + 0.01 EF / LF.
            final double e7 =
                    -2.0 * am +
                    so +
                    ss +
                    0.01 * ef / lf;

            // E8: AM - 0.25 SO CB^2 / CA.
            final double e8 =
                    am -
                    0.25 * so * cb * cb / ca;

            // E9 and E10.
            final double e9 =
                    impder -
                    27.5 * sd -
                    1.3 * sd * sd;

            final double e10 =
                    impk -
                    70.0 * sk +
                    8.6 * sk * sk;

            // E13 and E14.
            final double e13 =
                    qi -
                    1000.0 +
                    mv * mv / 24000.0;

            final double e14 =
                    1000.0 * pt -
                    md * pd;

            // E16. The SIF constant is -2, hence "- gconst" adds 2.
            final double e16 =
                    vn +
                    vs +
                    qf / 790.0 +
                    2.0 -
                    mz / cr +
                    dv * pt;

            // E18.
            final double e18 =
                    impmot -
                    1000.0 * pt / (pm + 20.0) -
                    12.0 * FastMath.sqrt(pt);

            // E26 and E27.
            final double e26 =
                    st -
                    1.25 * sr * nm;

            final double e27 =
                    sr -
                    md / ms;

            // E28: QV - 2.4 SX^(3/2) EL / sqrt(LX).
            final double e28 =
                    qv -
                    2.4 *
                    sx *
                    FastMath.sqrt(sx) *
                    el /
                    FastMath.sqrt(lx);

            // E29 and E30.
            final double e29 =
                    so -
                    0.785 * de * de * pt;

            final double e30 =
                    ss -
                    0.785 * ds * ds * pt;

            // E31.
            final double lfCubed =
                    lfSquared * lf;

            final double e31 =
                    cb -
                    2.0 *
                    (vn - ca * lfCubed) /
                    (lfSquared * (3.0 - so * lf));

            // E32.
            final double ratio =
                    mc * lx /
                    (50.0 * sr * el);

            final double ratioPowerThreeHalves =
                    ratio *
                    FastMath.sqrt(ratio);

            final double e32 =
                    impvoil -
                    1.15 *
                    sx *
                    (15.0 + 0.15 * sx) *
                    (ratioPowerThreeHalves + 8.0);

            return square(e4) +
                   square(e6) +
                   square(e7) +
                   square(e8) +
                   square(e9) +
                   square(e10) +
                   square(e13) +
                   square(e14) +
                   square(e16) +
                   square(e18) +
                   square(e26) +
                   square(e27) +
                   square(e28) +
                   square(e29) +
                   square(e30) +
                   square(e31) +
                   square(e32);
        }

        @Override
        public RealVector gradient(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Gradient is evaluated by finite differences.");
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Hessian is not required by the finite-difference option.");
        }
    }

    /**
     * The 15 linear equalities E1, E2, E3, E5, E11, E12, E15, E17,
     * E19, E20, E21, E22, E23, E24 and E25.
     */
    private static final class AvionEqualities
            extends EqualityConstraint {

        AvionEqualities() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {

            final double[] constraints = new double[M];

            constraints[0] =
                    point.getEntry(SD) -
                    0.13 * point.getEntry(SR);

            constraints[1] =
                    point.getEntry(SX) -
                    0.7 * point.getEntry(SR);

            constraints[2] =
                    point.getEntry(LX) -
                    point.getEntry(LR);

            constraints[3] =
                    point.getEntry(SF) -
                    point.getEntry(ST) -
                    2.0 * point.getEntry(SD) -
                    2.0 * point.getEntry(SX) -
                    2.0 * point.getEntry(SK);

            constraints[4] =
                    point.getEntry(IMPFUS) -
                    20.0 * point.getEntry(SF);

            constraints[5] =
                    point.getEntry(MD) -
                    2.0 * point.getEntry(MV);

            constraints[6] =
                    point.getEntry(QF) -
                    point.getEntry(QI) -
                    point.getEntry(QV);

            constraints[7] =
                    point.getEntry(IMPTRAIN) -
                    0.137 * point.getEntry(MV);

            constraints[8] =
                    point.getEntry(IMPNMOT) -
                    35.0 * point.getEntry(NM);

            constraints[9] =
                    point.getEntry(IMPPET) -
                    0.043 * point.getEntry(QI);

            constraints[10] =
                    point.getEntry(IMPPIL) -
                    200.0 * point.getEntry(NP);

            constraints[11] =
                    point.getEntry(IMPCAN) -
                    120.0 * point.getEntry(NG);

            constraints[12] =
                    point.getEntry(IMPSNA) -
                    300.0 * point.getEntry(NS) -
                    400.0;

            constraints[13] =
                    point.getEntry(MC) -
                    point.getEntry(MV) +
                    95.0 * point.getEntry(NP) +
                    70.0 * point.getEntry(NG) +
                    660.0 * point.getEntry(NM) +
                    0.5 * point.getEntry(QI) -
                    380.0;

            // The SIF constant is -290, hence "- gconst" adds 290.
            constraints[14] =
                    point.getEntry(MZ) -
                    point.getEntry(IMPTRAIN) +
                    point.getEntry(IMPNMOT) +
                    point.getEntry(IMPPET) +
                    point.getEntry(IMPPIL) +
                    point.getEntry(IMPCAN) +
                    point.getEntry(IMPSNA) +
                    290.0;

            return new ArrayRealVector(
                    constraints,
                    false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Jacobian is evaluated by finite differences.");
        }
    }

    /** Square helper. */
    private static double square(final double value) {
        return value * value;
    }

    /** Construct bounds without exposing the static arrays. */
    private static SimpleBounds bounds() {
        return new SimpleBounds(
                LOWER.clone(),
                UPPER.clone());
    }

    @Test
    public void testAVION2() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(10000),
                        new InitialGuess(START.clone()),
                        new ObjectiveFunction(new Objective()),
                        new AvionEqualities(),
                        bounds(),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}