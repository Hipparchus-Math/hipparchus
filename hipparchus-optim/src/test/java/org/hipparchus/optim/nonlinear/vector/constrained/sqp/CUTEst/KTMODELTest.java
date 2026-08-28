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

import java.util.Arrays;

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
 * Corrected economic formulation of CUTEst problem {@code KTMODEL}.
 *
 * <p>The source represents a four-sector dynamic macro-economic model over
 * thirty periods. The SIF declares 726 scalar variables. Fifteen variables
 * are fixed and are eliminated here, yielding the 711 free variables reported
 * by decoded CUTEst tables. The model has 450 nonlinear equalities.</p>
 *
 * <p>The historic SIF cards are internally inconsistent: the classification
 * says 720 variables although expansion produces 726, the file is classified
 * as having no objective, and several dependent-variable terms and matrix
 * entries are corrupted. This class restores the published Kendrick-Taylor
 * optimization model, using the 30-period SIF data and its corrected numerical
 * starting point. The CUTEst test is treated strictly as a feasibility
 * problem with a constant zero objective:</p>
 * <ul>
 *   <li>the capital transition is {@code K(t+1)=K(t)+mu*KAG};</li>
 *   <li>the output equation contains {@code Q(t,i)} and the complete CES
 *       outer exponent;</li>
 *   <li>the resource equation contains {@code C(t,i)};</li>
 *   <li>the initial labor supply is 8.6, consistent with the four fixed
 *       first-period labor values.</li>
 * </ul>
 *
 * <p>The supplied starting point is the corrected KTMODEL/KTMODEL2 numerical
 * point, not an interpolated or synthetic point. Derivatives are evaluated
 * through forward finite differences and the SQP BFGS approximation.</p>
 */
public class KTMODELTest {

    /** Number of planning periods. */
    private static final int PERIODS = 30;

    /** Number of sectors. */
    private static final int SECTORS = 4;

    /** Number of variables declared by the model before fixed-variable elimination. */
    private static final int FULL_N = 726;

    /** Number of fixed variables eliminated from the optimization vector. */
    private static final int FIXED_VARIABLES = 15;

    /** Number of free variables exposed to SQPOptimizerS2. */
    private static final int N = FULL_N - FIXED_VARIABLES;

    /** Number of nonlinear equality constraints. */
    private static final int M = 450;

    /** KTMODEL is a feasibility problem: the objective is identically zero. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    private static final double DISCOUNT_RATE = 0.03;
    private static final double DEBT_RATE = 0.05;
    private static final double EPSILON = 0.5;
    private static final double CONSUMPTION_LOWER = 1.0e-6;
    private static final double LABOR_INITIAL = 8.6;
    private static final double LABOR_GROWTH = 0.02;

    private static final double[] A = { 0.48, 0.33, 0.345, 0.3925 };
    private static final double[] B = { 0.85, 0.90, 0.91, 0.87 };
    private static final double[] MU = { 0.275, 0.35, 0.30, 0.35 };
    private static final double[] TAU = { 0.41, 1.26, 1.89, 0.47 };
    private static final double[] NU = { 0.03, 0.035, 0.025, 0.025 };
    private static final double[] BETA = { 0.35, 0.30, 0.25, 0.20 };
    private static final double[] RHO = { -0.166, 0.111, 0.111, 0.666 };
    private static final double[] D = { 0.0008, 0.090, 0.030, 0.004 };
    private static final double[] PI = { 0.63, 0.98, 0.10, 0.10 };

    private static final double[] K_INITIAL = { 2.02, 2.13, 1.26, 1.27 };
    private static final double[] K_TERMINAL = { 14.2, 20.0, 10.2, 10.3 };
    private static final double[] L_INITIAL = { 5.1, 0.84, 0.36, 2.3 };

    /** Exogenous exports by period and sector. */
    private static final double[][] EXPORTS = {
        { 0.064, 0.035, 0.106, 0.131 },
        { 0.070, 0.043, 0.128, 0.148 },
        { 0.076, 0.051, 0.154, 0.166 },
        { 0.082, 0.062, 0.185, 0.185 },
        { 0.088, 0.073, 0.220, 0.205 },
        { 0.094, 0.087, 0.261, 0.277 },
        { 0.099, 0.103, 0.308, 0.251 },
        { 0.103, 0.121, 0.362, 0.276 },
        { 0.107, 0.141, 0.424, 0.302 },
        { 0.110, 0.165, 0.495, 0.330 },
        { 0.112, 0.189, 0.530, 0.348 },
        { 0.113, 0.214, 0.567, 0.366 },
        { 0.115, 0.243, 0.607, 0.334 },
        { 0.115, 0.274, 0.649, 0.404 },
        { 0.116, 0.309, 0.695, 0.425 },
        { 0.116, 0.347, 0.744, 0.446 },
        { 0.115, 0.389, 0.796, 0.469 },
        { 0.114, 0.435, 0.851, 0.492 },
        { 0.110, 0.486, 0.911, 0.516 },
        { 0.108, 0.541, 0.975, 0.541 },
        { 0.114, 0.606, 1.006, 0.560 },
        { 0.121, 0.676, 1.038, 0.579 },
        { 0.127, 0.751, 1.070, 0.599 },
        { 0.134, 0.833, 1.102, 0.618 },
        { 0.142, 0.921, 1.134, 0.638 },
        { 0.150, 1.017, 1.198, 0.658 },
        { 0.158, 1.119, 1.198, 0.678 },
        { 0.166, 1.230, 1.230, 0.698 },
        { 0.175, 1.349, 1.261, 0.718 },
        { 0.185, 1.477, 1.292, 0.739 }
    };

    /** Leontief input-output matrix. */
    private static final double[][] AA = {
        { 0.10, 0.09, 0.17, 0.01 },
        { 0.09, 0.33, 0.24, 0.12 },
        { 0.04, 0.02, 0.12, 0.05 },
        { 0.03, 0.09, 0.09, 0.08 }
    };

    /** Capital-goods coefficient matrix. */
    private static final double[][] BB = {
        { 0.0,    0.0,    0.0,    0.0 },
        { 0.6908, 1.3109, 0.1769, 0.15 },
        { 0.001,  0.0199, 0.0022, 0.0 },
        { 0.0,    0.0,    0.0,    0.0 }
    };

    /**
     * Corrected numerical starting point.
     *
     * <p>Row layout:
     * ksi, K(1..4), gam, Q(1..4), C(1..4), del(1..4),
     * L(1..4), M2, M3.</p>
     */
    private static final double[][] PERIOD_START = {
        { 0.0, 2.02, 2.13, 1.26, 1.27, 0.25,
          1.5824706, 1.4338352, 0.93881068, 0.96418087, 1.063367, 0.30658512,
          0.60609873, 0.49885553, 0.1546822, 0.085635062, 0.1696786, 0.10188439,
          5.1, 0.84, 0.36, 2.3, 0.0, 0.0 },
        { 1.0220136, 2.1474812, 2.2087845, 1.3858756, 1.3567091, 0.29737259,
          1.7428778, 1.7779406, 0.98301901, 0.86430611, 1.1642574, 0.53705577,
          0.61546349, 0.35016122, 0.15273789, 0.091335202, 0.18578568, 0.11213136,
          5.4679576, 1.0647922, 0.36, 1.8792502, 0.0, 0.0 },
        { 2.0722324, 2.2749857, 2.29262, 1.5238659, 1.4517037, 0.33347851,
          1.9345542, 1.8338029, 1.0290609, 0.79543358, 1.3186216, 0.54921599,
          0.62581567, 0.25368323, 0.15004432, 0.088947588, 0.2330375, 0.099934103,
          5.9390836, 1.0442219, 0.36, 1.6041345, 0.0, 0.0 },
        { 3.1190326, 2.4017909, 2.3746838, 1.6907124, 1.5386167, 0.31821387,
          2.1003007, 2.335946, 1.079197, 0.62252534, 1.4098752, 0.5607578,
          0.63287538, 0.02, 0.17033776, 0.086698026, 0.29976728, 2.3057285,
          6.2555346, 1.3764943, 0.36, 1.1343599, 0.0, 0.0 },
        { 4.0780633, 2.5440561, 2.4550418, 1.8957129, 2.0225426, 0.52028373,
          1.8657186, 1.3815054, 1.0645461, 1.6234601, 1.2710468, 0.14272656,
          0.56412797, 1.0185121, 0.14647303, 0.02, 0.25442974, 0.2880789,
          4.897166, 0.62583741, 0.33205642, 3.4538568, 0.0, 0.0 },
        { 5.1558883, 2.6703667, 2.474698, 2.0846294, 2.2416789, 0.28968834,
          0.5432637, 1.7376584, 1.1146323, 2.9831498, 0.02, 0.02,
          0.5419143, 2.2559574, 0.17248857, 0.26267773, 0.21988009, 0.24578245,
          0.65197036, 0.81558091, 0.33243443, 7.6951092, 0.0, 0.0 },
        { 6.0111074, 2.8166152, 2.6877725, 2.2576552, 2.439789, 0.25,
          0.6385515, 2.5381646, 1.1616932, 3.0429391, 0.02, 0.50404675,
          0.51383192, 2.2085921, 0.17572694, 0.31078498, 0.26223977, 0.30454648,
          0.78449273, 1.2736726, 0.33229782, 7.2945337, 0.0, 0.0 },
        { 6.9588775, 2.9663999, 2.935592, 2.4595439, 2.6786901, 0.25,
          0.68561105, 2.8567701, 1.2039529, 3.2357489, 0.02, 0.59700994,
          0.47905577, 2.3279062, 0.19132812, 0.36898536, 0.27560213, 0.33241659,
          0.81219777, 1.3827015, 0.32933338, 7.3544641, 0.0, 0.0 },
        { 7.921917, 3.1286589, 3.2245978, 2.6735024, 2.9397901, 0.25,
          0.72504136, 3.1015554, 1.2472769, 3.4655554, 0.02, 0.60243108,
          0.43666074, 2.4866795, 0.20261128, 0.45058177, 0.30979842, 0.31619128,
          0.81899618, 1.424393, 0.32649971, 7.5063817, 0.0, 0.0 },
        { 8.8843305, 3.3003793, 3.5689516, 2.9121399, 3.1956171, 0.25,
          1.8260276, 3.356243, 1.2886361, 2.7912943, 0.98573746, 0.60582006,
          0.38561037, 1.7760203, 0.20527273, 0.55295396, 0.31760931, 0.32227669,
          3.3213387, 1.4565503, 0.32246508, 5.177442, 0.0, 0.0 },
        { 9.9828005, 3.4754291, 3.9806729, 3.1602272, 3.4595735, 0.25,
          1.4258132, 3.5710812, 1.3384312, 3.4982758, 0.58853717, 0.61062661,
          0.37149722, 2.4002291, 0.19720725, 0.60751629, 0.33190394, 0.36144987,
          2.0936888, 1.4510375, 0.3209893, 6.6176364, 0.0, 0.0 },
        { 11.065293, 3.6458342, 4.4348056, 3.4216498, 3.7536589, 0.25,
          2.1011243, 3.7519575, 1.3765265, 3.1454406, 1.1760991, 0.61403409,
          0.35537079, 2.0162173, 0.25459103, 0.63604791, 0.36632176, 0.34508681,
          3.5668903, 1.4217616, 0.31582976, 5.3885374, 0.0, 0.0 },
        { 12.202019, 3.8590172, 4.9178119, 3.7089706, 4.0409806, 0.25,
          2.0401711, 3.9875881, 1.4275351, 3.5690527, 1.0858806, 0.61585105,
          0.33706697, 2.3648201, 0.22626623, 0.71370006, 0.38672726, 0.38072648,
          3.1555493, 1.4184103, 0.31398397, 6.0189359, 0.0, 0.0 },
        { 13.350325, 4.0536371, 5.4581992, 4.0140528, 4.3566387, 0.25,
          1.672634, 4.2334771, 1.4815378, 4.288835, 0.71548558, 0.61586746,
          0.31640979, 2.9942417, 0.24252923, 0.7686972, 0.38427772, 0.42258141,
          2.1400039, 1.4130112, 0.31261687, 7.259385, 0.0, 0.0 },
        { 14.483882, 4.261642, 6.0445521, 4.3226533, 4.7052054, 0.25,
          2.4096844, 4.4716217, 1.5262078, 3.9326306, 1.3535949, 0.61385311,
          0.29321056, 2.5969776, 0.20510681, 0.86441012, 0.42648751, 0.46904628,
          3.5289963, 1.3994701, 0.30881275, 6.1102381, 0.0, 0.0 },
        { 15.669356, 4.4426568, 6.7014438, 4.6630523, 5.0901911, 0.25,
          2.6432141, 4.7268527, 1.5769644, 4.0858717, 1.530981, 0.60955576,
          0.26726697, 2.682562, 0.23082428, 0.92846898, 0.46589851, 0.46790452,
          3.8029607, 1.38748, 0.30612119, 6.0779058, 0.0, 0.0 },
        { 16.872523, 4.6444564, 7.4124788, 5.033937, 5.4797794, 0.25,
          2.6009726, 4.7686471, 1.478876, 4.6400999, 1.5010447, 0.60269928,
          0.23836218, 3.1786128, 0.29163801, 0.86251017, 0.48557089, 0.45764113,
          3.4547435, 1.2928171, 0.26901591, 6.7893805, 0.0, 0.13805102 },
        { 18.115324, 4.8927994, 8.0993444, 5.4232879, 5.8668506, 0.25,
          2.3862378, 4.6855447, 1.2845076, 5.4115919, 1.3418329, 0.59298155,
          0.20626372, 3.8995476, 0.30946569, 0.75583528, 0.49236328, 0.4585241,
          2.7926691, 1.1685015, 0.21535689, 7.8655487, 0.0, 0.3647531 },
        { 19.399892, 5.1560346, 8.7270425, 5.8228399, 6.2586086, 0.25,
          2.9201456, 4.8427455, 1.2462945, 5.3632252, 1.8177844, 0.58007218,
          0.17072245, 3.8038208, 0.3260612, 0.78143007, 0.50312127, 0.46975977,
          3.5577888, 1.1367583, 0.19800682, 7.3903638, 0.0, 0.44620737 },
        { 20.725169, 5.4333931, 9.3805658, 6.2349478, 6.6623325, 0.25,
          2.8125594, 4.5260909, 0.8900854, 6.1192754, 1.8054011, 0.56361011,
          0.13147137, 4.5408934, 0.42966282, 0.52545641, 0.48373848, 0.41643659,
          3.1034954, 0.96789336, 0.1253872, 8.3318001, 0.0, 0.81700713 },
        { 22.114094, 5.7854346, 9.8490656, 6.6387737, 7.0289624, 0.25,
          2.4714909, 4.673846, 0.5513209, 6.9963336, 1.5276523, 0.53980187,
          0.12499618, 5.3599726, 0.43386511, 0.58791767, 0.36452018, 0.51432036,
          2.313599, 0.94948227, 0.067347986, 9.4487184, 0.0, 1.1844369 },
        { 23.522926, 6.144334, 10.369505, 6.9587327, 7.4725745, 0.5648112,
          2.1194508, 4.9766009, 1.2206889, 7.7431245, 1.0556761, 0.512132,
          0.1187052, 5.9540304, 0.55403317, 0.49518402, 0.51447856, 0.38997648,
          1.6451019, 0.96937115, 0.16936869, 10.250889, 0.0, 0.62865308 },
        { 24.872617, 6.5868073, 10.818263, 7.3916469, 7.822904, 0.25,
          2.4624262, 6.4993134, 2.0020839, 7.7327219, 1.0881722, 0.48028903,
          0.11267162, 5.707437, 0.64095445, 1.1420748, 0.44979329, 0.38816623,
          1.9195659, 1.3034088, 0.29741902, 9.7750314, 0.0, 8.8633528e-18 },
        { 26.149412, 7.0904737, 11.745643, 7.7812212, 8.17338, 5.4210109e-18,
          4.7779224, 6.5783734, 1.4311656, 6.1626888, 3.2725896, 0.44394334,
          0.10697615, 4.211961, 0.71625934, 1.1903553, 0.511794, 0.42438596,
          5.0964954, 1.2306186, 0.18891916, 7.0453004, 0.0, 0.5628125 },
        { 27.571234, 7.6486448, 12.719704, 8.2196711, 8.5548198, 0.25,
          2.6044332, 6.863966, 1.7267834, 8.8695591, 1.2043019, 0.40274665,
          0.10170756, 6.7062258, 1.2125356, 0.91149536, 0.50374164, 0.38402991,
          1.7265353, 1.2115009, 0.22779777, 10.666726, 0.0, 0.36975497 },
        { 28.911459, 8.4845841, 13.508161, 8.6555493, 8.904832, 0.25,
          5.5459798, 7.107383, 1.730251, 6.850539, 3.8439836, 0.35633129,
          0.096963237, 4.7102326, 1.2346703, 1.0008649, 0.33096457, 0.45714631,
          5.2960824, 1.1884944, 0.21902994, 7.4056048, 0.0, 0.41703779 },
        { 30.390985, 9.3587943, 14.369912, 8.9573653, 9.316192, 0.25,
          3.8418761, 7.7786312, 2.1036351, 9.2024799, 2.153372, 0.30430981,
          0.092849648, 6.8204485, 2.7082306, 0.39299959, 0.5170604, 0.37036217,
          2.6177667, 1.2542909, 0.26785234, 10.251486, 0.0, 0.15777132 },
        { 31.813488, 10.827462, 14.741024, 9.408348, 9.6572339, 0.25,
          5.7084903, 7.9391579, 2.0987232, 8.4452598, 3.8202219, 0.24627437,
          0.089482908, 6.0306055, 3.4697818, 0.030985729, 0.16604599, 0.39735859,
          4.3544834, 1.2210708, 0.25628775, 8.8473818, 0.0, 0.22352392 },
        { 33.345674, 12.616277, 14.771871, 9.5673449, 10.022115, 0.25,
          7.4665302, 8.4015011, 1.9112471, 7.7480777, 5.3921149, 0.18179635,
          0.086989383, 5.2887296, 0.536157, 1.808448, 0.24220802, 0.065812597,
          5.8321176, 1.2566843, 0.22099455, 7.6630117, 0.0, 0.50035247 },
        { 34.95412, 13.096128, 16.196551, 9.7950426, 10.087013, 0.25,
          4.7930677, 20.901936, 2.8405831, 8.5438954, 1.6834514, 0.11042591,
          0.085506224, 4.8754142, 1.4465126, 8.4269429, 0.45228158, 0.22313589,
          2.6001895, 3.9627465, 0.34650299, 8.3628253, 0.0, 0.0 }
    };

    /** Fixed full-variable values; NaN marks a free variable. */
    private static final double[] FIXED = createFixedValues();

    /** Map from the 726 declared variables to the 711 free variables. */
    private static final int[] FULL_TO_FREE = createFullToFreeMap();

    /** Constant null objective used for the CUTEst feasibility problem. */
    private static final class Objective extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            return 0.0;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Gradient is evaluated by forward finite differences.");
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Hessian is supplied by the SQP BFGS approximation.");
        }
    }

    /** The 450 dynamic model equations. */
    private static final class ModelEqualities extends EqualityConstraint {

        ModelEqualities() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {

            final double[] residuals = new double[M];
            int row = 0;

            // Welfare accumulation: 30 equations.
            double discount = 1.0 + DISCOUNT_RATE;
            for (int t = 0; t < PERIODS; ++t) {
                double utility = 0.0;
                for (int sector = 0; sector < SECTORS; ++sector) {
                    utility += A[sector] *
                               FastMath.pow(KTMODELTest.value(point, cIndex(t, sector)),
                                            B[sector]);
                }
                residuals[row++] =
                        KTMODELTest.value(point, ksiIndex(t + 1)) -
                        KTMODELTest.value(point, ksiIndex(t)) -
                        utility / discount;
                discount *= 1.0 + DISCOUNT_RATE;
            }

            // Capital transition: 120 equations.
            for (int t = 0; t < PERIODS; ++t) {
                for (int sector = 0; sector < SECTORS; ++sector) {
                    final double capital =
                            KTMODELTest.value(point, kIndex(t, sector));
                    final double investment =
                            KTMODELTest.value(point, delIndex(t, sector));
                    final double effectiveInvestment =
                            absorptiveCapacity(capital, investment, sector);

                    residuals[row++] =
                            KTMODELTest.value(point, kIndex(t + 1, sector)) -
                            capital -
                            MU[sector] * effectiveInvestment;
                }
            }

            // Foreign-debt transition: 30 equations.
            for (int t = 0; t < PERIODS; ++t) {
                double imports = KTMODELTest.value(point, m2Index(t)) +
                                 KTMODELTest.value(point, m3Index(t));
                double exports = 0.0;

                for (int sector = 0; sector < SECTORS; ++sector) {
                    imports += D[sector] *
                               KTMODELTest.value(point, qIndex(t, sector));
                    imports += PI[sector] *
                               KTMODELTest.value(point, delIndex(t, sector));
                    exports += EXPORTS[t][sector];
                }

                residuals[row++] =
                        KTMODELTest.value(point, gamIndex(t + 1)) -
                        (1.0 + DEBT_RATE) *
                        KTMODELTest.value(point, gamIndex(t)) -
                        imports +
                        exports;
            }

            // CES production: 120 equations.
            for (int t = 0; t < PERIODS; ++t) {
                for (int sector = 0; sector < SECTORS; ++sector) {
                    final double capital =
                            KTMODELTest.value(point, kIndex(t, sector));
                    final double labor =
                            KTMODELTest.value(point, lIndex(t, sector));
                    final double inner =
                            BETA[sector] *
                            FastMath.pow(capital, -RHO[sector]) +
                            (1.0 - BETA[sector]) *
                            FastMath.pow(labor, -RHO[sector]);
                    final double ces =
                            FastMath.pow(inner, -1.0 / RHO[sector]);
                    final double productivity =
                            TAU[sector] *
                            FastMath.pow(1.0 + NU[sector], t + 1);

                    residuals[row++] =
                            KTMODELTest.value(point, qIndex(t, sector)) -
                            productivity * ces;
                }
            }

            // Sector resource balances: 120 equations.
            for (int t = 0; t < PERIODS; ++t) {
                for (int sector = 0; sector < SECTORS; ++sector) {
                    double productionAvailable = 0.0;
                    double capitalDemand = 0.0;

                    for (int j = 0; j < SECTORS; ++j) {
                        final double identity =
                                sector == j ? 1.0 : 0.0;
                        final double importCoefficient =
                                sector == j ? D[sector] : 0.0;
                        final double p =
                                identity -
                                AA[sector][j] +
                                importCoefficient;

                        productionAvailable +=
                                p *
                                KTMODELTest.value(point, qIndex(t, j));
                        capitalDemand +=
                                BB[sector][j] *
                                KTMODELTest.value(point, delIndex(t, j));
                    }

                    final double directImports;
                    if (sector == 1) {
                        directImports = KTMODELTest.value(point, m2Index(t));
                    } else if (sector == 2) {
                        directImports = KTMODELTest.value(point, m3Index(t));
                    } else {
                        directImports = 0.0;
                    }

                    residuals[row++] =
                            KTMODELTest.value(point, cIndex(t, sector)) -
                            productionAvailable +
                            capitalDemand +
                            EXPORTS[t][sector] -
                            directImports;
                }
            }

            // Aggregate labor supply: 30 equations.
            double laborSupply = LABOR_INITIAL;
            for (int t = 0; t < PERIODS; ++t) {
                double usedLabor = 0.0;
                for (int sector = 0; sector < SECTORS; ++sector) {
                    usedLabor += KTMODELTest.value(point, lIndex(t, sector));
                }
                residuals[row++] = usedLabor - laborSupply;
                laborSupply *= 1.0 + LABOR_GROWTH;
            }

            if (row != M) {
                throw new IllegalStateException(
                        "KTMODEL residual count is " + row + ", expected " + M);
            }

            return new ArrayRealVector(residuals, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Jacobian is evaluated by forward finite differences.");
        }
    }

    /** Absorptive-capacity investment function verified from the SIF start. */
    private static double absorptiveCapacity(final double capital,
                                             final double investment,
                                             final int sector) {
        final double ratio =
                1.0 +
                (EPSILON / MU[sector]) *
                investment / capital;
        return capital *
               (1.0 - FastMath.pow(ratio, -1.0 / EPSILON));
    }

    /** Construct the 711-entry corrected starting point. */
    private static double[] initialPoint() {
        final double[] full = new double[FULL_N];

        for (int t = 0; t < PERIODS; ++t) {
            final double[] row = PERIOD_START[t];

            full[ksiIndex(t)] = row[0];
            for (int sector = 0; sector < SECTORS; ++sector) {
                full[kIndex(t, sector)] = row[1 + sector];
            }
            full[gamIndex(t)] = row[5];

            for (int sector = 0; sector < SECTORS; ++sector) {
                full[qIndex(t, sector)] = row[6 + sector];
                full[cIndex(t, sector)] = row[10 + sector];
                full[delIndex(t, sector)] = row[14 + sector];
                full[lIndex(t, sector)] = row[18 + sector];
            }

            full[m2Index(t)] = row[22];
            full[m3Index(t)] = row[23];
        }

        final double[] reduced = new double[N];
        for (int fullIndex = 0; fullIndex < FULL_N; ++fullIndex) {
            final int freeIndex = FULL_TO_FREE[fullIndex];
            if (freeIndex >= 0) {
                reduced[freeIndex] = full[fullIndex];
            }
        }
        return reduced;
    }

    /** CUTEst bounds after eliminating the fifteen fixed variables. */
    private static SimpleBounds bounds() {
        final double[] lower = new double[N];
        final double[] upper = new double[N];

        Arrays.fill(lower, 0.0);
        Arrays.fill(upper, Double.POSITIVE_INFINITY);

        // gam(2),...,gam(29) are explicitly free in the SIF.
        for (int t = 1; t <= 28; ++t) {
            lower[freeIndex(gamIndex(t))] = Double.NEGATIVE_INFINITY;
        }

        // Every consumption variable has the positive COLO lower bound.
        for (int t = 0; t < PERIODS; ++t) {
            for (int sector = 0; sector < SECTORS; ++sector) {
                lower[freeIndex(cIndex(t, sector))] =
                        CONSUMPTION_LOWER;
            }
        }

        return new SimpleBounds(lower, upper);
    }

    /** Build fixed values in the original 726-variable declaration space. */
    private static double[] createFixedValues() {
        final double[] fixed = new double[FULL_N];
        Arrays.fill(fixed, Double.NaN);

        fixed[ksiIndex(0)] = 0.0;
        fixed[gamIndex(0)] = 0.0;
        fixed[gamIndex(PERIODS - 1)] = 8.0;

        for (int sector = 0; sector < SECTORS; ++sector) {
            fixed[kIndex(0, sector)] = K_INITIAL[sector];
            fixed[kIndex(PERIODS - 1, sector)] = K_TERMINAL[sector];
            fixed[lIndex(0, sector)] = L_INITIAL[sector];
        }

        return fixed;
    }

    /** Create the reduction map and verify the dimension exactly. */
    private static int[] createFullToFreeMap() {
        final int[] map = new int[FULL_N];
        Arrays.fill(map, -1);

        int free = 0;
        for (int full = 0; full < FULL_N; ++full) {
            if (Double.isNaN(FIXED[full])) {
                map[full] = free++;
            }
        }

        if (free != N) {
            throw new IllegalStateException(
                    "KTMODEL free-variable count is " + free +
                    ", expected " + N);
        }

        return map;
    }

    /** Read a declared variable, substituting fixed values where required. */
    private static double value(final RealVector point,
                                final int fullIndex) {
        final double fixed = FIXED[fullIndex];
        return Double.isNaN(fixed) ?
               point.getEntry(FULL_TO_FREE[fullIndex]) :
               fixed;
    }

    /** Return the reduced index of a variable known to be free. */
    private static int freeIndex(final int fullIndex) {
        final int index = FULL_TO_FREE[fullIndex];
        if (index < 0) {
            throw new IllegalArgumentException(
                    "Variable " + fullIndex + " is fixed.");
        }
        return index;
    }

    // Original declaration-space indices.
    private static int ksiIndex(final int t) { return t; }
    private static int gamIndex(final int t) { return 31 + t; }
    private static int kIndex(final int t, final int sector) {
        return 62 + 4 * t + sector;
    }
    private static int qIndex(final int t, final int sector) {
        return 186 + 4 * t + sector;
    }
    private static int cIndex(final int t, final int sector) {
        return 306 + 4 * t + sector;
    }
    private static int delIndex(final int t, final int sector) {
        return 426 + 4 * t + sector;
    }
    private static int lIndex(final int t, final int sector) {
        return 546 + 4 * t + sector;
    }
    private static int m2Index(final int t) { return 666 + t; }
    private static int m3Index(final int t) { return 696 + t; }

    @Test
    public void testKTMODEL() {
        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();
        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(10000),
                        new InitialGuess(initialPoint()),
                        new ObjectiveFunction(new Objective()),
                        new ModelEqualities(),
                        bounds(),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}