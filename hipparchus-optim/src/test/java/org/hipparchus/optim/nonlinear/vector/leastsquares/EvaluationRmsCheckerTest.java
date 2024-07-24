/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.optim.nonlinear.vector.leastsquares;

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem.Evaluation;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/** Unit tests for {@link EvaluationRmsChecker}. */
class EvaluationRmsCheckerTest {

    /** check {@link ConvergenceChecker#converged(int, Object, Object)}. */
    @Test
    void testConverged() {
        //setup
        ConvergenceChecker<Evaluation> checker = new EvaluationRmsChecker(0.1, 1);
        Evaluation e200 = mockEvaluation(200);
        Evaluation e1 = mockEvaluation(1);

        //action + verify
        //just matches rel tol
        assertTrue(checker.converged(0, e200, mockEvaluation(210)));
        //just matches abs tol
        assertTrue(checker.converged(0, e1, mockEvaluation(1.9)));
        //matches both
        assertTrue(checker.converged(0, e1, mockEvaluation(1.01)));
        //matches neither
        assertFalse(checker.converged(0, e200, mockEvaluation(300)));
    }

    /**
     * Create a mock {@link Evaluation}.
     *
     * @param rms the evaluation's rms.
     * @return a new mock evaluation.
     */
    private static Evaluation mockEvaluation(final double rms) {
        return new Evaluation() {
            public RealMatrix getCovariances(double threshold) {
                return null;
            }

            public RealVector getSigma(double covarianceSingularityThreshold) {
                return null;
            }

            public double getRMS() {
                return rms;
            }

            public RealMatrix getJacobian() {
                return null;
            }

            public double getCost() {
                return 0;
            }

            public double getChiSquare() {
                return 0;
            }

            public double getReducedChiSquare(int n) {
                return 0;
            }

            public RealVector getResiduals() {
                return null;
            }

            public RealVector getPoint() {
                return null;
            }
        };
    }

}
