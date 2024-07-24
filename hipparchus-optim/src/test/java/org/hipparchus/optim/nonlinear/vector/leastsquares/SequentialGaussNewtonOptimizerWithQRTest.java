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

package org.hipparchus.optim.nonlinear.vector.leastsquares;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.linear.QRDecomposer;
import org.hipparchus.optim.LocalizedOptimFormats;
import org.hipparchus.optim.SimpleVectorValueChecker;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem.Evaluation;
import org.junit.jupiter.api.Test;

import java.io.IOException;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * <p>Some of the unit tests are re-implementations of the MINPACK <a
 * href="http://www.netlib.org/minpack/ex/file17">file17</a> and <a
 * href="http://www.netlib.org/minpack/ex/file22">file22</a> test files.
 * The redistribution policy for MINPACK is available <a
 * href="http://www.netlib.org/minpack/disclaimer">here</a>/
 *
 */
public class SequentialGaussNewtonOptimizerWithQRTest
    extends AbstractSequentialLeastSquaresOptimizerAbstractTest {

    @Override
    public int getMaxIterations() {
        return 1000;
    }

    @Override
    public void defineOptimizer(Evaluation evaluation) {
        this.optimizer = new SequentialGaussNewtonOptimizer().
                         withDecomposer(new QRDecomposer(1.0e-11)).
                         withFormNormalEquations(false).
                         withEvaluation(evaluation);
    }

    @Override
    @Test
    public void testMoreEstimatedParametersUnsorted() {
        /*
         * Exception is expected with this optimizer
         */
        try {
            super.testMoreEstimatedParametersUnsorted();
            customFail(optimizer);
        } catch (MathIllegalStateException mise) {
            assertEquals(LocalizedOptimFormats.UNABLE_TO_SOLVE_SINGULAR_PROBLEM, mise.getSpecifier());
        }
    }

    @Test
    void testMaxEvaluations() throws Exception {
        try {
            CircleVectorial circle = new CircleVectorial();
            circle.addPoint( 30.0,  68.0);
            circle.addPoint( 50.0,  -6.0);
            circle.addPoint(110.0, -20.0);
            circle.addPoint( 35.0,  15.0);
            circle.addPoint( 45.0,  97.0);

            LeastSquaresProblem lsp = builder(circle)
                            .checkerPair(new SimpleVectorValueChecker(1e-30, 1e-30))
                            .maxIterations(Integer.MAX_VALUE)
                            .start(new double[]{98.680, 47.345})
                            .build();

            defineOptimizer(null);
            optimizer.optimize(lsp);
            customFail(optimizer);
        } catch (MathIllegalStateException e) {
            assertEquals(LocalizedCoreFormats.MAX_COUNT_EXCEEDED,
                                e.getSpecifier());
        }
    }


    @Override
    @Test
    public void testHahn1() throws IOException {
        try {
            /*
             * TODO This test leads to a singular problem with the Gauss-Newton
             * optimizer. This should be inquired.
             */
            super.testHahn1();
        } catch (MathIllegalStateException mise) {
            assertEquals(LocalizedOptimFormats.UNABLE_TO_SOLVE_SINGULAR_PROBLEM, mise.getSpecifier());
        }
    }

}
