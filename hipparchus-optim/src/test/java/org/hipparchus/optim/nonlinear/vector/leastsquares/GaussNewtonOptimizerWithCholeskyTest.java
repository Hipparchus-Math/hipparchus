/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.optim.nonlinear.vector.leastsquares;

import java.io.IOException;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.optim.SimpleVectorValueChecker;
import org.hipparchus.optim.nonlinear.vector.leastsquares.GaussNewtonOptimizer;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresOptimizer;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem;
import org.hipparchus.optim.nonlinear.vector.leastsquares.GaussNewtonOptimizer.Decomposition;
import org.junit.Test;

/**
 * <p>Some of the unit tests are re-implementations of the MINPACK <a
 * href="http://www.netlib.org/minpack/ex/file17">file17</a> and <a
 * href="http://www.netlib.org/minpack/ex/file22">file22</a> test files.
 * The redistribution policy for MINPACK is available <a
 * href="http://www.netlib.org/minpack/disclaimer">here</a>/
 *
 */
public class GaussNewtonOptimizerWithCholeskyTest
    extends AbstractLeastSquaresOptimizerAbstractTest {

    @Override
    public int getMaxIterations() {
        return 1000;
    }

    @Override
    public LeastSquaresOptimizer getOptimizer() {
        return new GaussNewtonOptimizer(Decomposition.CHOLESKY);
    }

    @Override
    @Test(expected=MathIllegalStateException.class)
    public void testMoreEstimatedParametersSimple() {
        /*
         * Exception is expected with this optimizer
         */
        super.testMoreEstimatedParametersSimple();
    }

    @Override
    @Test(expected=MathIllegalStateException.class)
    public void testMoreEstimatedParametersUnsorted() {
        /*
         * Exception is expected with this optimizer
         */
        super.testMoreEstimatedParametersUnsorted();
    }

    @Test(expected=MathIllegalStateException.class)
    public void testMaxEvaluations() throws Exception {
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

        optimizer.optimize(lsp);
    }

    @Override
    @Test(expected=MathIllegalStateException.class)
    public void testCircleFittingBadInit() {
        /*
         * This test does not converge with this optimizer.
         */
        super.testCircleFittingBadInit();
        fail(optimizer);
    }

    @Override
    @Test(expected=MathIllegalStateException.class)
    public void testHahn1()
        throws IOException {
        /*
         * TODO This test leads to a singular problem with the Gauss-Newton
         * optimizer. This should be inquired.
         */
        super.testHahn1();
        fail(optimizer);
    }

}
