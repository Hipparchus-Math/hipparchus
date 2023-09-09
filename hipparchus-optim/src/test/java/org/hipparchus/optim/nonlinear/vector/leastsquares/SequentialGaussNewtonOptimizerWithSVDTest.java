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

import java.io.IOException;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.euclidean.threed.Plane;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.SingularValueDecomposer;
import org.hipparchus.optim.SimpleVectorValueChecker;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresOptimizer.Optimum;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem.Evaluation;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * <p>Some of the unit tests are re-implementations of the MINPACK <a
 * href="http://www.netlib.org/minpack/ex/file17">file17</a> and <a
 * href="http://www.netlib.org/minpack/ex/file22">file22</a> test files.
 * The redistribution policy for MINPACK is available <a
 * href="http://www.netlib.org/minpack/disclaimer">here</a>/
 *
 */
public class SequentialGaussNewtonOptimizerWithSVDTest
    extends AbstractSequentialLeastSquaresOptimizerAbstractTest {

    @Override
    public int getMaxIterations() {
        return 1000;
    }

    @Override
    public void defineOptimizer(Evaluation evaluation) {
        this.optimizer = new SequentialGaussNewtonOptimizer().
                         withDecomposer(new SingularValueDecomposer()).
                         withFormNormalEquations(false).
                         withEvaluation(evaluation);
    }

    @Test
    public void testMaxEvaluations() throws Exception {
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

            fail(optimizer);
        } catch (MathIllegalStateException e) {
            Assert.assertEquals(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, e.getSpecifier());
        }
    }


    @Override
    @Test
    public void testHahn1() throws IOException {
        try {
            /*
             * When NOT FORMING normal equations, the optimizer diverges and hit max evaluations.
             * When FORMING normal equations, the optimizer converges,
             * but the results are very bad
             */
            super.testHahn1();
            fail(optimizer);
        } catch (MathIllegalStateException e) {
            Assert.assertEquals(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, e.getSpecifier());
        }
    }

    @Test
    @Override
    public void testGetIterations() {
        /* this diverges with SVD and no normal equations */
        try {
            super.testGetIterations();
            fail(optimizer);
        } catch (MathIllegalStateException e) {
            Assert.assertEquals(LocalizedCoreFormats.MAX_COUNT_EXCEEDED,
                                e.getSpecifier());
        }
    }

    @Test
    @Override
    public void testNonInvertible() throws Exception {
        /*  SVD can compute a solution to singular problems.
         *  In this case the target vector, b, is not in the
         *  span of the jacobian matrix, A. The closest point
         *  to b on the plane spanned by A is computed.
         */
        LinearProblem problem = new LinearProblem(new double[][]{
                {1, 2, -3},
                {2, 1, 3},
                {-3, 0, -9}
        }, new double[]{1, 1, 1});

        defineOptimizer(null);
        Optimum optimum = optimizer.optimize(problem.getBuilder().build());

        Plane span = new Plane(Vector3D.ZERO, new Vector3D(1, 2, -3), new Vector3D(2, 1, 0), TOl);
        double expected = FastMath.abs(span.getOffset(new Vector3D(1, 1, 1)));
        double actual = optimum.getResiduals().getNorm();

        //verify
        Assert.assertEquals(expected, actual, TOl);
    }

}
