/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements. See the NOTICE file distributed with this
 * work for additional information regarding copyright ownership. The ASF
 * licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * https://www.apache.org/licenses/LICENSE-2.0 Unless required by applicable law
 * or agreed to in writing, software distributed under the License is
 * distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for the specific language
 * governing permissions and limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.leastsquares;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.MultivariateMatrixFunction;
import org.hipparchus.analysis.MultivariateVectorFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.DiagonalMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem.Evaluation;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Pair;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * The only features tested here are utility methods defined
 * in {@link LeastSquaresProblem.Evaluation} that compute the
 * chi-square and parameters standard-deviations.
 */
public class EvaluationTest {

    /**
     * Create a {@link LeastSquaresBuilder} from a {@link StatisticalReferenceDataset}.
     *
     * @param dataset the source data
     * @return a builder for further customization.
     */
    public LeastSquaresBuilder builder(StatisticalReferenceDataset dataset) {
        StatisticalReferenceDataset.LeastSquaresProblem problem
                = dataset.getLeastSquaresProblem();
        final double[] start = dataset.getParameters();
        final double[] observed = dataset.getData()[1];
        final double[] weights = new double[observed.length];
        Arrays.fill(weights, 1d);

        return new LeastSquaresBuilder()
                .model(problem.getModelFunction(), problem.getModelFunctionJacobian())
                .target(observed)
                .weight(new DiagonalMatrix(weights))
                .start(start);
    }

    @Test
    void testComputeResiduals() {
        //setup
        RealVector point = new ArrayRealVector(2);
        Evaluation evaluation = new LeastSquaresBuilder()
                .target(new ArrayRealVector(new double[]{3,-1}))
                .model(point1 -> new Pair<RealVector, RealMatrix>(
                        new ArrayRealVector(new double[]{1, 2}),
                        MatrixUtils.createRealIdentityMatrix(2)
                ))
                .weight(MatrixUtils.createRealIdentityMatrix(2))
                .build()
                .evaluate(point);

        //action + verify
        assertArrayEquals(
            new double[]{2, -3},
            evaluation.getResiduals().toArray(),
            Precision.EPSILON);
    }

    @Test
    void testComputeCovariance() throws IOException {
        //setup
        RealVector point = new ArrayRealVector(2);
        Evaluation evaluation = new LeastSquaresBuilder()
                .model(point1 -> new Pair<RealVector, RealMatrix>(
                        new ArrayRealVector(2),
                        MatrixUtils.createRealDiagonalMatrix(new double[]{1, 1e-2})
                ))
                .weight(MatrixUtils.createRealDiagonalMatrix(new double[]{1, 1}))
                .target(new ArrayRealVector(2))
                .build()
                .evaluate(point);

        //action
        UnitTestUtils.customAssertEquals(
                "covariance",
                evaluation.getCovariances(FastMath.nextAfter(1e-4, 0.0)),
                MatrixUtils.createRealMatrix(new double[][]{{1, 0}, {0, 1e4}}),
                Precision.EPSILON
        );

        //singularity fail
        try {
            evaluation.getCovariances(FastMath.nextAfter(1e-4, 1.0));
            fail("Expected Exception");
        } catch (MathIllegalArgumentException e) {
            //expected
        }
    }

    @Test
    void testComputeValueAndJacobian() {
        //setup
        final RealVector point = new ArrayRealVector(new double[]{1, 2});
        Evaluation evaluation = new LeastSquaresBuilder()
                .weight(new DiagonalMatrix(new double[]{16, 4}))
                .model(actualPoint -> {
                    //verify correct values passed in
                    assertArrayEquals(
                            point.toArray(), actualPoint.toArray(), Precision.EPSILON);
                    //return values
                    return new Pair<RealVector, RealMatrix>(
                            new ArrayRealVector(new double[]{3, 4}),
                            MatrixUtils.createRealMatrix(new double[][]{{5, 6}, {7, 8}})
                    );
                })
                .target(new double[2])
                .build()
                .evaluate(point);

        //action
        RealVector residuals = evaluation.getResiduals();
        RealMatrix jacobian = evaluation.getJacobian();

        //verify
        assertArrayEquals(evaluation.getPoint().toArray(), point.toArray(), 0);
        assertArrayEquals(new double[]{-12, -8}, residuals.toArray(), Precision.EPSILON);
        UnitTestUtils.customAssertEquals(
                "jacobian",
                jacobian,
                MatrixUtils.createRealMatrix(new double[][]{{20, 24},{14, 16}}),
                Precision.EPSILON);
    }

    @Test
    void testComputeCost() throws IOException {
        final StatisticalReferenceDataset dataset
            = StatisticalReferenceDatasetFactory.createKirby2();

        final LeastSquaresProblem lsp = builder(dataset).build();

        final double expected = dataset.getResidualSumOfSquares();
        final double cost = lsp.evaluate(lsp.getStart()).getCost();
        final double actual = cost * cost;
        assertEquals(expected, actual, 1e-11 * expected, dataset.getName());
    }

    @Test
    void testComputeRMS() throws IOException {
        final StatisticalReferenceDataset dataset
            = StatisticalReferenceDatasetFactory.createKirby2();

        final LeastSquaresProblem lsp = builder(dataset).build();

        final double expected = FastMath.sqrt(dataset.getResidualSumOfSquares() /
                                              dataset.getNumObservations());
        final double actual = lsp.evaluate(lsp.getStart()).getRMS();
        assertEquals(expected, actual, 1e-11 * expected, dataset.getName());
    }

    @Test
    void testComputeSigma() throws IOException {
        final StatisticalReferenceDataset dataset
            = StatisticalReferenceDatasetFactory.createKirby2();

        final LeastSquaresProblem lsp = builder(dataset).build();

        final double[] expected = dataset.getParametersStandardDeviations();

        final Evaluation evaluation = lsp.evaluate(lsp.getStart());
        final double cost = evaluation.getCost();
        final RealVector sig = evaluation.getSigma(1e-14);
        final int dof = lsp.getObservationSize() - lsp.getParameterSize();
        for (int i = 0; i < sig.getDimension(); i++) {
            final double actual = FastMath.sqrt(cost * cost / dof) * sig.getEntry(i);
            assertEquals(expected[i], actual, 1e-6 * expected[i], dataset.getName() + ", parameter #" + i);
        }
    }

    @Test
    void testEvaluateCopiesPoint() throws IOException {
        //setup
        StatisticalReferenceDataset dataset
                = StatisticalReferenceDatasetFactory.createKirby2();
        LeastSquaresProblem lsp = builder(dataset).build();
        RealVector point = new ArrayRealVector(lsp.getParameterSize());

        //action
        Evaluation evaluation = lsp.evaluate(point);

        //verify
        assertNotSame(point, evaluation.getPoint());
        point.setEntry(0, 1);
        assertEquals(0, evaluation.getPoint().getEntry(0), 0);
    }

    @Test
    void testLazyEvaluation() {
        final RealVector dummy = new ArrayRealVector(new double[] { 0 });

        final LeastSquaresProblem p
            = LeastSquaresFactory.create(LeastSquaresFactory.model(dummyModel(), dummyJacobian()),
                                         dummy, dummy, null, null, 0, 0, true, null);

        // Should not throw because actual evaluation is deferred.
        final Evaluation eval = p.evaluate(dummy);

        try {
            eval.getResiduals();
            fail("Exception expected");
        } catch (RuntimeException e) {
            // Expecting exception.
            assertEquals("dummyModel", e.getMessage());
        }

        try {
            eval.getJacobian();
            fail("Exception expected");
        } catch (RuntimeException e) {
            // Expecting exception.
            assertEquals("dummyJacobian", e.getMessage());
        }
    }

    // MATH-1151
    @Test
    void testLazyEvaluationPrecondition() {
        final RealVector dummy = new ArrayRealVector(new double[] { 0 });

        // "ValueAndJacobianFunction" is required but we implement only
        // "MultivariateJacobianFunction".
        final MultivariateJacobianFunction m1 = notUsed -> new Pair<RealVector, RealMatrix>(null, null);

        try {
            // Should throw.
            LeastSquaresFactory.create(m1, dummy, dummy, null, null, 0, 0, true, null);
            fail("Expecting MathIllegalStateException");
        } catch (MathIllegalStateException e) {
            // Expected.
        }

        final MultivariateJacobianFunction m2 = new ValueAndJacobianFunction() {
                public Pair<RealVector, RealMatrix> value(RealVector notUsed) {
                    return new Pair<RealVector, RealMatrix>(null, null);
                }
                public RealVector computeValue(final double[] params) {
                    return null;
                }
                public RealMatrix computeJacobian(final double[] params) {
                    return null;
                }
            };

        // Should pass.
        LeastSquaresFactory.create(m2, dummy, dummy, null, null, 0, 0, true, null);
    }

    @Test
    void testDirectEvaluation() {
        final RealVector dummy = new ArrayRealVector(new double[] { 0 });

        final LeastSquaresProblem p
            = LeastSquaresFactory.create(LeastSquaresFactory.model(dummyModel(), dummyJacobian()),
                                         dummy, dummy, null, null, 0, 0, false, null);

        try {
            // Should throw.
            p.evaluate(dummy);
            fail("Exception expected");
        } catch (RuntimeException e) {
            // Expecting exception.
            // Whether it is model or Jacobian that caused it is not significant.
            final String msg = e.getMessage();
            assertTrue(msg.equals("dummyModel") ||
                              msg.equals("dummyJacobian"));
        }
    }

    /** Used for testing direct vs lazy evaluation. */
    private MultivariateVectorFunction dummyModel() {
        return p -> {
            throw new RuntimeException("dummyModel");
        };
    }

    /** Used for testing direct vs lazy evaluation. */
    private MultivariateMatrixFunction dummyJacobian() {
        return p -> {
            throw new RuntimeException("dummyJacobian");
        };
    }
}
