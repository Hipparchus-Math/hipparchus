/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.fitting;

import java.util.Collection;

import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.linear.DiagonalMatrix;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresBuilder;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem;

/**
 * Fits points to a {@link
 * org.hipparchus.analysis.polynomials.PolynomialFunction.Parametric polynomial}
 * function.
 * <br>
 * The size of the {@link #withStartPoint(double[]) initial guess} array defines the
 * degree of the polynomial to be fitted.
 * They must be sorted in increasing order of the polynomial's degree.
 * The optimal values of the coefficients will be returned in the same order.
 *
 */
public class PolynomialCurveFitter extends AbstractCurveFitter {
    /** Parametric function to be fitted. */
    private static final PolynomialFunction.Parametric FUNCTION = new PolynomialFunction.Parametric();
    /** Initial guess. */
    private final double[] initialGuess;
    /** Maximum number of iterations of the optimization algorithm. */
    private final int maxIter;

    /**
     * Constructor used by the factory methods.
     *
     * @param initialGuess Initial guess.
     * @param maxIter Maximum number of iterations of the optimization algorithm.
     * @throws MathRuntimeException if {@code initialGuess} is {@code null}.
     */
    private PolynomialCurveFitter(double[] initialGuess, int maxIter) {
        this.initialGuess = initialGuess.clone();
        this.maxIter = maxIter;
    }

    /**
     * Creates a default curve fitter.
     * Zero will be used as initial guess for the coefficients, and the maximum
     * number of iterations of the optimization algorithm is set to
     * {@link Integer#MAX_VALUE}.
     *
     * @param degree Degree of the polynomial to be fitted.
     * @return a curve fitter.
     *
     * @see #withStartPoint(double[])
     * @see #withMaxIterations(int)
     */
    public static PolynomialCurveFitter create(int degree) {
        return new PolynomialCurveFitter(new double[degree + 1], Integer.MAX_VALUE);
    }

    /**
     * Configure the start point (initial guess).
     * @param newStart new start point (initial guess)
     * @return a new instance.
     */
    public PolynomialCurveFitter withStartPoint(double[] newStart) {
        return new PolynomialCurveFitter(newStart.clone(),
                                         maxIter);
    }

    /**
     * Configure the maximum number of iterations.
     * @param newMaxIter maximum number of iterations
     * @return a new instance.
     */
    public PolynomialCurveFitter withMaxIterations(int newMaxIter) {
        return new PolynomialCurveFitter(initialGuess,
                                         newMaxIter);
    }

    /** {@inheritDoc} */
    @Override
    protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> observations) {
        // Prepare least-squares problem.
        final int len = observations.size();
        final double[] target  = new double[len];
        final double[] weights = new double[len];

        int i = 0;
        for (WeightedObservedPoint obs : observations) {
            target[i]  = obs.getY();
            weights[i] = obs.getWeight();
            ++i;
        }

        final AbstractCurveFitter.TheoreticalValuesFunction model =
                new AbstractCurveFitter.TheoreticalValuesFunction(FUNCTION, observations);

        if (initialGuess == null) {
            throw MathRuntimeException.createInternalError();
        }

        // Return a new least squares problem set up to fit a polynomial curve to the
        // observed points.
        return new LeastSquaresBuilder().
                maxEvaluations(Integer.MAX_VALUE).
                maxIterations(maxIter).
                start(initialGuess).
                target(target).
                weight(new DiagonalMatrix(weights)).
                model(model.getModelFunction(), model.getModelFunctionJacobian()).
                build();

    }

}
