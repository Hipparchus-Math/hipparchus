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

package org.hipparchus.analysis.solvers;

import org.hipparchus.analysis.polynomials.PolynomialFunction;

/**
 * Base class for solvers.
 *
 */
public abstract class AbstractPolynomialSolver
    extends BaseAbstractUnivariateSolver<PolynomialFunction>
    implements PolynomialSolver {
    /** Function. */
    private PolynomialFunction polynomialFunction;

    /**
     * Construct a solver with given absolute accuracy.
     *
     * @param absoluteAccuracy Maximum absolute error.
     */
    protected AbstractPolynomialSolver(final double absoluteAccuracy) {
        super(absoluteAccuracy);
    }
    /**
     * Construct a solver with given accuracies.
     *
     * @param relativeAccuracy Maximum relative error.
     * @param absoluteAccuracy Maximum absolute error.
     */
    protected AbstractPolynomialSolver(final double relativeAccuracy,
                                       final double absoluteAccuracy) {
        super(relativeAccuracy, absoluteAccuracy);
    }
    /**
     * Construct a solver with given accuracies.
     *
     * @param relativeAccuracy Maximum relative error.
     * @param absoluteAccuracy Maximum absolute error.
     * @param functionValueAccuracy Maximum function value error.
     */
    protected AbstractPolynomialSolver(final double relativeAccuracy,
                                       final double absoluteAccuracy,
                                       final double functionValueAccuracy) {
        super(relativeAccuracy, absoluteAccuracy, functionValueAccuracy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void setup(int maxEval, PolynomialFunction f,
                             double min, double max, double startValue) {
        super.setup(maxEval, f, min, max, startValue);
        polynomialFunction = f;
    }

    /** Get the coefficients of the polynomial function.
     * @return the coefficients of the polynomial function
     */
    protected double[] getCoefficients() {
        return polynomialFunction.getCoefficients();
    }
}
