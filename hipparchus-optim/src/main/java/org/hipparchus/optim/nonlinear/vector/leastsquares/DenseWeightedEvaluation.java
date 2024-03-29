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
package org.hipparchus.optim.nonlinear.vector.leastsquares;

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem.Evaluation;

/**
 * Applies a dense weight matrix to an evaluation.
 *
 */
class DenseWeightedEvaluation extends AbstractEvaluation {

    /** the unweighted evaluation */
    private final Evaluation unweighted;
    /** reference to the weight square root matrix */
    private final RealMatrix weightSqrt;

    /**
     * Create a weighted evaluation from an unweighted one.
     *
     * @param unweighted the evalutation before weights are applied
     * @param weightSqrt the matrix square root of the weight matrix
     */
    DenseWeightedEvaluation(final Evaluation unweighted,
                            final RealMatrix weightSqrt) {
        // weight square root is square, nR=nC=number of observations
        super(weightSqrt.getColumnDimension());
        this.unweighted = unweighted;
        this.weightSqrt = weightSqrt;
    }

    /* apply weights */

    /** {@inheritDoc} */
    @Override
    public RealMatrix getJacobian() {
        return weightSqrt.multiply(this.unweighted.getJacobian());
    }

    /** {@inheritDoc} */
    @Override
    public RealVector getResiduals() {
        return this.weightSqrt.operate(this.unweighted.getResiduals());
    }

    /* delegate */

    /** {@inheritDoc} */
    @Override
    public RealVector getPoint() {
        return unweighted.getPoint();
    }

}
