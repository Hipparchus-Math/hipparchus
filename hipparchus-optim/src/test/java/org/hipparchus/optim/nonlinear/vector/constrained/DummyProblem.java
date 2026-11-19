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
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

class DummyProblem implements SQPProblem {

    private final boolean hasBounds;
    private final boolean hasEquality;
    private final boolean hasInequality;

    public DummyProblem(final boolean hasBounds,
                        final boolean hasEquality,
                        final boolean hasInequality) {
        this.hasBounds     = hasBounds;
        this.hasEquality   = hasEquality;
        this.hasInequality = hasInequality;
    }

    @Override
    public int getdim() {
        return 3;
    }

    @Override
    public boolean hasBounds() {
        return hasBounds;
    }

    @Override
    public double[] getBoundsLB() {
        return hasBounds ? new double[] { 0.0, 0.0, 0.0 } : null;
    }

    @Override
    public double[] getBoundsUB() {
        return hasBounds ? new double[] { 100.0, 200.0, 300.0 } : null;
    }

    @Override
    public double getObjectiveEvaluation(final RealVector rv) {
        final double x = rv.getEntry(0);
        final double y = rv.getEntry(1);
        final double z = rv.getEntry(2);
        return x * x + y * y + z * z;
    }

    @Override
    public RealVector getObjectiveGradient(final RealVector rv) {
        final double x = rv.getEntry(0);
        final double y = rv.getEntry(1);
        final double z = rv.getEntry(2);
        return MatrixUtils.createRealVector(new double[] { 2 * x, 2 * y, 2 * z });
    }

    @Override
    public boolean hasEquality() {
        return hasEquality;
    }

    @Override
    public RealVector getEqConstraintEvaluation(final RealVector rv) {
        if (hasEquality) {
           final double x = rv.getEntry(0);
           final double y = rv.getEntry(1);
           return MatrixUtils.createRealVector(new double[] { x + y });
        } else {
            return null;
        }
    }

    @Override
    public RealMatrix getEqConstraintJacobian(final RealVector rv) {
        return hasEquality ?
               MatrixUtils.createRealMatrix(new double[][] { { 1.0, 1.0, 0.0 } }) :
               null;
    }

    @Override
    public RealVector getEqConstraintLB() {
        return hasEquality ?
               MatrixUtils.createRealVector(new double[] { 1.0 }) :
               null;
    }

    @Override
    public boolean hasInequality() {
        return hasInequality;
    }

    @Override
    public RealVector getIneqConstraintEvaluation(final RealVector rv) {
        if (hasInequality) {
            final double x = rv.getEntry(0);
            final double z = rv.getEntry(2);
            return MatrixUtils.createRealVector(new double[] { x + z });
        } else {
            return null;
        }
    }

    @Override
    public RealMatrix getIneqConstraintJacobian(final RealVector rv) {
        return hasInequality ?
               MatrixUtils.createRealMatrix(new double[][] { { 1.0, 0.0, 1.0 } }) :
               null;
    }

    @Override
    public RealVector getIneqConstraintLB() {
        return hasInequality ?
               MatrixUtils.createRealVector(new double[] { 2.0 }) :
               null;
    }

}
