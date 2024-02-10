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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.optim.LocalizedOptimFormats;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.MathUtils;

/**
 * Abstract class for Sequential Quadratic Programming solvers
 * @since 3.1
 */
public abstract class AbstractSQPOptimizer extends ConstraintOptimizer {

    /** Algorithm settings. */
    private SQPOption settings;

    /** Objective function. */
    private TwiceDifferentiableFunction obj;

    /** Equality constraint (may be null). */
    private EqualityConstraint eqConstraint;

    /** Inequality constraint (may be null). */
    private InequalityConstraint iqConstraint;

    /** Simple constructor.
     */
    protected AbstractSQPOptimizer() {
        this.settings = new SQPOption();
    }

    /** Getter for settings. */
    public SQPOption getSettings() {
        return settings;
    }

    /** Getter for objective function. */
    public TwiceDifferentiableFunction getObj() {
        return obj;
    }

    /** Getter for equality constraint. */
    public EqualityConstraint getEqConstraint() {
        return eqConstraint;
    }

    /** Getter for inequality constraint. */
    public InequalityConstraint getIqConstraint() {
        return iqConstraint;
    }

    @Override
    public LagrangeSolution optimize(OptimizationData... optData) {
        return super.optimize(optData);
    }

    @Override
    protected void parseOptimizationData(OptimizationData... optData) {
        super.parseOptimizationData(optData);
        for (OptimizationData data : optData) {

            if (data instanceof ObjectiveFunction) {
                obj = (TwiceDifferentiableFunction) ((ObjectiveFunction) data).getObjectiveFunction();
                continue;
            }

            if (data instanceof EqualityConstraint) {
                eqConstraint = (EqualityConstraint) data;
                continue;
            }
            if (data instanceof InequalityConstraint) {
                iqConstraint = (InequalityConstraint) data;
                continue;
            }

            if (data instanceof SQPOption) {
                settings = (SQPOption) data;
                continue;
            }

        }

        // if we got here, convexObjective exists
        int n = obj.dim();
        if (eqConstraint != null) {
            int nDual = eqConstraint.dimY();
            if (nDual >= n) {
                throw new MathIllegalArgumentException(LocalizedOptimFormats.CONSTRAINTS_RANK, nDual, n);
            }
            int nTest = eqConstraint.dim();
            if (nDual == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_NOT_ALLOWED);
            }
            MathUtils.checkDimension(nTest, n);
        }

    }


}
