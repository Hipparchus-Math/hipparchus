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
import org.hipparchus.linear.EigenDecompositionSymmetric;
import org.hipparchus.optim.LocalizedOptimFormats;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.MathUtils;

/**
 * Abstract class for Sequential Quadratic Programming solvers
 * @since 3.1
 */
public abstract class AbstractSQPOptimizer2 extends ConstraintOptimizer {

    /** Algorithm settings. */
    private SQPOption settings;

    /** Tolerance for symmetric matrix decomposition.
     * @since 4.1
     */
    private MatrixDecompositionTolerance matrixDecompositionTolerance;

    /** Objective function. */
    private TwiceDifferentiableFunction OBJ;

    /** Equality constraint (may be null). */
    private EqualityConstraint EQ;

    /** Inequality constraint (may be null). */
    private InequalityConstraint IQ;

    /** Inequality constraint (may be null). */
    private BoundedConstraint BOX;
    
    /** Simple Bounds (may be null). */
    private SimpleBounds SB;
    
    /** Default QPSolver. */
    private QPOptimizer QPSolver = new QPDualActiveSolver();

    /** Simple constructor.
     */
    protected AbstractSQPOptimizer2() {
        this.settings                     = new SQPOption();
        this.matrixDecompositionTolerance = new MatrixDecompositionTolerance(EigenDecompositionSymmetric.DEFAULT_EPSILON);
    }

    /** Getter for settings.
     * @return settings
     */
    public SQPOption getSettings() {
        return settings;
    }

    /** Getter for matrix decomposition tolerance.
     * @return matrix decomposition tolerance
     * @since 4.1
     */
    public MatrixDecompositionTolerance getMatrixDecompositionTolerance() {
        return matrixDecompositionTolerance;
    }

    /** Getter for objective function.
     * @return objective function
     */
    public TwiceDifferentiableFunction getObj() {
        return OBJ;
    }

    /** Getter for equality constraint.
     * @return equality constraint
     */
    public EqualityConstraint getEqConstraint() {
        return EQ;
    }

    /** Getter for inequality constraint.
     * @return inequality constraint
     */
    public InequalityConstraint getIqConstraint() {
        return IQ;
    }

     /** Getter for box constraint.
     * @return inequality constraint
     */
    public BoundedConstraint getBoxConstraint() {
        return BOX;
    }
    
     /** Getter for simple bounds.
     * @return simple bounds
     */
    public SimpleBounds getSimpleBounds() {
        return SB;
    }

    /** Getter for QP Solver.
     * @return QP Solver
     */
    public QPOptimizer getQPSolver() {
        return QPSolver;
    }

    @Override
    public LagrangeSolution optimize(OptimizationData... optData) {
        return super.optimize(optData);
    }

    @Override
    public void parseOptimizationData(OptimizationData... optData) {
        super.parseOptimizationData(optData);
        for (OptimizationData data : optData) {

            if (data instanceof SQPProblem) {
                SQPProblem problem = (SQPProblem) data;
                OBJ = new SQPObj(problem);
                
                EQ= (problem.hasEquality())?new SQPEq((SQPProblem) data):null;
                IQ= (problem.hasInequality())?new SQPIneq((SQPProblem) data):null;
                
                double[] lb = ((SQPProblem) data).getBoxConstraintLB();
                double[]ub = ((SQPProblem) data).getBoxConstraintUB();
                SB = (problem.hasBounds())?new SimpleBounds(lb,ub):null;
               
                continue;
            }
            
            if (data instanceof ObjectiveFunction) {
                OBJ = (TwiceDifferentiableFunction) ((ObjectiveFunction) data).getObjectiveFunction();
                continue;
            }

            if (data instanceof EqualityConstraint) {
                EQ = (EqualityConstraint) data;
                continue;
            }
            if (data instanceof InequalityConstraint) {
                IQ = (InequalityConstraint) data;
                continue;
            }

            if (data instanceof BoundedConstraint) {
                BOX = (BoundedConstraint) data;
                continue;
            }
            
            if (data instanceof SimpleBounds) {
                SB = (SimpleBounds) data;
                continue;
            }

            if (data instanceof SQPOption) {
                settings = (SQPOption) data;
            }

            if (data instanceof QPOptimizer) {
                QPSolver = (QPOptimizer) data;
            }

            if (data instanceof MatrixDecompositionTolerance) {
                matrixDecompositionTolerance = (MatrixDecompositionTolerance) data;
            }

        }

        // if we got here, convexObjective exists
        int n = OBJ.dim();
        if (EQ != null) {
            int nDual = EQ.dimY();
            if (nDual > n) {
                throw new MathIllegalArgumentException(LocalizedOptimFormats.CONSTRAINTS_RANK, nDual, n);
            }
            int nTest = EQ.dim();
            if (nDual == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_NOT_ALLOWED);
            }
            MathUtils.checkDimension(nTest, n);
        }

    }

}
