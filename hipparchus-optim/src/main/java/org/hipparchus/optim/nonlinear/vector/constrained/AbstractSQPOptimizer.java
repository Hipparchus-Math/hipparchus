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
