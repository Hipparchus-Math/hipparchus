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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;


/** Augmented Penalty Function.
 * <p>
 * This class computes the penalty function and its gradient, combining:
 * </p>
 * <ul>
 *   <li>The objective function</li>
 *   <li>Equality constraints</li>
 *   <li>Inequality constraints</li>
 * </ul>
 * <p>
 * Typical usage:
 * </p>
 * <ul>
 *   <li>Call update(...) before line search or Hessian update</li>
 *   <li>Use value(alpha) to evaluate penalty at x + alpha * dx(this store also evalutation of objective and constraints)</li>
 *   <li>Use gradient() to retrieve penalty gradient at current point</li>
 * </ul>
 */
public class MeritFunctionL2 {

    /** Objective function. */
    private final TwiceDifferentiableFunction objective;

    /** Equality constraints (may be null). */
    private final Constraint eqConstraint;

    /** Inequality constraints (may be null). */
    private final Constraint iqConstraint;

    /** Current point. */
    private RealVector x;

    /** Lagrange multipliers. */
    private RealVector y;

    /** Search direction. */
    private RealVector       dx;

    /** Multipliers for QP. */
    private RealVector u;

    /** R vector. */
    private final RealVector r;

    /** Gradient of the objective at current point. */
    private RealVector J;

    /** penalty gradient. */
    private double penaltyGradient;

    /** Objective function evaluation. */
    private double objEval;

    /** Equality evaluation. */
    private RealVector eqEval;

    /** Inequality evaluation. */
    private RealVector iqEval;

    /** Penalty evaluation. */
    private double pEval;

    /** Gradient of equality. */
    private RealMatrix JE;

    /** Gradient of inequality. */
    private RealMatrix JI;

    /**
     * Constructor.
     *
     * @param objective Objective function
     * @param eqConstraint Equality constraint (may be null)
     * @param iqConstraint Inequality constraint (may be null)
     * @param x current point
     */
    public MeritFunctionL2(final TwiceDifferentiableFunction objective,
                           final Constraint eqConstraint,
                           final Constraint iqConstraint,
                           final RealVector x) {
        this.objective = objective;
        this.eqConstraint = eqConstraint;
        this.iqConstraint = iqConstraint;
        this.x = new ArrayRealVector(x);
        int me = 0;
        int mi = 0;
        if (this.eqConstraint != null) {
            me = this.eqConstraint.dimY();
        }
        if (this.iqConstraint != null) {
            mi = this.iqConstraint.dimY();
        }
        final int m = me + mi;
        this.dx = new ArrayRealVector(x.getDimension());
        this.y = new ArrayRealVector(m);
        this.u = new ArrayRealVector(m);
        this.r = new ArrayRealVector(m,1.0);
        this.J = new ArrayRealVector(x.getDimension());
        this.eqEval = new ArrayRealVector(me);
        this.iqEval = new ArrayRealVector(mi);
        //this evaluate objective function contraints function and penaly
        this.value(0);
    }

    /** Update internal parameters for next penalty computation.
     * @param newJ Gradient of objective at current x
     * @param newJE Gradient of equality x
     * @param newJI Gradient of inequality at current x
     * @param newX Current iterate
     * @param newY Lagrange multipliers
     * @param newDx Search direction from QP
     * @param newU multiplier from QP
     */
    public void update(final RealVector newJ, final RealMatrix newJE, final RealMatrix newJI,
                       final RealVector newX, final RealVector newY,
                       final RealVector newDx, final RealVector newU) {
        this.J               = newJ;
        this.JE              = newJE;
        this.JI              = newJI;
        this.x               = newX;
        this.y               = newY;
        this.dx              = newDx;
        this.u               = newU;
        this.penaltyGradient = gradient();
    }

    /**
     * Penalty Gradient
     * @return penalty gradient
     */
    double getGradient() {
        return this.penaltyGradient;
    }

    /**
     * Get Last Objective Evaluation;
     * @return penalty gradient
     */
    double getObjEval() {
        return this.objEval;
    }

    /**
     * Get Last Inequality Constraints Evaluation;
     * @return penalty gradient
     */
    RealVector getIqEval() {
        return this.iqEval;
    }

     /**
     * Get Last Equality Constraints Evaluation;
     * @return penalty gradient
     */
    RealVector getEqEval() {
        return this.eqEval;
    }

    /**
     * Get last penalty evaluation.
     * @return lat penalty evaluation
     */
    double getPenaltyEval() {
        return this.pEval;
    }

    /**
     * Evaluate penalty function at x + alpha * dx.
     *
     * @param alpha Step length
     * @return penalty value
     */
    public double value(double alpha) {
        RealVector xAlpha = x.add(dx.mapMultiply(alpha));
        RealVector yAlpha = null;
        if (y.getDimension() > 0) {
            yAlpha = y.add(u.subtract(y).mapMultiply(alpha));
        }

        objEval = this.objective.value(xAlpha);
        double penalty =objEval;

        int me = 0;
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
            RealVector re = r.getSubVector(0, me);
            RealVector ye = yAlpha.getSubVector(0, me);
            eqEval = eqConstraint.value(xAlpha);
            RealVector g = eqEval.subtract(eqConstraint.getLowerBound());

            RealVector g2 = g.ebeMultiply(g);
            penalty -= ye.dotProduct(g) - 0.5 * re.dotProduct(g2);
        }

        int mi = 0;
        if (iqConstraint != null) {
            mi = iqConstraint.dimY();
            RealVector ri = r.getSubVector(me, mi);
            RealVector yi = yAlpha.getSubVector(me, mi);

            RealVector yk = yAlpha.getSubVector(me, mi);

            iqEval = iqConstraint.value(xAlpha);
            RealVector gk = iqEval.subtract(iqConstraint.getLowerBound());

            RealVector g = new ArrayRealVector(gk);
            RealVector mask = new ArrayRealVector(g.getDimension(), 1.0);

            for (int i = 0; i < gk.getDimension(); i++) {
                if (gk.getEntry(i) > (yk.getEntry(i) / ri.getEntry(i))) {
                    mask.setEntry(i, 0.0);
                    penalty -= 0.5 * yi.getEntry(i) * yi.getEntry(i) / ri.getEntry(i);
                }
            }

            RealVector g2 = g.ebeMultiply(g.ebeMultiply(mask));
            penalty -= yi.dotProduct(g.ebeMultiply(mask)) - 0.5 * ri.dotProduct(g2);
        }
        pEval = penalty;
        return penalty;
    }

    /**
     * Get penalty gradient at current x.
     *
     * @return penalty gradient
     */
    private double gradient() {
        if (y.getDimension() > 0) {
            return gradX().dotProduct(dx) + gradY().dotProduct(u.subtract(y));
        } else {
            return gradX().dotProduct(dx);
        }
    }


     public RealVector gradX() {
        RealVector partial = J;
        int me = 0;
        int mi;

        if (eqConstraint != null) {
            me = eqConstraint.dimY();
            RealVector re = r.getSubVector(0, me);
            RealVector ye = y.getSubVector(0, me);

            RealVector ge = this.eqEval.subtract(eqConstraint.getLowerBound());
            RealMatrix jacob = JE;
            RealVector term = jacob.preMultiply(ye.subtract(ge.ebeMultiply(re)));
            partial = partial.subtract(term);
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();

            RealVector ri = r.getSubVector(me, mi);
            RealVector yi = y.getSubVector(me, mi);
            RealVector gi = this.iqEval.subtract(iqConstraint.getLowerBound());
            RealMatrix jacob = JI;
            RealVector mask = new ArrayRealVector(mi, 1.0);

            for (int i = 0; i < gi.getDimension(); i++) {
                if (gi.getEntry(i) > yi.getEntry(i) / ri.getEntry(i)) {
                    mask.setEntry(i, 0.0);
                }
            }

            RealVector term=jacob.preMultiply((yi.subtract(gi.ebeMultiply(ri))).ebeMultiply(mask));
            partial=partial.subtract(term);
        }

        return partial;
    }


     public RealVector gradY() {

        int me = 0;
        int mi;
        RealVector partial = new ArrayRealVector(y.getDimension());
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
            RealVector g = this.eqEval.subtract(eqConstraint.getLowerBound());
            partial.setSubVector(0, g.mapMultiply(-1.0));
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();

            RealVector ri = r.getSubVector(me, mi);
            RealVector yi = y.getSubVector(me, mi);
            RealVector gi = this.iqEval.subtract(iqConstraint.getLowerBound());
            RealVector viri = new ArrayRealVector(mi, 0.0);
            for (int i = 0; i < gi.getDimension(); i++) {
                viri.setEntry(i,
                              gi.getEntry(i) > yi.getEntry(i) / ri.getEntry(i) ?
                              -yi.getEntry(i)/ri.getEntry(i) :
                              -gi.getEntry(i));
            }

            partial.setSubVector(me, viri);
        }

        return partial;
    }

    /**
     * Update weight vector Rj.
     * called after QP solution before update the penalty function
     * @param H hessina Matrix (updated after line search)
     * @param newY last estimate of multiplier (updated after line search)
     * @param newDx direction of x provided by QP solution
     * @param newU direction of y provided by QP solution
     * @param sigmaValue value of the additional variable of QP solution
     * @param iterations current iteration
     */
     public void updateRj(RealMatrix H, RealVector newY, RealVector newDx, RealVector newU, double sigmaValue, int iterations) {
        // calculate sigma vector that depends on iterations
        if (newY.getDimension() == 0) {
            return;
        }
        RealVector sigma = new ArrayRealVector(r.getDimension());
        for (int i = 0; i < sigma.getDimension(); i++) {
            final double appoggio = iterations / FastMath.sqrt(r.getEntry(i));
            sigma.setEntry(i, FastMath.min(1.0, appoggio));
        }

        int me = 0;
        int mi = 0;
        if (this.eqConstraint != null) {
            me = this.eqConstraint .dimY();
        }
        if (this.iqConstraint != null) {
            mi = this.iqConstraint.dimY();
        }

        RealVector sigmar = sigma.ebeMultiply(r);
        //(u-v)^2 or (ru-v)
        RealVector numerator = ((newU.subtract(newY)).ebeMultiply(newU.subtract(newY))).mapMultiply(2.0 * (mi + me));

        double denominator = newDx.dotProduct(H.operate(newDx)) * (1.0 - sigmaValue);
        RealVector r1 = new ArrayRealVector(r);
        if (this.eqConstraint != null) {
            for (int i = 0; i < me; i++) {
                r1.setEntry(i,
                            FastMath.max(sigmar.getEntry(i),
                                         numerator.getEntry(i) / denominator));
            }
        }
        if (this.iqConstraint != null) {
            for (int i = 0; i < mi; i++) {
                r1.setEntry(me + i,
                            FastMath.max(sigmar.getEntry(me + i),
                                         numerator.getEntry(me + i) / denominator));
            }
        }

        r.setSubVector(0, r1);
    }

    /** Reset R vector to unity.
     */
    public void resetRj() {
        if (y.getDimension() > 0) {
            this.r.set(1.0);
        }
     }

    /** Get search direction.
     * @return search direction
     */
    RealVector getDx() {
       return  this.dx;
    }

}
