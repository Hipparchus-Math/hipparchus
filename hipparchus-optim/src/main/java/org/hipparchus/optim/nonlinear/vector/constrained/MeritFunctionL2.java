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
import org.hipparchus.util.Precision;


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
    
    /** bounds constraints (may be null). */
    private final LinearInequalityConstraint bounds;

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
    
    /** R vector old. */
    private final RealVector rOld;

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
    
    /** bounds evaluation. */
    private RealVector bEval;

    /** Penalty evaluation. */
    private double pEval;

    /** Gradient of equality. */
    private RealMatrix JE;

    /** Gradient of inequality. */
    private RealMatrix JI;
    
    /** Gradient of bounds. */
    private RealMatrix JB;
    
    /** initial Sigma . */
    private double sigmaInit=1.0e-2;
    
    /** rMax . */
    private double rMax=1.0e9;

    /**
     * Constructor.
     *
     * @param objective Objective function
     * @param eqConstraint Equality constraint (may be null)
     * @param iqConstraint Inequality constraint (may be null)
     * @param x current point
     * @deprecated as of 4.1, replaced by {@link
     * #MeritFunctionL2(TwiceDifferentiableFunction, Constraint, Constraint,
     * LinearInequalityConstraint, RealVector)}
     */
    @Deprecated
    public MeritFunctionL2(final TwiceDifferentiableFunction objective,
                           final Constraint eqConstraint,
                           final Constraint iqConstraint,
                           final RealVector x) {
        this(objective, eqConstraint, iqConstraint, null, x);
    }

    /**
     * Constructor.
     *
     * @param objective Objective function
     * @param eqConstraint Equality constraint (may be null)
     * @param iqConstraint Inequality constraint (may be null)
     * @param bounds Linear inequality constraint (may be null)
     * @param x current point
     * @since 4.1
     */
    public MeritFunctionL2(final TwiceDifferentiableFunction objective,
                           final Constraint eqConstraint,
                           final Constraint iqConstraint,
                           final LinearInequalityConstraint bounds,
                           final RealVector x) {
        this.objective = objective;
        this.eqConstraint = eqConstraint;
        this.iqConstraint = iqConstraint;
        this.bounds=bounds;
        this.x = new ArrayRealVector(x);
        int me = 0;
        int mi = 0;
        int mb = 0;
        if (this.eqConstraint != null) {
            me = this.eqConstraint.dimY();
        }
        if (this.iqConstraint != null) {
            mi = this.iqConstraint.dimY();
        }
        
        if (this.bounds != null) {
            mb = this.bounds.dimY();
            JB=this.bounds.jacobian(null);
        }
        final int m = me + mi + mb;
        this.dx = new ArrayRealVector(x.getDimension());
        this.y = new ArrayRealVector(m);
        this.u = new ArrayRealVector(m);
        this.r = new ArrayRealVector(m,sigmaInit);
        this.rOld=new ArrayRealVector(m,sigmaInit);
        this.J = new ArrayRealVector(x.getDimension());
        this.eqEval = new ArrayRealVector(me);
        this.iqEval = new ArrayRealVector(mi);
        this.bEval = new ArrayRealVector(mb);
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
     * Get Last X;
     * @return x
     */
    RealVector getX() {
        return this.x;
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
     * Get Last Bounds Evaluation;
     * @return penalty gradient
     */
    RealVector getBEval() {
        return this.bEval;
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
        
        int mb = 0;
        if (bounds != null) {
            mb = bounds.dimY();
            RealVector ri = r.getSubVector(me+mi, mb);
            RealVector yi = yAlpha.getSubVector(me+mi, mb);

            RealVector yk = yAlpha.getSubVector(me+mi, mb);

            bEval = bounds.value(xAlpha);
            RealVector gk = bEval.subtract(bounds.getLowerBound());

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
        int mi = 0;
        int mb = 0;
      
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
       
        if (bounds != null) {
            mb = bounds.dimY();

            RealVector rb = r.getSubVector(me+mi, mb);
           
            RealVector yb = y.getSubVector(me+mi, mb);
            
            RealVector gb = this.bEval.subtract(bounds.getLowerBound());
           
            RealMatrix jacob = JB;
            RealVector maskb = new ArrayRealVector(mb, 1.0);

            for (int i = 0; i < gb.getDimension(); i++) {
                if (gb.getEntry(i) > yb.getEntry(i) / rb.getEntry(i)) {
                    maskb.setEntry(i, 0.0);
                }
            }

            RealVector termb=jacob.preMultiply((yb.subtract(gb.ebeMultiply(rb))).ebeMultiply(maskb));
            partial=partial.subtract(termb);
        }
       
        return partial;
    }


     public RealVector gradY() {

        int me = 0;
        int mi = 0;
        int mb = 0;
        RealVector partial = new ArrayRealVector(y.getDimension());
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
            RealVector g = this.eqEval.subtract(eqConstraint.getLowerBound());
            partial.setSubVector(0, g.mapMultiply(-1.0));
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();

            RealVector ri = r.getSubVector(me, mi);
            RealVector riOld = rOld.getSubVector(me, mi);
            RealVector yi = y.getSubVector(me, mi);
            RealVector gi = this.iqEval.subtract(iqConstraint.getLowerBound());
            RealVector viri = new ArrayRealVector(mi, 0.0);
            for (int i = 0; i < gi.getDimension(); i++) {
                viri.setEntry(i,
                              gi.getEntry(i) > yi.getEntry(i) / ri.getEntry(i) ?
                              -yi.getEntry(i)/riOld.getEntry(i) :
                              -gi.getEntry(i));
            }

            partial.setSubVector(me, viri);
        }
        if (bounds != null) {
            mb = bounds.dimY();

            RealVector rb = r.getSubVector(me+mi, mb);
            RealVector rbOld = rOld.getSubVector(me+mi, mb);
            RealVector yb = y.getSubVector(me+mi, mb);
            RealVector gb = this.bEval.subtract(bounds.getLowerBound());
            RealVector viri = new ArrayRealVector(mb, 0.0);
            for (int i = 0; i < gb.getDimension(); i++) {
                viri.setEntry(i,
                              gb.getEntry(i) > yb.getEntry(i) / rb.getEntry(i) ?
                              -yb.getEntry(i)/rbOld.getEntry(i) :
                              -gb.getEntry(i));
            }

            partial.setSubVector(me+mi, viri);
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
     
     
     public void updateRj(RealMatrix H,
                     RealVector newY,   
                     RealVector newDx,  
                     RealVector newU,  
                     double sigmaValue, 
                     int iterations) {  
         
    if (newY == null || newY.getDimension() == 0) {
        return;
    }
    this.rOld.setSubVector(0, r.copy());
    final int m = r.getDimension();
    final double dbd = newDx.dotProduct(H.operate(newDx));                         
    final int iter=iterations+1;
   
   
    double dbdsigma = rMax ;
    if (dbd * (1.0 - sigmaValue)  >Precision.SAFE_MIN) 
       
        dbdsigma = (2.0 * m) / (dbd * (1.0 - sigmaValue));
       
    
    
    

    for (int j = 0; j < m; j++) {
        final double rj = r.getEntry(j);
        final double diff = newU.getEntry(j) - newY.getEntry(j);  
        final double rNew =  Math.min(rMax,dbdsigma * diff*diff);
        
        
        
        final double rjSigmaj  = Math.min(rj, iter * Math.sqrt(rj)); 
        r.setEntry(j, Math.max(rNew, rjSigmaj));                       
    }
}



    /** Reset R vector to initial value.
     */
    public void resetRj() {
        if (y.getDimension() > 0) {
            this.r.set(sigmaInit);
            this.rOld.set(sigmaInit);
        }
     }

    /** Get search direction.
     * @return search direction
     */
    RealVector getDx() {
       return  this.dx;
    }

    //
    public void rUp()
    {
        if(r!=null)
        {   //rOld.setSubVector(0, r.copy());
            for(int i=0;i<r.getDimension();i++)
                r.setEntry(i, FastMath.min(1.0e9, r.getEntry(i)*10.0));
        }
           
            
    }
    public RealVector getR()
    {
            return this.r;
}
    public double getRmax()
    {
            return this.rMax;
    }
}
