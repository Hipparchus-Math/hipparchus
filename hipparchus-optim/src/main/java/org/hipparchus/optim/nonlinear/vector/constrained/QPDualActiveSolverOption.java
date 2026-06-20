package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.optim.OptimizationData;
import org.hipparchus.util.Precision;

/**
 *
 * @author rocca
 */
public class QPDualActiveSolverOption implements OptimizationData  {
    
    /** Default tolerance for constraint satisfaction */
    public static final double DEFAULT_EPSILON = 1.0e-15;//>0
    
    private double eps;
    
    
    public QPDualActiveSolverOption() 
    {
        this.eps=DEFAULT_EPSILON;
    }
    
    /** Set tolerance for convergence 
     *  typical 1.0e-16 to 1.0e-9
     * @param eps tolerance for convergence 
     */
    public void setEps(final double eps) {
        this.eps = eps;
    }
    
    /** get tolerance for convergence and active constraint evaluation.
     */
    public double getEps() {
        return eps;
    }
}
