package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.optim.OptimizationData;
import org.hipparchus.util.Precision;

/**
 *
 * @author rocca
 */
public class QPDualActiveSolverOption implements OptimizationData  {
    
    /** Default tolerance for constraint satisfaction */
    public static final double DEFAULT_EPSILON = Precision.EPSILON;//>0
    
     /** Default Matrix Mode */
    public static final QPMatrixMode DEFAULT_MATRIX_MODE = QPMatrixMode.FULL;
    
    /**
    *  Default Relaxed feasibility tolerance used only when both the primal and dual
     * step are not available. In that numerically stalled case, a target
     * constraint whose violation does not exceed this value is accepted as
     * numerically feasible.
     */   
    public static final double DEFAULT_EPSILON_RELAXATION = 1.0e-12;//>0
     
     
    /** Tolerance for constraint satisfaction */
    private double eps;
    
    /**
     * Relaxed feasibility tolerance used only when both the primal and dual
     * step are not available for ill conditioning problems. In that numerically stalled case, a target
     * constraint whose violation does not exceed this value is accepted as
     * numerically feasible. For the major part of the problems choosing 10.0*eps is enough
     */
    private double epsRelaxation;
    
    /** Matrix Mode */
    private QPMatrixMode matrixMode;
    
    
    public QPDualActiveSolverOption() 
    {
        this.eps=DEFAULT_EPSILON;
        this.epsRelaxation=DEFAULT_EPSILON_RELAXATION;
        this.matrixMode=DEFAULT_MATRIX_MODE;
    }
    
    /** Set tolerance for convergence 
     *  typical 1.0e-16 to 1.0e-9
     * @param eps tolerance for convergence 
     */
    public void setEps(final double eps) {
        this.eps = eps;
    }
    
    /** get tolerance for convergence and active constraint evaluation
     * @return .*/
    public double getEps() {
        return eps;
    }
    
     /** Set tolerance convergence for constraint in case of degenerancy
     *  typical 10.0*eps
     * @param eps  tolerance for convergence 
     */
    public void setEpsRelaxation(final double eps) {
        this.epsRelaxation = eps;
    }
    
    /** get tolerance for convergence for constraint in case of degenerancy
     * @return .*/
    public double getEpsRelaxation() {
        return epsRelaxation;
    }
    
     /** Set MatrixMode 
      *
      *@param matrixMode  
      */
    public void setMatrixMode(final QPMatrixMode matrixMode)
    {
        this.matrixMode=matrixMode;
    }
    
    /** Get MatrixMode
     * @return  */
     public QPMatrixMode getMatrixMode()
    {
      return  this.matrixMode;
    }
}
