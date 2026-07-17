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
    
    /** Default tolerance for constraint satisfaction in case of degenerancy  */    
    public static final double DEFAULT_EPSILON_RELAXATION = 1.0e-9;//>0
     
     
    /** Tolerance for constraint satisfaction */
    private double eps;
    
    /** Tolerance for constraint satisfaction in case of degenerancy */
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
     *  typical 1.0e-3 to 1.0e-16
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
