package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;

public class BFGSUpdater1 {

    private final LDL ldl;
    private final RealMatrix initialH;
    private final double EPS;
    private static final double GAMMA = 0.2;
    private final double[][]lcopy;
    private final double[]dcopy;
    private double skipCount=0;

    public BFGSUpdater1(final RealMatrix initialHess, final double eps, final boolean autoScale, final double decompositionEpsilon) {
        this.initialH = initialHess.copy();
        int n=this.initialH.getColumnDimension();
        this.EPS = eps;
        this.ldl = new LDL(initialHess.getRowDimension());
        lcopy=new double[n][n];
        dcopy=new double[n];
        resetHessian();
    }

    public int update(RealVector s, RealVector y1) {
        RealVector Hs = ldl.operate(s);
        double sHs = s.dotProduct(Hs);
        double sty = s.dotProduct(y1);

        if (sHs <= 0.0) {
            resetHessian();
            return 1;
        }
        boolean DAMPED=false;
        RealVector y = y1.copy();
        if (sty < GAMMA * sHs) {
            DAMPED=true;
            double phi = (1.0 - GAMMA) * sHs / (sHs - sty);
            y = y1.mapMultiply(phi).add(Hs.mapMultiply(1.0 - phi));
            sty = s.dotProduct(y);
        }

        if (!(sty > 0)) return 2;
        
        ldl.copyTo(lcopy,dcopy);
        // Upgrade
        ldl.update(y, 1.0 / sty);
        
        // Downdate
        if (!ldl.update(Hs, -1.0 / sHs)) {
            double gamma=1;
             if (!DAMPED) gamma =y.dotProduct(y) /sty;
//                gamma = sty /sHs;
                double th=1.0e-3;
                gamma = FastMath.max(th, FastMath.min(1.0/th, gamma));
                this.resetHessian(gamma);
                  
                return 3;
        }
       double[] d=ldl.getD();
       double mind=Double.POSITIVE_INFINITY;
       double maxd=0.0;
       for(int i=0;i<d.length;i++)
       {
           if(d[i]<mind) mind=d[i];
           if(d[i]>maxd) maxd=d[i];
       }
      
       double condRatio = mind / maxd;
boolean badCond = condRatio < 1.0e-12;
boolean badMin  = mind < 1.0e-12;

// CASO 4: Collasso Totale (Sia malcondizionata che singolare)
if ( badMin) {
    System.out.println("BFGS : condition < 1e-12 OR mind < 1e-12");
    skipCount+=1;
    if(skipCount<2){
    // La matrice è distrutta. Evitiamo divisioni per zero facendo uno skip.
    // Se il solver si blocca qui spesso, potrebbe essere necessario un reset totale.
    ldl.restoreFrom(lcopy, dcopy);
    return 4; 
    }
    else
    {
             double gamma=1;
             if (!DAMPED) gamma =y.dotProduct(y) /sty;
//                gamma = sty /sHs;
                double th=1.0e-3;
                gamma = FastMath.max(th, FastMath.min(1.0/th, gamma));
                this.resetHessian(gamma);
                skipCount=0;  
                return 5;
        
    }
}

        skipCount=0;
        return 0;
    }

    public void resetHessian() {
        ldl.factorizeNoPivot(initialH);
    }

    public void resetHessian(double gamma) {
        ldl.reset(gamma);
    }

    public RealMatrix getL() {
        
        double[][] lData = ldl.getLData();
        double[] dData = ldl.getD();
        int n=dData.length;
        double[][] lChol = new double[n][n];
        for (int j = 0; j < n; j++) {
            double sqrtDj = FastMath.sqrt(FastMath.max(0.0, dData[j]));
            lChol[j][j] = sqrtDj;
            for (int i = j + 1; i < n; i++) {
                lChol[i][j] = lData[i][j] * sqrtDj;
            }
        }
        return new Array2DRowRealMatrix(lChol, false);
    }

    public RealMatrix getHessian() {
        
        double[][] lData = ldl.getLData();
        double[] dData = ldl.getD();
        int n=dData.length;
        double[][] h = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int k = 0; k <= i; k++) {
                double sum = 0;
                for (int j = 0; j <= k; j++) {
                    double lij = (i == j) ? 1.0 : lData[i][j];
                    double lkj = (k == j) ? 1.0 : lData[k][j];
                    sum += lij * dData[j] * lkj;
                }
                h[i][k] = sum; h[k][i] = sum;
            }
        }
        return new Array2DRowRealMatrix(h, false);
    }
}