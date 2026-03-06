package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.*;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Robust LDL-based BFGS updater for SQP methods.
 *
 * H = L D Lᵀ
 * L = lower unit triangular (diag = 1)
 * D = diagonal pivot vector
 *
 * Implements:
 *  • Powell damping
 *  • Kraft-style rank-two BFGS update   H ← H + (y yᵀ)/(sᵀy)  − (Bs Bsᵀ)/(sᵀBs)
 *  • LDL rank-1 update/downdate (LDL001, Schittkowski-style, FO = 4)
 *  • Optional global scaling reset with gamma = (yᵀy)/(sᵀy)
 *  • Optional spectral regularization (H ← H + τI) via regularizeLDL()
 *
 * This matches the spirit of NLPQL / ROBUST.SQP.
 */
public class BFGSUpdaterLDL {

    /** Initial Hessian (or its approximation) */
    private final RealMatrix initialH;

    /** lower unit-triangular factor (L) */
    private RealMatrix L;

    /** diagonal pivot vector (D) */
    private RealVector d;

    /** dimension */
    private final int n;

    /** Powell damping parameter */
    private static final double GAMMA = 0.2;

    /** target condition number (for regularization) */
    private static final double COND_TARGET = 1.0e8;

    /** initial regularization scale */
    private static final double TAU_INIT = 1e-12;

    private final double sqrtEPS = FastMath.sqrt(Precision.EPSILON);
    private final double decompositionEpsilon;
    private final double EPS;
    private final boolean SCALE;
    private boolean skipped=false;
    private double gammaOld=1.0;

    public BFGSUpdaterLDL(final RealMatrix Hini,
                          final double eps,
                          final boolean autoScale,
                          final double decompositionEpsilon) {

        this.initialH = Hini.copy();
        this.n = Hini.getRowDimension();
        this.decompositionEpsilon = decompositionEpsilon;

        L = MatrixUtils.createRealIdentityMatrix(n);
        d = new ArrayRealVector(n);
        EPS = eps;
        SCALE = autoScale;
        resetHessian();
    }

    // ==============================================================
    // PUBLIC API
    // ==============================================================

    /**
     * Reconstructs the full Hessian H = L D Lᵀ from LDL factors.
     */
    public RealMatrix getHessian() {

        // H = L diag(d) Lᵀ
        // costruiamo LD = L con ogni colonna j moltiplicata per d_j
        RealMatrix LD = L.copy();

        for (int j = 0; j < n; j++) {
            double dj = d.getEntry(j);
            for (int i = j; i < n; i++) {
                LD.setEntry(i, j, LD.getEntry(i, j) * FastMath.sqrt(dj));
            }
        }

        // H = LD * Lᵀ
        return LD.multiplyTransposed(LD);
        // Se multiplyTransposed non esistesse, usare:
        // return LD.multiply(L.transpose());
    }

    /**
     * Reset Hessian to initialH (factorized as LDL).
     */
    public void resetHessian() {
        factorLDL(initialH, L, d);
    }

    /**
     * Reset Hessian to gamma * initialH (factorized as LDL).
     */
    public void resetHessian(double gamma) {
        RealMatrix scaled = initialH.scalarMultiply(gamma);
        factorLDL(scaled, L, d);
    }

    /**
     * Main BFGS update.
     *
     * @param s  search step (x_{k+1} - x_k)
     * @param y1 gradient difference (∇L_{k+1} - ∇L_k) (raw)
     *
     * @return
     *  0 = ok
     *  1 = curvature fail (sᵀBs ≤ 0)
     *  3 = gamma scaling reset triggered (SCALE logic)
     *  4 = downdate failed → Hessian reset with gamma
     */
    public int update(RealVector s, RealVector y1) {

        // Hs = B s = L D Lᵀ s
        RealVector Bs = applyHessian(s);
        double sBs = s.dotProduct(Bs);

        // Basic curvature check on Hs
        if (sBs <= 0.0) {
            return 1;
        }

        // sᵀy and damped y (Powell)
        double sty = s.dotProduct(y1);
        RealVector y = y1.copy();
        

        // Powell damping: enforce sᵀy ≥ γ sᵀBs
        if (sty < GAMMA * sBs) {
            System.out.println("BEDORE DAMP"+sty);
            double phi = (sBs - GAMMA * sBs) / (sBs - sty); // = (1-γ) sBs / (sBs - sᵀy)
            
            
            
            // y* = φ y + (1-φ) Bs
            y = y1.mapMultiply(phi).add(Bs.mapMultiply(1.0 - phi));
            sty = GAMMA * sBs;  // per costruzione
//            if (sty < GAMMA * sBs) return 2;
             System.out.println(phi+"AFTER DAMP"+sty);
            
        }   
            double yy = y.dotProduct(y);
            double gamma = yy / sty;
            
//            System.out.println((ss/sBs)+";"+(yy/sty));
//                if(yy/sty>1.0e8){
////             double gamma = yy / sty;
////
////           
////                gamma = FastMath.max(this.sqrtEPS, FastMath.min(1.0 / this.sqrtEPS, gamma));
//                this.resetHessian(1.0e8);
////             
//         }

        // Optional auto scaling reset (Schittkowski-like)
        if (SCALE) {

           
            

            if (gamma < FastMath.sqrt(EPS)) {
                gamma = FastMath.max(this.sqrtEPS, FastMath.min(1.0 / this.sqrtEPS, gamma));
                this.resetHessian(gamma);
                return 3;
            }
        }
        
        // ==========================================================
        // Rank-two BFGS update in LDL:
        //
        // H ← H + (y yᵀ)/(sᵀy) − (Bs Bsᵀ)/(sᵀBs)
        //
        // cioè due rank-1:
        //  • update con u = y / sqrt(sty)
        //  • downdate con v = Bs / sqrt(sBs)
        // ==========================================================

        double invSqrtSty = 1.0 / FastMath.sqrt(sty);
        double invSqrtSBS = 1.0 / FastMath.sqrt(sBs);

        RealVector u = y.mapMultiply(invSqrtSty);   // u uᵀ = y yᵀ / sty
        RealVector v = Bs.mapMultiply(invSqrtSBS);  // v vᵀ = Bs Bsᵀ / sBs
        RealMatrix Lcopy = L.copy();
        RealVector Dcopy = d.copy();
        // 1) LDL update: H ← H + u uᵀ
        double minLambda=d.getMinValue();
        double cond = estimateCondLDL();
        ldlUpdate(u);
        boolean result=ldlDowndate(v );
        
        // 2) LDL downdate: H ← H − v vᵀ
//        if (!result  ) {
//
//            // robust fallback: reset with gamma scaling
//            double yy = y.dotProduct(y);
//            double gamma = yy / sty;
//
//            gamma = FastMath.max(this.sqrtEPS, FastMath.min(1.0 / this.sqrtEPS, gamma));
//            this.resetHessian(gamma);
//           
//            return 4;
//        }
//        if(!result)
//        {
//           L.setSubMatrix(Lcopy.getData(), 0 ,0);
//            d.setSubVector(0, d);  
//        }
//         if(estimateCondLDL()*Precision.EPSILON*100.0>EPS)
//         {
//            localDiagonalRegularization();
//         }
//           if(estimateCondLDL()*Precision.EPSILON*100.0>EPS)
//           {
//            double yy = y.dotProduct(y);
//            double gamma = yy / sty;
//
//            gamma = FastMath.max(this.sqrtEPS, FastMath.min(1.0 / this.sqrtEPS, gamma));
//            this.resetHessian(gamma);
//           }
//           localDiagonalRegularization();  
        // Se vuoi, qui potresti eventualmente fare:
        if (d.getMinValue()<this.sqrtEPS ) { 
           

            gamma = FastMath.max(this.sqrtEPS, FastMath.min(1.0 / this.sqrtEPS, gamma));
            this.resetHessian(gammaOld);
            gammaOld=1;
            skipped=false;
            return 5;
        }
            else if (d.getMinValue()<this.sqrtEPS&& !skipped)
                    {
                    L.setSubMatrix(Lcopy.getData(), 0, 0);
                    d.setSubVector(0, Dcopy);
                    skipped=true;
                    return 4;
                    
                    }
        
        gammaOld=gamma;
        return 0;
    }

    // ==============================================================
    // LDL Mechanics
    // ==============================================================

    /** LDL factorization: H = L D Lᵀ, L unit-lower, D diagonal. */
    private void factorLDL(final RealMatrix H, final RealMatrix L, final RealVector d) {

        double[][] A = H.getData();
        int n = H.getRowDimension();

        // reset L=I
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                L.setEntry(i, j, 0.0);
            }
            L.setEntry(i, i, 1.0);
        }

        double[] D = new double[n];

        // classic LDL without pivoting
        for (int j = 0; j < n; j++) {
            double djj = A[j][j];

            for (int k = 0; k < j; k++) {
                double Ljk = L.getEntry(j, k);
                djj -= Ljk * Ljk * D[k];
            }
            D[j] = djj;

            for (int i = j + 1; i < n; i++) {
                double Lij = A[i][j];
                for (int k = 0; k < j; k++) {
                    Lij -= L.getEntry(i, k) * L.getEntry(j, k) * D[k];
                }
                L.setEntry(i, j, Lij / D[j]);
            }
        }

        for (int i = 0; i < n; i++) {
            d.setEntry(i, D[i]);
        }
    }

    /** Bs = L D Lᵀ s */
    private RealVector applyHessian(RealVector x) {
        // z = Lᵀ x
        double[] z = L.preMultiply(x).toArray();
        // z = D z
        for (int i = 0; i < n; i++) {
            z[i] *= d.getEntry(i);
        }
        // return L z
        return L.operate(new ArrayRealVector(z, false));
    }

    // ==============================================================
    // LDL rank-1 update & downdate (Schittkowski / LDL001)
    // ==============================================================

    /**
     * Rank-one LDL update: H ← H + z zᵀ
     * Implementa la parte "update" di LDL001, con soglia FO = 4
     * per distinguere crescita moderata/forte del pivot.
     *
     * @param zIn vettore di aggiornamento (dimensione n)
     */
    private void ldlUpdate(final RealVector zIn) {

        final int n = this.n;
        final double FO = 4.0;

        // U = z (verrà modificato in posto)
        final double[] U = zIn.toArray().clone();

        double T = 1.0;

        for (int j = 0; j < n; ++j) {

            double PJ = U[j];
            double DJ = d.getEntry(j);

            double TNEW = T + (PJ * PJ) / DJ;
            double TT   = TNEW / T;
            double TTI  = 1.0 / TT;
            double DJNEW = DJ * TT;
            double BETA  = PJ / (DJ * TNEW);

            int j1 = j + 1;
            if (j1 < n) {

                if (DJNEW <= FO * DJ) {
                    // Branch "moderate growth"
                    for (int r = j1; r < n; ++r) {
                        // U(r) = U(r) - PJ * L(r,j)
                        double Ur = U[r] - PJ * L.getEntry(r, j);
                        U[r] = Ur;
                        // L(r,j) = L(r,j) + BETA * U(r)
                        double CLrj = L.getEntry(r, j) + BETA * Ur;
                        L.setEntry(r, j, CLrj);
                    }
                } else {
                    // Branch "large growth" (stabilized)
                    for (int r = j1; r < n; ++r) {
                        double CLRJ = L.getEntry(r, j);
                        // L(r,j) = TTI * CLRJ + BETA * U(r)
                        double CLrj = TTI * CLRJ + BETA * U[r];
                        L.setEntry(r, j, CLrj);
                        // U(r) = U(r) - PJ * CLRJ
                        U[r] = U[r] - PJ * CLRJ;
                    }
                }
            }

            T = TNEW;
            d.setEntry(j, DJNEW);
        }
    }

    /**
     * Rank-one LDL downdate: H ← H − z zᵀ
     * Implementa la parte "downdate" di LDL001, robusta:
     *
     *  • calcola U = L^{-1} z (forward)
     *  • PDP = Σ u_i² / D_i
     *  • T = 1 − PDP (clamp a decompositionEpsilon se troppo piccolo)
     *  • sweep all'indietro con TNEW, TT, BETA, aggiornando L e D.
     *
     * @param zIn vettore di downdate (dimensione n)
     * @return true se i pivot restano positivi, false se SPD distrutta
     */
    private boolean ldlDowndate(final RealVector zIn) {

        final int n = this.n;
        final double eps = decompositionEpsilon;

        final double[] V = zIn.toArray().clone();
        final double[] U = new double[n];

        // U(0) = V(0); PDP = U(0)^2 / D(0)
        U[0] = V[0];
        double PDP = (U[0] * U[0]) / d.getEntry(0);

        // Forward: U(i) = V(i) - Σ_{j<i} L(i,j)*U(j), PDP += U(i)^2 / D(i)
        for (int i = 1; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += L.getEntry(i, j) * U[j];
            }
            double UI = V[i] - sum;
            U[i] = UI;
            PDP += (UI * UI) / d.getEntry(i);
        }

        // T = 1 - PDP
        double T = 1.0 - PDP;
        
        if (T <= this.decompositionEpsilon) {
            
           T = this.decompositionEpsilon;   // regolarizzazione minima (RobustSQP-style)
        }

        // Ciclo all'indietro
        for (int idx = 0; idx < n; ++idx) {
            int j = n - 1 - idx;

            double PJ = U[j];
            double DJ = d.getEntry(j);

            double TNEW = T + (PJ * PJ) / DJ;
            double TT   = T / TNEW;
            double DJNEW = DJ * TT;
            double BETA  = -PJ / (DJ * T);

            V[j] = PJ;
            T = TNEW;

            
                d.setEntry(j, DJNEW);
            

            int j1 = j + 1;
            if (j1 < n) {
                for (int r = j1; r < n; ++r) {
                    double CLRJ = L.getEntry(r, j);
                    // L(r,j) = CLRJ + BETA * V(r)
                    double CLrj = CLRJ + BETA * V[r];
                    L.setEntry(r, j, CLrj);
                    // V(r) = V(r) + PJ * CLRJ
                    V[r] = V[r] + PJ * CLRJ;
                }
            }
        }

        return true;
    }

    // ==============================================================
    // Check & Regularization
    // ==============================================================

    private double estimateCond() {
        double min = Double.POSITIVE_INFINITY, max = 0.0;
        for (int i = 0; i < n; i++) {
            double pi = FastMath.abs(d.getEntry(i));
            min = FastMath.min(min, pi);
            max = FastMath.max(max, pi);
        }
        return max / min;
    }

    private void regularizeLDL() {

        double cond = estimateCond();
        if (cond < COND_TARGET) {
            return;
        }

        // Partiamo dall'H attuale
        RealMatrix H = getHessian();
        double tau = TAU_INIT * d.getMaxValue();
        int iter = 0;

        while (iter < 3) {

            RealMatrix Hp = H.copy();
            for (int i = 0; i < n; i++) {
                Hp.addToEntry(i, i, tau);
            }

            try {
                factorLDL(Hp, L, d);
                if (estimateCond() < COND_TARGET) {
                    return;
                }
            } catch (Exception ignore) {
                // se fallisce, aumentiamo tau
            }

            tau *= 10.0;
            iter++;
        }
    }

    /**
     * Solves L z = b for z, with L unit-lower.
     */
    private RealVector solveLower(final RealVector b) {
        final int n = this.n;
        double[] z = b.toArray();
        for (int i = 0; i < n; i++) {
            double sum = z[i];
            for (int j = 0; j < i; j++) {
                sum -= L.getEntry(i, j) * z[j];
            }
            z[i] = sum;
        }
        return new ArrayRealVector(z, false);
    }
    
    /** Stima deterministica del condizionamento di H = L D Lᵀ. */
private double estimateCondLDL() {

    // spread dei pivot D
    double minD = Double.POSITIVE_INFINITY;
    double maxD = 0.0;
    for (int i = 0; i < n; ++i) {
        double di = FastMath.abs(d.getEntry(i));
        minD = FastMath.min(minD, di);
        maxD = FastMath.max(maxD, di);
    }
    double condDiag = maxD / minD;

    // spread delle norme di riga di L (L unit-lower)
    double minRow = Double.POSITIVE_INFINITY;
    double maxRow = 0.0;
    for (int i = 0; i < n; ++i) {
        double rowNorm2 = 1.0; // L(i,i)=1
        for (int j = 0; j < i; ++j) {
            double lij = L.getEntry(i, j);
            rowNorm2 += lij * lij;
        }
        minRow = FastMath.min(minRow, rowNorm2);
        maxRow = FastMath.max(maxRow, rowNorm2);
    }
    double spreadL = maxRow / minRow;

    return condDiag * spreadL;
}

    
    /** Regolarizzazione locale di D per limitare il condizionamento. */
private void localDiagonalRegularization() {

    double minD = Double.POSITIVE_INFINITY;
    double maxD = 0.0;
    for (int i = 0; i < n; ++i) {
        double di = FastMath.abs(d.getEntry(i));
        minD = FastMath.min(minD, di);
        maxD = FastMath.max(maxD, di);
    }

    // spread L
    double minRow = Double.POSITIVE_INFINITY;
    double maxRow = 0.0;
    for (int i = 0; i < n; ++i) {
        double rowNorm2 = 1.0;
        for (int j = 0; j < i; ++j) {
            double lij = L.getEntry(i, j);
            rowNorm2 += lij * lij;
        }
        minRow = FastMath.min(minRow, rowNorm2);
        maxRow = FastMath.max(maxRow, rowNorm2);
    }
    double spreadL = maxRow / minRow;

    double condApprox = (maxD / minD) * spreadL;
    if (condApprox <= COND_TARGET) {
        return; // già ok
    }

    // vogliamo (maxD / dFloor) * spreadL ≈ COND_TARGET  ⇒ dFloor ≈ maxD * spreadL / COND_TARGET
    double diagTargetCond = COND_TARGET / spreadL;
    if (diagTargetCond < 1.0) {
        diagTargetCond = 1.0;
    }
    double dFloor = maxD / diagTargetCond;

    // "floor" dei pivot troppo piccoli
    for (int i = 0; i < n; ++i) {
        double di = d.getEntry(i);
        if (FastMath.abs(di) < dFloor) {
            d.setEntry(i, FastMath.copySign(dFloor, di));
        }
    }
}


}
