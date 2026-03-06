package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;

/**
 * IQPR-like refinement using the existing QPDualActiveSolver which natively supports:
 *   - LinearEqualityConstraint
 *   - LinearInequalityConstraint (Ax >= b)
 *
 * No active-set freezing: the oracle updates the active set internally and returns
 * a full lambda vector (eq then ineq, with 0 for inactive constraints).
 *
 * This implementation provides:
 *   - Δ with history: Δ_k <= alphaGrow * Δ_{k-1}
 *   - Backstepping on failure or non-improvement: Δ *= backstepFactor
 *   - Rollback: never returns worse KKT residual than the input
 */
public final class IQPRRefiner {

    private final QPDualActiveSolver oracle;
    private final double targetEpsilon;
    private final int maxIterations;

    // Scaling policy (paper/reference spirit)
    private final double alphaGrow;        // max growth per outer iteration
    private final int maxBacksteps;        // backstepping attempts
    private final double backstepFactor;   // scale reduction per backstep

    public IQPRRefiner(final QPDualActiveSolver oracle,
                       final double targetEpsilon,
                       final int maxIterations) {
        this(oracle, targetEpsilon, maxIterations,
             FastMath.pow(2, 24),  // very permissive growth cap
             10,
             1e-2);
    }

    public IQPRRefiner(final QPDualActiveSolver oracle,
                       final double targetEpsilon,
                       final int maxIterations,
                       final double alphaGrow,
                       final int maxBacksteps,
                       final double backstepFactor) {
        this.oracle = oracle;
        this.targetEpsilon = targetEpsilon;
        this.maxIterations = maxIterations;
        this.alphaGrow = alphaGrow;
        this.maxBacksteps = maxBacksteps;
        this.backstepFactor = backstepFactor;
    }

    public LagrangeSolution refine(final QuadraticFunction qf,
                                   final LinearInequalityConstraint iqc,
                                   final LinearEqualityConstraint eqc,
                                   final LagrangeSolution currentSol,
                                   final IsCholesky isChol) {

        RealVector bestX = currentSol.getX().copy();
        if (hasNaN(bestX)) {
            return currentSol;
        }

        final int nEq = (eqc != null) ? eqc.dimY() : 0;
        final int nIq = (iqc != null) ? iqc.dimY() : 0;

        RealVector bestY = normalizeLambda(currentSol.getLambda(), nEq, nIq);

        KKTResiduals bestRes = calculateKKTResiduals(qf, iqc, eqc, bestX, bestY, isChol);
        double bestErr = bestRes.combined();

        if (bestErr <= targetEpsilon) {
            return new LagrangeSolution(bestX, bestY, qf.value(bestX));
        }

        RealVector x = bestX.copy();
        RealVector y = (bestY != null) ? bestY.copy() : null;

        // Δ history
        double deltaPrev = 1.0;

        for (int k = 0; k < maxIterations; k++) {

            final KKTResiduals res = calculateKKTResiduals(qf, iqc, eqc, x, y, isChol);
            if (res.combined() <= targetEpsilon) {
                bestX = x; bestY = y;
                bestRes = res; bestErr = res.combined();
                break;
            }

            // Δ_k = min(1/deltaP, 1/deltaD, alphaGrow*Δ_{k-1})
            final double invP = (res.deltaP > 0.0) ? 1.0 / res.deltaP : 1e18;
            final double invD = (res.deltaD > 0.0) ? 1.0 / res.deltaD : 1e18;
            double deltaK = FastMath.min(FastMath.min(invP, invD), alphaGrow * deltaPrev);

            boolean accepted = false;
            double tryDelta = deltaK;

            for (int bt = 0; bt <= maxBacksteps && !accepted; bt++) {

                // Correction QP in variable d:
                // min 1/2 d^T H d + (Δ*rc)^T d
                final QuadraticFunction refinedQF =
                        new QuadraticFunction(qf.getP(), res.rc.mapMultiply(tryDelta), 0.0);

                // Constraints on d:
                // Equality: Aeq d = Δ*(beq - Aeq x) = Δ*rbEq
                final LinearEqualityConstraint corrEq =
                        (eqc != null) ? new LinearEqualityConstraint(eqc.getA(), res.rbEq.mapMultiply(tryDelta)) : null;

                // Inequality (Ax >= b): Aiq d >= Δ*(b - Aiq x) = Δ*rbIneq
                final LinearInequalityConstraint corrIq =
                        (iqc != null) ? new LinearInequalityConstraint(iqc.jacobian(null), res.rbIneq.mapMultiply(tryDelta)) : null;

                // Solve (oracle updates active set internally)
                final LagrangeSolution corr = oracle.optimize(
                        new ObjectiveFunction(refinedQF),
                        corrEq,
                        corrIq,
                        isChol
                );

                if (corr == null || corr.getX() == null || hasNaN(corr.getX())) {
                    tryDelta *= backstepFactor;
                    continue;
                }

                // Update primal
                final RealVector dx = corr.getX().mapDivide(tryDelta);
                final RealVector xCand = x.add(dx);

                // Update multipliers: y_{k+1} = y_k + (lambdaCorr)/Δ
                RealVector yCand = (y != null) ? y.copy() : normalizeLambda(null, nEq, nIq);
                if (corr.getLambda() != null) {
                    final RealVector dy = corr.getLambda().mapDivide(tryDelta);

                    // eq part
                    for (int i = 0; i < nEq && i < dy.getDimension(); i++) {
                        yCand.setEntry(i, yCand.getEntry(i) + dy.getEntry(i));
                    }
                    // ineq part (clamp >= 0)
                    for (int i = 0; i < nIq; i++) {
                        final int pos = nEq + i;
                        final int src = nEq + i;
                        if (src < dy.getDimension()) {
                            yCand.setEntry(pos, yCand.getEntry(pos) + dy.getEntry(src));
                        }
                        yCand.setEntry(pos, FastMath.max(0.0, yCand.getEntry(pos)));
                    }
                }

                // Check improvement
                final KKTResiduals candRes = calculateKKTResiduals(qf, iqc, eqc, xCand, yCand, isChol);
                final double candErr = candRes.combined();

                if (candErr < bestErr) {
                    bestX = xCand;
                    bestY = yCand;
                    bestRes = candRes;
                    bestErr = candErr;

                    x = xCand;
                    y = yCand;

                    deltaPrev = tryDelta;
                    accepted = true;

                    if (bestErr <= targetEpsilon) {
                        break;
                    }
                } else {
                    tryDelta *= backstepFactor;
                }
            }

            if (!accepted) {
                break;
            }
            System.out.println("QP REFIN:"+k+"X:"+x+";Y:"+y);
            
        }

        return new LagrangeSolution(bestX, bestY, qf.value(bestX));
    }

    /**
     * IMPORTANT: this matches your convention: lambda = [eq ; ineq].
     * Inactive constraints already have lambda=0 (your buildSolution does this),
     * but after refinement update we enforce lambda_ineq >= 0.
     */
    private RealVector normalizeLambda(final RealVector lambda, final int nEq, final int nIq) {
        final int dim = nEq + nIq;
        if (dim == 0) {
            return null;
        }
        final org.hipparchus.linear.ArrayRealVector y =
                (lambda == null || lambda.getDimension() != dim || hasNaN(lambda))
                ? new org.hipparchus.linear.ArrayRealVector(dim)
                : new org.hipparchus.linear.ArrayRealVector(lambda.toArray());

        for (int i = 0; i < nIq; i++) {
            final int pos = nEq + i;
            y.setEntry(pos, FastMath.max(0.0, y.getEntry(pos)));
        }
        return y;
    }

    private KKTResiduals calculateKKTResiduals(final QuadraticFunction qf,
                                               final LinearInequalityConstraint iqc,
                                               final LinearEqualityConstraint eqc,
                                               final RealVector x,
                                               final RealVector y,
                                               final IsCholesky isChol) {

        // Hx
        final RealVector hx;
        if (isChol != null && isChol.isCholesky()) {
            final RealMatrix L = qf.getP();
            hx = L.operate(L.transpose().operate(x));
        } else {
            hx = qf.getP().operate(x);
        }

        // stationarity residual: rc = Hx + q - Aeq^T*λeq - Aiq^T*λiq
        RealVector rc = hx.add(qf.getQ());

        final int nEq = (eqc != null) ? eqc.dimY() : 0;
        final int nIq = (iqc != null) ? iqc.dimY() : 0;

        if (y != null) {
            if (nEq > 0) {
                rc = rc.subtract(eqc.getA().transpose().operate(y.getSubVector(0, nEq)));
            }
            if (nIq > 0) {
                rc = rc.subtract(iqc.jacobian(null).transpose().operate(y.getSubVector(nEq, nIq)));
            }
        }

        // primal residuals:
        // eq: rbEq = b - Aeq x
        // ineq Ax >= b: rbIneq = b - Aiq x (positive => violation)
        final RealVector rbEq = (nEq > 0) ? eqc.getLowerBound().subtract(eqc.getA().operate(x)) : null;
        final RealVector rbIneq = (nIq > 0) ? iqc.getLowerBound().subtract(iqc.jacobian(null).operate(x)) : null;

        return new KKTResiduals(rc, rbEq, rbIneq, y, nEq);
    }

    private boolean hasNaN(final RealVector v) {
        if (v == null) return true;
        for (final double d : v.toArray()) {
            if (Double.isNaN(d) || Double.isInfinite(d)) return true;
        }
        return false;
    }

    private static final class KKTResiduals {
        final RealVector rc, rbEq, rbIneq;
        final double deltaP, deltaD, deltaS;

        KKTResiduals(final RealVector rc,
                     final RealVector rbEq,
                     final RealVector rbIneq,
                     final RealVector y,
                     final int ineqOffset) {
            this.rc = rc;
            this.rbEq = rbEq;
            this.rbIneq = rbIneq;

            // primal violation
            double dp = 0.0;
            if (rbEq != null) dp = FastMath.max(dp, rbEq.getLInfNorm());
            if (rbIneq != null) {
                for (double v : rbIneq.toArray()) dp = FastMath.max(dp, FastMath.max(0.0, v));
            }
            this.deltaP = dp;

            // stationarity
            this.deltaD = (rc != null) ? rc.getLInfNorm() : Double.POSITIVE_INFINITY;

            // complementarity (max, not sum)
            double ds = 0.0;
            if (rbIneq != null && y != null && y.getDimension() >= ineqOffset + rbIneq.getDimension()) {
                for (int i = 0; i < rbIneq.getDimension(); i++) {
                    ds = FastMath.max(ds, FastMath.abs(rbIneq.getEntry(i) * y.getEntry(ineqOffset + i)));
                }
            }
            this.deltaS = ds;
        }

        double combined() {
            return FastMath.max(FastMath.max(deltaP, deltaD), deltaS);
        }
    }
}