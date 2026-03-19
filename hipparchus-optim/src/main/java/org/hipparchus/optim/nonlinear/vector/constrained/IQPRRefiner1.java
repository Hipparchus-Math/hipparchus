package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;

/**
 * IQPR-like refinement for QP problems of the form:
 *
 *   minimize   1/2 x' H x + q' x
 *   subject to Aeq x = beq
 *              Aiq x >= biq
 *
 * Correction QP in variable d:
 *
 *   minimize   1/2 d' H d + (Δ rD)' d
 *   subject to Aeq d = Δ rEq
 *              Aiq d >= Δ rIq
 *
 * with:
 *
 *   rEq = beq - Aeq x
 *   rIq = biq - Aiq x      (positive => violated inequality)
 *   rD  = Hx + q - Aeq' λeq - Aiq' λiq
 *
 * Update:
 *
 *   x_{k+1} = x_k + d / Δ
 *   y_{k+1} = y_k + λcorr / Δ
 *
 * Acceptance policy:
 *   - Iterate until residuals are below targetEpsilon or maxIterations is reached.
 *   - At the end, accept the refined solution only if its combined residual is
 *     smaller than the initial combined residual; otherwise return the initial
 *     solution.
 */
public final class IQPRRefiner1 {

    private final QPDualActiveSolverR oracle;
    private final double targetEpsilon;
    private final int maxIterations;
    private final double alphaGrow;
    private final int maxBacksteps;
    private final double backstepFactor;
    private final boolean verbose;

    public IQPRRefiner1(final QPDualActiveSolverR oracle,
                        final double targetEpsilon,
                        final int maxIterations) {
        this(oracle, targetEpsilon, maxIterations, FastMath.pow(2.0, 24.0), 8, 1.0e-2, true);
    }

    public IQPRRefiner1(final QPDualActiveSolverR oracle,
                        final double targetEpsilon,
                        final int maxIterations,
                        final double alphaGrow,
                        final int maxBacksteps,
                        final double backstepFactor,
                        final boolean verbose) {
        this.oracle = oracle;
        this.targetEpsilon = targetEpsilon;
        this.maxIterations = maxIterations;
        this.alphaGrow = alphaGrow;
        this.maxBacksteps = maxBacksteps;
        this.backstepFactor = backstepFactor;
        this.verbose = verbose;
    }

    public LagrangeSolution refine(final QuadraticFunction qf,
                                   final LinearInequalityConstraint iqc,
                                   final LinearEqualityConstraint eqc,
                                   final LagrangeSolution currentSol,
                                   final IsCholesky isChol) {

        final RealVector initialX = currentSol.getX().copy();
        if (hasNaN(initialX)) {
            return currentSol;
        }

        final int nEq = (eqc != null) ? eqc.dimY() : 0;
        final int nIq = (iqc != null) ? iqc.dimY() : 0;

        final RealVector initialY = normalizeLambda(currentSol.getLambda(), nEq, nIq);
        final KKTResiduals initialRes = calculateKKTResiduals(qf, iqc, eqc, initialX, initialY, isChol);
        final double initialErr = initialRes.combined();

        RealVector x = initialX.copy();
        RealVector y = (initialY != null) ? initialY.copy() : null;
        KKTResiduals currentRes = initialRes;
        double deltaPrev = 1.0;

        for (int k = 0; k < maxIterations; k++) {

            logIteration(k, x, y, currentRes);

            if (currentRes.deltaP <= targetEpsilon &&
                currentRes.deltaD <= targetEpsilon &&
                currentRes.deltaS <= targetEpsilon) {
                break;
            }

            final double invP = (currentRes.deltaP > 0.0) ? 1.0 / currentRes.deltaP : 1.0e18;
            final double invD = (currentRes.deltaD > 0.0) ? 1.0 / currentRes.deltaD : 1.0e18;
            final double deltaBase = FastMath.min(FastMath.min(invP, invD), alphaGrow * deltaPrev);

            boolean stepDone = false;

            for (int bt = 0; bt <= maxBacksteps; bt++) {

                final double deltaK = deltaBase * FastMath.pow(backstepFactor, bt);

                final QuadraticFunction refinedQF =
                        new QuadraticFunction(qf.getP(), currentRes.rD.mapMultiply(deltaK), 0.0);

                final LinearEqualityConstraint corrEq =
                        (eqc != null)
                                ? new LinearEqualityConstraint(eqc.getA(), currentRes.rEq.mapMultiply(deltaK))
                                : null;

                final LinearInequalityConstraint corrIq =
                        (iqc != null)
                                ? new LinearInequalityConstraint(iqc.jacobian(null), currentRes.rIq.mapMultiply(deltaK))
                                : null;

                final LagrangeSolution corr = oracle.optimize(
                        new ObjectiveFunction(refinedQF),
                        corrEq,
                        corrIq,
                        isChol
                );

                if (corr == null || corr.getX() == null || hasNaN(corr.getX())) {
                    continue;
                }

                x = x.add(corr.getX().mapDivide(deltaK));

                if (y == null && corr.getLambda() != null) {
                    y = new ArrayRealVector(corr.getLambda().getDimension());
                }

                if (y != null && corr.getLambda() != null) {
                    final RealVector dy = corr.getLambda().mapDivide(deltaK);
                    final int dim = FastMath.min(y.getDimension(), dy.getDimension());
                    for (int i = 0; i < dim; i++) {
                        y.addToEntry(i, dy.getEntry(i));
                    }
                }

                currentRes = calculateKKTResiduals(qf, iqc, eqc, x, y, isChol);
                deltaPrev = deltaK;
                stepDone = true;
                break;
            }

            if (!stepDone) {
                break;
            }
        }

        final KKTResiduals finalRes = calculateKKTResiduals(qf, iqc, eqc, x, y, isChol);
        final double finalErr = finalRes.combined();

        if (finalErr < initialErr) {
            logDecision(true);
            return new LagrangeSolution(x, y, qf.value(x));
        } else {
            logDecision(false);
            return new LagrangeSolution(initialX, initialY, qf.value(initialX));
        }
    }

    private RealVector normalizeLambda(final RealVector lambda, final int nEq, final int nIq) {
        final int dim = nEq + nIq;
        if (dim == 0) {
            return null;
        }
        if (lambda == null || lambda.getDimension() != dim || hasNaN(lambda)) {
            return new ArrayRealVector(dim);
        }
        return new ArrayRealVector(lambda.toArray());
    }

    private KKTResiduals calculateKKTResiduals(final QuadraticFunction qf,
                                               final LinearInequalityConstraint iqc,
                                               final LinearEqualityConstraint eqc,
                                               final RealVector x,
                                               final RealVector y,
                                               final IsCholesky isChol) {

        final RealVector hx;
        if (isChol != null && isChol.isCholesky()) {
            final RealMatrix L = qf.getP();
            hx = L.operate(L.transpose().operate(x));
        } else {
            hx = qf.getP().operate(x);
        }

        final int nEq = (eqc != null) ? eqc.dimY() : 0;
        final int nIq = (iqc != null) ? iqc.dimY() : 0;

        // rD = Hx + q - Aeq' λeq - Aiq' λiq
        RealVector rD = hx.add(qf.getQ());
        if (y != null) {
            if (nEq > 0) {
                rD = rD.subtract(eqc.getA().transpose().operate(y.getSubVector(0, nEq)));
            }
            if (nIq > 0) {
                rD = rD.subtract(iqc.jacobian(null).transpose().operate(y.getSubVector(nEq, nIq)));
            }
        }

        // rEq = beq - Aeq x
        final RealVector rEq =
                (nEq > 0)
                        ? eqc.getLowerBound().subtract(eqc.getA().operate(x))
                        : null;

        // rIq = biq - Aiq x   (positive => violated inequality)
        final RealVector rIq =
                (nIq > 0)
                        ? iqc.getLowerBound().subtract(iqc.jacobian(null).operate(x))
                        : null;

        return new KKTResiduals(rD, rEq, rIq, y, nEq);
    }

    private boolean hasNaN(final RealVector v) {
        if (v == null) {
            return true;
        }
        for (double d : v.toArray()) {
            if (Double.isNaN(d) || Double.isInfinite(d)) {
                return true;
            }
        }
        return false;
    }

    private void logIteration(final int k,
                              final RealVector x,
                              final RealVector y,
                              final KKTResiduals res) {
        if (!verbose) {
            return;
        }
        System.out.println("iter=" + k + " x=" + x + " y=" + y);
        System.out.println("  dP=" + res.deltaP +
                           " dD=" + res.deltaD +
                           " dS=" + res.deltaS +
                           " comb=" + res.combined());
    }

    private void logDecision(final boolean accepted) {
        if (!verbose) {
            return;
        }
        System.out.println(accepted ? "  ACCEPTED" : "  ROLLBACK");
    }

    private static final class KKTResiduals {
        final RealVector rD;
        final RealVector rEq;
        final RealVector rIq;
        final double deltaP;
        final double deltaD;
        final double deltaS;

        KKTResiduals(final RealVector rD,
                     final RealVector rEq,
                     final RealVector rIq,
                     final RealVector y,
                     final int ineqOffset) {
            this.rD = rD;
            this.rEq = rEq;
            this.rIq = rIq;

            // deltaP = max( ||rEq||∞ , max_i max(0, rIq_i) )
            double dp = 0.0;
            if (rEq != null) {
                dp = FastMath.max(dp, rEq.getLInfNorm());
            }
            if (rIq != null) {
                for (double v : rIq.toArray()) {
                    dp = FastMath.max(dp, FastMath.max(0.0, v));
                }
            }
            this.deltaP = dp;

            // deltaD = max( ||rD||∞ , max_i max(0, -yIq_i) )
            double dd = (rD != null) ? rD.getLInfNorm() : Double.POSITIVE_INFINITY;
            if (y != null) {
                for (int i = ineqOffset; i < y.getDimension(); i++) {
                    dd = FastMath.max(dd, FastMath.max(0.0, -y.getEntry(i)));
                }
            }
            this.deltaD = dd;

            // slack s = Aiq x - biq = -rIq
            // deltaS = max_i | yIq_i * s_i |
            double ds = 0.0;
            if (rIq != null && y != null && y.getDimension() >= ineqOffset + rIq.getDimension()) {
                for (int i = 0; i < rIq.getDimension(); i++) {
                    final double slack = -rIq.getEntry(i);
                    final double yi = y.getEntry(ineqOffset + i);
                    ds = FastMath.max(ds, FastMath.abs(yi * slack));
                }
            }
            this.deltaS = ds;
        }

        double combined() {
            return FastMath.max(deltaP, FastMath.max(deltaD, deltaS));
        }
    }
}