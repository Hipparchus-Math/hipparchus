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

import java.util.ArrayList;
import java.util.Collections;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.EigenDecompositionSymmetric;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;

/**
 * Sequential Quadratic Programming Optimizer.
 * <p>
 * min f(x)
 * <p>
 * q(x) = b1
 * <p>
 * h(x)&gt;=b2
 * <p>
 * Algorithm based on paper:"On the convergence of a sequential quadratic
 * programming method(Klaus Shittkowki,January 1982)"
 * @since 3.1
 * @deprecated as of 4.1, replaced by {@link SQPOptimizerS2}
 */
@Deprecated
public class SQPOptimizerS extends AbstractSQPOptimizer {

    /** Forgetting factor. */
    private static final int FORGETTING_FACTOR = 10;

    /** Jacobian constraint. */
    private RealMatrix constraintJacob;

    /** Value of the equality constraint. */
    private RealVector equalityEval;

    /** Value of the inequality constraint. */
    private RealVector inequalityEval;

    /** Gradient of the objective function. */
    private RealVector functionGradient;

    /** Hessian of the objective function. */
    private RealMatrix hessian;

    /** {@inheritDoc} */
    @Override
    public LagrangeSolution doOptimize() {
        int me = 0;
        int mi = 0;

        //EQUALITY CONSTRAINT
        if (this.getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        //INEQUALITY CONSTRAINT
        if (this.getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
        }

        double alpha;
        double rho = 100.0;
        int failedSearch = 0;
        RealVector x;
        if (this.getStartPoint() != null) {
            x = new ArrayRealVector(this.getStartPoint());
        } else {
            x = new ArrayRealVector(this.getObj().dim());
        }

        RealVector y = new ArrayRealVector(me + mi, 0.0);
        RealVector r = new ArrayRealVector(me + mi, 1.0);
        ArrayList<Double> oldPenalty = new ArrayList<>();
        //INITIAL VALUES
        double functionEval = this.getObj().value(x);
        functionGradient = this.getObj().gradient(x);
        double maxGrad = functionGradient.getLInfNorm();


        if (this.getEqConstraint() != null) {

            equalityEval = this.getEqConstraint().value(x);
        }
        if (this.getIqConstraint() != null) {

            inequalityEval = this.getIqConstraint().value(x);
        }
        constraintJacob = computeJacobianConstraint(x);

        if (this.getEqConstraint() != null) {
            maxGrad = FastMath.max(maxGrad, equalityEval.getLInfNorm());
        }
        if (this.getIqConstraint() != null) {
            maxGrad = FastMath.max(maxGrad, inequalityEval.getLInfNorm());
        }

        if (!getSettings().useFunHessian()) {
            hessian = MatrixUtils.createRealIdentityMatrix(x.getDimension());
        } else {
            hessian = this.getObj().hessian(x);
        }

        for (int i = 0; i < this.getMaxIterations(); i++) {
            iterations.increment();

            alpha = 1.0;
            LagrangeSolution sol1 = null;

            int qpLoop = 0;
            double sigma = maxGrad;
            double currentPenaltyGrad;

            //LOOP TO FIND SOLUTION WITH SIGMA<SIGMA THRESHOLD
            while (sigma > getSettings().getSigmaMax() && qpLoop < getSettings().getQpMaxLoop()) {
                sol1 = solveQP(x, y, rho);
                sigma = sol1.getValue();

                if (sigma > getSettings().getSigmaMax()) {
                    rho = getSettings().getRhoCons() * rho;
                    qpLoop += 1;

                }

            }
            //IF SIGMA>SIGMA THRESHOLD ASSIGN DIRECTION FROM PENALTY GRADIENT
            final RealVector dx;
            final RealVector dy;
            if (qpLoop == getSettings().getQpMaxLoop()) {

                dx = (MatrixUtils.inverse(hessian).operate(penaltyFunctionGradX(functionGradient, x, y, r))).mapMultiply(-1.0);
                dy = y.subtract(penaltyFunctionGradY(x, y, r));
                sigma = 0;

            } else {
                dx = sol1.getX();
                sigma = sol1.getValue();
                if (sigma < 0) {
                    sigma = 0;
                }
                dy = sol1.getLambda();
                r = updateRj(hessian, y, dx, dy, r, sigma);
            }

            currentPenaltyGrad = penaltyFunctionGrad(functionGradient, dx, dy, x, y, r);
            int search = 0;

            rho = updateRho(dx, dy, hessian, constraintJacob, sigma);

            double currentPenalty = penaltyFunction(functionEval, 0, x, y, dx, dy.subtract(y), r);

            double alphaF = this.getObj().value(x.add(dx.mapMultiply(alpha)));
            double alfaPenalty = penaltyFunction(alphaF, alpha, x, y, dx, dy.subtract(y), r);


            //LINE SEARCH

            while ((alfaPenalty - currentPenalty) >= getSettings().getMu() * alpha * currentPenaltyGrad &&
                    search < getSettings().getMaxLineSearchIteration()) {
                double alfaStar = -0.5 * alpha * alpha * currentPenaltyGrad / (-alpha * currentPenaltyGrad + alfaPenalty - currentPenalty);


                alpha = FastMath.max(getSettings().getB() * alpha, FastMath.min(1.0, alfaStar));
                alphaF = this.getObj().value(x.add(dx.mapMultiply(alpha)));
                alfaPenalty = penaltyFunction(alphaF, alpha, x, y, dx, dy.subtract(y), r);
                search = search + 1;

            }



            if (getSettings().getConvCriteria() == 0) {
                if (dx.mapMultiply(alpha).dotProduct(hessian.operate(dx.mapMultiply(alpha))) < getSettings().getEps() * getSettings().getEps()) {
//                    x = x.add(dx.mapMultiply(alfa));
//                    y = y.add((dy.subtract(y)).mapMultiply(alfa));
                    break;
                }
            } else {
                if (alpha * dx.getNorm() < getSettings().getEps() * (1 + x.getNorm())) {
//                    x = x.add(dx.mapMultiply(alfa));
//                    y = y.add((dy.subtract(y)).mapMultiply(alfa));
                    break;
                }

            }

            if (search == getSettings().getMaxLineSearchIteration()) {
                failedSearch += 1;
            }

            boolean notMonotone = false;
            if (search == getSettings().getMaxLineSearchIteration()) {

                alpha = 1.0;
                search = 0;
                Double max = Collections.max(oldPenalty);
                if (oldPenalty.size() == 1) {
                    max = max * 1.3;
                }
                while ((alfaPenalty - max) >= getSettings().getMu() * alpha * currentPenaltyGrad &&
                       search < getSettings().getMaxLineSearchIteration()) {

                    double alfaStar = -0.5 * alpha * alpha * currentPenaltyGrad / (-alpha * currentPenaltyGrad + alfaPenalty - currentPenalty);


                    alpha = FastMath.max(getSettings().getB() * alpha, FastMath.min(1.0, alfaStar));
                    // alfa = FastMath.min(1.0, FastMath.max(this.b * alfa, alfaStar));
                    // alfa = FastMath.max(this.b * alfa, alfaStar);
                    alphaF = this.getObj().value(x.add(dx.mapMultiply(alpha)));
                    alfaPenalty = penaltyFunction(alphaF, alpha, x, y, dx, dy.subtract(y), r);
                    search = search + 1;

                }
                notMonotone = true;
            }

            //UPDATE ALL FUNCTION
            RealVector oldGradient = functionGradient;
            RealMatrix oldJacob = constraintJacob;
            // RealVector old1 = penaltyFunctionGradX(oldGradient,x, y.add((dy.subtract(y)).mapMultiply(alfa)),r);
            RealVector old1 = lagrangianGradX(oldGradient, oldJacob, x, y.add((dy.subtract(y)).mapMultiply(alpha)));
            functionEval = alphaF;
            functionGradient = this.getObj().gradient(x.add(dx.mapMultiply(alpha)));
            if (this.getEqConstraint() != null) {

                equalityEval = this.getEqConstraint().value(x.add(dx.mapMultiply(alpha)));
            }
            if (this.getIqConstraint() != null) {

                inequalityEval = this.getIqConstraint().value(x.add(dx.mapMultiply(alpha)));
            }

            constraintJacob = computeJacobianConstraint(x.add(dx.mapMultiply(alpha)));
            // RealVector new1 = penaltyFunctionGradX(functionGradient,x.add(dx.mapMultiply(alfa)), y.add((dy.subtract(y)).mapMultiply(alfa)),r);
            RealVector new1 = lagrangianGradX(functionGradient, constraintJacob, x.add(dx.mapMultiply(alpha)), y.add((dy.subtract(y)).mapMultiply(alpha)));

            if (failedSearch == 2) {
                hessian = MatrixUtils.createRealIdentityMatrix(x.getDimension());
                failedSearch = 0;
            }

            if (!notMonotone) {
                hessian = BFGSFormula(hessian, dx, alpha, new1, old1);
            }

            //STORE PENALTY IN QUEQUE TO PERFORM NOT MONOTONE LINE SEARCH IF NECESSARY
            // oldPenalty.add(alfaPenalty);
            oldPenalty.add(currentPenalty);
            if (oldPenalty.size() > FORGETTING_FACTOR) {
                oldPenalty.remove(0);
            }

            x = x.add(dx.mapMultiply(alpha));
            y = y.add((dy.subtract(y)).mapMultiply(alpha));
        }

        return new LagrangeSolution(x, y, functionEval);
    }

    private double penaltyFunction(double currentF, double alfa, RealVector x, RealVector y, RealVector dx, RealVector uv, RealVector r) {
        // the set of constraints is the same as the previous one but they must be evaluated with the increment

        int me = 0;
        int mi;
        double partial = currentF;
        RealVector yalfa = y.add(uv.mapMultiply(alfa));
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
            RealVector re = r.getSubVector(0, me);
            RealVector ye = yalfa.getSubVector(0, me);
            RealVector g = this.getEqConstraint().value(x.add(dx.mapMultiply(alfa))).subtract(getEqConstraint().getLowerBound());
            RealVector g2 = g.ebeMultiply(g);
            partial -= ye.dotProduct(g) - 0.5 * re.dotProduct(g2);
        }

        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
            RealVector ri = r.getSubVector(me, mi);
            RealVector yi = yalfa.getSubVector(me, mi);
            RealVector yk = y.getSubVector(me, mi);
            RealVector gk = this.inequalityEval.subtract(getIqConstraint().getLowerBound());
            RealVector g = this.getIqConstraint().value(x.add(dx.mapMultiply(alfa))).subtract(getIqConstraint().getLowerBound());
            RealVector mask = new ArrayRealVector(g.getDimension(), 1.0);

            for (int i = 0; i < gk.getDimension(); i++) {

                if (gk.getEntry(i) > (yk.getEntry(i) / ri.getEntry(i))) {

                    mask.setEntry(i, 0.0);

                    partial -= 0.5 * yi.getEntry(i) * yi.getEntry(i) / ri.getEntry(i);
                }

            }

            RealVector g2 = g.ebeMultiply(g.ebeMultiply(mask));
            partial -= yi.dotProduct(g.ebeMultiply(mask)) - 0.5 * ri.dotProduct(g2);

        }

        return partial;
    }

    private double penaltyFunctionGrad(RealVector currentGrad, RealVector dx, RealVector dy, RealVector x, RealVector y, RealVector r) {

        int me = 0;
        int mi;
        double partial = currentGrad.dotProduct(dx);
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
            RealVector re = r.getSubVector(0, me);
            RealVector ye = y.getSubVector(0, me);
            RealVector dye = dy.getSubVector(0, me);
            RealVector ge = this.getEqConstraint().value(x).subtract(getEqConstraint().getLowerBound());
            RealMatrix jacob = this.getEqConstraint().jacobian(x);
            RealVector firstTerm = jacob.transpose().operate(ye);
            RealVector secondTerm = jacob.transpose().operate(ge.ebeMultiply(re));
//partial -= firstTerm.dotProduct(dx) - secondTerm.dotProduct(dx) + g.dotProduct(dye);
            partial += -firstTerm.dotProduct(dx) + secondTerm.dotProduct(dx) - ge.dotProduct(dye.subtract(ye));
        }

        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();

            RealVector ri = r.getSubVector(me, mi);
            RealVector dyi = dy.getSubVector(me, mi);
            RealVector yi = y.getSubVector(me, mi);

            RealVector gi = this.getIqConstraint().value(x).subtract(getIqConstraint().getLowerBound());
            RealMatrix jacob = this.getIqConstraint().jacobian(x);

            RealVector mask = new ArrayRealVector(mi, 1.0);
            RealVector viri = new ArrayRealVector(mi, 0.0);

            for (int i = 0; i < gi.getDimension(); i++) {

                if (gi.getEntry(i) > yi.getEntry(i) / ri.getEntry(i)) {
                    mask.setEntry(i, 0.0);

                } else {
                    viri.setEntry(i, yi.getEntry(i) / ri.getEntry(i));
                }

            }
            RealVector firstTerm = jacob.transpose().operate(yi.ebeMultiply(mask));
            RealVector secondTerm = jacob.transpose().operate(gi.ebeMultiply(ri.ebeMultiply(mask)));

            partial -= firstTerm.dotProduct(dx) - secondTerm.dotProduct(dx) + (gi.ebeMultiply(mask)).dotProduct(dyi.subtract(yi)) + viri.dotProduct(dyi.subtract(yi));
            // partial -= firstTerm.dotProduct(dx) - secondTerm.dotProduct(dx) + g.dotProduct(dyi.ebeMultiply(mask));
        }

        return partial;
    }

    private RealVector penaltyFunctionGradX(RealVector currentGrad, RealVector x, RealVector y, RealVector r) {

        int me = 0;
        int mi;
        RealVector partial = currentGrad.copy();
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
            RealVector re = r.getSubVector(0, me);
            RealVector ye = y.getSubVector(0, me);

            RealVector ge = this.getEqConstraint().value(x).subtract(getEqConstraint().getLowerBound());
            RealMatrix jacob = this.getEqConstraint().jacobian(x);

            RealVector firstTerm = jacob.transpose().operate(ye);
            RealVector secondTerm = jacob.transpose().operate(ge.ebeMultiply(re));

            partial = partial.subtract(firstTerm.subtract(secondTerm));
        }

        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();

            RealVector ri = r.getSubVector(me, mi);

            RealVector yi = y.getSubVector(me, mi);
            RealVector gi = this.getIqConstraint().value(x).subtract(getIqConstraint().getLowerBound());
            RealMatrix jacob = this.getIqConstraint().jacobian(x);

            RealVector mask = new ArrayRealVector(mi, 1.0);

            for (int i = 0; i < gi.getDimension(); i++) {

                if (gi.getEntry(i) > yi.getEntry(i) / ri.getEntry(i)) {
                    mask.setEntry(i, 0.0);

                }

            }
            RealVector firstTerm = jacob.transpose().operate(yi.ebeMultiply(mask));
            RealVector secondTerm = jacob.transpose().operate(gi.ebeMultiply(ri.ebeMultiply(mask)));
            partial = partial.subtract(firstTerm.subtract(secondTerm));
        }

        return partial;
    }

    private RealVector penaltyFunctionGradY(RealVector x, RealVector y, RealVector r) {

        int me = 0;
        int mi;
        RealVector partial = new ArrayRealVector(y.getDimension());
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
            RealVector g = this.getEqConstraint().value(x).subtract(getEqConstraint().getLowerBound());
            partial.setSubVector(0, g.mapMultiply(-1.0));
        }

        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();

            RealVector ri = r.getSubVector(me, mi);

            RealVector yi = y.getSubVector(me, mi);
            RealVector gi = this.getIqConstraint().value(x).subtract(getIqConstraint().getLowerBound());

            RealVector mask = new ArrayRealVector(mi, 1.0);

            RealVector viri = new ArrayRealVector(mi, 0.0);

            for (int i = 0; i < gi.getDimension(); i++) {

                if (gi.getEntry(i) > yi.getEntry(i) / ri.getEntry(i)) {
                    mask.setEntry(i, 0.0);

                } else {
                    viri.setEntry(i, yi.getEntry(i) / ri.getEntry(i));
                }

            }

            partial.setSubVector(me, gi.ebeMultiply(mask).add(viri).mapMultiply(-1.0));
        }

        return partial;
    }

    private RealVector updateRj(RealMatrix H, RealVector y, RealVector dx, RealVector dy, RealVector r, double additional) { //r = updateRj(currentH,dx,y,u,r,sigm);
        //CALCULATE SIGMA VECTOR THAT DEPEND FROM ITERATION
        RealVector sigma = new ArrayRealVector(r.getDimension());
        for (int i = 0; i < sigma.getDimension(); i++) {
            final double appoggio = iterations.getCount() / FastMath.sqrt(r.getEntry(i));
            sigma.setEntry(i, FastMath.min(1.0, appoggio));
        }

        int me = 0;
        int mi = 0;
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
        }

        RealVector sigmar = sigma.ebeMultiply(r);
        //(u-v)^2 or (ru-v)
        RealVector numerator = ((dy.subtract(y)).ebeMultiply(dy.subtract(y))).mapMultiply(2.0 * (mi + me));

        double denominator = dx.dotProduct(H.operate(dx)) * (1.0 - additional);
        RealVector r1 = r.copy();
        if (getEqConstraint() != null) {
            for (int i = 0; i < me; i++) {
                r1.setEntry(i, FastMath.max(sigmar.getEntry(i), numerator.getEntry(i) / denominator));
            }
        }
        if (getIqConstraint() != null) {
            for (int i = 0; i < mi; i++) {
                r1.setEntry(me + i, FastMath.max(sigmar.getEntry(me + i), numerator.getEntry(me + i) / denominator));
            }
        }


        return r1;
    }

    private double updateRho(RealVector dx, RealVector dy, RealMatrix H, RealMatrix jacobianG, double additionalVariable) {
        double num = 10.0 * FastMath.pow(dx.dotProduct(jacobianG.transpose().operate(dy)), 2);
        double den = (1.0 - additionalVariable) * (1.0 - additionalVariable) * dx.dotProduct(H.operate(dx));
        return FastMath.max(10.0, num / den);
    }

    private LagrangeSolution solveQP(RealVector x, RealVector y, double rho) {

        RealMatrix H = hessian;
        RealVector g = functionGradient;

        int me = 0;
        int mi = 0;
        int add = 0;
        boolean violated = false;
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        if (getIqConstraint() != null) {

            mi = getIqConstraint().dimY();
            violated = inequalityEval.subtract(getIqConstraint().getLowerBound()).getMinValue() <= getSettings().getEps() ||
                       y.getMaxValue() >= 0;

        }
        // violated = true;
        if (me > 0 || violated) {
            add = 1;
        }

        RealMatrix H1 = new Array2DRowRealMatrix(H.getRowDimension() + add, H.getRowDimension() + add);
        H1.setSubMatrix(H.getData(), 0, 0);
        if (add == 1) {
            H1.setEntry(H.getRowDimension(), H.getRowDimension(), rho);
        }

        RealVector g1 = new ArrayRealVector(g.getDimension() + add);
        g1.setSubVector(0, g);

        LinearEqualityConstraint eqc = null;
        RealVector conditioneq;
        if (getEqConstraint() != null) {
            RealMatrix eqJacob = constraintJacob.getSubMatrix(0, me - 1, 0, x.getDimension() - 1);

            RealMatrix Ae = new Array2DRowRealMatrix(me, x.getDimension() + add);
            RealVector be = new ArrayRealVector(me);
            Ae.setSubMatrix(eqJacob.getData(), 0, 0);
            conditioneq = this.equalityEval.subtract(getEqConstraint().getLowerBound());
            Ae.setColumnVector(x.getDimension(), conditioneq.mapMultiply(-1.0));

            be.setSubVector(0, getEqConstraint().getLowerBound().subtract(this.equalityEval));
            eqc = new LinearEqualityConstraint(Ae, be);

        }
        LinearInequalityConstraint iqc = null;

        if (getIqConstraint() != null) {
//
            RealMatrix iqJacob = constraintJacob.getSubMatrix(me, me + mi - 1, 0, x.getDimension() - 1);

            RealMatrix Ai = new Array2DRowRealMatrix(mi, x.getDimension() + add);
            RealVector bi = new ArrayRealVector(mi);
            Ai.setSubMatrix(iqJacob.getData(), 0, 0);

            RealVector conditioniq = this.inequalityEval.subtract(getIqConstraint().getLowerBound());

            if (add == 1) {

                for (int i = 0; i < conditioniq.getDimension(); i++) {

                    if (!(conditioniq.getEntry(i) <= getSettings().getEps() || y.getEntry(me + i) > 0)) {
                        conditioniq.setEntry(i, 0);
                    }
                }

                Ai.setColumnVector(x.getDimension(), conditioniq.mapMultiply(-1.0));

            }
            bi.setSubVector(0, getIqConstraint().getLowerBound().subtract(this.inequalityEval));
            iqc = new LinearInequalityConstraint(Ai, bi);

        }
        LinearBoundedConstraint bc = null;
        if (add == 1) {

            RealMatrix sigmaA = new Array2DRowRealMatrix(1, x.getDimension() + 1);
            sigmaA.setEntry(0, x.getDimension(), 1.0);
            ArrayRealVector lb = new ArrayRealVector(1, 0.0);
            ArrayRealVector ub = new ArrayRealVector(1, 1.0);

            bc = new LinearBoundedConstraint(sigmaA, lb, ub);

        }

        QuadraticFunction q = new QuadraticFunction(H1, g1, 0);

        ADMMQPOptimizer solver = new ADMMQPOptimizer();
        LagrangeSolution sol = solver.optimize(new ObjectiveFunction(q), iqc, eqc, bc, getMatrixDecompositionTolerance());

        final double sigma;
        if (add == 1) {
            sigma = sol.getX().getEntry(x.getDimension());
        } else {
            sigma = 0;
        }

        return new LagrangeSolution(sol.getX().getSubVector(0, x.getDimension()),
                                    sol.getLambda().getSubVector(0, me + mi),
                                    sigma);

    }

    private RealMatrix computeJacobianConstraint(RealVector x) {
        int me = 0;
        int mi = 0;
        RealMatrix je = null;
        RealMatrix ji = null;
        if (this.getEqConstraint() != null) {
            me = this.getEqConstraint().dimY();
            je = this.getEqConstraint().jacobian(x);
        }

        if (this.getIqConstraint() != null) {
            mi = this.getIqConstraint().dimY();
            ji = this.getIqConstraint().jacobian(x);
        }

        RealMatrix jacobian = new Array2DRowRealMatrix(me + mi, x.getDimension());
        if (me > 0) {
            jacobian.setSubMatrix(je.getData(), 0, 0);
        }
        if (mi > 0) {
            jacobian.setSubMatrix(ji.getData(), me, 0);
        }

        return jacobian;
    }
    /*
     *DAMPED BFGS FORMULA
     */

    private RealMatrix BFGSFormula(RealMatrix oldH, RealVector dx, double alfa, RealVector newGrad, RealVector oldGrad) {

        RealMatrix oldH1 = oldH.copy();
        RealVector y1 = newGrad.subtract(oldGrad);
        RealVector s = dx.mapMultiply(alfa);
        double theta = 1.0;
        if (s.dotProduct(y1) < 0.2 * s.dotProduct(oldH.operate(s))) {
            theta = 0.8 * s.dotProduct(oldH.operate(s)) / (s.dotProduct(oldH.operate(s)) - s.dotProduct(y1));
        }

        RealVector y = y1.mapMultiply(theta).add((oldH.operate(s)).mapMultiply(1.0 - theta));

        RealMatrix firstTerm = y.outerProduct(y).scalarMultiply(1.0 / s.dotProduct(y));
        RealMatrix secondTerm = oldH1.multiply(s.outerProduct(s)).multiply(oldH);
        double thirtTerm = s.dotProduct(oldH1.operate(s));
        RealMatrix Hnew = oldH1.add(firstTerm).subtract(secondTerm.scalarMultiply(1.0 / thirtTerm));
        //RESET HESSIAN IF NOT POSITIVE DEFINITE
        EigenDecompositionSymmetric dsX = new EigenDecompositionSymmetric(Hnew,
                                                                          getMatrixDecompositionTolerance().getEpsMatrixDecomposition(),
                                                                          true);
        double min = new ArrayRealVector(dsX.getEigenvalues()).getMinValue();
        if (min < 0) {
            Hnew = MatrixUtils.createRealIdentityMatrix(oldH.getRowDimension());
        }
        return Hnew;

    }

}
