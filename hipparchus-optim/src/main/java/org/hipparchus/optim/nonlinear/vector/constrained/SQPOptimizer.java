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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.MathUtils;

/**
 * Sequential Quadratic Programming Optimizer.
 * <br/>
 * min f(x)
 * <br/>
 * q(x)=b1
 * <br/>
 * h(x)>=b2
 * <br/>
 * Algorithm based on paper:"On the convergence of a sequential quadratic
 * programming method(Klaus Shittkowki,January 1982)"
 * @since 3.1
 */
public class SQPOptimizer extends ConstraintOptimizer {

    int forgetFactor = 20;
    //Convergence and Constraint tolerance
    double epsilon = SQPOption.SQP_eps;//>0
    //QP SUBPROBLEM PARAMETER
    private double rhoConstant = SQPOption.SQP_rhoCons;//rho>1
    private double sigmaMax = SQPOption.SQP_sigmaMax;//0<sigma<1
    private int qpMaxLoop = SQPOption.SQP_qpMaxLoop;
    //ARMIJO CONDITION PARAMETER
    private double mu = SQPOption.SQP_mu;//[0,0.5]
    //LINE SEARCH PARAMETER
    private double b = SQPOption.SQP_b;//[0;1]
    private boolean useFunHessian = SQPOption.SQP_useFunHessian;
    private int maxLineSearchIteration = 10;

    private TwiceDifferentiableFunction obj;
    private LinearEqualityConstraint eqConstraint;
    private LinearInequalityConstraint iqConstraint;
    private LinearBoundedConstraint bqConstraint;
    private ArrayRealVector xStart;
    //QP SOLVER
    QPOptimizer qpSolver = new ADMMQPOptimizer();

    @Override
    public LagrangeSolution optimize(OptimizationData... optData) {
        return super.optimize(optData);
    }

    @Override
    protected void parseOptimizationData(OptimizationData... optData) {
        super.parseOptimizationData(optData);
        for (OptimizationData data : optData) {

            if (data instanceof ObjectiveFunction) {
                obj = (TwiceDifferentiableFunction) ((ObjectiveFunction) data).getObjectiveFunction();
                continue;
            }

            if (data instanceof LinearEqualityConstraint) {
                eqConstraint = (LinearEqualityConstraint) data;
                continue;
            }
            if (data instanceof LinearInequalityConstraint) {
                iqConstraint = (LinearInequalityConstraint) data;
                continue;
            }

            if (data instanceof LinearBoundedConstraint) {
                bqConstraint = (LinearBoundedConstraint) data;
                continue;
            }

            if (data instanceof QPOptimizer) {
                qpSolver = (QPOptimizer) data;
                continue;

            }

            if (data instanceof SQPOption) {
                epsilon = ((SQPOption) data).eps;
                sigmaMax = ((SQPOption) data).sigmaMax;
                qpMaxLoop = ((SQPOption) data).qpMaxLoop;
                rhoConstant = ((SQPOption) data).rhoCons;
                mu = ((SQPOption) data).mu;
                b = ((SQPOption) data).b;
                maxLineSearchIteration = ((SQPOption) data).maxLineSearchIteration;
                useFunHessian = ((SQPOption) data).useFunHessian;
                continue;

            }

        }
        // if we got here, convexObjective exists
        int n = obj.dim();
        if (eqConstraint != null) {
            int nDual = eqConstraint.dimY();
            if (nDual >= n) {
                throw new IllegalArgumentException("Rank of constraints must be < domain dimension");
            }
            int nTest = eqConstraint.dim();
            if (nDual == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_NOT_ALLOWED);
            }
            MathUtils.checkDimension(nTest, n);
        }

    }

    @Override
    public LagrangeSolution doOptimize() {
        int me = 0;
        int mi = 0;
        int mb = 0;

        //EQUALITY CONSTRAINT
        if (this.eqConstraint != null) {
            me = eqConstraint.dimY();
        }
        //INEQUALITY CONSTRAINT
        if (this.iqConstraint != null) {
            mi = iqConstraint.dimY();
        }
        //BOUNDED CONSTRAINT
        if (this.bqConstraint != null) {
            mb = bqConstraint.dimY();
        }

        double alfa = 1.0;
        double rho = 1.0;

        RealVector x = null;
        if (this.getStartPoint() != null) {
            x = new ArrayRealVector(this.getStartPoint());
        } else {
            x = new ArrayRealVector(this.obj.dim());
        }

        RealVector y = new ArrayRealVector(me + mi, 0.0);
        RealVector dx = new ArrayRealVector(x.getDimension());
        RealVector dy = new ArrayRealVector(y.getDimension());

        RealVector r = new ArrayRealVector(me + mi, 1.0);
        RealVector u = new ArrayRealVector(me + mi, 0);
        ArrayList<Double> oldPenalty = new ArrayList<Double>();
        //INITIAL VALUES
        double currentF = this.obj.value(x);

        RealMatrix currentH = null;

        if (useFunHessian == false) {
            currentH = MatrixUtils.createRealIdentityMatrix(x.getDimension());
        } else {
            currentH = this.obj.hessian(x);
        }

        RealVector currentGrad = this.obj.gradient(x);
        RealMatrix currentjacobianGrad = computeJacobianGrad(x);
        int penaltyPositive = 0;
        for (int i = 0; i < this.getMaxIterations(); i++) {
            RealVector dx1 = new ArrayRealVector(x.getDimension());

            iterations.increment();



            alfa = 1.0;
            LagrangeSolution sol1 = null;



            int qpLoop = 0;
            double sigma = 10;
            double currentPenaltyGrad = 0;


            //LOOP TO FIND SOLUTION WITH SIGMA<SIGMA TRESHOLD
            while (sigma > sigmaMax && qpLoop < this.qpMaxLoop) {
                sol1 = solveQP(currentH, currentGrad, x, y, rho,currentF);
                sigma = sol1.getValue();
                if (sigma > sigmaMax) {
                    rho = this.rhoConstant * rho;
                    qpLoop += 1;
                }

            }
            //IF SIGMA>SIGMA TRESHOLD ASSIGN DIRECTION FROM PENALTY GRADIENT
            if (qpLoop == this.qpMaxLoop) {
                penaltyPositive = 0;
                dx = (MatrixUtils.inverse(currentH).operate(penaltyFunctionGradX(currentGrad, x, y, r))).mapMultiply(-1.0);
                dy = y.subtract(penaltyFunctionGradY(currentGrad, x, y, r));
                sigma = 0;

            } else {
                dx = sol1.getX();
                sigma = sol1.getValue();
                dy = sol1.getLambda();
                r = updateRj(currentH, x, y, dx, dy, r, sigma);
            }

            currentPenaltyGrad = penaltyFunctionGrad(currentGrad, dx, dy, x, y, r);
            rho = updateRho(dx, dy, currentH, currentjacobianGrad, sigma, rho);
            if ((dx.getNorm() < (Math.sqrt(epsilon) * (1 + x.getNorm()))))
            {

                    break;
            }


            double currentPenalty = penaltyFunction(currentF, x, y, r, y);
            //STORE PENALTY IN QUEQUE TO PERFORM NOT MONOTONE LINE SEARCH IF NECESSARY
            oldPenalty.add(currentPenalty);
            if (oldPenalty.size() > this.forgetFactor)
            {
                oldPenalty.remove(0);
            }
            double alfaF = this.obj.value(x.add(dx.mapMultiply(alfa)));
            double alfaPenalty = penaltyFunction(alfaF, x.add(dx.mapMultiply(alfa)), y.add((dy.subtract(y)).mapMultiply(alfa)), r, y.add((dy.subtract(y)).mapMultiply(alfa)));


            //LINE SEARCH
            int search = 0;
            while ((alfaPenalty - currentPenalty) >= this.mu * alfa * currentPenaltyGrad && search < maxLineSearchIteration) {
                double alfaStar = -0.5 * alfa * alfa * currentPenaltyGrad / (-alfa * currentPenaltyGrad + alfaPenalty - currentPenalty);
                 alfaStar = Math.min(1, alfaStar);

                alfa = Math.max(this.b * alfa, alfaStar);
                alfaF = this.obj.value(x.add(dx.mapMultiply(alfa)));
                alfaPenalty = penaltyFunction(alfaF, x.add(dx.mapMultiply(alfa)), y.add((dy.subtract(y)).mapMultiply(alfa)), r, y.add((dy.subtract(y)).mapMultiply(alfa)));
                search = search + 1;

            }

           //MONOTONE LINE SEARCH
            if (search == maxLineSearchIteration  || currentPenaltyGrad <=0  ) {

                alfa = 1.0;
                search = 0;
                Double max = Collections.min(oldPenalty);
                if (oldPenalty.size() == 1) {
                    max = max * 1.3;
                }
                while ((alfaPenalty - max) >= this.mu * alfa * currentPenaltyGrad && search < maxLineSearchIteration) {

                    double alfaStar = -0.5 * alfa * alfa * currentPenaltyGrad / (-alfa * currentPenaltyGrad + alfaPenalty - currentPenalty);
                    alfaStar = Math.min(1, alfaStar);

                    alfa = Math.max(this.b * alfa, alfaStar);
                    alfaF = this.obj.value(x.add(dx.mapMultiply(alfa)));
                    alfaPenalty = penaltyFunction(alfaF, x.add(dx.mapMultiply(alfa)), y.add((dy.subtract(y)).mapMultiply(alfa)), r, y.add((dy.subtract(y)).mapMultiply(alfa)));
                    search = search + 1;

                }
            }
            RealVector oldGrad = currentGrad.copy();
            currentGrad = this.obj.gradient(x.add(dx.mapMultiply(alfa)));
            currentF = alfaF;

            if (search == maxLineSearchIteration) {
                currentH = MatrixUtils.createRealIdentityMatrix(x.getDimension());

            } else {
                    if (currentPenaltyGrad <= 0) {

                        RealMatrix oldH = currentH.copy();
                        currentH = BFGSFormula(oldH, x, dx, alfa, currentGrad, oldGrad);
                    }
                }


                currentjacobianGrad = computeJacobianGrad(x.add(dx.mapMultiply(alfa)));
                x = x.add(dx.mapMultiply(alfa));
                y = y.add((dy.subtract(y)).mapMultiply(alfa));
                u = dy.copy();



        }




        return new LagrangeSolution(x, y, this.obj.value(x));
    }

    private double penaltyFunction(double currentF, RealVector x, RealVector y, RealVector r, RealVector yk) {
        // the set of constraints is the same as the previous one but they must be evaluated with the increment

        int me = 0;
        int mi = 0;
        double partial = currentF;

        if (eqConstraint != null) {
            me = eqConstraint.dimY();
            RealVector re = r.getSubVector(0, me);
            RealVector ye = y.getSubVector(0, me);
            RealVector g = this.eqConstraint.value(x).subtract(eqConstraint.getLowerBound());
            RealVector g2 = g.ebeMultiply(g);
            partial -= ye.dotProduct(g) - 0.5 * re.dotProduct(g2);
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();
            RealVector ri = r.getSubVector(me, mi);
            RealVector yi = y.getSubVector(me, mi);
            RealVector yik = yk.getSubVector(me, mi);
            RealVector g = this.iqConstraint.value(x).subtract(iqConstraint.getLowerBound());
            RealVector mask = new ArrayRealVector(g.getDimension(), 1.0);

            for (int i = 0; i < g.getDimension(); i++) {

                if (g.getEntry(i) > yik.getEntry(i) / ri.getEntry(i)) {

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
        int mi = 0;
        double partial = currentGrad.dotProduct(dx);
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
            RealVector re = r.getSubVector(0, me);
            RealVector ye = y.getSubVector(0, me);
            RealVector dye = dy.getSubVector(0, me);
            RealVector ge = this.eqConstraint.value(x).subtract(eqConstraint.getLowerBound());
            RealMatrix jacob = this.eqConstraint.jacobian(x);
            RealVector firstTerm = (jacob.transpose().operate(ye));
            RealVector secondTerm = (jacob.transpose().operate(ge.ebeMultiply(re)));
//partial -= firstTerm.dotProduct(dx) - secondTerm.dotProduct(dx) + g.dotProduct(dye);
            partial -= firstTerm.dotProduct(dx) - secondTerm.dotProduct(dx) + ge.dotProduct(dye.subtract(ye));
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();

            RealVector ri = r.getSubVector(me, mi);
            RealVector dyi = dy.getSubVector(me, mi);
            RealVector yi = y.getSubVector(me, mi);

            RealVector gi = this.iqConstraint.value(x).subtract(iqConstraint.getLowerBound());
            RealMatrix jacob = this.iqConstraint.jacobian(x);

            RealVector mask = new ArrayRealVector(mi, 1.0);
            RealVector viri = new ArrayRealVector(mi, 0.0);

            for (int i = 0; i < gi.getDimension(); i++) {

                if (gi.getEntry(i) > yi.getEntry(i) / ri.getEntry(i)) {
                    mask.setEntry(i, 0.0);

                } else {
                    viri.setEntry(i, yi.getEntry(i) / ri.getEntry(i));
                }

            }
            RealVector firstTerm = (jacob.transpose().operate(yi.ebeMultiply(mask)));
            RealVector secondTerm = (jacob.transpose().operate(gi.ebeMultiply(ri.ebeMultiply(mask))));

            partial -= firstTerm.dotProduct(dx) - secondTerm.dotProduct(dx) + gi.dotProduct(dyi.subtract(yi).ebeMultiply(mask)) + viri.dotProduct(dyi.subtract(yi));
            // partial -= firstTerm.dotProduct(dx) - secondTerm.dotProduct(dx) + g.dotProduct(dyi.ebeMultiply(mask));
        }

        return partial;
    }

    private RealVector LagrangianeGradX(RealVector currentGrad, RealVector x, RealVector y) {

        int me = 0;
        int mi = 0;
        RealVector partial = currentGrad.copy();
        if (eqConstraint != null) {
            me = eqConstraint.dimY();

            RealVector ye = y.getSubVector(0, me);

            RealMatrix jacob = this.eqConstraint.jacobian(x);

            RealVector firstTerm = (jacob.transpose().operate(ye));

            partial = partial.subtract(firstTerm);
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();

            RealVector yi = y.getSubVector(me, mi);

            RealMatrix jacob = this.iqConstraint.jacobian(x);

            RealVector firstTerm = (jacob.transpose().operate(yi));

            partial = partial.subtract(firstTerm);
        }
        return partial;
    }

    private RealVector penaltyFunctionGradX(RealVector currentGrad, RealVector x, RealVector y, RealVector r) {

        int me = 0;
        int mi = 0;
        RealVector partial = currentGrad.copy();
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
            RealVector re = r.getSubVector(0, me);
            RealVector ye = y.getSubVector(0, me);

            RealVector ge = this.eqConstraint.value(x).subtract(eqConstraint.getLowerBound());
            RealMatrix jacob = this.eqConstraint.jacobian(x);

            RealVector firstTerm = (jacob.transpose().operate(ye));
            RealVector secondTerm = (jacob.transpose().operate(ge.ebeMultiply(re)));

            partial = partial.subtract(firstTerm.subtract(secondTerm));
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();

            RealVector ri = r.getSubVector(me, mi);

            RealVector yi = y.getSubVector(me, mi);
            RealVector gi = this.iqConstraint.value(x).subtract(iqConstraint.getLowerBound());
            RealMatrix jacob = this.iqConstraint.jacobian(x);

            RealVector mask = new ArrayRealVector(mi, 1.0);

            for (int i = 0; i < gi.getDimension(); i++) {

                if (gi.getEntry(i) > yi.getEntry(i) / ri.getEntry(i)) {
                    mask.setEntry(i, 0.0);

                }

            }
            RealVector firstTerm = (jacob.transpose().operate(yi.ebeMultiply(mask)));
            RealVector secondTerm = (jacob.transpose().operate(gi.ebeMultiply(ri.ebeMultiply(mask))));
            partial = partial.subtract(firstTerm.subtract(secondTerm));
        }

        return partial;
    }

    private RealVector penaltyFunctionGradY(RealVector currentGrad, RealVector x, RealVector y, RealVector r) {

        int me = 0;
        int mi = 0;
        RealVector partial = new ArrayRealVector(y.getDimension());
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
            RealVector g = this.eqConstraint.value(x).subtract(eqConstraint.getLowerBound());
            partial.setSubVector(0, g.mapMultiply(-1.0));
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();

            RealVector ri = r.getSubVector(me, mi);

            RealVector yi = y.getSubVector(me, mi);
            RealVector gi = this.iqConstraint.value(x).subtract(iqConstraint.getLowerBound());

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

    private RealVector updateRj(RealMatrix H, RealVector x, RealVector y, RealVector dx, RealVector dy, RealVector r, double additional) { //r = updateRj(currentH,dx,y,u,r,sigm);
        //CALCULATE SIGMA VECTOR THAT DEPEND FROM ITERATION
        RealVector sigma = new ArrayRealVector(r.getDimension());
        for (int i = 0; i < sigma.getDimension(); i++) {
            double appoggio = 0;
            // if (r.getEntry(i)>epsilon)
            appoggio = ((this.iterations.getCount()) / (Math.sqrt(r.getEntry(i))));
            sigma.setEntry(i, Math.min(1.0, appoggio));
        }

        int me = 0;
        int mi = 0;
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
        }
        if (iqConstraint != null) {
            mi = iqConstraint.dimY();
        }

        RealVector sigmar = sigma.ebeMultiply(r);
        //(u-v)^2 or (ru-v)
        RealVector numeratore = (dy.subtract(y).ebeMultiply(dy.subtract(y))).mapMultiply(2.0 * (mi + me));
       // RealVector numeratore = (dy.subtract(y)).mapMultiply(2.0 * (mi + me));

        //RealVector numeratore = dy.subtract(y).ebeMultiply(dy.subtract(y)).mapMultiply(2.0*(mi + me));

        double denominatore = dx.dotProduct(H.operate(dx)) * (1.0 - additional);
        RealVector r1 = r.copy();
        if (eqConstraint != null) {
            RealVector ye = y.getSubVector(0, me);
            RealVector g = eqConstraint.value(x).subtract(eqConstraint.getLowerBound());
            for (int i = 0; i < me; i++) {

                // if (g.getEntry(i)<=epsilon || ye.getEntry(i)>0 )
                r1.setEntry(i, Math.max(sigmar.getEntry(i), numeratore.getEntry(i) / denominatore));

            }
        }
        if (iqConstraint != null) {
            RealVector yi = y.getSubVector(me, mi);
            RealVector g = iqConstraint.value(x).subtract(iqConstraint.getLowerBound());
            for (int i = 0; i < mi; i++) {
                //if (g.getEntry(i)<=epsilon || yi.getEntry(i)>0 )
                r1.setEntry(me + i, Math.max(sigmar.getEntry(me + i), numeratore.getEntry(me + i) / denominatore));

            }
        }


        return r1;
    }

    private double updateRho(RealVector dx, RealVector dy, RealMatrix H, RealMatrix jacobianG, double additionalVariable, double rho) {
        double num = 10.0 * Math.pow(dy.dotProduct(jacobianG.operate(dx)), 2);
        double den = (1.0 - additionalVariable) * (1.0 - additionalVariable) * dx.dotProduct(H.operate(dx));

        return Math.max(1.0, num / den);
    }

    private LagrangeSolution solveQP(RealMatrix H, RealVector g, RealVector x, RealVector y, double rho,double evaluateFunction) {

        int me = 0;
        int mi = 0;
        int mb = 0;
        int add = 0;
        boolean violated = false;
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
        }
        if (iqConstraint != null) {

            mi = iqConstraint.dimY();
            violated = (iqConstraint.value(x).subtract(iqConstraint.getLowerBound()).getMaxValue() <= epsilon) || (y.getMinValue() > 0);

        }

        if (me > 0 || violated == true) {
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
        RealVector conditioneq = null;
        if (eqConstraint != null) {

            RealMatrix Ae = new Array2DRowRealMatrix(me, x.getDimension() + add);
            RealVector be = new ArrayRealVector(me);
            Ae.setSubMatrix(eqConstraint.jacobian(x).getData(), 0, 0);
            conditioneq = eqConstraint.value(x).subtract(eqConstraint.getLowerBound());
            Ae.setColumnVector(x.getDimension(), conditioneq.mapMultiply(-1.0));

            be.setSubVector(0, eqConstraint.getLowerBound().subtract(eqConstraint.value(x)));
            eqc = new LinearEqualityConstraint(Ae, be);

        }
        LinearInequalityConstraint iqc = null;

        if (iqConstraint != null) {

            RealMatrix Ai = new Array2DRowRealMatrix(mi, x.getDimension() + add);
            RealVector bi = new ArrayRealVector(mi);
            Ai.setSubMatrix(iqConstraint.jacobian(x).getData(), 0, 0);

            RealVector conditioniq = iqConstraint.value(x).subtract(iqConstraint.getLowerBound());

            if (add == 1) {

                for (int i = 0; i < conditioniq.getDimension(); i++) {

                    if (conditioniq.getEntry(i) > epsilon && y.getEntry(me + i) <= 0) {
                        conditioniq.setEntry(i, 0);
                    }
                }

                Ai.setColumnVector(x.getDimension(), conditioniq.mapMultiply(-1.0));

            }
            bi.setSubVector(0, iqConstraint.getLowerBound().subtract(iqConstraint.value(x)));
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

        QuadraticFunction q = new QuadraticFunction(H1, g1, evaluateFunction);
        LagrangeSolution sol = null;
        double sigma = 0;
        ADMMQPOptimizer solver = new ADMMQPOptimizer();
        sol = solver.optimize(new ObjectiveFunction(q), iqc, eqc, bc);

        if (add == 1) {
            sigma = sol.getX().getEntry(x.getDimension());
        } else {
            sigma = 0;
        }

        LagrangeSolution solnew = new LagrangeSolution(sol.getX().getSubVector(0, x.getDimension()), sol.getLambda().getSubVector(0, me + mi), sigma);
        return solnew;

    }

    private LagrangeSolution solveQP1(RealMatrix H, RealVector g, RealVector x, RealVector y, double rho) {

        RealMatrix H1 = new Array2DRowRealMatrix(H.getRowDimension(), H.getRowDimension());
        H1.setSubMatrix(H.getData(), 0, 0);

        RealVector g1 = new ArrayRealVector(g.getDimension());
        g1.setSubVector(0, g);
        QuadraticFunction q = new QuadraticFunction(H1, g1, 0);

        int me = 0;
        int mi = 0;

        LinearEqualityConstraint eqc = null;
        if (eqConstraint != null) {
            me = eqConstraint.dimY();
            RealMatrix Ae = new Array2DRowRealMatrix(me, x.getDimension());
            RealVector be = new ArrayRealVector(me);
            Ae.setSubMatrix(eqConstraint.jacobian(x).getData(), 0, 0);

            be.setSubVector(0, eqConstraint.getLowerBound().subtract(eqConstraint.value(x)));
            eqc = new LinearEqualityConstraint(Ae, be);
        }
        LinearInequalityConstraint iqc = null;
        if (iqConstraint != null) {
            mi = iqConstraint.dimY();
            RealMatrix Ai = new Array2DRowRealMatrix(mi, x.getDimension());
            RealVector bi = new ArrayRealVector(mi);
            Ai.setSubMatrix(iqConstraint.jacobian(x).getData(), 0, 0);

            bi.setSubVector(0, iqConstraint.getLowerBound().subtract(iqConstraint.value(x)));

            iqc = new LinearInequalityConstraint(Ai, bi);

        }

        LagrangeSolution sol = null;
        double sigma = 0;
        sol = qpSolver.optimize(new ObjectiveFunction(q), iqc, eqc);

        sigma = 0;

        LagrangeSolution solnew = new LagrangeSolution(sol.getX().getSubVector(0, x.getDimension()), sol.getLambda().getSubVector(0, me + mi), sigma);
        return solnew;

    }

    private RealMatrix computeJacobianGrad(RealVector x) {
        int me = 0;
        int mi = 0;
        RealMatrix je = null;
        RealMatrix ji = null;
        if (this.eqConstraint != null) {
            me = this.eqConstraint.dimY();
            je = this.eqConstraint.jacobian(x);
        }

        if (this.iqConstraint != null) {
            mi = this.iqConstraint.dimY();
            ji = this.iqConstraint.jacobian(x);
        }

        RealMatrix giacobian = new Array2DRowRealMatrix(me + mi, x.getDimension());
        if (me > 0) {
            giacobian.setSubMatrix(je.getData(), 0, 0);
        }
        if (mi > 0) {
            giacobian.setSubMatrix(ji.getData(), me, 0);
        }

        return giacobian;
    }

    RealMatrix BFGSFormula(RealMatrix oldH, RealVector x, RealVector dir, double alfa, RealVector newGrad, RealVector oldGrad) {

        RealMatrix oldH1 = oldH;
        RealVector y1 = newGrad.subtract(oldGrad);
        RealVector s = dir.mapMultiply(alfa);
        double theta = 1.0;
        if (s.dotProduct(y1) < 0.1 * s.dotProduct(oldH.operate(s))) {
            theta = 0.9 * s.dotProduct(oldH.operate(s)) / (s.dotProduct(oldH.operate(s)) - s.dotProduct(y1));
        }
        RealVector y = y1.mapMultiply(theta).add(oldH.operate(s).mapMultiply(1.0 - theta));

        RealMatrix firstTerm = y.outerProduct(y).scalarMultiply(1.0 / s.dotProduct(y));
        RealMatrix secondTerm = oldH1.operate(s).outerProduct(oldH1.operate(s));
        double thirtTerm = s.dotProduct(oldH1.operate(s));
        RealMatrix Hnew = oldH1.add(firstTerm).subtract(secondTerm.scalarMultiply(1.0 / thirtTerm));
        return Hnew;
    }

}
