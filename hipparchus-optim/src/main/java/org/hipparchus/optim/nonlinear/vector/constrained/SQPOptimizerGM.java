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
 * <br/>
 * min f(x)
 * <br/>
 * q(x)=b1
 * <br/>
 * h(x)>=b2
 * <br/>
 * Algorithm based on paper:"Some Theoretical properties of an augmented lagrangian merit function
 * (Gill,Murray,Sauders,Wriht,April 1986)"
 * @since 3.1
 */
public class SQPOptimizerGM extends AbstractSQPOptimizer {

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

    /** Simple constructor.
     */
    public SQPOptimizerGM() {
        super();
    }

    @Override
    public LagrangeSolution doOptimize() {
        int me = 0;
        int mi = 0;

        //EQUALITY CONSTRAINT
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        //INEQUALITY CONSTRAINT
        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
        }

        double alfa = 1.0;
        double rho = 100.0;
        RealVector x = null;
        if (getStartPoint() != null) {
            x = new ArrayRealVector(getStartPoint());
        } else {
            x = new ArrayRealVector(getObj().dim());
        }

        RealVector y = new ArrayRealVector(me + mi, 0.0);
         //INITIAL VALUES
        double functionEval = getObj().value(x);
        functionGradient = getObj().gradient(x);
        double maxGrad = functionGradient.getLInfNorm();

        if (getEqConstraint() != null) {
          equalityEval = getEqConstraint().value(x);
        }
        if (getIqConstraint() != null) {
            inequalityEval = getIqConstraint().value(x);
        }
        constraintJacob = computeJacobianConstraint(x);

        if (getEqConstraint() != null) {
            maxGrad = FastMath.max(maxGrad, equalityEval.getLInfNorm());
        }
        if (getIqConstraint() != null) {
            maxGrad = FastMath.max(maxGrad, inequalityEval.getLInfNorm());
        }

        if (!getSettings().getUseFunHessian()) {
            hessian= MatrixUtils.createRealIdentityMatrix(x.getDimension()).scalarMultiply(maxGrad);
        } else {
            hessian = getObj().hessian(x);
        }


        rho = 0;

        for (int i = 0; i < getMaxIterations(); i++) {

            iterations.increment();

            alfa = 1.0;

            LagrangeSolution sol1 = null;
            //SOLVE QP
            sol1 = solveQP(x,y);
            RealVector p = sol1.getX();
            RealVector e = sol1.getLambda().subtract(y);

            RealVector s = calculateSvector(y,rho);
            RealVector se = null;
            RealMatrix jacobi = null;
            RealVector q = null;
            RealVector qe = null;

            //TEST CON SI SOLO PER INEQUALITY AND Q FOR ALL
            if (getEqConstraint() != null) {
                se = new ArrayRealVector(me);
                jacobi = constraintJacob.getSubMatrix(0, me - 1, 0, x.getDimension() - 1);
                qe = new ArrayRealVector(me);
            }
            if (getIqConstraint() != null) {
                jacobi = constraintJacob.getSubMatrix(me, me + mi - 1, 0, x.getDimension() - 1);
                q = jacobi.operate(p).add(inequalityEval.subtract(getIqConstraint().getLowerBound())).subtract(s);
            }



            //CALCULATE PENALTY GRADIENT
            //
            double penaltyGradient = penaltyFunctionGrad(p, y, s, e, qe, q, rho);
            ArrayRealVector g = new ArrayRealVector(me + mi);
            if (me > 0) {
                g.setSubVector(0,equalityEval.subtract(getEqConstraint().getLowerBound()));
            }
            if (mi > 0) {
                g.setSubVector(me,inequalityEval.subtract(getIqConstraint().getLowerBound()).subtract(s));
            }

            double rhoSegnato = 2.0 * e.getNorm()/g.getNorm();

           // rho = rhoSegnato;
            //if (!(penaltyGradient<=-0.5 * p.dotProduct(hessian.operate(p))))
           while (!(penaltyGradient <= -0.5 * p.dotProduct(hessian.operate(p)))) {

                 rho = FastMath.max(rhoSegnato,2.0 * rho);

                penaltyGradient = penaltyFunctionGrad(p, y, s, e, qe, q, rho);
            }
            //LINE SEARCH
            double alfaEval = getObj().value(x.add(p.mapMultiply(alfa)));

            double alfaPenalty = 0;
            RealVector sineq = null;
            RealVector seq = null;
            if (se != null) {
                seq = se.add(qe.mapMultiply(alfa));
            }
            if (s != null) {
                sineq = s.add(q.mapMultiply(alfa));
            }

            double currentPenalty = penaltyFunction(functionEval,x,y,se,s,rho);
            alfaPenalty = penaltyFunction(alfaEval,x.add(p.mapMultiply(alfa)),y.add(e.mapMultiply(alfa)),seq,sineq,rho);



            int search = 0;
            while ((alfaPenalty -currentPenalty) >= getSettings().getMu() * alfa *  penaltyGradient &&
                   search < getSettings().getMaxLineSearchIteration()) {

                double alfaStar = -0.5 * alfa * alfa *  penaltyGradient/ (-alfa *  penaltyGradient + alfaPenalty - currentPenalty);

                alfa =  FastMath.max(getSettings().getB() * alfa, FastMath.min(1.0,alfaStar));
                alfaEval = getObj().value(x.add(p.mapMultiply(alfa)));
                if (se != null) {
                    seq = se.add(qe.mapMultiply(alfa));
                }
                if (s != null) {
                    sineq = s.add(q.mapMultiply(alfa));
                }
                alfaPenalty = penaltyFunction(alfaEval,x.add(p.mapMultiply(alfa)),y.add(e.mapMultiply(alfa)),seq,sineq,rho);

                search = search + 1;

            }


            if (getSettings().getConvCriteria() == 0) {
                if (p.mapMultiply(alfa).dotProduct(hessian.operate(p.mapMultiply(alfa))) < getSettings().getEps() * getSettings().getEps()) {
                    break;
                }
            } else {
                if (alfa * p.getNorm() <FastMath.sqrt(getSettings().getEps()) * (1 + x.getNorm())) {
                    break;
                }

            }


            //UPDATE ALL FUNCTION
            RealVector oldGradient = functionGradient;
            RealMatrix oldJacob = constraintJacob;
            RealVector old1 = lagrangianGradX(oldGradient,oldJacob,x,y.add(e.mapMultiply(alfa)));
            functionEval = alfaEval;
            functionGradient = getObj().gradient(x.add(p.mapMultiply(alfa)));
            constraintJacob = computeJacobianConstraint(x.add(p.mapMultiply(alfa)));
            RealVector new1 = lagrangianGradX(functionGradient,constraintJacob,x.add(p.mapMultiply(alfa)),y.add(e.mapMultiply(alfa)));
            hessian = BFGSFormula(hessian, p, alfa, new1, old1);

            if (getEqConstraint() != null) {
                equalityEval = getEqConstraint().value(x.add(p.mapMultiply(alfa)));
            }
            if (getIqConstraint() != null) {
                inequalityEval = getIqConstraint().value(x.add(p.mapMultiply(alfa)));
            }

            x = x.add(p.mapMultiply(alfa));
            y = y.add(e.mapMultiply(alfa));
        }

        functionEval = getObj().value(x);
        functionGradient = getObj().gradient(x);

        constraintJacob = computeJacobianConstraint(x);
        return new LagrangeSolution(x, y, functionEval);
    }

    private RealVector calculateSvector(RealVector y, double rho) {
        int me = 0;
        int mi = 0;
        RealVector si = null;
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
            si = new ArrayRealVector(mi);
            RealVector yi = y.getSubVector(me,mi);
            for (int i = 0; i < inequalityEval.getDimension(); i++) {
                if (rho == 0) {
                    si.setEntry(i, FastMath.max(0, inequalityEval.getEntry(i)));
                } else {
                    si.setEntry(i, FastMath.max(0, inequalityEval.getEntry(i) - yi.getEntry(i) / rho));
                }
            }
        }
        return si;
    }

    private double penaltyFunction(double currentF, RealVector x, RealVector y, RealVector se, RealVector s, double rho) {
        // the set of constraints is the same as the previous one but they must be evaluated with the increment

        int me = 0;
        int mi = 0;
        double partial = currentF;

        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();

            RealVector ye = y.getSubVector(0, me);
            RealVector g = getEqConstraint().value(x).subtract(getEqConstraint().getLowerBound());
            partial += -ye.dotProduct(g.subtract(se)) + 0.5 * rho*(g.subtract(se)).dotProduct(g.subtract(se));
        }

        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();

            RealVector yi = y.getSubVector(me, mi);


            RealVector g = getIqConstraint().value(x).subtract(getIqConstraint().getLowerBound());





            partial+= -yi.dotProduct(g.subtract(s)) +0.5 * rho*(g.subtract(s)).dotProduct(g.subtract(s));;

        }

        return partial;
    }

    private double penaltyFunctionGrad(RealVector p, RealVector y,RealVector s,RealVector e,RealVector qe, RealVector q,double rho) {

        int me = 0;
        int mi = 0;
        double partial = functionGradient.dotProduct(p);
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();

            RealVector ye = y.getSubVector(0, me);
            RealVector ee = e.getSubVector(0, me);
          //  RealVector qe = q.getSubVector(0, me);


            RealVector ge = equalityEval.subtract(getEqConstraint().getLowerBound());
            RealMatrix jacob = constraintJacob.getSubMatrix(0, me-1,0,p.getDimension()-1);
            double term2 = p.dotProduct(jacob.transpose().operate(ye));
            double term3 = p.dotProduct(jacob.transpose().operate(ge))*rho;
            double term4 = ge.dotProduct(ee);
            double term5 = ye.dotProduct(qe);
            double term6 = qe.dotProduct(ge)*rho;
            partial +=-term2 + term3-term4 + term5-term6;
        }

        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();

            RealVector yi = y.getSubVector(me, mi);
            RealVector ei = e.getSubVector(me, mi);
            RealVector qi = q;
            RealVector si = s;

            RealVector gi = inequalityEval.subtract(getIqConstraint().getLowerBound());
            RealMatrix jacob = constraintJacob.getSubMatrix(me,+me + mi-1,0,p.getDimension()-1);
            double term2 = p.dotProduct(jacob.transpose().operate(yi));
            double term3 = p.dotProduct(jacob.transpose().operate(gi.subtract(si)))*rho;
            double term4 = (gi.subtract(si)).dotProduct(ei);
            double term5 = yi.dotProduct(qi);
            double term6 = qi.dotProduct(gi.subtract(si))*rho;

            partial +=-term2 + term3-term4 + term5-term6;
        }

        return partial;
    }

    private LagrangeSolution solveQP(RealVector x, RealVector y) {

        RealMatrix H = hessian;
        RealVector g = functionGradient;

        int me = 0;
        int mi = 0;
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
        }
        if (getIqConstraint() != null) {

            mi = getIqConstraint().dimY();

        }

        RealMatrix H1 = new Array2DRowRealMatrix(H.getRowDimension() , H.getRowDimension() );
        H1.setSubMatrix(H.getData(), 0, 0);
        RealVector g1 = new ArrayRealVector(g.getDimension() );
        g1.setSubVector(0, g);

        LinearEqualityConstraint eqc = null;
        if (getEqConstraint() != null) {
            RealMatrix eqJacob = constraintJacob.getSubMatrix(0,me-1,0,x.getDimension()-1);

            RealMatrix Ae = new Array2DRowRealMatrix(me, x.getDimension() );
            RealVector be = new ArrayRealVector(me);
            Ae.setSubMatrix(eqJacob.getData(), 0, 0);



            be.setSubVector(0, getEqConstraint().getLowerBound().subtract(equalityEval));
            eqc = new LinearEqualityConstraint(Ae, be);

        }
        LinearInequalityConstraint iqc = null;

        if (getIqConstraint() != null) {
//
            RealMatrix iqJacob = constraintJacob.getSubMatrix(me,me + mi-1,0,x.getDimension()-1);

            RealMatrix Ai = new Array2DRowRealMatrix(mi, x.getDimension());
            RealVector bi = new ArrayRealVector(mi);
            Ai.setSubMatrix(iqJacob.getData(), 0, 0);

            bi.setSubVector(0, getIqConstraint().getLowerBound().subtract(inequalityEval));



            iqc = new LinearInequalityConstraint(Ai, bi);

        }


        QuadraticFunction q = new QuadraticFunction(H1, g1, 0);
        LagrangeSolution sol = null;

        ADMMQPOptimizer solver = new ADMMQPOptimizer();
        sol = solver.optimize(new ObjectiveFunction(q), iqc, eqc);

        LagrangeSolution solnew = new LagrangeSolution(sol.getX(), sol.getLambda(),0.0);
        return solnew;

    }

    private RealMatrix computeJacobianConstraint(RealVector x) {
        int me = 0;
        int mi = 0;
        RealMatrix je = null;
        RealMatrix ji = null;
        if (getEqConstraint() != null) {
            me = getEqConstraint().dimY();
            je = getEqConstraint().jacobian(x);
        }

        if (getIqConstraint() != null) {
            mi = getIqConstraint().dimY();
            ji = getIqConstraint().jacobian(x);
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

        RealMatrix oldH1 = oldH;
        RealVector y1 = newGrad.subtract(oldGrad);
        RealVector s = dx.mapMultiply(alfa);

        double theta = 1.0;
        if (s.dotProduct(y1) <= 0.2 * s.dotProduct(oldH.operate(s))) {
            theta = 0.8 * s.dotProduct(oldH.operate(s)) / (s.dotProduct(oldH.operate(s)) - s.dotProduct(y1));
        }

        RealVector y = y1.mapMultiply(theta).add(oldH.operate(s).mapMultiply(1.0 - theta));

        RealMatrix firstTerm = y.outerProduct(y).scalarMultiply(1.0 / s.dotProduct(y));
        RealMatrix secondTerm = oldH1.multiply(s.outerProduct(s)).multiply(oldH);
        double thirtTerm = s.dotProduct(oldH1.operate(s));
        RealMatrix Hnew = oldH1.add(firstTerm).subtract(secondTerm.scalarMultiply(1.0 / thirtTerm));
        //RESET HESSIAN IF NOT POSITIVE DEFINITE
        EigenDecompositionSymmetric dsX = new EigenDecompositionSymmetric(Hnew);
        double min = new ArrayRealVector(dsX.getEigenvalues()).getMinValue();
        if (new ArrayRealVector(dsX.getEigenvalues()).getMinValue() < 0) {

            RealMatrix diag = dsX.getD().
                              add(MatrixUtils.createRealIdentityMatrix(dx.getDimension()).
                                  scalarMultiply(getSettings().getEps() - min));
            Hnew = dsX.getV().multiply(diag.multiply(dsX.getVT()));

        }
        return Hnew;

    }

}
