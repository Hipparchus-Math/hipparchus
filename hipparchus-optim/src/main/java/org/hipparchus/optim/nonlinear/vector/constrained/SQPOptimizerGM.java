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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.EigenDecompositionSymmetric;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
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
 * Algorithm based on paper:"Some Theoretical properties of an augmented lagrangiane merit function
 * (Gill,Murray,Sauders,Wriht,April 1986)"
 * @since 3.1
 */
public class SQPOptimizerGM extends ConstraintOptimizer {

    int forgetFactor = 10;

    int convCriteria = SQPOption.SQP_ConvCriteria;
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
    private int maxLineSearchIteration = SQPOption.SQP_maxLineSearchIteration;

    private TwiceDifferentiableFunction obj;
    private EqualityConstraint eqConstraint;
    private InequalityConstraint iqConstraint;
    private BoundedConstraint bqConstraint;
    private ArrayRealVector xStart;

    private RealMatrix constraintJacob;
    private RealVector equalityEval;
    private RealVector inequalityEval;
    private double functionEval;
    private RealVector functionGradient;
    private RealMatrix hessian;


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

            if (data instanceof EqualityConstraint) {
                eqConstraint = (EqualityConstraint) data;
                continue;
            }
            if (data instanceof InequalityConstraint) {
                iqConstraint = (InequalityConstraint) data;
                continue;
            }

            if (data instanceof BoundedConstraint) {
                bqConstraint = (BoundedConstraint) data;
                continue;
            }

            if (data instanceof QPOptimizer) {
                qpSolver = (QPOptimizer) data;
                continue;

            }

            if (data instanceof SQPOption) {
                convCriteria = ((SQPOption) data).convCriteria;
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
        double rho = 100.0;
        int failedSearch = 0;
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
         //INITIAL VALUES
        functionEval = this.obj.value(x);
        functionGradient = this.obj.gradient(x);
        double maxGrad = functionGradient.getLInfNorm();





        if (this.eqConstraint!=null)
        {

          equalityEval = this.eqConstraint.value(x);
        }
        if (this.iqConstraint!=null)
        {

            inequalityEval = this.iqConstraint.value(x);
        }
        constraintJacob = computeJacobianConstraint(x);

         if (this.eqConstraint != null) {
            maxGrad = FastMath.max(maxGrad, equalityEval.getLInfNorm());
        }
        if (this.iqConstraint != null) {
            maxGrad = FastMath.max(maxGrad, inequalityEval.getLInfNorm());
        }

          if (useFunHessian == false) {
            hessian= MatrixUtils.createRealIdentityMatrix(x.getDimension()).scalarMultiply(maxGrad);
        } else {
            hessian = this.obj.hessian(x);
        }


        rho = 0;

        for (int i = 0; i < this.getMaxIterations(); i++) {
            RealVector dx1 = new ArrayRealVector(x.getDimension());

            iterations.increment();


            alfa = 1.0;

            LagrangeSolution sol1 = null;
            //SOLVE QP
            sol1 = solveQP(x,y);
            RealVector p = sol1.getX();
            RealVector e = sol1.getLambda().subtract(y);
//           if (p.getNorm()<epsilon) break;

            // RealVector q = sol1.getX().getSubVector(x.getDimension()+2 * me, mi);
             RealVector s = calculateSvector(y,rho);
             RealVector se = null;
             RealMatrix jacobi = null;
             RealVector q = null;
             RealVector qe = null;

            //TEST CON SI SOLO PER INEQUALITY AND Q FOR ALL
            RealVector appoggio = new ArrayRealVector(me + mi);
             if (eqConstraint!=null)
            {
            se = new ArrayRealVector(me);
            jacobi = constraintJacob.getSubMatrix(0,me-1,0,x.getDimension()-1);
            qe = new ArrayRealVector(me);
            //qe = jacobi.operate(p).add(equalityEval.subtract(eqConstraint.getLowerBound()));

            }
            if (iqConstraint!=null)
            {
            jacobi = constraintJacob.getSubMatrix(me,me + mi-1,0,x.getDimension()-1);
            q = jacobi.operate(p).add(inequalityEval.subtract(iqConstraint.getLowerBound())).subtract(s);
            }



            //CALCULATE PENALTY GRADIENT
            //
            double penaltyGradient = penaltyFunctionGrad(p,y,s,e,qe,q,rho);
            ArrayRealVector g = new ArrayRealVector(me + mi);
            if (me > 0) {
                g.setSubVector(0,equalityEval.subtract(this.eqConstraint.getLowerBound()));
            }
            if (mi > 0) {
                g.setSubVector(me,inequalityEval.subtract(this.iqConstraint.getLowerBound()).subtract(s));
            }

            double rhoSegnato = 2.0 * e.getNorm()/g.getNorm();

           // rho = rhoSegnato;
            //if (!(penaltyGradient<=-0.5 * p.dotProduct(hessian.operate(p))))
           while(!(penaltyGradient<=-0.5 * p.dotProduct(hessian.operate(p))))
            {

                 rho = FastMath.max(rhoSegnato,2.0 * rho);

                penaltyGradient = penaltyFunctionGrad(p, y,s,e,qe, q,rho);
            }
            //LINE SEARCH
            double alfaEval = this.obj.value(x.add(p.mapMultiply(alfa)));

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
               while ((alfaPenalty -currentPenalty) >= this.mu * alfa *  penaltyGradient && search < maxLineSearchIteration) {

                    double alfaStar = -0.5 * alfa * alfa *  penaltyGradient/ (-alfa *  penaltyGradient + alfaPenalty - currentPenalty);
//

                   alfa =  FastMath.max(this.b * alfa, FastMath.min(1.0,alfaStar));
                 // alfa = alfa * 0.5;
                    // alfa = FastMath.min(1.0, FastMath.max(this.b * alfa, alfaStar));
                   // alfa = FastMath.max(this.b * alfa, alfaStar);
                    alfaEval = this.obj.value(x.add(p.mapMultiply(alfa)));
                    if (se != null) {
                        seq = se.add(qe.mapMultiply(alfa));
                    }
                    if (s != null) {
                        sineq = s.add(q.mapMultiply(alfa));
                    }
                    alfaPenalty = penaltyFunction(alfaEval,x.add(p.mapMultiply(alfa)),y.add(e.mapMultiply(alfa)),seq,sineq,rho);

                    search = search + 1;

                }


            if (convCriteria == 0) {
                if (p.mapMultiply(alfa).dotProduct(hessian.operate(p.mapMultiply(alfa))) < epsilon * epsilon)
                {
//                    x = x.add(dx.mapMultiply(alfa));
//                    y = y.add((dy.subtract(y)).mapMultiply(alfa));
                    break;
                }
            } else {
                if (alfa * p.getNorm() <FastMath.sqrt(epsilon) * (1 + x.getNorm()))
                {
//                    x = x.add(dx.mapMultiply(alfa));
//                    y = y.add((dy.subtract(y)).mapMultiply(alfa));
                    break;
                }

            }


            //UPDATE ALL FUNCTION
            double oldFunction = functionEval;
            RealVector oldGradient = functionGradient;
            RealMatrix oldJacob = constraintJacob;
            RealVector old1 = LagrangianeGradX(oldGradient,oldJacob,x,y.add(e.mapMultiply(alfa)),rho);
            functionEval = alfaEval;
            functionGradient = this.obj.gradient(x.add(p.mapMultiply(alfa)));
            constraintJacob = computeJacobianConstraint(x.add(p.mapMultiply(alfa)));
            RealVector new1 = LagrangianeGradX(functionGradient,constraintJacob,x.add(p.mapMultiply(alfa)),y.add(e.mapMultiply(alfa)),rho);
            hessian = BFGSFormula(hessian, p, alfa, new1, old1);

        if (this.eqConstraint!=null)
        {

          equalityEval = this.eqConstraint.value(x.add(p.mapMultiply(alfa)));
        }
        if (this.iqConstraint!=null)
        {

            inequalityEval = this.iqConstraint.value(x.add(p.mapMultiply(alfa)));
        }


        x = x.add(p.mapMultiply(alfa));
        y = y.add(e.mapMultiply(alfa));
        }




        functionEval = this.obj.value(x);
        functionGradient = this.obj.gradient(x);

        constraintJacob = computeJacobianConstraint(x);
        double constraintCheck = constraintCheck(x);

        double dlagrange = LagrangianeGradX(functionGradient,constraintJacob,x, y,rho).getNorm();
        return new LagrangeSolution(x, y, functionEval);
    }

    private RealVector calculateSvector(RealVector y, double rho) {
        int me = 0;
        int mi = 0;
        RealVector si = null;
        if (this.eqConstraint != null) {
            me = this.eqConstraint.dimY();
        }
        if (this.iqConstraint!=null) {
            mi = this.iqConstraint.dimY();
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

        if (eqConstraint != null) {
            me = eqConstraint.dimY();

            RealVector ye = y.getSubVector(0, me);
            RealVector g = this.eqConstraint.value(x).subtract(eqConstraint.getLowerBound());
            RealVector g2 = g.ebeMultiply(g);
            partial += -ye.dotProduct(g.subtract(se)) + 0.5 * rho*(g.subtract(se)).dotProduct(g.subtract(se));
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();

            RealVector yi = y.getSubVector(me, mi);


            RealVector g = this.iqConstraint.value(x).subtract(iqConstraint.getLowerBound());





            partial+= -yi.dotProduct(g.subtract(s)) +0.5 * rho*(g.subtract(s)).dotProduct(g.subtract(s));;

        }

        return partial;
    }

    private double penaltyFunctionGrad(RealVector p, RealVector y,RealVector s,RealVector e,RealVector qe, RealVector q,double rho) {

        int me = 0;
        int mi = 0;
        double partial = functionGradient.dotProduct(p);
        if (eqConstraint != null) {
            me = eqConstraint.dimY();

            RealVector ye = y.getSubVector(0, me);
            RealVector ee = e.getSubVector(0, me);
          //  RealVector qe = q.getSubVector(0, me);


            RealVector ge = this.equalityEval.subtract(eqConstraint.getLowerBound());
            RealMatrix jacob = this.constraintJacob.getSubMatrix(0, me-1,0,p.getDimension()-1);
            double term2 = p.dotProduct(jacob.transpose().operate(ye));
            double term3 = p.dotProduct(jacob.transpose().operate(ge))*rho;
            double term4 = ge.dotProduct(ee);
            double term5 = ye.dotProduct(qe);
            double term6 = qe.dotProduct(ge)*rho;




            partial +=-term2 + term3-term4 + term5-term6;
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();


            RealVector yi = y.getSubVector(me, mi);
            RealVector ei = e.getSubVector(me, mi);
            RealVector qi = q;
            RealVector si = s;

            RealVector gi = this.inequalityEval.subtract(iqConstraint.getLowerBound());
            RealMatrix jacob = this.constraintJacob.getSubMatrix(me,+me + mi-1,0,p.getDimension()-1);
            double term2 = p.dotProduct(jacob.transpose().operate(yi));
            double term3 = p.dotProduct(jacob.transpose().operate(gi.subtract(si)))*rho;
            double term4 = (gi.subtract(si)).dotProduct(ei);
            double term5 = yi.dotProduct(qi);
            double term6 = qi.dotProduct(gi.subtract(si))*rho;




            partial +=-term2 + term3-term4 + term5-term6;
        }

        return partial;
    }

    private RealVector LagrangianeGradX(RealVector currentGrad,RealMatrix jacobConstraint,RealVector x, RealVector y,double rho) {

        int me = 0;
        int mi = 0;
        RealVector partial = currentGrad.copy();
        if (eqConstraint != null) {
            me = eqConstraint.dimY();

            RealVector ye = y.getSubVector(0, me);
            RealVector ge = this.equalityEval.subtract(eqConstraint.getLowerBound());
            RealMatrix jacobe = jacobConstraint.getSubMatrix(0, me-1,0,x.getDimension()-1);

            RealVector firstTerm = (jacobe.transpose().operate(ye));

           // partial = partial.subtract(firstTerm).add(jacobe.transpose().operate(ge).mapMultiply(rho));
            partial = partial.subtract(firstTerm);
        }

        if (iqConstraint != null) {
            mi = iqConstraint.dimY();

            RealVector yi = y.getSubVector(me, mi);
           RealMatrix jacobi = jacobConstraint.getSubMatrix(me,me + mi-1,0,x.getDimension()-1);
           RealVector gi = this.inequalityEval.subtract(iqConstraint.getLowerBound());

            RealVector firstTerm = (jacobi.transpose().operate(yi));

           // partial = partial.subtract(firstTerm).add(jacobi.transpose().operate(gi).mapMultiply(rho));;
             partial = partial.subtract(firstTerm);
        }
        return partial;
    }










    private LagrangeSolution solveQP(RealVector x, RealVector y) {

        RealMatrix H = hessian;
        RealVector g = functionGradient;

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

        }



        RealMatrix H1 = new Array2DRowRealMatrix(H.getRowDimension() , H.getRowDimension() );
        H1.setSubMatrix(H.getData(), 0, 0);
        RealVector g1 = new ArrayRealVector(g.getDimension() );
        g1.setSubVector(0, g);

        LinearEqualityConstraint eqc = null;
        RealVector conditioneq = null;
        if (eqConstraint != null) {
            RealMatrix eqJacob = constraintJacob.getSubMatrix(0,me-1,0,x.getDimension()-1);

            RealMatrix Ae = new Array2DRowRealMatrix(me, x.getDimension() );
            RealVector be = new ArrayRealVector(me);
            Ae.setSubMatrix(eqJacob.getData(), 0, 0);



            be.setSubVector(0, eqConstraint.getLowerBound().subtract(this.equalityEval));
            eqc = new LinearEqualityConstraint(Ae, be);

        }
        LinearInequalityConstraint iqc = null;

        if (iqConstraint != null) {
//
            RealMatrix iqJacob = constraintJacob.getSubMatrix(me,me + mi-1,0,x.getDimension()-1);

            RealMatrix Ai = new Array2DRowRealMatrix(mi, x.getDimension());
            RealVector bi = new ArrayRealVector(mi);
            Ai.setSubMatrix(iqJacob.getData(), 0, 0);

            bi.setSubVector(0, iqConstraint.getLowerBound().subtract(this.inequalityEval));



            iqc = new LinearInequalityConstraint(Ai, bi);

        }


        QuadraticFunction q = new QuadraticFunction(H1, g1, 0);
        LagrangeSolution sol = null;
        double sigma = 0;

        ADMMQPOptimizer solver = new ADMMQPOptimizer();
        sol = solver.optimize(new ObjectiveFunction(q), iqc, eqc);

        LagrangeSolution solnew = new LagrangeSolution(sol.getX(), sol.getLambda(),0.0);
        return solnew;

    }

    private LagrangeSolution solveQP1(RealVector x, RealVector y) {

        RealMatrix H = hessian;
        RealVector g = functionGradient;

        int me = 0;
        int mi = 0;
        int mb = 0;
        int add = 0;

        if (eqConstraint != null) {
            me = eqConstraint.dimY();
        }
        if (iqConstraint != null) {

            mi = iqConstraint.dimY();

        }



        RealMatrix H1 = new Array2DRowRealMatrix(H.getRowDimension() + mi + me * 2, H.getRowDimension() + mi + me * 2);
        H1.setSubMatrix(H.getData(), 0, 0);
        RealVector g1 = new ArrayRealVector(g.getDimension() + mi + me * 2);
        g1.setSubVector(0, g);

        LinearEqualityConstraint eqc = null;
        RealVector conditioneq = null;
        if (eqConstraint != null) {
            RealMatrix eqJacob = constraintJacob.getSubMatrix(0,me-1,0,x.getDimension()-1);

            RealMatrix Ae = new Array2DRowRealMatrix(me, x.getDimension() +mi + me * 2);
            RealVector be = new ArrayRealVector(me);
            Ae.setSubMatrix(eqJacob.getData(), 0, 0);
            Ae.setSubMatrix(MatrixUtils.createRealIdentityMatrix(me).getData(), 0, x.getDimension());
            Ae.setSubMatrix(MatrixUtils.createRealIdentityMatrix(me).scalarMultiply(-1.0).getData(), 0, x.getDimension()+me);
            conditioneq = this.equalityEval.subtract(eqConstraint.getLowerBound());


            be.setSubVector(0, eqConstraint.getLowerBound().subtract(this.equalityEval));
            eqc = new LinearEqualityConstraint(Ae, be);

        }
        LinearInequalityConstraint iqc = null;

        if (iqConstraint != null && me>0) {
//
            RealMatrix iqJacob = constraintJacob.getSubMatrix(me,me + mi-1,0,x.getDimension()-1);

            RealMatrix Ai = new Array2DRowRealMatrix(2 * me + 2*mi, x.getDimension()+2 * me + mi);
            RealVector bi = new ArrayRealVector(2 * mi + 2*me);
            Ai.setSubMatrix(iqJacob.getData(), 0, 0);
            Ai.setSubMatrix(MatrixUtils.createRealIdentityMatrix(mi).scalarMultiply(-1).getData(),0, x.getDimension()+2 * me);
            Ai.setSubMatrix(MatrixUtils.createRealIdentityMatrix(me).getData(),mi, x.getDimension());
            Ai.setSubMatrix(MatrixUtils.createRealIdentityMatrix(me).getData(),mi + me, x.getDimension()+me);
             Ai.setSubMatrix(MatrixUtils.createRealIdentityMatrix(mi).getData(),mi + 2*me, x.getDimension()+2 * me);
//            Ai.setSubMatrix(MatrixUtils.createRealIdentityMatrix(mi).getData(),mi, x.getDimension());
            bi.setSubVector(0, iqConstraint.getLowerBound().subtract(this.inequalityEval));



            iqc = new LinearInequalityConstraint(Ai, bi);

        }


        QuadraticFunction q = new QuadraticFunction(H1, g1, 0);
        LagrangeSolution sol = null;
        double sigma = 0;

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
    /*
     *DAMPED BFGS FORMULA
     */

    private RealMatrix BFGSFormula(RealMatrix oldH, RealVector dx, double alfa, RealVector newGrad, RealVector oldGrad) {

        RealMatrix oldH1 = oldH;
        RealVector y1 = newGrad.subtract(oldGrad);
        RealVector s = dx.mapMultiply(alfa);
//        if (s.dotProduct(y1)>=0)
//        {
        double theta = 1.0;
        if (s.dotProduct(y1) <= 0.2 * s.dotProduct(oldH.operate(s)) ) {

            theta =(0.8 * s.dotProduct(oldH.operate(s)) / (s.dotProduct(oldH.operate(s)) - s.dotProduct(y1)));
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

            RealMatrix diag = dsX.getD().add(MatrixUtils.createRealIdentityMatrix(dx.getDimension()).scalarMultiply((-1.0 * min + epsilon)));
            Hnew = dsX.getV().multiply(diag.multiply(dsX.getVT()));

        }
        return Hnew;
//        }
//        else return oldH;

    }

    private double constraintCheck(RealVector x) {
        // the set of constraints is the same as the previous one but they must be evaluated with the increment

        int me = 0;
        int mi = 0;
        double partial = 0;

        if (eqConstraint != null) {
            me = eqConstraint.dimY();

            RealVector g = this.eqConstraint.value(x).subtract(eqConstraint.getLowerBound());
            for (int i = 0; i < g.getDimension(); i++) {

                partial += FastMath.abs(g.getEntry(i));
            }
        }

        if (iqConstraint != null) {

            RealVector g = this.iqConstraint.value(x).subtract(iqConstraint.getLowerBound());

            for (int i = 0; i < g.getDimension(); i++) {

                partial += FastMath.abs(FastMath.min(g.getEntry(i), 0));
            }

        }

        return partial;
    }

}
