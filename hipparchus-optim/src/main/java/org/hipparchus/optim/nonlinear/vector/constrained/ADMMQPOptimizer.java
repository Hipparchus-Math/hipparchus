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
import java.util.List;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.ConjugateGradient;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Alternating Direction Method of Multipliers Quadratic Programming Optimizer.
 * <br/>
 * min 0.5XT * Q*X + G*X * a
 * <br/>
 * A * X = b1
 * <br/>
 * B * X>=b2
 * <br/>
 * lb<=C * X<=ub
 * <br/>
 * Algorithm based on paper:"An Operator Splitting Solver for Quadratic Programs(Bartolomeo Stellato, Goran Banjac, Paul Goulart, Alberto Bemporad, Stephen Boyd,February 13 2020)"
 * @since 3.1
 */

public class ADMMQPOptimizer extends QPOptimizer {


    private double eps =ADMMQPOption.ADMM_eps;
    private double epsInf = ADMMQPOption.ADMM_epsInf;
    private double sigma =ADMMQPOption.ADMM_sigma;
    private double alfa =ADMMQPOption.ADMM_alfa;
    private boolean scale = ADMMQPOption.ADMM_scale;
    private int scale_maxIteration =ADMMQPOption.ADMM_scale_maxIteration;
    private boolean rho_update = ADMMQPOption.ADMM_rho_update;
    private double rho_max = ADMMQPOption.ADMM_rho_max;
    private double rho_min = ADMMQPOption.ADMM_rho_min;
    private int rho_maxIteration = ADMMQPOption.ADMM_rho_maxIteration;
    private boolean polish = ADMMQPOption.ADMM_polish;
    private int polish_iteration =ADMMQPOption.ADMM_polish_iteration;


    private LinearEqualityConstraint eqConstraint = null;
    private LinearInequalityConstraint iqConstraint = null;
    private LinearBoundedConstraint bqConstraint = null;
    private QuadraticFunction function;

    private RealVector xStart;
    //DEFAULT SOLVER
    private ADMMQPKKT1 ADMMSolver = new ADMMQPKKT1();
    ConjugateGradient indirect = new ConjugateGradient(100000,eps,true);
    private ADMMQPConvergenceChecker checker;
    private ADMMQPConvergenceChecker checkerRho;
    private boolean converged = false;
    private double rho = 0.1;

    @Override
   public ConvergenceChecker<LagrangeSolution> getConvergenceChecker() {
       return this.checker;
   }


    @Override
    public LagrangeSolution optimize(OptimizationData... optData) {
        return super.optimize(optData);
    }

    @Override
    protected void parseOptimizationData(OptimizationData... optData) {
        super.parseOptimizationData(optData);
        for (OptimizationData data: optData) {

             if (data instanceof ObjectiveFunction) {
                function = (QuadraticFunction) ((ObjectiveFunction)data).getObjectiveFunction();
                continue;
            }

            if (data instanceof LinearEqualityConstraint) {
                eqConstraint = (LinearEqualityConstraint)data;
                continue;
            }
            if (data instanceof LinearInequalityConstraint) {
                iqConstraint = (LinearInequalityConstraint)data;
                continue;
            }

            if (data instanceof LinearBoundedConstraint) {

                bqConstraint = (LinearBoundedConstraint)data;
                continue;
            }
             if (data instanceof ADMMQPOption) {
                alfa = ((ADMMQPOption) data).alfa;
                eps =((ADMMQPOption) data).eps;
                epsInf=((ADMMQPOption) data).epsInf;
                polish=((ADMMQPOption) data).polish;
                polish_iteration=((ADMMQPOption) data).polish_iteration;
                rho_max=((ADMMQPOption) data).rho_max;
                rho_maxIteration=((ADMMQPOption) data).rho_maxIteration;
                rho_min=((ADMMQPOption) data).rho_min;
                rho_update=((ADMMQPOption) data).rho_update;
                scale=((ADMMQPOption) data).scale;
                scale_maxIteration=((ADMMQPOption) data).scale_maxIteration;
                sigma=((ADMMQPOption) data).sigma;
                continue;

            }


           if (data instanceof InitialGuess) {
               this.xStart = new ArrayRealVector (((InitialGuess)data).getInitialGuess());

                continue;
            }


        }
        // if we got here, convexObjective exists
        int n = function.dim();
        if (eqConstraint != null) {
            int nDual = eqConstraint.dimY();
            if (nDual >= n) {
                throw new IllegalArgumentException("Rank of constraints must be < domain dimension");
            }
            int nTest = eqConstraint.A.getColumnDimension();
            if (nDual == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_NOT_ALLOWED);
            }
            MathUtils.checkDimension(nTest, n);
        }

    }

    @Override
    public LagrangeSolution doOptimize() {
        final int n = function.dim();
        int me = 0;
        int mi = 0;
        int mb = 0;
        int rhoupdateCount = 0;

        //PHASE 1 First Solution


       //QUADRATIC TERM
        RealMatrix H = this.function.getH();
       //GRADIENT
        RealVector q = this.function.getC();


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


        RealVector lb = new ArrayRealVector(me + mi + mb);
        RealVector ub = new ArrayRealVector(me + mi + mb);

      //COMPOSE A MATRIX AND LOWER AND UPPER BOUND
        RealMatrix A = new Array2DRowRealMatrix(me + mi + mb,n);
        if (this.eqConstraint!=null)
        {
            A.setSubMatrix(eqConstraint.jacobian(null).getData(), 0, 0);
            lb.setSubVector(0,eqConstraint.getLowerBound());
            ub.setSubVector(0,eqConstraint.getUpperBound());
        }
        if (this.iqConstraint!=null)
        {
            A.setSubMatrix(iqConstraint.jacobian(null).getData(), me, 0);
            ub.setSubVector(me,iqConstraint.getUpperBound());
            lb.setSubVector(me,iqConstraint.getLowerBound());
        }

         if (mb>0)
        {


            A.setSubMatrix(bqConstraint.jacobian(null).getData(), me + mi, 0);
            ub.setSubVector(me + mi,bqConstraint.getUpperBound());
            lb.setSubVector(me + mi,bqConstraint.getLowerBound());
        }
       this.getStartPoint();

        checker = new ADMMQPConvergenceChecker(H, A, q, eps, eps);







        //SETUP WORKING MATRIX
        RealMatrix Hw = H.copy();
        RealMatrix Aw = A.copy();
        RealVector qw = q.copy();
        RealVector ubw = ub.copy();
        RealVector lbw = lb.copy();
        RealVector x =null;
        if (this.getStartPoint() != null) {
            x = new ArrayRealVector(this.getStartPoint());
        } else {
            x = new ArrayRealVector(this.function.dim());
        }

        ADMMQPModifiedRuizEquilibrium dec = new ADMMQPModifiedRuizEquilibrium(H, A,q,lb,ub);

        if (scale) {
           //
            dec.normalize(eps,scale_maxIteration);
            Hw = dec.getScaledH();
            Aw = dec.getScaledA();
            qw = dec.getScaledq();
            lbw = dec.getScaledLUb(lb);
            ubw = dec.getScaledLUb(ub);

            x = dec.scaleX(x.copy());

        }

        checkerRho = new ADMMQPConvergenceChecker(Hw, Aw, qw, eps, eps);
        //SETUP VECTOR SOLUTION

        RealVector z = Aw.operate(x);
        RealVector y = new ArrayRealVector(me + mi + mb);

        this.ADMMSolver.initialize(Hw,Aw,qw,me,lbw,ubw,rho,sigma,alfa);
        RealVector xstar = null;
        RealVector ystar = null;
        RealVector zstar = null;

        while (this.iterations.getCount() <= this.iterations.getMaximalCount()) {
            ADMMQPSolution sol = this.ADMMSolver.iterate(x, y, z);
            x = sol.getX();
            y = sol.getLambda();
            z = sol.getZ();
            //new ArrayRealVector(me + mi + mb);
            if (rhoupdateCount < this.rho_maxIteration) {
                double rp = checkerRho.residualPrime(x, z);
                double rd=  checkerRho.residualDual(x, y);
                double maxP = checkerRho.maxPrimal(x, z);
                double maxD = checkerRho.maxDual(x, y);
                boolean updated = manageRho(me,rp,rd,maxP,maxD);

                if (updated) {
                    ++rhoupdateCount;
                }
            }


            if (scale) {

                xstar = dec.unscaleX(x);
                ystar = dec.unscaleY(y);
                zstar = dec.unscaleZ(z);

            } else {

                xstar = x.copy();
                ystar = y.copy();
                zstar = z.copy();

            }

            double rp        = ((ADMMQPConvergenceChecker) this.getConvergenceChecker()).residualPrime(xstar, zstar);
            double rd        = ((ADMMQPConvergenceChecker) this.getConvergenceChecker()).residualDual(xstar, ystar);
            double maxPrimal = ((ADMMQPConvergenceChecker) this.getConvergenceChecker()).maxPrimal(xstar, zstar);
            double maxDual   = ((ADMMQPConvergenceChecker) this.getConvergenceChecker()).maxDual(xstar, ystar);





            if (((ADMMQPConvergenceChecker)this.getConvergenceChecker()).converged(this.iterations.getCount(), rp, rd, maxPrimal, maxDual)) {
                converged = true;
                break;
            }
            this.iterations.increment();

        }



        //SOLUTION POLISHING

        ADMMQPSolution finalSol = null;
        if (polish) {
            finalSol = polish(Hw, Aw, qw, lbw, ubw, x, y, z);
            if (scale) {
                xstar = dec.unscaleX(finalSol.getX());
                ystar = dec.unscaleY(finalSol.getLambda());
                zstar = dec.unscaleZ(finalSol.getZ());
            } else {
                xstar = finalSol.getX();
                ystar = finalSol.getLambda();
                zstar = finalSol.getZ();
            }
        }
        for (int i = 0; i < me + mi; i++) {
            // ystar.setEntry(i,FastMath.abs(ystar.getEntry(i)));
            ystar.setEntry(i,-ystar.getEntry(i));
        }

        return new LagrangeSolution(xstar,ystar,function.value(xstar));

    }


   public boolean isConverged() {
        return this.converged;
   }



   private ADMMQPSolution polish(RealMatrix H,RealMatrix A,RealVector q,RealVector lb,RealVector ub, RealVector x, RealVector y, RealVector z) {

       List<double[]> Aentry = new ArrayList<double[]>();
       List<Double>  lubEntry = new ArrayList<Double>();

       List<Double>  yEntry = new ArrayList<Double>();
       List<Integer> indexLowerBound = new ArrayList<Integer>();
       List<Integer> indexBound = new ArrayList<Integer>();

       // FIND ACTIVE ON LOWER BAND
       for (int j = 0; j <A.getRowDimension(); j++) {
           if (z.getEntry(j) - lb.getEntry(j) < -y.getEntry(j)) {  // lower-active

               indexBound.add(j);
               Aentry.add(A.getRow(j));
               lubEntry.add(lb.getEntry(j));
               yEntry.add(y.getEntry(j));

           }
       }
       //FIND ACTIVE ON UPPER BAND
       for (int j = 0; j <A.getRowDimension(); j++) {
           if (-z.getEntry(j) + ub.getEntry(j) < y.getEntry(j)) { // lower-active

               Aentry.add(A.getRow(j));
               lubEntry.add(ub.getEntry(j));
               yEntry.add(y.getEntry(j));
               indexBound.add(j);

           }

       }
       RealMatrix Aactive = null;
       RealVector lub = null;

       RealVector ystar;
       RealVector xstar = x.copy();
       //!Aentry.isEmpty()
       if (!Aentry.isEmpty()) {

           Aactive = new Array2DRowRealMatrix(Aentry.toArray(new double[Aentry.size()][]));
           lub = new ArrayRealVector(lubEntry.toArray(new Double[lubEntry.size()]));
           ystar = new ArrayRealVector(yEntry.toArray(new Double[yEntry.size()]));
           this.ADMMSolver.initialize(H, Aactive,q,0, lub,lub,sigma,sigma,alfa);
           // this.ADMMSolver.initialize(H, Aactive, MatrixUtils.createRealIdentityMatrix(Aactive.getRowDimension()).scalarMultiply(this.sigmaPolishing), this.sigmaPolishing);

           for (int i = 0; i < polish_iteration; i++) {
               RealVector kttx = (H.operate(xstar)).add(Aactive.transpose().operate(ystar));
               RealVector ktty = Aactive.operate(xstar);
               RealVector b1 = q.mapMultiply(-1.0).subtract(kttx);
               RealVector b2 = lub.mapMultiply(1.0).subtract(ktty);
               ADMMQPSolution dxz = this.ADMMSolver.solve(b1,b2);
               xstar = xstar.add(dxz.getX());
               ystar = ystar.add(dxz.getV());


           }
           //         for(int i = 0;i<indexBound.size();i++)
           //             y.setEntry(indexBound.get(i), ystar.getEntry(i));

           return new ADMMQPSolution(xstar,null,y,A.operate(xstar));
       } else {
           return new ADMMQPSolution(x,null,y,z);
       }
   }

   private double rhoEstimation(double rp,double rd,double rpMax,double rdMax) {
       // Return rho estimate
       double rhoEstimate = rho * FastMath.sqrt((rp * rdMax)/(rd * rpMax));
       rhoEstimate = FastMath.min(FastMath.max(rhoEstimate, rho_min),rho_max);
       return rhoEstimate;
   }// Constrain

   private boolean manageRho(int me,double rp, double rd, double maxPrimal, double maxDual) {
       boolean updated = false;
       if (rho_update)
       {

           double rhonew = rhoEstimation(rp,rd,maxPrimal,maxDual);
           if ((rhonew>rho * 5.0) || (rhonew<rho/5.0))
           {

               this.rho = rhonew;
               updated = true;

               this.ADMMSolver.updateSigmaRho(this.sigma,me,this.rho);
           }
       }
       return updated;
   }

   private RealMatrix createPenaltyMatrix(RealMatrix A, int me, double rho) {
       RealMatrix R;
       R = MatrixUtils.createRealIdentityMatrix(A.getRowDimension());

       for (int i = 0; i < R.getRowDimension(); i++) {
           if (i < me) {
               R.setEntry(i, i, rho * 1000.0);
           } else {
               R.setEntry(i, i, rho);
           }
       }
       return R;
   }

}
