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
import org.hipparchus.linear.DecompositionSolver;
import org.hipparchus.linear.EigenDecompositionSymmetric;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.SingularValueDecomposition;

/** Alternative Direction Method of Multipliers Solver.
 * @since 3.1
 */

public class ADMMQPKKT1 extends KKTSolver<ADMMQPSolution, RealMatrix, RealVector> {

    private DecompositionSolver dsX;
    private RealMatrix H;
    private RealMatrix A;
    private RealMatrix R;
    private double sigma;
    private RealVector q;
    private RealMatrix Rinv;
    private RealVector lb;
    private RealVector ub;
    private double rho;
    private double alfa;
    private Array2DRowRealMatrix M;

    public ADMMQPKKT1(RealMatrix H, RealMatrix A, RealMatrix R, double sigma) {
        this.H = H.add(MatrixUtils.createRealIdentityMatrix(H.getColumnDimension()).scalarMultiply(sigma));
        this.A = A.copy();
        this.R = R.copy();
        this.sigma = sigma;
        M = new Array2DRowRealMatrix(H.getRowDimension()+A.transpose().getRowDimension(),H.getRowDimension()+A.getRowDimension());
        M.setSubMatrix(H.getData(),0,0);
        M.setSubMatrix(A.getData(),H.getRowDimension(),0);
        M.setSubMatrix(A.transpose().getData(),0,H.getRowDimension());
        M.setSubMatrix(MatrixUtils.inverse(R).scalarMultiply(-1.0).getData(),H.getRowDimension(),H.getRowDimension());
       // dsX=(new EigenDecomposition(M)).getSolver();
        dsX=(new SingularValueDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
       // dsX=(new CholeskyDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
       // dsX=(new EigenDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
    }

    ADMMQPKKT1() {

    }

    @Override
    public ADMMQPSolution solve(RealVector b1, final RealVector b2) {
        RealVector z = dsX.solve(new ArrayRealVector((ArrayRealVector) b1,b2));
//        RealVector z1 = dsX.solve(b1.add((A.transpose().operate(R.operate(b2)))));
//        RealVector z2 = R.operate((A.operate(z1)).subtract(b2));
        return new ADMMQPSolution(z.getSubVector(0,b1.getDimension()), z.getSubVector(b1.getDimension(), b2.getDimension()));

    }

    @Override
    public RealMatrix getKKTMatrix(RealMatrix H, RealMatrix A, RealMatrix R) {

        RealMatrix KKT = MatrixUtils.createRealMatrix(H.getRowDimension() + A.getRowDimension(), H.getRowDimension() + A.getRowDimension());
        KKT.setSubMatrix(H.getData(), 0, 0);
        KKT.setSubMatrix(A.getData(), H.getRowDimension(), 0);
        KKT.setSubMatrix(A.transpose().getData(), 0, H.getColumnDimension());
        KKT.setSubMatrix(R.getData(), H.getRowDimension(), H.getColumnDimension());
        return KKT;
    }


    public void updateSigmaRho(double sigma, int me, double rho) {
        this.sigma = sigma;
        this.H = H.add(MatrixUtils.createRealIdentityMatrix(H.getColumnDimension()).scalarMultiply(sigma));
        createPenaltyMatrix(me, rho);
        M = new Array2DRowRealMatrix(H.getRowDimension()+A.getRowDimension(),H.getRowDimension()+A.getRowDimension());
        M.setSubMatrix(H.getData(),0,0);
        M.setSubMatrix(A.getData(),H.getRowDimension(),0);
        M.setSubMatrix(A.transpose().getData(),0,H.getRowDimension());
        M.setSubMatrix(Rinv.scalarMultiply(-1.0).getData(),H.getRowDimension(),H.getRowDimension());
        dsX=(new EigenDecompositionSymmetric(M)).getSolver();
        // dsX=(new CholeskyDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
      // dsX=(new SingularValueDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
         //dsX=(new CholeskyDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
      // dsX=(new EigenDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();

    }


    public void initialize(RealMatrix H, RealMatrix A, RealVector q, int me, RealVector lb, RealVector ub, double rho, double sigma, double alfa) {
        this.lb = lb;
        this.ub = ub;
        this.rho = rho;
        this.alfa = alfa;
         this.sigma = sigma;
        this.H = H.add(MatrixUtils.createRealIdentityMatrix(H.getColumnDimension()).scalarMultiply(sigma));
        this.A = A.copy();
        this.q = q.copy();
        createPenaltyMatrix(me, rho);

        M = new Array2DRowRealMatrix(H.getRowDimension()+A.getRowDimension(),H.getRowDimension()+A.getRowDimension());
        M.setSubMatrix(H.getData(),0,0);
        M.setSubMatrix(A.getData(),H.getRowDimension(),0);
        M.setSubMatrix(A.transpose().getData(),0,H.getRowDimension());
        M.setSubMatrix(Rinv.scalarMultiply(-1.0).getData(),H.getRowDimension(),H.getRowDimension());
      dsX=(new EigenDecompositionSymmetric(M)).getSolver();

     //dsX=(new QRDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
    // dsX=(new SingularValueDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
    //    dsX=(new CholeskyDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
       // dsX=(new EigenDecomposition(H.add(A.transpose().multiply(R.multiply(A))))).getSolver();
    }


    private void createPenaltyMatrix(int me, double rho) {
        this.R = MatrixUtils.createRealIdentityMatrix(A.getRowDimension());

        for (int i = 0; i < R.getRowDimension(); i++) {
            if (i < me) {
                R.setEntry(i, i, rho * 1000.0);

            } else {
                R.setEntry(i, i, rho);

            }
        }
        this.Rinv = MatrixUtils.inverse(R);
    }

    @Override
    public ADMMQPSolution iterate(RealVector... previousSol) {
        double onealfa = 1.0 - alfa;
        //SAVE OLD VALUE
        RealVector xold = previousSol[0].copy();
        RealVector yold = previousSol[1].copy();
        RealVector zold = previousSol[2].copy();

        //UPDATE RIGHT VECTOR
        RealVector b1 = (previousSol[0].mapMultiply(sigma)).subtract(q);
        RealVector b2 = previousSol[2].subtract(Rinv.operate(previousSol[1]));

        //SOLVE KKT SYSYEM
        ADMMQPSolution sol = this.solve(b1, b2);
        RealVector xtilde = sol.getX();
        RealVector vtilde = sol.getV();

        //UPDATE ZTILDE
        RealVector ztilde = zold.add(Rinv.operate(vtilde.subtract(yold)));
        //UPDATE X
        previousSol[0] = (xtilde.mapMultiply(alfa)).add((xold.mapMultiply(onealfa)));

        //UPDATE Z PARTIAL
        RealVector zpartial = (ztilde.mapMultiply(alfa)).add(zold.mapMultiply(onealfa)).add(Rinv.operate(yold));

        //PROJECT ZPARTIAL AND UPDATE Z
        for (int j = 0; j < previousSol[2].getDimension(); j++) {
            previousSol[2].setEntry(j, Math.min(Math.max(zpartial.getEntry(j), lb.getEntry(j)), ub.getEntry(j)));
        }

        //UPDATE Y
        RealVector ytilde = ztilde.mapMultiply(alfa).add(zold.mapMultiply(onealfa).subtract(previousSol[2]));
        previousSol[1] = yold.add(R.operate(ytilde));

        return new ADMMQPSolution(previousSol[0], vtilde, previousSol[1], previousSol[2]);
    }

}
