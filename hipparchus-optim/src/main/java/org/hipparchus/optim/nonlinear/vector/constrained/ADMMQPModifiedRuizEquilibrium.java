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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;

public class ADMMQPModifiedRuizEquilibrium {
    private double MIN_SCALING = 1.0e6; ///< minimum scaling value
    private double MAX_SCALING = 1.0e+6; ///< maximum scaling value
    final RealMatrix H;
    final RealMatrix A;
    final RealVector q;
    final RealVector lb;
    final RealVector ub;

    public RealMatrix D;
    public RealMatrix E;
    public double c;
    public RealMatrix Dinv;
    public RealMatrix Einv;
    public double cinv;

    public ADMMQPModifiedRuizEquilibrium(RealMatrix H, RealMatrix A, RealVector q, RealVector lb, RealVector ub) {
        this.H  = H;
        this.A  = A;
        this.q  = q;
        this.lb = lb;
        this.ub = ub;
    }

    public void normalize(double epsilon, int maxIteration) {

        int iteration = 0;
        double c = 1.0;
        RealVector gamma = new ArrayRealVector(H.getRowDimension() + A.getRowDimension());
        RealVector gammaD = new ArrayRealVector(H.getRowDimension());
        RealVector gammaE = new ArrayRealVector(A.getRowDimension());

        double lambda;
        RealVector q1 = q.copy();
        RealMatrix H1 = H.copy();
        RealMatrix A1 = A.copy();
        RealVector lb1 = lb.copy();
        RealMatrix D = null;
        RealMatrix E = null;
        RealMatrix SD = MatrixUtils.createRealIdentityMatrix(H.getRowDimension());
        RealMatrix SE = MatrixUtils.createRealIdentityMatrix(A.getRowDimension());
        RealVector H1norm = new ArrayRealVector(H1.getColumnDimension());
        while (iteration < maxIteration) {
       // while (1.0 - gamma.getLInfNorm() > epsilon && iteration < maxIteration) {

            for (int i = 0; i < H1.getColumnDimension(); i++) {
                double norma = (new ArrayRealVector(H1.getColumn(i), A1.getColumn(i))).getLInfNorm();
                gammaD.setEntry(i, norma);

            }


            for (int i = 0; i < A1.getRowDimension(); i++) {
                double norma = A1.getRowVector(i).getLInfNorm();
                gammaE.setEntry(i, norma);
            }

            gammaD = limit_scaling(gammaD);
            gammaE = limit_scaling(gammaE);

            for (int i = 0; i < gammaD.getDimension(); i++) {

                gammaD.setEntry(i, 1.0 / FastMath.sqrt(gammaD.getEntry(i)));

            }
            for (int i = 0; i < gammaE.getDimension(); i++) {

              gammaE.setEntry(i, 1.0 / FastMath.sqrt(gammaE.getEntry(i)));
            }



            D = MatrixUtils.createRealDiagonalMatrix(gammaD.toArray());
            E = MatrixUtils.createRealDiagonalMatrix(gammaE.toArray());



            H1 = D.multiply(H1.copy()).multiply(D.copy());
            q1 = D.operate(q1.copy());
            A1 = E.multiply(A1.copy()).multiply(D.copy());
            lb1 = E.operate(lb1.copy());

            //
            for (int i = 0; i < H1.getRowDimension(); i++) {
              //  double norma = H1.getColumnVector(i).getLInfNorm();
                double norma = (new ArrayRealVector(H1.getRow(i))).getLInfNorm();
                H1norm.setEntry(i, norma);
            }
//            for (int i = 0; i < H1.getColumnDimension(); i++) {
//              //  double norma = H1.getColumnVector(i).getLInfNorm();
//                double norma = (new ArrayRealVector(H1.getColumn(i),A1.getColumn(i))).getLInfNorm();
//                H1norm.setEntry(i, norma);
//            }
            double qnorm = q1.getLInfNorm();
            double lbnorm = lb1.getLInfNorm();
            double qnorm1 = FastMath.max(qnorm,lbnorm);
            if (qnorm==0) {
                qnorm = 1.0;
            }
            qnorm = limit_scaling(new ArrayRealVector(1,qnorm)).getEntry(0);
            double mean = 0;
            for (int i = 0;i<H1norm.getDimension();i++) {
                mean+=H1norm.getEntry(i) / H1norm.getDimension();
            }

            lambda = 1.0/limit_scaling(new ArrayRealVector(1,FastMath.max(mean, qnorm))).getEntry(0);

            H1 = H1.copy().scalarMultiply(lambda);
           // A1 = A1.copy().scalarMultiply(lambda);
            q1 = q1.copy().mapMultiply(lambda);
            c = lambda * c;
            gamma = new ArrayRealVector(gammaD.toArray(), gammaE.toArray());
            SD = D.multiply(SD.copy());
            SE = E.multiply(SE.copy());
            iteration += 1;
        }
        this.E = SE.copy();
        this.D = SD.copy();
        this.c = c;
        this.Einv = MatrixUtils.inverse(SE);
        this.Dinv = MatrixUtils.inverse(SD);
        this.cinv = 1.0 / c;




    }

    public RealMatrix getScaledH() {


        return (D.multiply(H).multiply(D)).scalarMultiply(c);
    }

    public RealMatrix getScaledA() {

        return E.multiply(A).multiply(D);
    }

    public RealVector getScaledq() {

        return (D.operate(q.mapMultiply(c)));
    }

    public RealVector getScaledLUb(RealVector lb1) {
        RealVector lb = new ArrayRealVector(lb1.getDimension());
        for (int i = 0; i < lb1.getDimension(); i++) {
            if (lb1.getEntry(i) != Double.POSITIVE_INFINITY) {
                lb.setEntry(i, E.getEntry(i, i) * lb1.getEntry(i));
            } else {
                lb.setEntry(i, Double.POSITIVE_INFINITY);
            }
        }
        return lb;
    }

    public RealVector unscaleX(RealVector x) {
        return D.operate(x);
    }

    public RealVector unscaleY(RealVector y) {

        return (E.operate(y)).mapMultiply(cinv);
    }

    public RealVector unscaleZ(RealVector z) {
        return Einv.operate(z);
    }

    RealVector scaleX(RealVector x) {
        return Dinv.operate(x);
    }

    private RealVector limit_scaling(RealVector v) {

        RealVector result = new ArrayRealVector(v.getDimension());
        for (int i = 0; i <v.getDimension() ; i++) {


            result.setEntry(i,  v.getEntry(i) < MIN_SCALING ? 1.0 : v.getEntry(i));
            result.setEntry(i,  v.getEntry(i) > MAX_SCALING ? MAX_SCALING : v.getEntry(i));

        }
        return result;
    }
}
