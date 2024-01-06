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

    /** Minimum scaling value. */
    private static final double MIN_SCALING = 1.0e6;

    /** Maximum scaling value. */
    private static final double MAX_SCALING = 1.0e+6;

    /** Square matrix of weights for quadratic terms. */
    private final RealMatrix H;

    /** Vector of weights for linear terms. */
    private final RealVector q;

    /** constraints coefficients matrix. */
    private final RealMatrix A;

    /** TBC. */
    private RealMatrix D;

    /** TBC. */
    private RealMatrix E;

    /** TBC. */
    private double c;

    /** Inverse of D. */
    private RealMatrix Dinv;

    /** Inverse of E. */
    private RealMatrix Einv;

    /** Inverse of c. */
    private double cinv;

    public ADMMQPModifiedRuizEquilibrium(RealMatrix H, RealMatrix A, RealVector q) {
        this.H  = H;
        this.A  = A;
        this.q  = q;
    }

    public void normalize(double epsilon, int maxIteration) {

        int iteration = 0;
        this.c = 1.0;
        RealVector gammaD = new ArrayRealVector(H.getRowDimension());
        RealVector gammaE = new ArrayRealVector(A.getRowDimension());

        double lambda;
        RealVector q1 = q.copy();
        RealMatrix H1 = H.copy();
        RealMatrix A1 = A.copy();
        RealMatrix diagD = null;
        RealMatrix diagE = null;
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

            gammaD = limitScaling(gammaD);
            gammaE = limitScaling(gammaE);

            for (int i = 0; i < gammaD.getDimension(); i++) {
                gammaD.setEntry(i, 1.0 / FastMath.sqrt(gammaD.getEntry(i)));
            }

            for (int i = 0; i < gammaE.getDimension(); i++) {
                gammaE.setEntry(i, 1.0 / FastMath.sqrt(gammaE.getEntry(i)));
            }

            diagD = MatrixUtils.createRealDiagonalMatrix(gammaD.toArray());
            diagE = MatrixUtils.createRealDiagonalMatrix(gammaE.toArray());

            H1 = diagD.multiply(H1.copy()).multiply(diagD.copy());
            q1 = diagD.operate(q1.copy());
            A1 = diagE.multiply(A1.copy()).multiply(diagD.copy());

            for (int i = 0; i < H1.getRowDimension(); i++) {
                double norma = (new ArrayRealVector(H1.getRow(i))).getLInfNorm();
                H1norm.setEntry(i, norma);
            }
            double qnorm = q1.getLInfNorm();
            if (qnorm == 0) {
                qnorm = 1.0;
            }
            qnorm = limitScaling(new ArrayRealVector(1, qnorm)).getEntry(0);
            double mean = 0;
            for (int i = 0; i < H1norm.getDimension(); i++) {
                mean+=H1norm.getEntry(i) / H1norm.getDimension();
            }

            lambda = 1.0 / limitScaling(new ArrayRealVector(1, FastMath.max(mean, qnorm))).getEntry(0);

            H1 = H1.copy().scalarMultiply(lambda);
            q1 = q1.copy().mapMultiply(lambda);
            c *= lambda;
            SD = diagD.multiply(SD.copy());
            SE = diagE.multiply(SE.copy());
            iteration += 1;
        }
        this.E    = SE.copy();
        this.D    = SD.copy();
        this.Einv = MatrixUtils.inverse(SE);
        this.Dinv = MatrixUtils.inverse(SD);
        this.cinv = 1.0 / c;

    }

    public RealMatrix getScaledH() {
        return D.multiply(H).multiply(D).scalarMultiply(c);
    }

    public RealMatrix getScaledA() {
        return E.multiply(A).multiply(D);
    }

    public RealVector getScaledq() {
        return D.operate(q.mapMultiply(c));
    }

    public RealVector getScaledLUb(RealVector lb1) {
        RealVector scaledLUb = new ArrayRealVector(lb1.getDimension());
        for (int i = 0; i < lb1.getDimension(); i++) {
            if (lb1.getEntry(i) != Double.POSITIVE_INFINITY) {
                scaledLUb.setEntry(i, E.getEntry(i, i) * lb1.getEntry(i));
            } else {
                scaledLUb.setEntry(i, Double.POSITIVE_INFINITY);
            }
        }
        return scaledLUb;
    }

    public RealVector unscaleX(RealVector x) {
        return D.operate(x);
    }

    public RealVector unscaleY(RealVector y) {
        return E.operate(y).mapMultiply(cinv);
    }

    public RealVector unscaleZ(RealVector z) {
        return Einv.operate(z);
    }

    RealVector scaleX(RealVector x) {
        return Dinv.operate(x);
    }

    private RealVector limitScaling(RealVector v) {

        RealVector result = new ArrayRealVector(v.getDimension());
        for (int i = 0; i < v.getDimension(); i++) {
            result.setEntry(i,  v.getEntry(i) < MIN_SCALING ? 1.0 : v.getEntry(i));
            result.setEntry(i,  v.getEntry(i) > MAX_SCALING ? MAX_SCALING : v.getEntry(i));
        }

        return result;

    }

}
