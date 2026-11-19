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
import org.hipparchus.linear.DecompositionSolver;
import org.hipparchus.linear.EigenDecompositionSymmetric;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;

/** Alternative Direction Method of Multipliers Solver.
 * @since 3.1
 */
public class ADMMQPKKT implements KarushKuhnTuckerSolver<ADMMQPSolution> {

    /** Tolerance for symmetric matrices decomposition.
     * @since 4.1
     */
    private double decompositionEpsilon;

    /** Square matrix of weights for quadratic terms. */
    private RealMatrix H;

    /** Vector of weights for linear terms. */
    private RealVector q;

    /** Constraints coefficients matrix. */
    private RealMatrix A;

    /** Regularization term sigma for Karush–Kuhn–Tucker solver. */
    private double sigma;

    /** TBC. */
    private RealMatrix R;

    /** Inverse of R. */
    private RealMatrix Rinv;

    /** Lower bound. */
    private RealVector lb;

    /** Upper bound. */
    private RealVector ub;

    /** Alpha filter for ADMM iteration. */
    private double alpha;

    /** Constrained problem KKT matrix. */
    private RealMatrix M; // NOPMD

    /** Solver for M. */
    private DecompositionSolver dsX;

    /** Simple constructor.
     * <p>
     * BEWARE, nothing is initialized here, it is {@link #initialize(RealMatrix, RealMatrix,
     * RealVector, int, RealVector, RealVector, double, double, double) initialize} <em>must</em>
     * be called before using the instance.
     * </p>
     */
    ADMMQPKKT() {
        decompositionEpsilon = EigenDecompositionSymmetric.DEFAULT_EPSILON;
    }

    /** {@inheritDoc} */
    @Override
    public ADMMQPSolution solve(RealVector b1, final RealVector b2) {
        RealVector z = dsX.solve(new ArrayRealVector((ArrayRealVector) b1, b2));
        return new ADMMQPSolution(z.getSubVector(0,b1.getDimension()), z.getSubVector(b1.getDimension(), b2.getDimension()));
    }

    /** Update tolerance for matrix decomposition
     * @param newDecompositionEpsilon tolerance for symmetric matrix decomposition
     * @since 4.1
     */
    public void updateDecompositionEpsilon(final double newDecompositionEpsilon) {
        this.decompositionEpsilon = newDecompositionEpsilon;
    }

    /** Update steps
     * @param newSigma new regularization term sigma for Karush–Kuhn–Tucker solver
     * @param me number of equality constraints
     * @param rho new step size
     */
    public void updateSigmaRho(double newSigma, int me, double rho) {
        this.sigma = newSigma;
        this.H = H.add(MatrixUtils.createRealIdentityMatrix(H.getColumnDimension()).scalarMultiply(newSigma));
        createPenaltyMatrix(me, rho);
        M =  MatrixUtils.createRealMatrix(H.getRowDimension() + A.getRowDimension(),
                                          H.getRowDimension() + A.getRowDimension());
        M.setSubMatrix(H.getData(), 0,0);
        M.setSubMatrix(A.getData(), H.getRowDimension(),0);
        M.setSubMatrix(A.transpose().getData(), 0, H.getRowDimension());
        M.setSubMatrix(Rinv.scalarMultiply(-1.0).getData(), H.getRowDimension(),H.getRowDimension());
        dsX = new EigenDecompositionSymmetric(M, decompositionEpsilon, true).getSolver();
    }

    /** Initialize problem
     * @param newH square matrix of weights for quadratic term
     * @param newA constraints coefficients matrix
     * @param newQ TBD
     * @param me number of equality constraints
     * @param newLb lower bound
     * @param newUb upper bound
     * @param rho step size
     * @param newSigma regularization term sigma for Karush–Kuhn–Tucker solver
     * @param newAlpha alpha filter for ADMM iteration
     */
    public void initialize(RealMatrix newH, RealMatrix newA, RealVector newQ,
                           int me, RealVector newLb, RealVector newUb,
                           double rho, double newSigma, double newAlpha) {
        this.lb = newLb;
        this.ub = newUb;
        this.alpha = newAlpha;
        this.sigma = newSigma;
        this.H = newH.add(MatrixUtils.createRealIdentityMatrix(newH.getColumnDimension()).scalarMultiply(newSigma));
        this.A = newA.copy();
        this.q = newQ.copy();
        createPenaltyMatrix(me, rho);

        M =  MatrixUtils.createRealMatrix(newH.getRowDimension() + newA.getRowDimension(),
                                          newH.getRowDimension() + newA.getRowDimension());
        M.setSubMatrix(newH.getData(),0,0);
        M.setSubMatrix(newA.getData(),newH.getRowDimension(),0);
        M.setSubMatrix(newA.transpose().getData(),0,newH.getRowDimension());
        M.setSubMatrix(Rinv.scalarMultiply(-1.0).getData(),newH.getRowDimension(),newH.getRowDimension());
        dsX = new EigenDecompositionSymmetric(M, decompositionEpsilon, true).getSolver();
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

    /** {@inheritDoc} */
    @Override
    public ADMMQPSolution iterate(RealVector... previousSol) {
        double onealfa = 1.0 - alpha;
        //SAVE OLD VALUE
        RealVector xold = previousSol[0].copy();
        RealVector yold = previousSol[1].copy();
        RealVector zold = previousSol[2].copy();

        //UPDATE RIGHT VECTOR
        RealVector b1 = previousSol[0].mapMultiply(sigma).subtract(q);
        RealVector b2 = previousSol[2].subtract(Rinv.operate(previousSol[1]));

        //SOLVE KKT SYSYEM
        ADMMQPSolution sol = solve(b1, b2);
        RealVector xtilde = sol.getX();
        RealVector vtilde = sol.getV();

        //UPDATE ZTILDE
        RealVector ztilde = zold.add(Rinv.operate(vtilde.subtract(yold)));
        //UPDATE X
        previousSol[0] = xtilde.mapMultiply(alpha).add(xold.mapMultiply(onealfa));

        //UPDATE Z PARTIAL
        RealVector zpartial = ztilde.mapMultiply(alpha).add(zold.mapMultiply(onealfa)).add(Rinv.operate(yold));

        //PROJECT ZPARTIAL AND UPDATE Z
        for (int j = 0; j < previousSol[2].getDimension(); j++) {
            previousSol[2].setEntry(j, FastMath.min(FastMath.max(zpartial.getEntry(j), lb.getEntry(j)), ub.getEntry(j)));
        }

        //UPDATE Y
        RealVector ytilde = ztilde.mapMultiply(alpha).add(zold.mapMultiply(onealfa).subtract(previousSol[2]));
        previousSol[1] = yold.add(R.operate(ytilde));

        return new ADMMQPSolution(previousSol[0], vtilde, previousSol[1], previousSol[2]);
    }

}
