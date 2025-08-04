package org.hipparchus.filtering.kalman;

import org.hipparchus.filtering.kalman.extended.NonLinearEvolution;
import org.hipparchus.filtering.kalman.extended.NonLinearProcess;
import org.hipparchus.linear.DecompositionSolver;
import org.hipparchus.linear.QRDecomposition;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.SingularValueDecomposition;

public class CramerRaoBound<T extends Measurement> {

    /** Simple constructor.
     * @param process non-linear process to estimate
     * @param R measurement noise covariance
     */
    public CramerRaoBound(final NonLinearProcess<T> process, final RealMatrix P0, final RealMatrix R) {
        this.process = process;
        this.P = P0;
        this.R = R;
    }

    /** Compute Cramér-Rao bound for a given step and a given measurement.
     * @param trueState true state
     * @param measurement measurement
     * @return Cramér-Rao bound
     */
    RealMatrix computeStepBound(final RealVector trueState, final T measurement) {
        final NonLinearEvolution evolution = process.getEvolution(measurement.getTime(), trueState, measurement);
        final RealMatrix F = evolution.getStateTransitionMatrix();
        final RealMatrix H = evolution.getMeasurementJacobian();
        final RealMatrix Q = evolution.getProcessNoiseMatrix();

        // J_meas = H.T * (R^-1 * H)
        DecompositionSolver solver = new QRDecomposition(R).getSolver();
        RealMatrix Rm1_H = solver.solve(H);
        RealMatrix J_meas = H.transpose().multiply(Rm1_H);

        // J_dyn = inv(F * P * F.T + Qk)
        RealMatrix F_P = F.multiply(P);
        RealMatrix F_P_F_transpose = F_P.multiply(F.transpose());
        RealMatrix F_P_F_transpose_Q = F_P_F_transpose.add(Q);
        SingularValueDecomposition svd = new SingularValueDecomposition(F_P_F_transpose_Q);
        RealMatrix J_dyn = svd.getSolver().getInverse();

        // J = J_dyn + J_meas
        RealMatrix J = J_meas.add(J_dyn);

        // P = inv(J) (BCR)
        SingularValueDecomposition svd_Jpost = new SingularValueDecomposition(J);
        P = svd_Jpost.getSolver().getInverse();

        return P;
    }

    private final NonLinearProcess<T> process;
    private final RealMatrix R;
    private RealMatrix P;
}
