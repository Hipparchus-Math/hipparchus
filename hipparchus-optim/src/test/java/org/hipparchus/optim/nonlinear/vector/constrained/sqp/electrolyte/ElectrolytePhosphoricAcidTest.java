package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import static org.junit.jupiter.api.Assertions.assertTrue;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/**
 * ripopt electrolyte problem 5: phosphoric-acid speciation.
 */
public class ElectrolytePhosphoricAcidTest {

    /** Reference Gibbs-energy value reported by the benchmark. */
    private static final double EXPECTED_OBJECTIVE =
            -0.05531209113684935;

    /**
     * The Gibbs-energy value is used to reject a different equilibrium basin.
     * Chemical validity is additionally checked through pH and trace PO4^3-.
     */
    private static final double OBJECTIVE_TOLERANCE =
            1.0e-4;

    /** Maximum accepted equality-constraint residual. */
    private static final double FEASIBILITY_TOLERANCE =
            1.0e-4;

    /** Total phosphate molality. */
    private static final double P_TOTAL =
            0.01;

    /** Number of optimization variables. */
    private static final int N = 6;

    /**
     * Species charges:
     * H3PO4, H2PO4-, HPO4^2-, PO4^3-, H+, OH-.
     */
    private static final double[] CHARGES = {
        0.0, -1.0, -2.0, -3.0, 1.0, -1.0
    };

    /** Extended Debye-Huckel ion-size parameters. */
    private static final double[] DH_A = {
        0.0, 4.5, 4.0, 4.0, 9.0, 3.5
    };

    /** Extended Debye-Huckel b-dot parameters. */
    private static final double[] DH_B = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };

    /** Standard-state chemical potentials. */
    private static final double[] MU0 = {
        0.0,
        2.148 * ElectrolyteTestSupport.LN10,
        (2.148 + 7.199) * ElectrolyteTestSupport.LN10,
        (2.148 + 7.199 + 12.35) *
                ElectrolyteTestSupport.LN10,
        0.0,
        14.0 * ElectrolyteTestSupport.LN10
    };

    /** Phosphate balance and electroneutrality equations. */
    private static final class Equality extends EqualityConstraint {

        Equality() {
            super(new ArrayRealVector(2));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector x) {

            final double h3po4 =
                    FastMath.exp(x.getEntry(0));
            final double h2po4 =
                    FastMath.exp(x.getEntry(1));
            final double hpo4 =
                    FastMath.exp(x.getEntry(2));
            final double po4 =
                    FastMath.exp(x.getEntry(3));
            final double hydrogen =
                    FastMath.exp(x.getEntry(4));
            final double hydroxide =
                    FastMath.exp(x.getEntry(5));

            return new ArrayRealVector(new double[] {
                h3po4 + h2po4 + hpo4 + po4 - P_TOTAL,
                -h2po4 -
                2.0 * hpo4 -
                3.0 * po4 +
                hydrogen -
                hydroxide
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {

            final double[] molality =
                    new double[N];

            for (int i = 0; i < N; ++i) {
                molality[i] =
                        FastMath.exp(x.getEntry(i));
            }

            final RealMatrix jacobian =
                    MatrixUtils.createRealMatrix(2, N);

            for (int i = 0; i < 4; ++i) {
                jacobian.setEntry(0, i, molality[i]);
            }

            jacobian.setEntry(1, 1, -molality[1]);
            jacobian.setEntry(1, 2, -2.0 * molality[2]);
            jacobian.setEntry(1, 3, -3.0 * molality[3]);
            jacobian.setEntry(1, 4, molality[4]);
            jacobian.setEntry(1, 5, -molality[5]);

            return jacobian;
        }
    }

    @Test
    public void testPhosphoricAcid() {

        final LagrangeSolution solution =
                ElectrolyteTestSupport.solve(
                        "Phosphoric acid",
                        new ElectrolyteTestSupport.GibbsObjective(
                                MU0,
                                CHARGES,
                                DH_A,
                                DH_B),
                        new Equality(),
                        new double[] {
                            FastMath.log(0.005),
                            FastMath.log(0.004),
                            FastMath.log(1.0e-5),
                            FastMath.log(1.0e-15),
                            FastMath.log(0.005),
                            FastMath.log(1.0e-10)
                        },
                        new SimpleBounds(
                                ElectrolyteTestSupport.filled(N, -55.0),
                                ElectrolyteTestSupport.filled(N, 0.0)),
                        EXPECTED_OBJECTIVE,
                        OBJECTIVE_TOLERANCE,
                        FEASIBILITY_TOLERANCE);

        final RealVector x =
                solution.getX();

        final double hydrogenMolality =
                FastMath.exp(x.getEntry(4));
        final double pH =
                -FastMath.log10(hydrogenMolality);
        final double phosphate =
                FastMath.exp(x.getEntry(3));

        assertTrue(
                pH > 1.0 && pH < 4.0,
                "Chemically invalid phosphoric-acid pH: " + pH);

        assertTrue(
                phosphate < 1.0e-6,
                "PO4^3- is not a trace species: " + phosphate);
    }
}
