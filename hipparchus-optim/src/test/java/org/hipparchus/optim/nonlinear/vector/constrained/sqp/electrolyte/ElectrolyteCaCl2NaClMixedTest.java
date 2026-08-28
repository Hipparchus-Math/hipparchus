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
 * ripopt electrolyte problem 4: mixed CaCl2 and NaCl speciation.
 */
public class ElectrolyteCaCl2NaClMixedTest {

    /** Reference Gibbs-energy value reported by the benchmark. */
    private static final double EXPECTED_OBJECTIVE =
            -0.7723717286693826;

    /**
     * The original electrolyte regression checks chemical consistency rather
     * than requiring bitwise agreement of the Gibbs energy. This tolerance
     * accommodates converged solutions produced by different globalization
     * strategies while still rejecting a different equilibrium basin.
     */
    private static final double OBJECTIVE_TOLERANCE =
            1.0e-4;

    /** Maximum accepted equality-constraint residual. */
    private static final double FEASIBILITY_TOLERANCE =
            1.0e-4;

    /** Tolerance used for the conserved analytical totals. */
    private static final double COMPOSITION_TOLERANCE =
            1.0e-3;

    /** Number of optimization variables. */
    private static final int N = 6;

    /** Species charges: Ca2+, Na+, Cl-, H+, OH-, CaOH+. */
    private static final double[] CHARGES = {
        2.0, 1.0, -1.0, 1.0, -1.0, 1.0
    };

    /** Extended Debye-Huckel ion-size parameters. */
    private static final double[] DH_A = {
        6.0, 4.0, 3.0, 9.0, 3.5, 4.0
    };

    /** Extended Debye-Huckel b-dot parameters. */
    private static final double[] DH_B = {
        0.165, 0.075, 0.015, 0.0, 0.0, 0.0
    };

    /** Standard-state chemical potentials. */
    private static final double[] MU0 = {
        0.0,
        0.0,
        0.0,
        0.0,
        14.0 * ElectrolyteTestSupport.LN10,
        (14.0 - 1.3) * ElectrolyteTestSupport.LN10
    };

    /** Mass-balance and electroneutrality equations. */
    private static final class Equality extends EqualityConstraint {

        Equality() {
            super(new ArrayRealVector(4));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector x) {

            final double calcium =
                    FastMath.exp(x.getEntry(0));
            final double sodium =
                    FastMath.exp(x.getEntry(1));
            final double chloride =
                    FastMath.exp(x.getEntry(2));
            final double hydrogen =
                    FastMath.exp(x.getEntry(3));
            final double hydroxide =
                    FastMath.exp(x.getEntry(4));
            final double calciumHydroxide =
                    FastMath.exp(x.getEntry(5));

            return new ArrayRealVector(new double[] {
                calcium + calciumHydroxide - 0.05,
                sodium - 0.1,
                chloride - 0.2,
                2.0 * calcium +
                sodium +
                hydrogen +
                calciumHydroxide -
                chloride -
                hydroxide
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {

            final double calcium =
                    FastMath.exp(x.getEntry(0));
            final double sodium =
                    FastMath.exp(x.getEntry(1));
            final double chloride =
                    FastMath.exp(x.getEntry(2));
            final double hydrogen =
                    FastMath.exp(x.getEntry(3));
            final double hydroxide =
                    FastMath.exp(x.getEntry(4));
            final double calciumHydroxide =
                    FastMath.exp(x.getEntry(5));

            final RealMatrix jacobian =
                    MatrixUtils.createRealMatrix(4, N);

            jacobian.setEntry(0, 0, calcium);
            jacobian.setEntry(0, 5, calciumHydroxide);

            jacobian.setEntry(1, 1, sodium);

            jacobian.setEntry(2, 2, chloride);

            jacobian.setEntry(3, 0, 2.0 * calcium);
            jacobian.setEntry(3, 1, sodium);
            jacobian.setEntry(3, 2, -chloride);
            jacobian.setEntry(3, 3, hydrogen);
            jacobian.setEntry(3, 4, -hydroxide);
            jacobian.setEntry(3, 5, calciumHydroxide);

            return jacobian;
        }
    }

    @Test
    public void testCaCl2NaClMixed() {

        final LagrangeSolution solution =
                ElectrolyteTestSupport.solve(
                        "CaCl2+NaCl mixed",
                        new ElectrolyteTestSupport.GibbsObjective(
                                MU0,
                                CHARGES,
                                DH_A,
                                DH_B),
                        new Equality(),
                        new double[] {
                            FastMath.log(0.05),
                            FastMath.log(0.1),
                            FastMath.log(0.2),
                            FastMath.log(1.0e-7),
                            FastMath.log(1.0e-7),
                            FastMath.log(1.0e-6)
                        },
                        new SimpleBounds(
                                ElectrolyteTestSupport.filled(N, -35.0),
                                ElectrolyteTestSupport.filled(N, 1.0)),
                        EXPECTED_OBJECTIVE,
                        OBJECTIVE_TOLERANCE,
                        FEASIBILITY_TOLERANCE);

        final RealVector x = solution.getX();

        final double calcium =
                FastMath.exp(x.getEntry(0));
        final double sodium =
                FastMath.exp(x.getEntry(1));
        final double chloride =
                FastMath.exp(x.getEntry(2));
        final double calciumHydroxide =
                FastMath.exp(x.getEntry(5));

        assertTrue(
                FastMath.abs(calcium + calciumHydroxide - 0.05) <=
                COMPOSITION_TOLERANCE,
                "Unexpected total calcium concentration.");

        assertTrue(
                FastMath.abs(sodium - 0.1) <=
                COMPOSITION_TOLERANCE,
                "Unexpected sodium concentration: " + sodium);

        assertTrue(
                FastMath.abs(chloride - 0.2) <=
                COMPOSITION_TOLERANCE,
                "Unexpected chloride concentration: " + chloride);
    }
}
