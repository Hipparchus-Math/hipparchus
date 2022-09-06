package org.hipparchus.ode;

import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.RotationConvention;
import org.hipparchus.geometry.euclidean.threed.RotationOrder;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.special.elliptic.jacobi.CopolarN;
import org.hipparchus.special.elliptic.jacobi.JacobiElliptic;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;


public class TestProblem8 extends TestProblemAbstract {

    /** Moments of inertia. */
    final double i1;
    final double i2;
    final double i3;

    /** Moments of inertia converted. */
    final double i1C;
    final double i2C;
    final double i3C;

    /**Substraction of inertias. */
    final double i32;
    final double i31;
    final double i21;

    /** Initial state. */
    final double[] y0;
    
    /** Converted initial state. */
    final double[] y0C;

    /** Converted axes. */
    final Vector3D[] axes;

    /** Twice the angular kinetic energy. */
    final double twoE;

    /** Square of kinetic momentum. */
    final double m2;

    /** State scaling factor. */
    final double o1Scale;

    /** State scaling factor. */
    final double o2Scale;

    /** State scaling factor. */
    final double o3Scale;

    /** Jacobi elliptic function. */
    final JacobiElliptic jacobi;

    /** Elliptic modulus k2 (k2 = m). */
    final double k2;

    /** Time scaling factor. */
    public final double tScale;

    /** Time reference for rotation rate. */
    final double tRef;

    /**Offset rotation  between initial inertial frame and the frame with moment vector and Z axis aligned. */
    Rotation mAlignedToInert;

    /** Rotation to switch to the converted axes frame. */
    final Rotation sortedToBody;

    /** Period of rotation rate. */
    final double period;

    /** Slope of the linear part of the phi model. */
    final double phiSlope;

    /** DenseOutputModel of phi. */
    final DenseOutputModel phiQuadratureModel;

    /** Integral part of quadrature model over one period. */
    final double integOnePeriod;
    private static int count = 0;

    /**
     * Simple constructor.
     * @param t0 initial time
     * @param t1 final time
     * @param omega0 initial rotation rate
     * @param r0 initial rotation
     * @param i1 inertia along first axis
     * @param i2 inertia along second axis
     * @param i3 inertia along third axis
     */
    public TestProblem8(final double t0, final double t1, final Vector3D omega0, final Rotation r0,
                        final double i1, final double i2, final double i3) {
        //Arguments in the super constructor :
        //Initial time, Primary state (o1, o2, o3, q0, q1, q2, q3), Final time, Error scale
        super(t0,
              new double[] {
                  omega0.getX(), omega0.getY(), omega0.getZ(),
                  r0.getQ0(), r0.getQ1(), r0.getQ2(), r0.getQ3()
              },
              t1,
              new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });
        this.i1 = i1;
        this.i2 = i2;
        this.i3 = i3;

        y0 = getInitialState().getPrimaryState();

        final double o12 = y0[0] * y0[0];
        final double o22 = y0[1] * y0[1];
        final double o32 = y0[2] * y0[2];

        twoE    =  i1 * o12 + i2 * o22 + i3 * o32;
        m2      =  i1 * i1 * o12 + i2 * i2 * o22 + i3 * i3 * o32;

        Vector3D[] axesP = { Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_K };

        double[] i = { i1, i2, i3 };
        y0C = y0.clone();

        System.out.print(++count + " " + i1 + " " + i2 + " " + i3 + " " + (m2 / twoE));
        if (i[0] > i[1]) {
            System.out.print(" 1↔2");
            Vector3D z = axesP[0];
            axesP[0] = axesP[1];
            axesP[1] = z;
            axesP[2] = axesP[2].negate();

            final double y = i[0];
            i[0] = i[1];
            i[1] = y;

            final double v = y0C[0];
            y0C[0] = y0C[1];
            y0C[1] = v;
            y0C[2] = -y0C[2];
        } else {
            System.out.print("  ∅ ");
        }

        if (i[1] > i[2]) {
            System.out.print(" 2↔3");
            Vector3D z = axesP[1];
            axesP[1] = axesP[2];
            axesP[2] = z;
            axesP[0] = axesP[0].negate();

            final double y = i[1];
            i[1] = i[2];
            i[2] = y;

            final double v = y0C[1];
            y0C[1] = y0C[2];
            y0C[2] = v;
            y0C[0] = -y0C[0];
        } else {
            System.out.print("  ∅ ");
        }

        if (i[0] > i[1]) {
            System.out.print(" 1↔2");
            Vector3D z = axesP[0];
            axesP[0] = axesP[1];
            axesP[1] = z;
            axesP[2] = axesP[2].negate();

            final double y = i[0];
            i[0] = i[1];
            i[1] = y;

            final double v = y0C[0];
            y0C[0] = y0C[1];
            y0C[1] = v;
            y0C[2] = -y0C[2];
        } else {
            System.out.print("  ∅ ");
        }

        final double condition = (twoE == 0) ? 0.0 : m2 / twoE;

        final boolean clockwise;
        if (condition < i[1]) {
            System.out.print(" 1↔3");
            clockwise = true;
            Vector3D z = axesP[0];
            axesP[0] = axesP[2];
            axesP[2] = z;
            axesP[1] = axesP[1].negate();

            final double y = i[0];
            i[0] = i[2];
            i[2] = y;

            final double v = y0C[0];
            y0C[0] = y0C[2];
            y0C[2] = v;
            y0C[1] = -y0C[1];

        } else {
            System.out.print("  ∅ ");
            clockwise = false;
        }

        i1C = i[0];
        i2C = i[1];
        i3C = i[2];

        axes = axesP;
        System.out.print("   " + axes[0].getX() + " " + axes[0].getY() + " " + axes[0].getZ());
        System.out.print("   " + axes[1].getX() + " " + axes[1].getY() + " " + axes[1].getZ());
        System.out.print("   " + axes[2].getX() + " " + axes[2].getY() + " " + axes[2].getZ());

        // convert initial conditions to Euler angles such the M is aligned with Z in sorted computation frame
        final Vector3D omega0Sorted = new Vector3D(y0C[0], y0C[1], y0C[2]);
        final Vector3D m0Sorted     = new Vector3D(i1C * omega0Sorted.getX(), i2C * omega0Sorted.getY(), i3C * omega0Sorted.getZ());
        final double   phi0         = 0; // this angle can be set arbitrarily, so 0 is a fair value (Eq. 37.13 - 37.14)
        final double   theta0       = FastMath.acos(m0Sorted.getZ() / m0Sorted.getNorm());
        final double   psi0         = FastMath.atan2(m0Sorted.getX(), m0Sorted.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!

        // compute offset rotation between inertial frame aligned with momentum and regular inertial frame
        final Rotation mAlignedToSorted = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                                       phi0, theta0, psi0);

        sortedToBody = new Rotation(Vector3D.PLUS_I, Vector3D.PLUS_J, axesP[0], axesP[1]);

        //Est-il nécessaire de garder le r0COnvertedAxis qui n'est peut être pas adapté et ne règle peut être pas le problème
        mAlignedToInert = sortedToBody.applyTo(r0).applyInverseTo(mAlignedToSorted);
        //mAlignedToInert = r0.applyInverseTo(mAlignedToBody);        

        i32  = i3C - i2C;
        i31  = i3C - i1C;
        i21  = i2C - i1C;

        // Ω is always o1Scale * cn((t-tref) * tScale), o2Scale * sn((t-tref) * tScale), o3Scale * dn((t-tref) * tScale)
        tScale  = FastMath.copySign(FastMath.sqrt(i32 * (m2 - twoE * i1C) / (i1C * i2C * i3C)),
                                    clockwise ? -omega0Sorted.getZ() : omega0Sorted.getZ());
        System.out.print(" " + tScale);
        o1Scale = FastMath.sqrt((twoE * i3C - m2) / (i1C * i31));
        o2Scale = FastMath.sqrt((twoE * i3C - m2) / (i2C * i32));
        o3Scale = FastMath.copySign(FastMath.sqrt((m2 - twoE * i1C) / (i3C * i31)), omega0Sorted.getZ());

        k2     = (twoE == 0) ? 0.0 : i21 * (twoE * i3C - m2) / (i32 * (m2 - twoE * i1C));
        jacobi = JacobiEllipticBuilder.build(k2);
        period = 4 * LegendreEllipticIntegral.bigK(k2) / tScale;

        if (o1Scale == 0) {
            // special case where twoE * i3C = m2, then o2Scale is also zero
            // motion is exactly along one axis
            tRef = t0;
        } else {
            if (FastMath.abs(omega0Sorted.getX()) >= FastMath.abs(omega0Sorted.getY())) {
                if (omega0Sorted.getX() >= 0) {
                    // omega is roughly towards +I
                    tRef = t0 - jacobi.arcsn(omega0Sorted.getY() / o2Scale) / tScale;
                    System.out.print(" t₀-arcsn");
                } else {
                    // omega is roughly towards -I
                    tRef = t0 + jacobi.arcsn(omega0Sorted.getY() / o2Scale) / tScale - 0.5 * period;
                    System.out.print(" t₀+arcsn");
                }
            } else {
                if (omega0Sorted.getY() >= 0) {
                    // omega is roughly towards +J
                    tRef = t0 - jacobi.arccn(omega0Sorted.getX() / o1Scale) / tScale;
                    System.out.print(" t₀-arccn");
                } else {
                    // omega is roughly towards -J
                    tRef = t0 + jacobi.arccn(omega0Sorted.getX() / o1Scale) / tScale;
                    System.out.print(" t₀+arccn");
                }
            }
        }

        phiSlope           = FastMath.sqrt(m2) / i3C;
        phiQuadratureModel = computePhiQuadratureModel(t0);
        integOnePeriod     = phiQuadratureModel.getInterpolatedState(phiQuadratureModel.getFinalTime()).getPrimaryState()[0];

    }

    private DenseOutputModel computePhiQuadratureModel(final double t0) {

        final double i32 = i3C - i2C;
        final double i31 = i3C - i1C;
        final double i21 = i2C - i1C;

        // coefficients for φ model
        final double b = phiSlope * i32 * i31;
        final double c = i1C * i32;
        final double d = i3C * i21;

        // integrate the quadrature phi term on one period
        final DormandPrince853Integrator integ = new DormandPrince853Integrator(1.0e-6 * period, 1.0e-2 * period,
                                                                                phiSlope * period * 1.0e-13, 1.0e-13);
        final DenseOutputModel model = new DenseOutputModel();
        integ.addStepHandler(model);

        integ.integrate(new OrdinaryDifferentialEquation() {

            /** {@inheritDoc} */
            @Override
            public int getDimension() {
                return 1;
            }

            /** {@inheritDoc} */
           @Override
            public double[] computeDerivatives(final double t, final double[] y) {
                final double sn = jacobi.valuesN((t - tRef) * tScale).sn();
                return new double[] {
                        b / (c + d * sn * sn)
                };
            }

        }, new ODEState(t0, new double[1]), t0 + period);

        return model;

    }

    public TfmState computeTorqueFreeMotion(double t) {

        // Computation of omega
        final CopolarN valuesN = jacobi.valuesN((t - tRef) * tScale);
        final Vector3D omegaP  = new Vector3D(o1Scale * valuesN.cn(), o2Scale * valuesN.sn(), o3Scale * valuesN.dn());
        final Vector3D omega   = sortedToBody.applyTo(omegaP);

        // Computation of the Euler angles
        final double   psi       = FastMath.atan2(i1C * omegaP.getX(), i2C * omegaP.getY());
        final double   theta     = FastMath.acos(omegaP.getZ() / phiSlope);
        final double   phiLinear = phiSlope * t;

        // Integration for the computation of phi
        final double t0            = getInitialTime();
        final int    nbPeriods     = (int) FastMath.floor((t - t0) / period);
        final double tStartInteg   = t0 + nbPeriods * period;
        final double integPartial  = phiQuadratureModel.getInterpolatedState(t - tStartInteg).getPrimaryState()[0];
        final double phiQuadrature = nbPeriods * integOnePeriod + integPartial;

        final double phi = phiLinear + phiQuadrature;

        // Computation of the quaternion

        // Rotation between computation frame (aligned with momentum) and sorted computation frame
        //(It is simply the angles equations provided by L&L)
        final Rotation alignedToSorted = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                                      phi, theta, psi);

        // combine with offset rotation to get back to regular inertial frame
        // Inert -> aligned + aligned -> body = inert -> body (What the user wants)
        Rotation inertToSorted = alignedToSorted.applyTo(mAlignedToInert.revert());

        Rotation bodyToOriginalFrame = sortedToBody.applyTo(inertToSorted);

        return new TfmState(t, omega, bodyToOriginalFrame, phi, theta, psi, sortedToBody, mAlignedToInert);

    }

    public double[] doComputeDerivatives(double t, double[] y) {

        final  double[] yDot = new double[getDimension()];

        // compute the derivatives using Euler equations
        yDot[0] = y[1] * y[2] * (i2 - i3) / i1;
        yDot[1] = y[2] * y[0] * (i3 - i1) / i2;
        yDot[2] = y[0] * y[1] * (i1 - i2) / i3;

        // compute the derivatives using Qdot = 0.5 * Omega_inertialframe * Q
        yDot[3] = 0.5 * (-y[0] * y[4] -y[1] * y[5] -y[2] * y[6]);
        yDot[4] = 0.5 * (y[0] * y[3] +y[2] * y[5] -y[1] * y[6]);
        yDot[5] = 0.5 * (y[1] * y[3] -y[2] * y[4] +y[0] * y[6]);
        yDot[6] = 0.5 * (y[2] * y[3] +y[1] * y[4] -y[0] * y[5]);

        return yDot;

    }

    public double[] computeTheoreticalState(double t) {
        final TfmState tfm = computeTorqueFreeMotion(t);
        return new double[] {
                tfm.getOmega().getX(),
                tfm.getOmega().getY(),
                tfm.getOmega().getZ(),
                tfm.getRotation().getQ0(),
                tfm.getRotation().getQ1(),
                tfm.getRotation().getQ2(),
                tfm.getRotation().getQ3()
        };
    }

    public static class TfmState {
        private final double   t;
        private final Vector3D omega;
        private final Rotation rotation;
        private final double phi;
        private final double theta;
        private final double psi;
        private final Rotation sortedToBody;
        private final Rotation mAlignedToInert;
        private TfmState(final double t, final Vector3D omega, final Rotation rotation,
                         final double phi, final double theta, final double psi,
                         final Rotation sortedToBody, final Rotation mAlignedToInert) {
            this.t               = t;
            this.omega           = omega;
            this.rotation        = rotation;
            this.phi             = phi;
            this.theta           = theta;
            this.psi             = psi;
            this.sortedToBody     = sortedToBody;
            this.mAlignedToInert = mAlignedToInert;
        }
        public double getT() {
            return t;
        }
        public Vector3D getOmega() {
            return omega;
        }
        public Rotation getRotation() {
            return rotation;
        }
        public double getPhi() {
            return phi;
        }
        public double getTheta() {
            return theta;
        }
        public double getPsi() {
            return psi;
        }
        public Rotation getSortedToBody() {
            return sortedToBody;
        }
        public Rotation getMAlignedToInert() {
            return mAlignedToInert;
        }
    }

}//Fin du programme
