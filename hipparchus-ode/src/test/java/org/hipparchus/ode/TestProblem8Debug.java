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


public class TestProblem8Debug extends TestProblemAbstract {

    /** Moments of inertia. */
    final double i1;
    final double i2;
    final double i3;

    /** Moments of inertia converted. */
    final double i1C;
    final double i2C;
    final double i3C;

    /** Subtraction of inertias. */
    final double i32;
    final double i31;
    final double i21;

    /**Initial state. */
    final double y0[];

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
    final double tScale;

    /** Time reference for rotation rate. */
    final double tRef;

    /**Offset rotation  between initial inertial frame and the frame with moment vector and Z axis aligned. */
    Rotation mAlignedToInert;

    /** Initial converted rotation. */
    final Rotation r0;

    /** Rotation to switch to the converted axes frame. */
    final Rotation convertAxes;

    /**DenseOutputModel of phi. */
    final double period;
    final double phiSlope;
    final DenseOutputModel phiQuadratureModel;
    final double integOnePeriod;

    //2ème état

    /** Moments of inertia. */
    final double i1DEUX;
    final double i2DEUX;
    final double i3DEUX;

    /** Moments of inertia converted. */
    final double i1CDEUX;
    final double i2CDEUX;
    final double i3CDEUX;

    /**Substraction of inertias. */
    final double i32DEUX;
    final double i31DEUX;
    final double i21DEUX;

    /** Twice the angular kinetic energy. */
    final double twoEDEUX;

    /** Square of kinetic momentum. */
    final double m2DEUX;

    /** State scaling factor. */
    final double o1ScaleDEUX;

    /** State scaling factor. */
    final double o2ScaleDEUX;

    /** State scaling factor. */
    final double o3ScaleDEUX;

    /** Jacobi elliptic function. */
    final JacobiElliptic jacobiDEUX;

    /** Elliptic modulus k2 (k2 = m). */
    final double k2DEUX;

    /** Time scaling factor. */
    final double tScaleDEUX;

    /** Time reference for rotation rate. */
    final double tRefDEUX;

    /**Offset rotation  between initial inertial frame and the frame with moment vector and Z axis aligned. */
    Rotation mAlignedToInertDEUX;

    final Rotation convertAxesDEUX;

    /**DenseOutputModel of phi. */
    final double periodDEUX;
    final double phiSlopeDEUX;
    final DenseOutputModel phiQuadratureModelDEUX;
    final double integOnePeriodDEUX;

    /**
     * Simple constructor.
     */
    public TestProblem8Debug() {
        //Arguments in the super constructor :
        //Initial time, Primary state (o1, o2, o3, q0, q1, q2, q3), Final time, Error scale
        super(0.0, new double[] {5.0, 0.0, 4.0, 0.9 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.437 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.0, 0.0,   0.0, 5.0, -4.0, -0.3088560588509295, 0.6360879930568342, 0.6360879930568341, 0.3088560588509294 }, 20.0, new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });
//        super(0.0, new double[] {5.0, 0.0, 4.0, 0.9 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.437 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.0, 0.0,   5.0, 0.0, 4.0, 0.9 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.437 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.0, 0.0}, 20.0, new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });

        i1 = 3.0 / 8.0;
        i2 = 1.0 / 2.0;
        i3 = 5.0 / 8.0;

        i2DEUX = 3.0 / 8.0;
        i1DEUX = 1.0 / 2.0;
        i3DEUX = 5.0 / 8.0;

        y0 = getInitialState().getPrimaryState();
        final double t0 = getInitialState().getTime();

        //1er état


        final double o12 = y0[0] * y0[0];
        final double o22 = y0[1] * y0[1];
        final double o32 = y0[2] * y0[2];

        twoE    =  i1 * o12 + i2 * o22 + i3 * o32;
        m2      =  i1 * i1 * o12 + i2 * i2 * o22 + i3 * i3 * o32;

        final Sorted sorted = sortInertiaAxis(i1, i2, i3, m2, twoE,
                                              new Vector3D(y0[0], y0[1], y0[2]),
                                              new Rotation(y0[3], y0[4], y0[5], y0[6], true));

        i1C = sorted.i1;
        i2C = sorted.i2;
        i3C = sorted.i3;

        // convert initial conditions to Euler angles such the M is aligned with Z in computation frame
        final Vector3D omega0Body = sorted.omega;
        r0         = sorted.rotation;
        final Vector3D m0Body     = new Vector3D(i1C * omega0Body.getX(), i2C * omega0Body.getY(), i3C * omega0Body.getZ());

        final double   phi0       = 0; // this angle can be set arbitrarily, so 0 is a fair value (Eq. 37.13 - 37.14)
        final double   theta0     = FastMath.acos(m0Body.getZ() / m0Body.getNorm());
        final double   psi0       = FastMath.atan2(m0Body.getX(), m0Body.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!

        //Compute offset rotation between inertial frame aligned with momentum and regular inertial frame
        final Rotation mAlignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                                     phi0, theta0, psi0);

        convertAxes = sorted.convertAxes;

        Rotation r0ConvertedAxis = convertAxes.applyTo(r0);

        mAlignedToInert = r0ConvertedAxis.applyInverseTo(mAlignedToBody);
        //mAlignedToInert = r0.applyInverseTo(mAlignedToBody);

        i32  = i3C - i2C;
        i31  = i3C - i1C;
        i21  = i2C - i1C;

        tScale  = FastMath.sqrt(i32 * (m2 - twoE * i1C) / (i1C * i2C * i3C));
        o1Scale = FastMath.sqrt((twoE * i3C - m2) / (i1C * i31));
        o2Scale = FastMath.sqrt((twoE * i3C - m2) / (i2C * i32));
        o3Scale = FastMath.sqrt((m2 - twoE * i1C) / (i3C * i31));

        k2     = i21 * (twoE * i3C - m2) / (i32 * (m2 - twoE * i1C));

        jacobi = JacobiEllipticBuilder.build(k2);

        if (omega0Body.getY() == 0) {
            tRef = t0;
        } else {
            tRef = t0 - jacobi.arcsn(omega0Body.getY() / o2Scale) / tScale;
        }

        period             = 4 * LegendreEllipticIntegral.bigK(k2) / tScale;
        phiSlope           = FastMath.sqrt(m2) / i3;
        phiQuadratureModel = computePhiQuadratureModel(i1C, i2C, i3C, k2, m2, jacobi, period, phiSlope, t0, tRef, tScale);
        integOnePeriod     = phiQuadratureModel.getInterpolatedState(phiQuadratureModel.getFinalTime()).getPrimaryState()[0];

        //2ème état

        final double o12DEUX = y0[7] * y0[7];
        final double o22DEUX = y0[8] * y0[8];
        final double o32DEUX = y0[9] * y0[9];

        twoEDEUX    =  i1DEUX * o12DEUX + i2DEUX * o22DEUX + i3DEUX * o32DEUX;
        m2DEUX      =  i1DEUX * i1DEUX * o12DEUX + i2DEUX * i2DEUX * o22DEUX + i3DEUX * i3DEUX * o32DEUX;

        final Sorted sortedDEUX = sortInertiaAxis(i1DEUX, i2DEUX, i3DEUX, m2DEUX, twoEDEUX,
                                                  new Vector3D(y0[7], y0[8], y0[9]),
                                                  new Rotation(y0[10], y0[11], y0[12], y0[13], true));
        i1CDEUX = sortedDEUX.i1;
        i2CDEUX = sortedDEUX.i2;
        i3CDEUX = sortedDEUX.i3;

        // convert initial conditions to Euler angles such the M is aligned with Z in computation frame
        final Vector3D omega0BodyDEUX = sortedDEUX.omega;
        final Vector3D m0BodyDEUX     = new Vector3D(i1CDEUX * omega0BodyDEUX.getX(), i2CDEUX * omega0BodyDEUX.getY(), i3CDEUX * omega0BodyDEUX.getZ());

        final double   phi0DEUX       = 0; // this angle can be set arbitrarily, so 0 is a fair value (Eq. 37.13 - 37.14)
        final double   theta0DEUX     = FastMath.acos(m0BodyDEUX.getZ() / m0BodyDEUX.getNorm());
        final double   psi0DEUX       = FastMath.atan2(m0BodyDEUX.getX(), m0BodyDEUX.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!

        //Compute offset rotation between inertial frame aligned with momentum and regular inertial frame
        final Rotation mAlignedToBodyDEUX = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                         phi0DEUX, theta0DEUX, psi0DEUX);

        convertAxesDEUX = sortedDEUX.convertAxes;;

        final Rotation r0ConvertedAxisDEUX = sortedDEUX.rotation;
        mAlignedToInertDEUX = r0ConvertedAxisDEUX.applyInverseTo(mAlignedToBodyDEUX);
        //mAlignedToInert = r0.applyInverseTo(mAlignedToBody);

        i32DEUX  = i3CDEUX - i2CDEUX;
        i31DEUX  = i3CDEUX - i1CDEUX;
        i21DEUX  = i2CDEUX - i1CDEUX;

        tScaleDEUX  = FastMath.sqrt(i32DEUX * (m2DEUX - twoEDEUX * i1CDEUX) / (i1CDEUX * i2CDEUX * i3CDEUX));
        o1ScaleDEUX = FastMath.sqrt((twoEDEUX * i3CDEUX - m2DEUX) / (i1CDEUX * i31DEUX));
        o2ScaleDEUX = FastMath.sqrt((twoEDEUX * i3CDEUX - m2DEUX) / (i2CDEUX * i32DEUX));
        o3ScaleDEUX = FastMath.sqrt((m2DEUX - twoEDEUX * i1CDEUX) / (i3CDEUX * i31DEUX));

        k2DEUX     = i21DEUX * (twoEDEUX * i3CDEUX - m2DEUX) / (i32DEUX * (m2DEUX - twoEDEUX * i1CDEUX));

        jacobiDEUX = JacobiEllipticBuilder.build(k2DEUX);

        if (omega0BodyDEUX.getY() == 0){
            tRefDEUX = t0;
        } else {
            tRefDEUX   = t0 - jacobiDEUX.arcsn(omega0BodyDEUX.getY() / o2ScaleDEUX) / tScaleDEUX;
        }

        periodDEUX             = 4 * LegendreEllipticIntegral.bigK(k2DEUX) / tScaleDEUX;
        phiSlopeDEUX           = FastMath.sqrt(m2DEUX) / i3DEUX;
        phiQuadratureModelDEUX = computePhiQuadratureModel(i1CDEUX, i2CDEUX, i3CDEUX,
                                                           k2DEUX, m2DEUX, jacobiDEUX, periodDEUX, phiSlopeDEUX,
                                                           t0, tRefDEUX, tScaleDEUX);
        integOnePeriodDEUX     = phiQuadratureModelDEUX.getInterpolatedState(phiQuadratureModelDEUX.getFinalTime()).getPrimaryState()[0];

    }//Fin du constructeur

    public static Sorted sortInertiaAxis(final double i1Init, final double i2Init, final double i3Init,
                                         final double m2Init, final double twoEInit,
                                         final Vector3D omegaInit, final Rotation rInit) {

        Sorted sorted = new Sorted(i1Init, i2Init, i3Init, omegaInit, rInit, Rotation.IDENTITY);

        if (sorted.i1 > sorted.i2) {
            final Rotation ij = new Rotation(Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_J, Vector3D.PLUS_I);
            sorted = new Sorted(sorted.i2, sorted.i1, sorted.i3, ij.applyTo(sorted.omega),
                                ij.applyTo(sorted.rotation), ij.applyTo(sorted.convertAxes));
        }

        if (sorted.i2 > sorted.i3) {
            final Rotation jk = new Rotation(Vector3D.PLUS_J, Vector3D.PLUS_K, Vector3D.PLUS_K, Vector3D.PLUS_J);
            sorted = new Sorted(sorted.i2, sorted.i1, sorted.i3, jk.applyTo(sorted.omega),
                                jk.applyTo(sorted.rotation), jk.applyTo(sorted.convertAxes));
        }

        if (sorted.i1 > sorted.i2) {
            final Rotation ij = new Rotation(Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_J, Vector3D.PLUS_I);
            sorted = new Sorted(sorted.i2, sorted.i1, sorted.i3, ij.applyTo(sorted.omega),
                                ij.applyTo(sorted.rotation), ij.applyTo(sorted.convertAxes));
        }

        if (m2Init / twoEInit < sorted.i2) {
            final Rotation ik = new Rotation(Vector3D.PLUS_I, Vector3D.PLUS_K, Vector3D.PLUS_K, Vector3D.PLUS_I);
            sorted = new Sorted(sorted.i2, sorted.i1, sorted.i3, ik.applyTo(sorted.omega),
                                ik.applyTo(sorted.rotation), ik.applyTo(sorted.convertAxes));
        }

        return sorted;

    }

    private static DenseOutputModel computePhiQuadratureModel(final double i1, final double i2, final double i3,
                                                              final double k2, final double m2, final JacobiElliptic jacobi,
                                                              final double period, final double phiSlope,
                                                              final double t0, final double tRef, final double tScale) {

        final double i32  = i3 - i2;
        final double i31  = i3 - i1;
        final double i21  = i2 - i1;

        // coefficients for φ model
        final double b = phiSlope * i32 * i31;
        final double c = i1 * i32;
        final double d = i3 * i21;

        //Integrate the quadrature phi term on one period
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

    public TfmState[] computeTorqueFreeMotion(double t) {

        // Computation of omega
        final CopolarN valuesN = jacobi.valuesN((t - tRef) * tScale);
        final Vector3D omegaP  = new Vector3D(o1Scale * valuesN.cn(), o2Scale * valuesN.sn(), o3Scale * valuesN.dn());
        final Vector3D omega   = convertAxes.applyInverseTo(omegaP);

        // Computation of the Euler angles
        // Compute rotation rate
        final double   psi       = FastMath.atan2(i1C * omegaP.getX(), i2C * omegaP.getY());
        final double   theta     = FastMath.acos(omegaP.getZ() / phiSlope);
        final double   phiLinear = phiSlope * t;

        // Integration for the computation of phi
        final double t0 = getInitialTime();
        final int nbPeriods = (int) FastMath.floor((t - t0) / period);//floor = entier inférieur = nb période entière
        final double tStartInteg = t0 + nbPeriods * period;//partie de période à la fin entre tau Integ et tau end
        final double integPartial = phiQuadratureModel.getInterpolatedState(t - tStartInteg).getPrimaryState()[0];// a vérifier, partie de l'intégrale apres le nb entier de période
        final double phiQuadrature = nbPeriods * integOnePeriod + integPartial;

        final double phi = phiLinear + phiQuadrature;

        // Computation of omega
        final CopolarN valuesNDEUX = jacobiDEUX.valuesN((t - tRefDEUX) * tScaleDEUX);
        final Vector3D omegaPDEUX = new Vector3D(o1ScaleDEUX * valuesNDEUX.cn(), o2ScaleDEUX * valuesNDEUX.sn(), o3ScaleDEUX * valuesNDEUX.dn());
        final Vector3D omegaDEUX  = convertAxesDEUX.applyInverseTo(omegaPDEUX);

        // Computation of the Euler angles
        final double   psiDEUX           = FastMath.atan2(i1CDEUX * omegaPDEUX.getX(), i2CDEUX * omegaPDEUX.getY());
        final double   thetaDEUX         = FastMath.acos(omegaPDEUX.getZ() / phiSlopeDEUX);
        final double   phiLinearDEUX     = phiSlopeDEUX * t;

        //Integration for the computation of phi
        final int nbPeriodsDEUX = (int) FastMath.floor((t - t0) / periodDEUX);//floor = entier inférieur = nb période entière
        final double tStartIntegDEUX = t0 + nbPeriodsDEUX * periodDEUX;//partie de période à la fin entre tau Integ et tau end
        final double integPartialDEUX = phiQuadratureModelDEUX.getInterpolatedState(t - tStartIntegDEUX).getPrimaryState()[0];// a vérifier, partie de l'intégrale apres le nb entier de période
        final double phiQuadratureDEUX = nbPeriodsDEUX * integOnePeriodDEUX + integPartialDEUX;

        final double phiDEUX = phiLinearDEUX + phiQuadratureDEUX;

        // Computation of the quaternion

        // Rotation between computation frame (aligned with momentum) and body
        //(It is simply the angles equations provided by L&L)
        final Rotation alignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                                    phi, theta, psi);

        // combine with offset rotation to get back to regular inertial frame
        //Inert -> aligned + aligned -> body = inert -> body (What the user wants)
        Rotation inertToBody = alignedToBody.applyTo(mAlignedToInert.revert());//alignedToBody.applyTo(mAlignedToInert.revert());
        Rotation bodyToOriginalFrame = convertAxes.applyInverseTo(inertToBody);//(inertToBody.applyInverseTo(convertAxes)).revert();

        //Computation of the quaternion

        // Rotation between computation frame (aligned with momentum) and body
        //(It is simply the angles equations provided by L&L)
        final Rotation alignedToBodyDEUX = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                                        phiDEUX, thetaDEUX, psiDEUX);

        // combine with offset rotation to get back to regular inertial frame
        //Inert -> aligned + aligned -> body = inert -> body (What the user wants)
        Rotation inertToBodyDEUX = alignedToBodyDEUX.applyTo(mAlignedToInertDEUX.revert());//alignedToBody.applyTo(mAlignedToInert.revert());
        Rotation bodyToOriginalFrameDEUX = convertAxesDEUX.applyInverseTo(inertToBodyDEUX);//(inertToBody.applyInverseTo(convertAxes)).revert();

        return new TfmState[] {
            new TfmState(t, omega,     bodyToOriginalFrame,     phi,     theta,     psi,     convertAxes,     mAlignedToInert),
            new TfmState(t, omegaDEUX, bodyToOriginalFrameDEUX, phiDEUX, thetaDEUX, psiDEUX, convertAxesDEUX, mAlignedToInertDEUX)
        };

    }

    public double[] doComputeDerivatives(double t, double[] y) {

        final  double[] yDot = new double[getDimension()];

        // compute the derivatives using Euler equations
        yDot[0] = y[1] * y[2] * (i2 - i3) / i1;
        yDot[1] = y[2] * y[0] * (i3 - i1) / i2;
        yDot[2] = y[0] * y[1] * (i1 - i2) / i3;

        // compute the derivatives using Qpoint = 0.5 * Omega_inertialframe * Q
        yDot[3] = 0.5 * (-y[0] * y[4] -y[1] * y[5] -y[2] * y[6]);
        yDot[4] = 0.5 * (y[0] * y[3] +y[2] * y[5] -y[1] * y[6]);
        yDot[5] = 0.5 * (y[1] * y[3] -y[2] * y[4] +y[0] * y[6]);
        yDot[6] = 0.5 * (y[2] * y[3] +y[1] * y[4] -y[0] * y[5]);

        // compute the derivatives using Euler equations
        yDot[7] = y[8] * y[9] * (i2DEUX - i3DEUX) / i1DEUX;
        yDot[8] = y[9] * y[7] * (i3DEUX - i1DEUX) / i2DEUX;
        yDot[9] = y[7] * y[8] * (i1DEUX - i2DEUX) / i3DEUX;

        // compute the derivatives using Qpoint = 0.5 * Omega_inertialframe * Q
        yDot[10] = 0.5 * (-y[7] * y[11] -y[8] * y[12] -y[9] * y[13]);
        yDot[11] = 0.5 * (y[7] * y[10] +y[9] * y[12] -y[8] * y[13]);
        yDot[12] = 0.5 * (y[8] * y[10] -y[9] * y[11] +y[7] * y[13]);
        yDot[13] = 0.5 * (y[9] * y[10] +y[8] * y[11] -y[7] * y[12]);
        return yDot;

    }

    public double[] computeTheoreticalState(double t) {

        final TfmState[] tfm = computeTorqueFreeMotion(t);

        return new double[] {
            tfm[0].getOmega().getX(),
            tfm[0].getOmega().getY(),
            tfm[0].getOmega().getZ(),
            tfm[0].getRotation().getQ0(),
            tfm[0].getRotation().getQ1(),
            tfm[0].getRotation().getQ2(),
            tfm[0].getRotation().getQ3(),

            tfm[1].getOmega().getX(),
            tfm[1].getOmega().getY(),
            tfm[1].getOmega().getZ(),
            tfm[1].getRotation().getQ0(),
            tfm[1].getRotation().getQ1(),
            tfm[1].getRotation().getQ2(),
            tfm[1].getRotation().getQ3()
        };
    }

    private static class Sorted {
        private final double   i1;
        private final double   i2;
        private final double   i3;
        private final Vector3D omega;
        private final Rotation rotation;
        private final Rotation convertAxes;
        private Sorted(final double i1, final double i2, final double i3,
                      final Vector3D omega, final Rotation rotation, final Rotation convertAxes) {
            this.i1              = i1;
            this.i2              = i2;
            this.i3              = i3;
            this.omega           = omega;
            this.rotation        = rotation;
            this.convertAxes     = convertAxes;
        }
    }

    public static class TfmState {
        private final double   t;
        private final Vector3D omega;
        private final Rotation rotation;
        private final double phi;
        private final double theta;
        private final double psi;
        private final Rotation convertAxes;
        private final Rotation mAlignedToInert;
        private TfmState(final double t, final Vector3D omega, final Rotation rotation,
                         final double phi, final double theta, final double psi,
                         final Rotation convertAxes, final Rotation mAlignedToInert) {
            this.t               = t;
            this.omega           = omega;
            this.rotation        = rotation;
            this.phi             = phi;
            this.theta           = theta;
            this.psi             = psi;
            this.convertAxes     = convertAxes;
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
        public Rotation getConvertAxes() {
            return convertAxes;
        }
        public Rotation getMAlignedToInert() {
            return mAlignedToInert;
        }
    }

}//Fin du programme
