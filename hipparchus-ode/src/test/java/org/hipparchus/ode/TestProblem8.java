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

	/**
	 * Simple constructor.
	 */
	public TestProblem8() {
		//Arguments in the super constructor :
		//Intital time, Primary state (t0, o1, o2, o3, q0, q1, q2, q3), Final time, Error scale

		super(0.0, new double[] {FastMath.sqrt(80/3), 0.0, 4.0, 0.9999999, 0.0, 0.0, 0.0}, 20.0, new double[] { 1.0, 1.0, 1.0 });
		i1 = 3.0 / 8.0;
		i2 = 1.0 / 2.0;
		i3 = 5.0 / 8.0;
	}

	public double[][] TorqueFreeMotion(double i1, double i2, double i3, double t0, double[] y0, double t) {

		/** Twice the angular kinetic energy. */
		final double twoE;

		/** Square of kinetic momentum. */
		final double m2;

		/** Time scaling factor. */
		final double tScale;

		/** Time reference for rotation rate. */
		final double tRef;

		/** Offset rotation from inertial computation frame aligned with momentum to regular inertial frame. */
		Rotation mAlignedToInert;

		/** Coefficients for φ model. */
		final double phiSlope;
		final double b;
		final double c;
		final double d;

		/**Period with respect to time. */
		final double period;

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


		/**DenseOutputModel of phi. */
		final DenseOutputModel phiQuadratureModel;
		final double integOnePeriod;

		final double i32  = i3 - i2;
		final double i31  = i3 - i1;
		final double i21  = i2 - i1;

		final double o12 = y0[0] * y0[0];
		final double o22 = y0[1] * y0[1];
		final double o32 = y0[2] * y0[2];

		twoE    =  i1 * o12 + i2 * o22 + i3 * o32;
		m2      =  i1 * i1 * o12 + i2 * i2 * o22 + i3 * i3 * o32;


		tScale  = FastMath.sqrt(i32 * (m2 - twoE * i1) / (i1 * i2 * i3));
		o1Scale = FastMath.sqrt((twoE * i3 - m2) / (i1 * i31));
		o2Scale = FastMath.sqrt((twoE * i3 - m2) / (i2 * i32));
		o3Scale = FastMath.sqrt((m2 - twoE * i1) / (i3 * i31));

		k2     = i21 * (twoE * i3 - m2) / (i32 * (m2 - twoE * i1));

		jacobi = JacobiEllipticBuilder.build(k2);
		tRef   = t0 - jacobi.arcsn(y0[1] / o2Scale) / tScale;


		// convert initial conditions to Euler angles such the M is aligned with Z in computation frame
		final Vector3D omega0Body = new Vector3D(y0[0], y0[1], y0[2]);
		final Rotation r0         = new Rotation(y0[3], y0[4], y0[5], y0[6], true);
		final Vector3D m0Body     = new Vector3D(i1 * omega0Body.getX(), i2 * omega0Body.getY(), i3 * omega0Body.getZ());
		final double   phi0       = 0; // this angle can be set arbitrarily, so 0 is a fair value
		final double   theta0     = FastMath.acos(m0Body.getZ() / m0Body.getNorm());
		final double   psi0       = FastMath.atan2(m0Body.getX(), m0Body.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!

		// compute offset rotation between inertial frame aligned with momentum and regular inertial frame
		final Rotation mAlignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
				phi0, theta0, psi0);
		mAlignedToInert = r0.applyInverseTo(mAlignedToBody);

		// coefficients for φ model
		period = 4 * LegendreEllipticIntegral.bigK(k2) / tScale;
		phiSlope = FastMath.sqrt(m2) / i3;// =a, pente de la partie linéaire
		b = phiSlope * i32 * i31;
		c = i1 * i32;
		d = i3 * i21;

		//Integrate the quadrature phi term on one period
		final DormandPrince853Integrator integ = new DormandPrince853Integrator(1.0e-6 * period, 1.0e-2 * period,
				phiSlope * period * 1.0e-13, 1.0e-13);
		phiQuadratureModel = new DenseOutputModel();
		integ.addStepHandler(phiQuadratureModel);
		integ.integrate(new OrdinaryDifferentialEquation() {

			/** {@inheritDoc} */
			@Override
			public int getDimension() {
				return 1;
			}
			@Override
			public double[] computeDerivatives(final double t, final double[] y) {
				final double sn = jacobi.valuesN((t - tRef) * tScale).sn();
				return new double[] {
						b / (c + d * sn * sn)
				};
			}

		}, new ODEState(t0, new double[1]), t0 + period);

		integOnePeriod = phiQuadratureModel.getInterpolatedState(phiQuadratureModel.getFinalTime()).getPrimaryState()[0];


		//Computation of omega
		final CopolarN valuesN = jacobi.valuesN((t - tRef) * tScale);

		final Vector3D omega = new Vector3D(
				o1Scale * valuesN.cn(),
				o2Scale * valuesN.sn(),
				o3Scale * valuesN.dn());

		//Computation of the euler angles
		// Compute rotation rate
		final double   o1            = o1Scale * valuesN.cn();
		final double   o2            = o2Scale * valuesN.sn();
		final double   o3            = o3Scale * valuesN.dn();

		// compute angles
		final double psi         = FastMath.atan2(i1 * o1, i2 * o2);
		final double   theta         = FastMath.acos(o3 / phiSlope);
		final double   phiLinear     = phiSlope * t;

		//Integration for the computation of phi
		final int nbPeriods = (int) FastMath.floor((t - t0) / period);//floor = entier inférieur = nb période entière
		final double tStartInteg = t0 + nbPeriods * period;//partie de période à la fin entre tau Integ et tau end
		final double integPartial = phiQuadratureModel.getInterpolatedState(t - tStartInteg).getPrimaryState()[0];// a vérifier, partie de l'intégrale apres le nb entier de période
		final double phiQuadrature = nbPeriods * integOnePeriod + integPartial;

		final double phi = phiLinear + phiQuadrature;


		//Computation of the quaternion
		// Rotation between computation frame (aligned with momentum) and body
		final Rotation alignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
				phi, theta, psi);

		// combine with offset rotation to get back to regular inertial frame
		final Rotation quaternion = alignedToBody.applyTo(mAlignedToInert.revert());


		double[][] data = {{omega.getX(), omega.getY(), omega.getZ()}, {phi, theta, psi}, {quaternion.getQ0(), quaternion.getQ1(), quaternion.getQ2(), quaternion.getQ3()}};
		return data;
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
		return yDot;

	}

	public double[] computeTheoreticalState(double t) {

		double[] y0 = getInitialState().getPrimaryState();
		double t0 = getInitialState().getTime();

		final double[][] TFM = TorqueFreeMotion(i1, i2, i3, t0, y0, t);
		final double[] omega = TFM[0];
		final double[] quaternion = TFM[2];
		final double[] angles = TFM[1];


		return new double[] {
				omega[0],
				omega[1],
				omega[2],
				quaternion[0],
				quaternion[1],
				quaternion[2],
				quaternion[3],
				angles[0],
				angles[1],
				angles[2]
		};
	}

}//Fin du programme






//    //Rotation in the body frame
//    final Rotation alignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
//    phi, theta, psi);
//
//    //Rotation in the user random frame ?
//    final Rotation inertToBody = alignedToBody.applyTo(mAlignedToInert.revert());
