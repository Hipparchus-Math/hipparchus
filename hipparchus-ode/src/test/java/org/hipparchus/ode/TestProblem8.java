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
import org.junit.Assert;


public class TestProblem8 extends TestProblemAbstract {

	/** Moments of inertia. */
	final double i1;
	final double i2;
	final double i3;

	/** Twice the angular kinetic energy. */
	final double twoE;

	/** Square of kinetic momentum. */
	final double m2;
	/**
	 * Simple constructor.
	 */
	public TestProblem8() {
		//Arguments in the super constructor :
		//Intital time, Primary state (o1, o2, o3, q0, q1, q2, q3), Final time, Error scale
		super(0.0, new double[] {5.0, 0.0, 4.0, 0.9, 0.437, 0.0, 0.0}, 20.0, new double[] { 1.0, 1.0, 1.0 });
		i2 = 3.0 / 8.0;
		i1 = 1.0 / 2.0;
		i3 = 5.0 / 8.0;


		final double[] y0 = getInitialState().getPrimaryState();

		final double o12 = y0[0] * y0[0];
		final double o22 = y0[1] * y0[1];
		final double o32 = y0[2] * y0[2];

		twoE    =  i1 * o12 + i2 * o22 + i3 * o32;
		m2      =  i1 * i1 * o12 + i2 * i2 * o22 + i3 * i3 * o32;

	}

	public double[][] computeTorqueFreeMotion(double i1P, double i2P, double i3P, double t0, double[] y0, double t) {



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

		double[] Y0 = {y0[0],y0[1],y0[2], y0[3],y0[4],y0[5],y0[6]};
		double[] I = {i1P, i2P, i3P};
		
		Vector3D[] axes = {Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_K};

		double[] i = {i1P, i2P, i3P};



		System.out.println("omega avant : "+y0[0]+" "+ y0[1]+" "+ y0[2]);
		System.out.println("Initial : iA1 = "+i1P+" "+axes[0]);
		System.out.println("iA2 = "+i2P + " "+axes[1]);
		System.out.println("iA3 = "+i3P + " "+axes[2]);
		if (i[0] > i[1]) {
			Vector3D z = axes[0];
			axes[0] = axes[1];
			axes[1] = z;
			axes[2] = axes[2].negate();

			double y = i[0];
			i[0] = i[1];
			i[1] = y;

			double v = y0[0];
			y0[0] = y0[1];
			y0[1] = v;
			y0[2] = -y0[2];
			System.out.println("1ere boucle : iA1 = "+i[0]+" "+axes[0]);
			System.out.println("iA2 = "+i[1] + " "+axes[1]);
			System.out.println("iA3 = "+i[2] + " "+axes[2]);
		}

		if (i[1] > i[2]) {
			Vector3D z = axes[1];
			axes[1] = axes[2];
			axes[2] = z;
			axes[0] = axes[0].negate();

			double y = i[1];
			i[1] = i[2];
			i[2] = y;

			double v = y0[1];
			y0[1] = y0[2];
			y0[2] = v;
			y0[0] = -y0[0];
			System.out.println("2ere boucle : iA1 = "+i[0]+" "+axes[0]);
			System.out.println("iA2 = "+i[1] + " "+axes[1]);
			System.out.println("iA3 = "+i[2] + " "+axes[2]);
		}

		if (i[0] > i[1]) {
			Vector3D z = axes[0];
			axes[0] = axes[1];
			axes[1] = z;
			axes[2] = axes[2].negate();

			double y = i[0];
			i[0] = i[1];
			i[1] = y;

			double v = y0[0];
			y0[0] = y0[1];
			y0[1] = v;
			y0[2] = -y0[2];
			System.out.println("3ere boucle : iA1 = "+i[0]+" "+axes[0]);
			System.out.println("iA2 = "+i[1] + " "+axes[1]);
			System.out.println("iA3 = "+i[2] + " "+axes[2]);
		}

		final double condition;
		if( y0[0] == 0 && y0[1] == 0 && y0[2] == 0) {
			condition = 0.0;
		}else {
			condition = m2/twoE;
		}

		System.out.println("Condition : "+condition);
		if(condition < i[1]) {
			Vector3D z = axes[0];
			axes[0] = axes[2];
			axes[2] = z;
			axes[1] = axes[1].negate();

			double y = i[0];
			i[0] = i[2];
			i[2] = y;

			double v = y0[0];
			y0[0] = y0[2];
			y0[2] = v;
			y0[1] = - y0[1];
			System.out.println("repère final  : iA1 = "+i[0]+" "+axes[0]);
			System.out.println("iA2 = "+i[1] + " "+axes[1]);
			System.out.println("iA3 = "+i[2] + " "+axes[2]);
		}
		System.out.println("omega après : "+y0[0]+" "+ y0[1]+" "+ y0[2]);
		System.out.println("axes après : "+axes[0]+" "+ axes[1]+" "+ axes[2]);
		System.out.println("inerties après : "+i[0]+" "+ i[1]+" "+ i[2]);

		final double i1 = i[0];
		final double i2 = i[1];
		final double i3 = i[2];

		final double i32  = i3 - i2;
		final double i31  = i3 - i1;
		final double i21  = i2 - i1;

		tScale  = FastMath.sqrt(i32 * (m2 - twoE * i1) / (i1 * i2 * i3));
		o1Scale = FastMath.sqrt((twoE * i3 - m2) / (i1 * i31));
		o2Scale = FastMath.sqrt((twoE * i3 - m2) / (i2 * i32));
		o3Scale = FastMath.sqrt((m2 - twoE * i1) / (i3 * i31));

		if( y0[0] == 0 && y0[1] == 0 && y0[2] == 0) {
			k2 = 0.0;
		} else {
			k2     = i21 * (twoE * i3 - m2) / (i32 * (m2 - twoE * i1));
		}

		jacobi = JacobiEllipticBuilder.build(k2);

		if( y0[0] == 0 && y0[1] == 0 && y0[2] == 0) {
			tRef = t0;
		}
		else if (y0[1] == 0){
			tRef = t0 - jacobi.arcsn(0) / tScale;
			System.out.println("Tref0 : "+tRef);
		} else {
			tRef   = t0 - jacobi.arcsn(y0[1] / o2Scale) / tScale;
			System.out.println("Tref : "+tRef);
		}

		// convert initial conditions to Euler angles such the M is aligned with Z in computation frame
		final Vector3D omega0Body = new Vector3D(y0[0], y0[1], y0[2]);
		final Rotation r0         = new Rotation(y0[3], y0[4], y0[5], y0[6], true);
		final Vector3D m0Body     = new Vector3D(i1 * omega0Body.getX(), i2 * omega0Body.getY(), i3 * omega0Body.getZ());

		
//		final Vector3D omega0Body = new Vector3D(Y0[0], Y0[1], Y0[2]);
//		final Rotation r0         = new Rotation(Y0[3], Y0[4], Y0[5], Y0[6], true);
//		final Vector3D m0Body     = new Vector3D(I[0] * omega0Body.getX(), I[1] * omega0Body.getY(), I[2] * omega0Body.getZ());

		final double   phi0       = 0; // this angle can be set arbitrarily, so 0 is a fair value
		final double   theta0;
		if( y0[0] == 0 && y0[1] == 0 && y0[2] == 0) {
			theta0 = 0.0;
		} else {
			theta0 =  FastMath.acos(m0Body.getZ() / m0Body.getNorm());
		}
		final double   psi0       = FastMath.atan2(m0Body.getX(), m0Body.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!

		
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
				if( y0[0] == 0 && y0[1] == 0 && y0[2] == 0) {
					return new double[] {
							0.0
					};
				} else {
					return new double[] {
							b / (c + d * sn * sn)
					};
				}
			}

		}, new ODEState(t0, new double[1]), t0 + period);

		integOnePeriod = phiQuadratureModel.getInterpolatedState(phiQuadratureModel.getFinalTime()).getPrimaryState()[0];


		//Compute offset rotation between inertial frame aligned with momentum and regular inertial frame
		final Rotation mAlignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
				phi0, theta0, psi0);

		Rotation convertAxes = new Rotation( Vector3D.PLUS_I, Vector3D.PLUS_J, axes[0], axes[1] );

		Vector3D axe1 = new Vector3D(1.0, 0.0, 0.0);
		Vector3D axe1c = convertAxes.applyTo(axe1);
		System.out.println("AAAAAA : "+axe1c.getX()+" "+axe1c.getY()+" "+axe1c.getZ());
		Vector3D axe2 = new Vector3D(0.0, 1.0, 0.0);
		Vector3D axe2c = convertAxes.applyTo(axe2);
		System.out.println("AAAAAA : "+axe2c.getX()+" "+axe2c.getY()+" "+axe2c.getZ());
		Vector3D axe3 = new Vector3D(0.0, 0.0,1.0);
		Vector3D axe3c = convertAxes.applyTo(axe3);
		System.out.println("AAAAAA : "+axe3c.getX()+" "+axe3c.getY()+" "+axe3c.getZ());

		System.out.println("Quaternion avant : "+r0.getQ0()+" "+r0.getQ1()+" "+r0.getQ2()+" "+r0.getQ3());
		
		Rotation r0ConvertedAxis = convertAxes.applyTo(r0);
		
		System.out.println("Quaternion converti : "+r0ConvertedAxis.getQ0()+" "+r0ConvertedAxis.getQ1()+" "+r0ConvertedAxis.getQ2()+" "+r0ConvertedAxis.getQ3());

		System.out.println("Axe : "+r0ConvertedAxis.getAxis(RotationConvention.FRAME_TRANSFORM));
		
		mAlignedToInert = r0ConvertedAxis.applyInverseTo(mAlignedToBody);

		
		//Computation of omega
		final CopolarN valuesN = jacobi.valuesN((t - tRef) * tScale);

		final Vector3D omegaP = new Vector3D(
				o1Scale * valuesN.cn(),
				o2Scale * valuesN.sn(),
				o3Scale * valuesN.dn());

		final Vector3D omega = convertAxes.applyTo(omegaP);

		System.out.println("omega fonction : "+o1Scale * valuesN.cn()+" "+o2Scale * valuesN.sn()+" "+o3Scale * valuesN.dn());
		System.out.println("omega fonction conversion : "+omega.getX()+" "+omega.getY()+" "+omega.getZ());

		//Computation of the euler angles
		//Compute rotation rate
		final double   o1            = omega.getX();//o1Scale * valuesN.cn();
		final double   o2            = omega.getY();//o2Scale * valuesN.sn();
		final double   o3            = omega.getZ();//o3Scale * valuesN.dn();

		//		double thetaroot = FastMath.sqrt((i3*(m2-twoE*i1)) / (m2*i31));
		//		double psiroot = FastMath.sqrt((i1*i32) / (i2*i31));
		//		psiroot * (valuesN.cn()/valuesN.sn());
		//		thetaroot * valuesN.dn();

		//Compute angles
		final double   psi           = FastMath.atan2(i1 * o1, i2 * o2);
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
		//(It is simply the angles equations provided by L&L)
		final Rotation alignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
				phi, theta, psi);

		// combine with offset rotation to get back to regular inertial frame
		//Inert -> aligned + aligned -> body = inert -> body (What the user wants)
		Rotation inertToBody = alignedToBody.applyTo(mAlignedToInert.revert());

		Rotation convertAxesReverse = new Rotation(axes[0], axes[1], Vector3D.PLUS_I, Vector3D.PLUS_J);

		Vector3D axe1cr = convertAxes.applyInverseTo(axes[0]);
		System.out.println("BBBBBB : "+axe1cr.getX()+" "+axe1cr.getY()+" "+axe1cr.getZ());
		Vector3D axe2cr = convertAxes.applyInverseTo(axes[1]);
		System.out.println("BBBBBB : "+axe2cr.getX()+" "+axe2cr.getY()+" "+axe2cr.getZ());
		Vector3D axe3cr = convertAxes.applyInverseTo(axes[2]);
		System.out.println("BBBBBB : "+axe3cr.getX()+" "+axe3cr.getY()+" "+axe3cr.getZ());

		Vector3D axe = inertToBody.getAxis(RotationConvention.FRAME_TRANSFORM).negate();
		double angle = inertToBody.getAngle();
		Rotation bodyToInert = new Rotation(axe, angle, RotationConvention.FRAME_TRANSFORM);
		
		Rotation originalFrame = convertAxes.applyInverseTo(inertToBody);

		double[] angles = originalFrame.getAngles(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM);
		double[][] data = {{omega.getX(), omega.getY(), omega.getZ()}, {angles[0], angles[1], angles[2]}, {originalFrame.getQ0(), originalFrame.getQ1(), originalFrame.getQ2(), originalFrame.getQ3()}, {axes[0].getX(), axes[0].getY(), axes[0].getZ()}, {axes[1].getX(), axes[1].getY(), axes[1].getZ()},{i1, i2, i3}};
		return data;
	}

	public double[] doComputeDerivatives(double t, double[] y) {

		final  double[] yDot = new double[getDimension()];
		double[] y0 = getInitialState().getPrimaryState();
		double t0 = getInitialState().getTime();

		//		System.out.println("Numerique omega : "+y[0]+" "+y[1]+" "+y[2]);
		//		final double[][] tfm = computeTorqueFreeMotion(i1, i2, i3, t0, y0, t);
		//		final double[] omega = tfm[0];
		//		final double[] I = tfm[5];
		//		final double[] q = tfm[2];
		//
		//		// compute the derivatives using Euler equations
		//		yDot[0] = omega[1] * omega[2] * (I[1] - I[2]) / I[0];
		//		yDot[1] = omega[2] * omega[0] * (I[2] - I[0]) / I[1];
		//		yDot[2] = omega[0] * omega[1] * (I[0] - I[1]) / I[2];
		//
		//		// compute the derivatives using Qpoint = 0.5 * Omega_inertialframe * Q
		//		yDot[3] = 0.5 * (-omega[0] * q[1] -omega[1] * q[2] -omega[2] * q[3]);
		//		yDot[4] = 0.5 * (omega[0] * q[0] +omega[2] * q[2] -omega[1] * q[3]);
		//		yDot[5] = 0.5 * (omega[1] * q[0] -omega[2] * q[1] +omega[0] * q[3]);
		//		yDot[6] = 0.5 * (omega[2] * q[0] +omega[1] * q[1] -omega[0] * q[2]);

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

		final double[][] tfm = computeTorqueFreeMotion(i1, i2, i3, t0, y0, t);
		final double[] omega = tfm[0];
		final double[] quaternion = tfm[2];
		final double[] angles = tfm[1];
		final double[] X = tfm[3];
		final double[] Y = tfm[4];

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
				angles[2],
				X[0],
				X[1],
				X[2],
				Y[0],
				Y[1],
				Y[2]
		};
	}

}//Fin du programme
