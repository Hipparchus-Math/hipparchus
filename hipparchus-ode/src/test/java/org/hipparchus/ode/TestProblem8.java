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
	final double tScale;

	/** Time reference for rotation rate. */
	final double tRef;

	/**Offset rotation  between initial inertial frame and the frame with moment vector and Z axis aligned. */
	Rotation mAlignedToInert;

	/** Initial converted rotation. */
	final Rotation r0;

	/** Rotation to switch to the converted axes frame. */
	final Rotation convertAxes;

	/**
	 * Simple constructor.
	 */
	public TestProblem8() {
		//Arguments in the super constructor :
		//Intital time, Primary state (o1, o2, o3, q0, q1, q2, q3), Final time, Error scale
		super(0.0, new double[] {5.0, 0.0, 4.0, 0.9, 0.437, 0.0, 0.0}, 20.0, new double[] { 1.0, 1.0, 1.0 });
		i1 = 3.0 / 8.0;
		i2 = 1.0 / 2.0;
		i3 = 5.0 / 8.0;


		y0 = getInitialState().getPrimaryState();
		final double t0 = getInitialState().getTime();

		final double o12 = y0[0] * y0[0];
		final double o22 = y0[1] * y0[1];
		final double o32 = y0[2] * y0[2];

		twoE    =  i1 * o12 + i2 * o22 + i3 * o32;
		m2      =  i1 * i1 * o12 + i2 * i2 * o22 + i3 * i3 * o32;

		Vector3D[] axesP = {Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_K};

		double[] i = {i1, i2, i3};
		double[] y0P = y0;

		System.out.println("omega avant : "+y0[0]+" "+ y0[1]+" "+ y0[2]);
		System.out.println("Initial : iA1 = "+i1+" "+axesP[0]);
		System.out.println("iA2 = "+i2 + " "+axesP[1]);
		System.out.println("iA3 = "+i3 + " "+axesP[2]);
		if (i[0] > i[1]) {
			Vector3D z = axesP[0];
			axesP[0] = axesP[1];
			axesP[1] = z;
			axesP[2] = axesP[2].negate();

			double y = i[0];
			i[0] = i[1];
			i[1] = y;

			double v = y0P[0];
			y0P[0] = y0P[1];
			y0P[1] = v;
			y0P[2] = -y0P[2];
			System.out.println("1ere boucle : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}

		if (i[1] > i[2]) {
			Vector3D z = axesP[1];
			axesP[1] = axesP[2];
			axesP[2] = z;
			axesP[0] = axesP[0].negate();

			double y = i[1];
			i[1] = i[2];
			i[2] = y;

			double v = y0P[1];
			y0P[1] = y0P[2];
			y0P[2] = v;
			y0P[0] = -y0P[0];
			System.out.println("2ere boucle : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}

		if (i[0] > i[1]) {
			Vector3D z = axesP[0];
			axesP[0] = axesP[1];
			axesP[1] = z;
			axesP[2] = axesP[2].negate();

			double y = i[0];
			i[0] = i[1];
			i[1] = y;

			double v = y0P[0];
			y0P[0] = y0P[1];
			y0P[1] = v;
			y0P[2] = -y0P[2];
			System.out.println("3ere boucle : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}

		final double condition;
		if( y0P[0] == 0 && y0P[1] == 0 && y0P[2] == 0) {
			condition = 0.0;
		}else {
			condition = m2/twoE;
		}

		System.out.println("Condition : "+condition);
		if(condition < i[1]) {
			Vector3D z = axesP[0];
			axesP[0] = axesP[2];
			axesP[2] = z;
			axesP[1] = axesP[1].negate();

			double y = i[0];
			i[0] = i[2];
			i[2] = y;

			double v = y0P[0];
			y0P[0] = y0P[2];
			y0P[2] = v;
			y0P[1] = - y0P[1];
			System.out.println("repère final  : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}
		System.out.println("omega après : "+y0P[0]+" "+ y0P[1]+" "+ y0P[2]);
		System.out.println("axes après : "+axesP[0]+" "+ axesP[1]+" "+ axesP[2]);
		System.out.println("inerties après : "+i[0]+" "+ i[1]+" "+ i[2]);

		i1C = i[0];
		i2C = i[1];
		i3C = i[2];

		axes = axesP;

		y0C = y0P;

		// convert initial conditions to Euler angles such the M is aligned with Z in computation frame
		final Vector3D omega0Body = new Vector3D(y0C[0], y0C[1], y0C[2]);
		r0         = new Rotation(y0C[3], y0C[4], y0C[5], y0C[6], true);
		final Vector3D m0Body     = new Vector3D(i1C * omega0Body.getX(), i2C * omega0Body.getY(), i3C * omega0Body.getZ());

		final double   phi0       = 0; // this angle can be set arbitrarily, so 0 is a fair value (Eq. 37.13 - 37.14)
		final double   theta0 =  FastMath.acos(m0Body.getZ() / m0Body.getNorm());
		final double   psi0       = FastMath.atan2(m0Body.getX(), m0Body.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!

		//Compute offset rotation between inertial frame aligned with momentum and regular inertial frame
		final Rotation mAlignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
				phi0, theta0, psi0);

		convertAxes = new Rotation( Vector3D.PLUS_I, Vector3D.PLUS_J, axesP[0], axesP[1] );

		Rotation r0ConvertedAxis = convertAxes.applyTo(r0);

		//Est-il nécéssaire de garder le r0COnvertedAxis qui n'est peut être pas adapté et ne règle peut être pas le problème
		mAlignedToInert = r0ConvertedAxis.applyInverseTo(mAlignedToBody);
		//mAlignedToInert = r0.applyInverseTo(mAlignedToBody);		



		i32  = i3C - i2C;
		i31  = i3C - i1C;
		i21  = i2C - i1C;

		tScale  = FastMath.sqrt(i32 * (m2 - twoE * i1C) / (i1C * i2C * i3C));
		o1Scale = FastMath.sqrt((twoE * i3C - m2) / (i1C * i31));
		o2Scale = FastMath.sqrt((twoE * i3C - m2) / (i2C * i32));
		o3Scale = FastMath.sqrt((m2 - twoE * i1C) / (i3C * i31));

		if( y0[0] == 0 && y0[1] == 0 && y0[2] == 0) {
			k2 = 0.0;
		} else {
			k2     = i21 * (twoE * i3C - m2) / (i32 * (m2 - twoE * i1C));
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
		System.out.println((y0[1] / o2Scale));


	}

	public double[][] computeTorqueFreeMotion(double i1, double i2, double i3, double t0, double[] y0, double t) {





		/** Coefficients for φ model. */
		final double phiSlope;
		final double b;
		final double c;
		final double d;

		/**Period with respect to time. */
		final double period;

		/**DenseOutputModel of phi. */
		final DenseOutputModel phiQuadratureModel;
		final double integOnePeriod;

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
		},new ODEState(t0, new double[1]), t0 + period);

		integOnePeriod = phiQuadratureModel.getInterpolatedState(phiQuadratureModel.getFinalTime()).getPrimaryState()[0];


		//		Rotation convertAxes = new Rotation( Vector3D.PLUS_I, Vector3D.PLUS_J, axes[0], axes[1] );
		//
		//		Vector3D axe1 = new Vector3D(1.0, 0.0, 0.0);
		//		Vector3D axe1c = convertAxes.applyTo(axe1);
		//		System.out.println("AAAAAA : "+axe1c.getX()+" "+axe1c.getY()+" "+axe1c.getZ());
		//		Vector3D axe2 = new Vector3D(0.0, 1.0, 0.0);
		//		Vector3D axe2c = convertAxes.applyTo(axe2);
		//		System.out.println("AAAAAA : "+axe2c.getX()+" "+axe2c.getY()+" "+axe2c.getZ());
		//		Vector3D axe3 = new Vector3D(0.0, 0.0,1.0);
		//		Vector3D axe3c = convertAxes.applyTo(axe3);
		//		System.out.println("AAAAAA : "+axe3c.getX()+" "+axe3c.getY()+" "+axe3c.getZ());
		//
		//		System.out.println("Quaternion avant : "+r0.getQ0()+" "+r0.getQ1()+" "+r0.getQ2()+" "+r0.getQ3());
		//
		//		Rotation r0ConvertedAxis = convertAxes.applyTo(r0);
		//
		//		System.out.println("Quaternion converti : "+r0ConvertedAxis.getQ0()+" "+r0ConvertedAxis.getQ1()+" "+r0ConvertedAxis.getQ2()+" "+r0ConvertedAxis.getQ3());
		//
		//		System.out.println("Axe : "+r0ConvertedAxis.getAxis(RotationConvention.FRAME_TRANSFORM));
		//
		//		//Rotation r0Aligned = mAlignedToBody.applyInverseTo(r0ConvertedAxis);
		//		
		//		mAlignedToInert = r0ConvertedAxis.applyInverseTo(mAlignedToBody);


		//Computation of omega
		final CopolarN valuesN = jacobi.valuesN((t - tRef) * tScale);

		final Vector3D omegaP = new Vector3D(
				o1Scale * valuesN.cn(),
				o2Scale * valuesN.sn(),
				o3Scale * valuesN.dn());

		//We need to convert again omega because the formula used isn't converted with the transformations before
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

		System.out.println("PSI = "+ psi+"\n"+"THETA = "+theta+"\n"+"PHI LIN = "+phiLinear);
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
		Rotation inertToBody = alignedToBody.applyTo(mAlignedToInert.revert());//alignedToBody.applyTo(mAlignedToInert.revert());

		Rotation convertAxesReverse = new Rotation(axes[0], axes[1], Vector3D.PLUS_I, Vector3D.PLUS_J);

		Vector3D axe1cr = convertAxesReverse.applyTo(axes[0]);
		System.out.println("BBBBBB : "+axe1cr.getX()+" "+axe1cr.getY()+" "+axe1cr.getZ());
		Vector3D axe2cr = convertAxesReverse.applyTo(axes[1]);
		System.out.println("BBBBBB : "+axe2cr.getX()+" "+axe2cr.getY()+" "+axe2cr.getZ());
		Vector3D axe3cr = convertAxesReverse.applyTo(axes[2]);
		System.out.println("BBBBBB : "+axe3cr.getX()+" "+axe3cr.getY()+" "+axe3cr.getZ());

		Rotation bodyToOriginalFrame = convertAxes.applyInverseTo(inertToBody);//(inertToBody.applyInverseTo(convertAxes)).revert();

		double[] angles = bodyToOriginalFrame.getAngles(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM);
		double[][] data = {{omega.getX(), omega.getY(), omega.getZ()}, {angles[0], angles[1], angles[2]}, {bodyToOriginalFrame.getQ0(), bodyToOriginalFrame.getQ1(), bodyToOriginalFrame.getQ2(), bodyToOriginalFrame.getQ3()}, {axes[0].getX(), axes[0].getY(), axes[0].getZ()}, {axes[1].getX(), axes[1].getY(), axes[1].getZ()},{i1, i2, i3}};
		return data;
	}

	public double[] doComputeDerivatives(double t, double[] y) {

		final  double[] yDot = new double[getDimension()];
		double[] y0 = getInitialState().getPrimaryState(); //Initial state before conversion
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

		double t0 = getInitialState().getTime();

		final double[][] tfm = computeTorqueFreeMotion(i1C, i2C, i3C, t0, y0C, t);
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
