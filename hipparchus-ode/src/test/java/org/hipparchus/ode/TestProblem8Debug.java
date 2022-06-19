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

	/**Substraction of inertias. */
	final double i32;
	final double i31;
	final double i21;

	/**Initial state. */
	final double y0[];

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

	/** Converted initial state. */
	final double[] y0CDEUX;

	/** Converted axes. */
	final Vector3D[] axesDEUX;

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
//		super(0.0, new double[] {5.0, 0.0, 4.0, 0.9 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.437 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.0, 0.0,   4.0, 0.0, 5.0, -0.3088560588509295, 0.6360879930568342, -0.3088560588509294, 0.6360879930568341 }, 20.0, new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });
		super(0.0, new double[] {5.0, 0.0, 4.0, 0.9 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.437 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.0, 0.0,   5.0, 0.0, 4.0, 0.9 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.437 / FastMath.sqrt(0.9 * 0.9 + 0.437 * 0.437), 0.0, 0.0}, 20.0, new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });

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

		final double[][] converted = sortInertiaAxis();
		final double[] i = converted[1];
		final Vector3D axeX = new Vector3D(converted[2][0], converted[2][1], converted[2][2]);
		final Vector3D axeY = new Vector3D(converted[3][0], converted[3][1], converted[3][2]);

		i1C = i[0];
		i2C = i[1];
		i3C = i[2];

		axes = new Vector3D[] {axeX, axeY, axeX.crossProduct(axeY)};

		y0C = converted[0].clone();

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

		convertAxes = new Rotation( Vector3D.PLUS_I, Vector3D.PLUS_J, axes[0], axes[1] );

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

		if (y0C[1] == 0){
			tRef = t0;
			System.out.println("Tref0 : "+tRef);
		} else {
			tRef   = t0 - jacobi.arcsn(y0C[1] / o2Scale) / tScale;
			System.out.println("Tref : "+tRef);
		}
		System.out.println((y0C[1] / o2Scale));

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

		final double[][] convertedDEUX = sortInertiaAxisDEUX();
		final double[] iDEUX = convertedDEUX[1];
		final Vector3D axeXDEUX = new Vector3D(convertedDEUX[2][0], convertedDEUX[2][1], convertedDEUX[2][2]);
		final Vector3D axeYDEUX = new Vector3D(convertedDEUX[3][0], convertedDEUX[3][1], convertedDEUX[3][2]);

		i1CDEUX = iDEUX[0];
		i2CDEUX = iDEUX[1];
		i3CDEUX = iDEUX[2];

		axesDEUX = new Vector3D[] {axeXDEUX, axeYDEUX, axeXDEUX.crossProduct(axeYDEUX)};

		y0CDEUX = convertedDEUX[0].clone();

		// convert initial conditions to Euler angles such the M is aligned with Z in computation frame
		final Vector3D omega0BodyDEUX = new Vector3D(y0CDEUX[7], y0CDEUX[8], y0CDEUX[9]);
		final Vector3D m0BodyDEUX     = new Vector3D(i1CDEUX * omega0BodyDEUX.getX(), i2CDEUX * omega0BodyDEUX.getY(), i3CDEUX * omega0BodyDEUX.getZ());

		final double   phi0DEUX       = 0; // this angle can be set arbitrarily, so 0 is a fair value (Eq. 37.13 - 37.14)
		final double   theta0DEUX =  FastMath.acos(m0BodyDEUX.getZ() / m0BodyDEUX.getNorm());
		final double   psi0DEUX       = FastMath.atan2(m0BodyDEUX.getX(), m0BodyDEUX.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!

		//Compute offset rotation between inertial frame aligned with momentum and regular inertial frame
		final Rotation mAlignedToBodyDEUX = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
				phi0DEUX, theta0DEUX, psi0DEUX);

		convertAxesDEUX = new Rotation(convertedDEUX[5][0], convertedDEUX[5][1], convertedDEUX[5][2], convertedDEUX[5][3], false);
		
		
		final Rotation r0ConvertedAxisDEUX = new Rotation(convertedDEUX[4][0], convertedDEUX[4][1], convertedDEUX[4][2], convertedDEUX[4][3], false);
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

		if (y0CDEUX[8] == 0){
			tRefDEUX = t0;
		} else {
			tRefDEUX   = t0 - jacobiDEUX.arcsn(y0CDEUX[8] / o2ScaleDEUX) / tScaleDEUX;
		}

        periodDEUX             = 4 * LegendreEllipticIntegral.bigK(k2DEUX) / tScaleDEUX;
        phiSlopeDEUX           = FastMath.sqrt(m2DEUX) / i3DEUX;
        phiQuadratureModelDEUX = computePhiQuadratureModel(i1CDEUX, i2CDEUX, i3CDEUX,
                                                           k2DEUX, m2DEUX, jacobiDEUX, periodDEUX, phiSlopeDEUX,
                                                           t0, tRefDEUX, tScaleDEUX);
        integOnePeriodDEUX     = phiQuadratureModelDEUX.getInterpolatedState(phiQuadratureModelDEUX.getFinalTime()).getPrimaryState()[0];

	}//Fin du constructeur

	public double[][] sortInertiaAxis() {

		Vector3D[] axesP = {Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_K};

		double[] i = {i1, i2, i3};
		double[] y0C = y0.clone();

		System.out.println("omega avant : "+y0C[0]+" "+ y0C[1]+" "+ y0C[2]);
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

			double v = y0C[0];
			y0C[0] = y0C[1];
			y0C[1] = v;
			y0C[2] = -y0C[2];
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

			double v = y0C[1];
			y0C[1] = y0C[2];
			y0C[2] = v;
			y0C[0] = -y0C[0];
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

			double v = y0C[0];
			y0C[0] = y0C[1];
			y0C[1] = v;
			y0C[2] = -y0C[2];
			System.out.println("3ere boucle : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}

		final double condition = m2/twoE;

		System.out.println("Condition : "+condition);
		if(condition < i[1]) {
			Vector3D z = axesP[0];
			axesP[0] = axesP[2];
			axesP[2] = z;
			axesP[1] = axesP[1].negate();

			double y = i[0];
			i[0] = i[2];
			i[2] = y;

			double v = y0C[0];
			y0C[0] = y0C[2];
			y0C[2] = v;
			y0C[1] = - y0C[1];
			System.out.println("repère final  : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}
		System.out.println("omega après : "+y0C[0]+" "+ y0C[1]+" "+ y0C[2]);
		System.out.println("axes après : "+axesP[0]+" "+ axesP[1]+" "+ axesP[2]);
		System.out.println("inerties après : "+i[0]+" "+ i[1]+" "+ i[2]);

        final Rotation r0DEUX              = new Rotation(y0C[3], y0C[4], y0C[5], y0C[6], true);
        final Rotation convertAxesDEUX     = new Rotation( Vector3D.PLUS_I, Vector3D.PLUS_J, axesP[0], axesP[1] );
        final Rotation r0ConvertedAxisDEUX = convertAxesDEUX.applyTo(r0DEUX);
        y0C[3] = r0ConvertedAxisDEUX.getQ0();
        y0C[4] = r0ConvertedAxisDEUX.getQ1();
        y0C[5] = r0ConvertedAxisDEUX.getQ2();
        y0C[6] = r0ConvertedAxisDEUX.getQ3();

        return new double[][] { y0C, i, {axesP[0].getX(), axesP[0].getY(), axesP[0].getZ()}, {axesP[1].getX(), axesP[1].getY(), axesP[1].getZ()} };
	}

	public double[][] sortInertiaAxisDEUX() {

		Vector3D[] axesP = {Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_K};
		double[] i = {i1DEUX, i2DEUX, i3DEUX};

		double[] y0C = y0.clone();

		System.out.println("omega avant : "+y0C[7]+" "+ y0C[8]+" "+ y0C[9]);
		System.out.println("Initial : iA1 = "+i1DEUX+" "+axesP[0]);
		System.out.println("iA2 = "+i2DEUX + " "+axesP[1]);
		System.out.println("iA3 = "+i3DEUX + " "+axesP[2]);
		if (i[0] > i[1]) {
			Vector3D z = axesP[0];
			axesP[0] = axesP[1];
			axesP[1] = z;
			axesP[2] = axesP[2].negate();

			double y = i[0];
			i[0] = i[1];
			i[1] = y;

			double v = y0C[7];
			y0C[7] = y0C[8];
			y0C[8] = v;
			y0C[9] = -y0C[9];
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

			double v = y0C[8];
			y0C[8] = y0C[9];
			y0C[9] = v;
			y0C[7] = -y0C[7];
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

			double v = y0C[7];
			y0C[7] = y0C[8];
			y0C[8] = v;
			y0C[9] = -y0C[9];
			System.out.println("3ere boucle : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}

		final double condition = m2DEUX/twoEDEUX;

		System.out.println("Condition : "+condition);
		if(condition < i[1]) {
			Vector3D z = axesP[0];
			axesP[0] = axesP[2];
			axesP[2] = z;
			axesP[1] = axesP[1].negate();

			double y = i[0];
			i[0] = i[2];
			i[2] = y;

			double v = y0C[7];
			y0C[7] = y0C[9];
			y0C[9] = v;
			y0C[8] = - y0C[8];
			System.out.println("repère final  : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}
		System.out.println("omega après : "+y0C[7]+" "+ y0C[8]+" "+ y0C[9]);
		System.out.println("axes après : "+axesP[0]+" "+ axesP[1]+" "+ axesP[2]);
		System.out.println("inerties après : "+i[0]+" "+ i[1]+" "+ i[2]);

		
		final Rotation r0DEUX         = new Rotation(y0C[10], y0C[11], y0C[12], y0C[13], true);
		
		final Rotation convertAxesDEUX = new Rotation( Vector3D.PLUS_I, Vector3D.PLUS_J, axesP[0], axesP[1] );

		final Rotation r0ConvertedAxisDEUX = convertAxesDEUX.applyTo(r0DEUX);
        y0C[10] = r0ConvertedAxisDEUX.getQ0();
        y0C[11] = r0ConvertedAxisDEUX.getQ1();
        y0C[12] = r0ConvertedAxisDEUX.getQ2();
        y0C[13] = r0ConvertedAxisDEUX.getQ3();
		
		return new double[][] { y0C, i, {axesP[0].getX(), axesP[0].getY(), axesP[0].getZ()}, {axesP[1].getX(), axesP[1].getY(), axesP[1].getZ()}, {r0ConvertedAxisDEUX.getQ0(), r0ConvertedAxisDEUX.getQ1(), r0ConvertedAxisDEUX.getQ2(), r0ConvertedAxisDEUX.getQ3()}
		, {convertAxesDEUX.getQ0(), convertAxesDEUX.getQ1(), convertAxesDEUX.getQ2(), convertAxesDEUX.getQ3()}};
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

	public double[][][] computeTorqueFreeMotion(double i1, double i2, double i3, double t0, double[] y0, double t, double i1DEUX, double i2DEUX, double i3DEUX) {

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


		//Computation of omega
		final CopolarN valuesNDEUX = jacobiDEUX.valuesN((t - tRefDEUX) * tScaleDEUX);

		final Vector3D omegaPDEUX = new Vector3D(
				o1ScaleDEUX * valuesNDEUX.cn(),
				o2ScaleDEUX * valuesNDEUX.sn(),
				o3ScaleDEUX * valuesNDEUX.dn());

		//We need to convert again omega because the formula used isn't converted with the transformations before
		final Vector3D omegaDEUX = convertAxesDEUX.applyTo(omegaPDEUX);


		//Computation of the euler angles
		//Compute rotation rate
		final double   o1            = omega.getX();//o1Scale * valuesN.cn();
		final double   o2            = omega.getY();//o2Scale * valuesN.sn();
		final double   o3            = omega.getZ();//o3Scale * valuesN.dn();

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


		//Computation of the euler angles
		//Compute rotation rate
		final double   o1DEUX            = omegaDEUX.getX();//o1Scale * valuesN.cn();
		final double   o2DEUX            = omegaDEUX.getY();//o2Scale * valuesN.sn();
		final double   o3DEUX            = omegaDEUX.getZ();//o3Scale * valuesN.dn();

		//Compute angles
		final double   psiDEUX           = FastMath.atan2(i1DEUX * o1DEUX, i2DEUX * o2DEUX);
		final double   thetaDEUX         = FastMath.acos(o3DEUX / phiSlopeDEUX);
		final double   phiLinearDEUX     = phiSlopeDEUX * t;

		//Integration for the computation of phi
		final int nbPeriodsDEUX = (int) FastMath.floor((t - t0) / periodDEUX);//floor = entier inférieur = nb période entière
		final double tStartIntegDEUX = t0 + nbPeriodsDEUX * periodDEUX;//partie de période à la fin entre tau Integ et tau end
		final double integPartialDEUX = phiQuadratureModelDEUX.getInterpolatedState(t - tStartIntegDEUX).getPrimaryState()[0];// a vérifier, partie de l'intégrale apres le nb entier de période
		final double phiQuadratureDEUX = nbPeriodsDEUX * integOnePeriodDEUX + integPartialDEUX;

		final double phiDEUX = phiLinearDEUX + phiQuadratureDEUX;

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


		//Computation of the quaternion

		// Rotation between computation frame (aligned with momentum) and body
		//(It is simply the angles equations provided by L&L)
		final Rotation alignedToBodyDEUX = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
				phiDEUX, thetaDEUX, psiDEUX);

		// combine with offset rotation to get back to regular inertial frame
		//Inert -> aligned + aligned -> body = inert -> body (What the user wants)
		Rotation inertToBodyDEUX = alignedToBodyDEUX.applyTo(mAlignedToInertDEUX.revert());//alignedToBody.applyTo(mAlignedToInert.revert());

		Rotation convertAxesReverseDEUX = new Rotation(axesDEUX[0], axesDEUX[1], Vector3D.PLUS_I, Vector3D.PLUS_J);

		Vector3D axe1crDEUX = convertAxesReverseDEUX.applyTo(axesDEUX[0]);
		System.out.println("BBBBBBDEUX : "+axe1crDEUX.getX()+" "+axe1crDEUX.getY()+" "+axe1crDEUX.getZ());
		Vector3D axe2crDEUX = convertAxesReverseDEUX.applyTo(axesDEUX[1]);
		System.out.println("BBBBBBDEUX : "+axe2crDEUX.getX()+" "+axe2crDEUX.getY()+" "+axe2crDEUX.getZ());
		Vector3D axe3crDEUX = convertAxesReverseDEUX.applyTo(axesDEUX[2]);
		System.out.println("BBBBBBDEUX : "+axe3crDEUX.getX()+" "+axe3crDEUX.getY()+" "+axe3crDEUX.getZ());

		Rotation bodyToOriginalFrameDEUX = convertAxesDEUX.applyInverseTo(inertToBodyDEUX);//(inertToBody.applyInverseTo(convertAxes)).revert();



		double[] angles = bodyToOriginalFrame.getAngles(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM);


		double[] anglesDEUX = bodyToOriginalFrameDEUX.getAngles(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM);


		double[][] data = {{omega.getX(), omega.getY(), omega.getZ()}, {angles[0], angles[1], angles[2]}, {bodyToOriginalFrame.getQ0(), bodyToOriginalFrame.getQ1(), bodyToOriginalFrame.getQ2(), bodyToOriginalFrame.getQ3()}, {axes[0].getX(), axes[0].getY(), axes[0].getZ()}, {axes[1].getX(), axes[1].getY(), axes[1].getZ()},{i1, i2, i3}};
		double[][] dataDEUX = {{omegaDEUX.getX(), omegaDEUX.getY(), omegaDEUX.getZ()}, {anglesDEUX[0], anglesDEUX[1], anglesDEUX[2]}, {bodyToOriginalFrameDEUX.getQ0(), bodyToOriginalFrameDEUX.getQ1(), bodyToOriginalFrameDEUX.getQ2(), bodyToOriginalFrameDEUX.getQ3()}, {axesDEUX[0].getX(), axesDEUX[0].getY(), axesDEUX[0].getZ()}, {axesDEUX[1].getX(), axesDEUX[1].getY(), axesDEUX[1].getZ()},{i1DEUX, i2DEUX, i3DEUX}};
		double[][][] allData = {data, dataDEUX};
		return allData;
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

		double t0 = getInitialState().getTime();

		final double[][][] tfm = computeTorqueFreeMotion(i1C, i2C, i3C, t0, y0C, t, i1CDEUX, i2CDEUX, i3CDEUX);
		final double[] omega = tfm[0][0];
		final double[] quaternion = tfm[0][2];
		final double[] angles = tfm[0][1];
		final double[] X = tfm[0][3];
		final double[] Y = tfm[0][4];


		final double[] omegaDEUX = tfm[1][0];
		final double[] quaternionDEUX = tfm[1][2];
		final double[] anglesDEUX = tfm[1][1];
		final double[] XDEUX = tfm[1][3];
		final double[] YDEUX = tfm[1][4];

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
				Y[2],


				omegaDEUX[0],
				omegaDEUX[1],
				omegaDEUX[2],
				quaternionDEUX[0],
				quaternionDEUX[1],
				quaternionDEUX[2],
				quaternionDEUX[3],
				anglesDEUX[0],
				anglesDEUX[1],
				anglesDEUX[2],
				XDEUX[0],
				XDEUX[1],
				XDEUX[2],
				YDEUX[0],
				YDEUX[1],
				YDEUX[2]
		};
	}

}//Fin du programme
