package org.hipparchus.ode;

import static org.junit.Assert.*;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import org.junit.Assert;
import org.junit.Test;
import org.orekit.torque_free.TorqueFreeMotion;
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.RotationConvention;
import org.hipparchus.geometry.euclidean.threed.RotationOrder;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.EmbeddedRungeKuttaIntegrator;
import org.hipparchus.ode.nonstiff.EmbeddedRungeKuttaIntegratorAbstractTest;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;


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
    	//Calls, Primary state (o1, o2, o3, q0, q1, q2, q3), Final time, Error scale
    	
        super(0.0, new double[] { FastMath.sqrt(80/3), 0.0, 4.0, 0.9999999, 0.0, 0.0, 0.0}, 10.0, new double[] { 1.0, 1.0, 1.0 });
        i1 = 3.0 / 8.0;
        i2 = 1.0 / 2.0;
        i3 = 5.0 / 8.0;
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
		TorqueFreeMotion TFM = new TorqueFreeMotion(i1, i2, i3, y0);
        final Vector3D omega = TFM.rotationVector(t).getRotationVector();
        final Rotation quaternion = TFM.quaternion(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM, t);
        ArrayList<Double> angles = new ArrayList<>();
        angles = TFM.angles(t);

        
        
		return new double[] {
				omega.getX(),
				omega.getY(),
				omega.getZ(),
				quaternion.getQ0(),
				quaternion.getQ1(),
				quaternion.getQ2(),
				quaternion.getQ3(),
				angles.get(0),
				angles.get(1),
				angles.get(2)
		};
    }
	
}//Fin du programme