/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.ode;

import java.util.Arrays;

import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.RotationConvention;
import org.hipparchus.geometry.euclidean.threed.RotationOrder;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.DecompositionSolver;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.QRDecomposer;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.special.elliptic.jacobi.CopolarN;
import org.hipparchus.special.elliptic.jacobi.JacobiElliptic;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;

/**
 * This class is used in the junit tests for the ODE integrators.

 * <p>This specific problem correspond to torque-free motion of a solid body
 * with moments of inertia I₁, I₂, and I₃ with respect to body axes x, y, and z.
 * We use here the notations from Landau and Lifchitz Course of Theoretical Physics,
 * Mechanics vol 1.
 * </p>
 * <p>
 * The equations of torque-free motion are given in the solid body frame by
 * equation 36.5:
 * <pre>
 *    I₁ dΩ₁/dt + (I₃ - I₂) Ω₂ Ω₃ = 0
 *    I₂ dΩ₂/dt + (I₁ - I₃) Ω₃ Ω₁ = 0
 *    I₃ dΩ₃/dt + (I₂ - I₁) Ω₁ Ω₂ = 0
 * </pre>
 * <p>
 * This problem solves the full motion (rotation rate and rotation), whereas
 * {@code TestProblem7} only solves the rotation rate part.
 * </p>
 * @param <T> the type of the field elements
 */
public class TestProblem8 extends TestProblemAbstract {

    /** Inertia tensor. */
    final RealMatrix inertiaTensor;

    /** Solver for inertia tensor. */
    final DecompositionSolver inertiaSolver;

    /** Inertia sorted to get a motion about axis 3. */
    final Inertia sortedInertia;

    /** State scaling factor. */
    final double o1Scale;

    /** State scaling factor. */
    final double o2Scale;

    /** State scaling factor. */
    final double o3Scale;

    /** Jacobi elliptic function. */
    final JacobiElliptic jacobi;

    /** Time scaling factor. */
    public final double tScale;

    /** Time reference for rotation rate. */
    final double tRef;

    /** Offset rotation  between initial inertial frame and the frame with moment vector and Z axis aligned. */
    Rotation inertToAligned;

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

    /**
     * Simple constructor.
     * @param t0 initial time
     * @param t1 final time
     * @param omega0 initial rotation rate
     * @param r0 initial rotation
     * @param i1 inertia along first axis
     * @param a1 first principal inertia axis
     * @param i2 inertia along second axis
     * @param a2 second principal inertia axis
     * @param i3 inertia along third axis
     * @param a3 third principal inertia axis
     */
    public TestProblem8(final double t0, final double t1, final Vector3D omega0, final Rotation r0,
                        final double i1, final Vector3D a1,
                        final double i2, final Vector3D a2,
                        final double i3, final Vector3D a3) {
        // Arguments in the super constructor :
        // Initial time, Primary state (o1, o2, o3, q0, q1, q2, q3), Final time, Error scale
        super(t0,
              new double[] {
                  omega0.getX(), omega0.getY(), omega0.getZ(),
                  r0.getQ0(), r0.getQ1(), r0.getQ2(), r0.getQ3()
              },
              t1,
              new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });

        // build inertia tensor
        final Vector3D n1  = a1.normalize();
        final Vector3D n2  = a2.normalize();
        final Vector3D n3  = (Vector3D.dotProduct(Vector3D.crossProduct(a1, a2), a3) > 0 ?
                              a3.normalize() : a3.normalize().negate());
        final RealMatrix q = MatrixUtils.createRealMatrix(3, 3);
        q.setEntry(0, 0, n1.getX());
        q.setEntry(0, 1, n1.getY());
        q.setEntry(0, 2, n1.getZ());
        q.setEntry(1, 0, n2.getX());
        q.setEntry(1, 1, n2.getY());
        q.setEntry(1, 2, n2.getZ());
        q.setEntry(2, 0, n3.getX());
        q.setEntry(2, 1, n3.getY());
        q.setEntry(2, 2, n3.getZ());
        final RealMatrix d = MatrixUtils.createRealDiagonalMatrix(new double[] { i1, i2, i3});
        this.inertiaTensor = q.multiply(d.multiplyTransposed(q));
        this.inertiaSolver = new QRDecomposer(1.0e-10).decompose(inertiaTensor);

        final Vector3D m0 = new Vector3D(i1 * Vector3D.dotProduct(omega0, n1), n1,
                                         i2 * Vector3D.dotProduct(omega0, n2), n2,
                                         i3 * Vector3D.dotProduct(omega0, n3), n3);

        // sort axes in increasing moments of inertia order
        Inertia inertia = new Inertia(new InertiaAxis(i1, n1), new InertiaAxis(i2, n2), new InertiaAxis(i3, n3));
        if (inertia.getInertiaAxis1().getI() > inertia.getInertiaAxis2().getI()) {
            inertia = inertia.swap12();
        }
        if (inertia.getInertiaAxis2().getI() > inertia.getInertiaAxis3().getI()) {
            inertia = inertia.swap23();
        }
        if (inertia.getInertiaAxis1().getI() > inertia.getInertiaAxis2().getI()) {
            inertia = inertia.swap12();
        }

        // in order to simplify implementation, we want the motion to be about axis 3
        // which is either the minimum or the maximum inertia axis
        final double  o1                = Vector3D.dotProduct(omega0, n1);
        final double  o2                = Vector3D.dotProduct(omega0, n2);
        final double  o3                = Vector3D.dotProduct(omega0, n3);
        final double  o12               = o1 * o1;
        final double  o22               = o2 * o2;
        final double  o32               = o3 * o3;
        final double  twoE              = i1 * o12 + i2 * o22 + i3 * o32;
        final double  m2                = i1 * i1 * o12 + i2 * i2 * o22 + i3 * i3 * o32;
        final double  separatrixInertia = (twoE == 0) ? 0.0 : m2 / twoE;
        final boolean clockwise;
        if (separatrixInertia < inertia.getInertiaAxis2().getI()) {
            // motion is about minimum inertia axis
            // we swap axes to put them in decreasing moments order
            // motion will be clockwise about axis 3
            clockwise = true;
            inertia   = inertia.swap13();
        } else {
            // motion is about maximum inertia axis
            // we keep axes in increasing moments order
            // motion will be counter-clockwise about axis 3
            clockwise = false;
        }
        sortedInertia = inertia;

        final double i1C = inertia.getInertiaAxis1().getI();
        final double i2C = inertia.getInertiaAxis2().getI();
        final double i3C = inertia.getInertiaAxis3().getI();
        final double i32 = i3C - i2C;
        final double i31 = i3C - i1C;
        final double i21 = i2C - i1C;

        // convert initial conditions to Euler angles such the M is aligned with Z in aligned computation frame
        sortedToBody   = new Rotation(Vector3D.PLUS_I, Vector3D.PLUS_J, inertia.getInertiaAxis1().getA(), inertia.getInertiaAxis2().getA());
        final Vector3D omega0Sorted = sortedToBody.applyInverseTo(omega0);
        final Vector3D m0Sorted     = sortedToBody.applyInverseTo(m0);
        final double   phi0         = 0; // this angle can be set arbitrarily, so 0 is a fair value (Eq. 37.13 - 37.14)
        final double   theta0       = FastMath.acos(m0Sorted.getZ() / m0Sorted.getNorm());
        final double   psi0         = FastMath.atan2(m0Sorted.getX(), m0Sorted.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!

        // compute offset rotation between inertial frame aligned with momentum and regular inertial frame
        final Rotation alignedToSorted0 = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                                       phi0, theta0, psi0);
        inertToAligned = alignedToSorted0.applyInverseTo(sortedToBody.applyInverseTo(r0));

        // Ω is always o1Scale * cn((t-tref) * tScale), o2Scale * sn((t-tref) * tScale), o3Scale * dn((t-tref) * tScale)
        tScale  = FastMath.copySign(FastMath.sqrt(i32 * (m2 - twoE * i1C) / (i1C * i2C * i3C)),
                                    clockwise ? -omega0Sorted.getZ() : omega0Sorted.getZ());
        o1Scale = FastMath.sqrt((twoE * i3C - m2) / (i1C * i31));
        o2Scale = FastMath.sqrt((twoE * i3C - m2) / (i2C * i32));
        o3Scale = FastMath.copySign(FastMath.sqrt((m2 - twoE * i1C) / (i3C * i31)), omega0Sorted.getZ());

        final double k2 = (twoE == 0) ? 0.0 : i21 * (twoE * i3C - m2) / (i32 * (m2 - twoE * i1C));
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
                } else {
                    // omega is roughly towards -I
                    tRef = t0 + jacobi.arcsn(omega0Sorted.getY() / o2Scale) / tScale - 0.5 * period;
                }
            } else {
                if (omega0Sorted.getY() >= 0) {
                    // omega is roughly towards +J
                    tRef = t0 - jacobi.arccn(omega0Sorted.getX() / o1Scale) / tScale;
                } else {
                    // omega is roughly towards -J
                    tRef = t0 + jacobi.arccn(omega0Sorted.getX() / o1Scale) / tScale;
                }
            }
        }

        phiSlope           = FastMath.sqrt(m2) / i3C;
        phiQuadratureModel = computePhiQuadratureModel(t0);
        integOnePeriod     = phiQuadratureModel.getInterpolatedState(phiQuadratureModel.getFinalTime()).getPrimaryState()[0];

    }

    private DenseOutputModel computePhiQuadratureModel(final double t0) {

        final double i1C = sortedInertia.getInertiaAxis1().getI();
        final double i2C = sortedInertia.getInertiaAxis2().getI();
        final double i3C = sortedInertia.getInertiaAxis3().getI();

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

    public double[] computeTheoreticalState(double t) {

        final double t0            = getInitialTime();

        // angular velocity
        final CopolarN valuesN     = jacobi.valuesN((t - tRef) * tScale);
        final Vector3D omegaSorted = new Vector3D(o1Scale * valuesN.cn(), o2Scale * valuesN.sn(), o3Scale * valuesN.dn());
        final Vector3D omegaBody   = sortedToBody.applyTo(omegaSorted);

        // first Euler angles are directly linked to angular velocity
        final double   psi         = FastMath.atan2(sortedInertia.getInertiaAxis1().getI() * omegaSorted.getX(),
                                                    sortedInertia.getInertiaAxis2().getI() * omegaSorted.getY());
        final double   theta       = FastMath.acos(omegaSorted.getZ() / phiSlope);
        final double   phiLinear   = phiSlope * (t - t0);

        // third Euler angle results from a quadrature
        final int    nbPeriods     = (int) FastMath.floor((t - t0) / period);
        final double tStartInteg   = t0 + nbPeriods * period;
        final double integPartial  = phiQuadratureModel.getInterpolatedState(t - tStartInteg).getPrimaryState()[0];
        final double phiQuadrature = nbPeriods * integOnePeriod + integPartial;
        final double phi           = phiLinear + phiQuadrature;

        // rotation between computation frame (aligned with momentum) and sorted computation frame
        // (it is simply the angles equations provided by Landau & Lifchitz)
        final Rotation alignedToSorted = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                                      phi, theta, psi);

        // combine with offset rotation to get back from regular inertial frame to body frame
        Rotation inertToBody = sortedToBody.applyTo(alignedToSorted.applyTo(inertToAligned));

        return new double[] {
            omegaBody.getX(), omegaBody.getY(), omegaBody.getZ(),
            inertToBody.getQ0(), inertToBody.getQ1(), inertToBody.getQ2(), inertToBody.getQ3()
        };

    }

    public double[] doComputeDerivatives(double t, double[] y) {

        final  double[] yDot = new double[getDimension()];

        // compute the derivatives using Euler equations
        final double[]   omega    = Arrays.copyOfRange(y, 0, 3);
        final double[]   minusOiO = Vector3D.crossProduct(new Vector3D(omega),
                                                          new Vector3D(inertiaTensor.operate(omega))).negate().toArray();
        final RealVector omegaDot = inertiaSolver.solve(MatrixUtils.createRealVector(minusOiO));
        yDot[0] = omegaDot.getEntry(0);
        yDot[1] = omegaDot.getEntry(1);
        yDot[2] = omegaDot.getEntry(2);

        // compute the derivatives using Qdot = 0.5 * Omega_inertialframe * Q
        yDot[3] = 0.5 * (-y[0] * y[4] - y[1] * y[5] - y[2] * y[6]);
        yDot[4] = 0.5 * ( y[0] * y[3] + y[2] * y[5] - y[1] * y[6]);
        yDot[5] = 0.5 * ( y[1] * y[3] - y[2] * y[4] + y[0] * y[6]);
        yDot[6] = 0.5 * ( y[2] * y[3] + y[1] * y[4] - y[0] * y[5]);

        return yDot;

    }

    /** Container for inertia of a 3D object.
     * <p>
     * Instances of this class are immutable
     * </p>
    */
   public static class Inertia {

       /** Inertia along first axis. */
       private final InertiaAxis iA1;

       /** Inertia along second axis. */
       private final InertiaAxis iA2;

       /** Inertia along third axis. */
       private final InertiaAxis iA3;

       /** Simple constructor from principal axes.
        * @param iA1 inertia along first axis
        * @param iA2 inertia along second axis
        * @param iA3 inertia along third axis
        */
       public Inertia(final InertiaAxis iA1, final InertiaAxis iA2, final InertiaAxis iA3) {
           this.iA1 = iA1;
           this.iA2 = iA2;
           this.iA3 = iA3;
       }

       /** Swap axes 1 and 2.
        * @return inertia with swapped axes
        */
       public Inertia swap12() {
           return new Inertia(iA2, iA1, iA3.negate());
       }

       /** Swap axes 1 and 3.
        * @return inertia with swapped axes
        */
       public Inertia swap13() {
           return new Inertia(iA3, iA2.negate(), iA1);
       }

       /** Swap axes 2 and 3.
        * @return inertia with swapped axes
        */
       public Inertia swap23() {
           return new Inertia(iA1.negate(), iA3, iA2);
       }

       /** Get inertia along first axis.
        * @return inertia along first axis
        */
       public InertiaAxis getInertiaAxis1() {
           return iA1;
       }

       /** Get inertia along second axis.
        * @return inertia along second axis
        */
       public InertiaAxis getInertiaAxis2() {
           return iA2;
       }

       /** Get inertia along third axis.
        * @return inertia along third axis
        */
       public InertiaAxis getInertiaAxis3() {
           return iA3;
       }

   }

   /** Container for moment of inertia and associated inertia axis.
     * <p>
     * Instances of this class are immutable
     * </p>
     */
    public static class InertiaAxis {

        /** Moment of inertia. */
        private final double i;

        /** Inertia axis. */
        private final Vector3D a;

        /** Simple constructor to pair a moment of inertia with its associated axis.
         * @param i moment of inertia
         * @param a inertia axis
         */
        public InertiaAxis(final double i, final Vector3D a) {
            this.i = i;
            this.a = a;
        }

        /** Reverse the inertia axis.
         * @return new container with reversed axis
         */
        public InertiaAxis negate() {
            return new InertiaAxis(i, a.negate());
        }

        /** Get the moment of inertia.
         * @return moment of inertia
         */
        public double getI() {
            return i;
        }

        /** Get the inertia axis.
         * @return inertia axis
         */
        public Vector3D getA() {
            return a;
        }

    }

}
