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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.geometry.euclidean.threed.FieldRotation;
import org.hipparchus.geometry.euclidean.threed.FieldVector3D;
import org.hipparchus.geometry.euclidean.threed.RotationConvention;
import org.hipparchus.geometry.euclidean.threed.RotationOrder;
import org.hipparchus.linear.FieldDecompositionSolver;
import org.hipparchus.linear.FieldMatrix;
import org.hipparchus.linear.FieldQRDecomposer;
import org.hipparchus.linear.FieldVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.ode.nonstiff.DormandPrince853FieldIntegrator;
import org.hipparchus.special.elliptic.jacobi.FieldCopolarN;
import org.hipparchus.special.elliptic.jacobi.FieldJacobiElliptic;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

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
 * {@code TestFieldProblem7} only solves the rotation rate part.
 * </p>
 * @param <T> the type of the field elements
 */
public class TestFieldProblem8<T extends CalculusFieldElement<T>>
    extends TestFieldProblemAbstract<T> {

    /** Inertia tensor. */
    final FieldMatrix<T> inertiaTensor;

    /** Solver for inertia tensor. */
    final FieldDecompositionSolver<T> inertiaSolver;

    /** Inertia sorted to get a motion about axis 3. */
    final Inertia<T> sortedInertia;

    /** State scaling factor. */
    final T o1Scale;

    /** State scaling factor. */
    final T o2Scale;

    /** State scaling factor. */
    final T o3Scale;

    /** Jacobi elliptic function. */
    final FieldJacobiElliptic<T> jacobi;

    /** Time scaling factor. */
    public final T tScale;

    /** Time reference for rotation rate. */
    final T tRef;

    /** Offset rotation  between initial inertial frame and the frame with moment vector and Z axis aligned. */
    FieldRotation<T> inertToAligned;

    /** Rotation to switch to the converted axes frame. */
    final FieldRotation<T> sortedToBody;

    /** Period of rotation rate. */
    final T period;

    /** Slope of the linear part of the phi model. */
    final T phiSlope;

    /** DenseOutputModel of phi. */
    final FieldDenseOutputModel<T> phiQuadratureModel;

    /** Integral part of quadrature model over one period. */
    final T integOnePeriod;

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
    public TestFieldProblem8(final T t0, final T t1, final FieldVector3D<T> omega0, final FieldRotation<T> r0,
                             final T i1, final FieldVector3D<T> a1,
                             final T i2, final FieldVector3D<T> a2,
                             final T i3, final FieldVector3D<T> a3) {
        // Arguments in the super constructor :
        // Initial time, Primary state (o1, o2, o3, q0, q1, q2, q3), Final time, Error scale
        super(t0,
              toArray(omega0.getX(), omega0.getY(), omega0.getZ(),
                      r0.getQ0(), r0.getQ1(), r0.getQ2(), r0.getQ3()),
              t1,
              toArray(t0.getField(), 7, 1.0));

        // build inertia tensor
        final FieldVector3D<T> n1  = a1.normalize();
        final FieldVector3D<T> n2  = a2.normalize();
        final FieldVector3D<T> n3  = (FieldVector3D.dotProduct(FieldVector3D.crossProduct(a1, a2), a3).getReal() > 0 ?
                                     a3.normalize() : a3.normalize().negate());
        final FieldMatrix<T> q = MatrixUtils.createFieldMatrix(t0.getField(), 3, 3);
        q.setEntry(0, 0, n1.getX());
        q.setEntry(0, 1, n1.getY());
        q.setEntry(0, 2, n1.getZ());
        q.setEntry(1, 0, n2.getX());
        q.setEntry(1, 1, n2.getY());
        q.setEntry(1, 2, n2.getZ());
        q.setEntry(2, 0, n3.getX());
        q.setEntry(2, 1, n3.getY());
        q.setEntry(2, 2, n3.getZ());
        final FieldMatrix<T> d = MatrixUtils.createFieldDiagonalMatrix(toArray(i1, i2, i3));
        this.inertiaTensor = q.multiply(d.multiplyTransposed(q));
        this.inertiaSolver = new FieldQRDecomposer<>(t0.getField().getZero().newInstance(1.0e-10)).decompose(inertiaTensor);

        final FieldVector3D<T> m0 = new FieldVector3D<>(i1.multiply(FieldVector3D.dotProduct(omega0, n1)), n1,
                                                        i2.multiply(FieldVector3D.dotProduct(omega0, n2)), n2,
                                                        i3.multiply(FieldVector3D.dotProduct(omega0, n3)), n3);

        // sort axes in increasing moments of inertia order
        Inertia<T> inertia = new Inertia<>(new InertiaAxis<>(i1, n1), new InertiaAxis<>(i2, n2), new InertiaAxis<>(i3, n3));
        if (inertia.getInertiaAxis1().getI().subtract(inertia.getInertiaAxis2().getI()).getReal() > 0) {
            inertia = inertia.swap12();
        }
        if (inertia.getInertiaAxis2().getI().subtract(inertia.getInertiaAxis3().getI()).getReal() > 0) {
            inertia = inertia.swap23();
        }
        if (inertia.getInertiaAxis1().getI().subtract(inertia.getInertiaAxis2().getI()).getReal() > 0) {
            inertia = inertia.swap12();
        }

        // in order to simplify implementation, we want the motion to be about axis 3
        // which is either the minimum or the maximum inertia axis
        final T  o1                = FieldVector3D.dotProduct(omega0, n1);
        final T  o2                = FieldVector3D.dotProduct(omega0, n2);
        final T  o3                = FieldVector3D.dotProduct(omega0, n3);
        final T  o12               = o1.multiply(o1);
        final T  o22               = o2.multiply(o2);
        final T  o32               = o3.multiply(o3);
        final T  twoE              = i1.multiply(o12).add(i2.multiply(o22)).add(i3.multiply(o32));
        final T  m2                = i1.multiply(i1).multiply(o12).add(i2.multiply(i2).multiply(o22)).add(i3.multiply(i3).multiply(o32));
        final T  separatrixInertia = (twoE.isZero()) ? t0.getField().getZero() : m2.divide(twoE);
        final boolean clockwise;
        if (separatrixInertia.subtract(inertia.getInertiaAxis2().getI()).getReal() < 0) {
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

        final T i1C = inertia.getInertiaAxis1().getI();
        final T i2C = inertia.getInertiaAxis2().getI();
        final T i3C = inertia.getInertiaAxis3().getI();
        final T i32 = i3C.subtract(i2C);
        final T i31 = i3C.subtract(i1C);
        final T i21 = i2C.subtract(i1C);

         // convert initial conditions to Euler angles such the M is aligned with Z in sorted computation frame
        sortedToBody   = new FieldRotation<>(FieldVector3D.getPlusI(t0.getField()),
                                             FieldVector3D.getPlusJ(t0.getField()),
                                             inertia.getInertiaAxis1().getA(),
                                             inertia.getInertiaAxis2().getA());
        final FieldVector3D<T> omega0Sorted = sortedToBody.applyInverseTo(omega0);
        final FieldVector3D<T> m0Sorted     = sortedToBody.applyInverseTo(m0);
        final T   phi0         = t0.getField().getZero(); // this angle can be set arbitrarily, so 0 is a fair value (Eq. 37.13 - 37.14)
        final T   theta0       = FastMath.acos(m0Sorted.getZ().divide(m0Sorted.getNorm()));
        final T   psi0         = FastMath.atan2(m0Sorted.getX(), m0Sorted.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!

        // compute offset rotation between inertial frame aligned with momentum and regular inertial frame
        final FieldRotation<T> alignedToSorted0 = new FieldRotation<>(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                                                      phi0, theta0, psi0);
        inertToAligned = alignedToSorted0.applyInverseTo(sortedToBody.applyInverseTo(r0));

        // Ω is always o1Scale * cn((t-tref) * tScale), o2Scale * sn((t-tref) * tScale), o3Scale * dn((t-tref) * tScale)
        tScale  = FastMath.copySign(FastMath.sqrt(i32.multiply(m2.subtract(twoE.multiply(i1C))).divide((i1C.multiply(i2C).multiply(i3C)))),
                                    clockwise ? omega0Sorted.getZ().negate() : omega0Sorted.getZ());
        o1Scale = FastMath.sqrt(twoE.multiply(i3C).subtract(m2).divide(i1C.multiply(i31)));
        o2Scale = FastMath.sqrt(twoE.multiply(i3C).subtract(m2).divide(i2C.multiply(i32)));
        o3Scale = FastMath.copySign(FastMath.sqrt(m2.subtract(twoE.multiply(i1C)).divide(i3C.multiply(i31))),
                                    omega0Sorted.getZ());

        final T k2 = (twoE.isZero()) ?
                     t0.getField().getZero() :
                     i21.multiply(twoE.multiply(i3C).subtract(m2)).
                         divide(i32.multiply(m2.subtract(twoE.multiply(i1C))));
        jacobi = JacobiEllipticBuilder.build(k2);
        period = LegendreEllipticIntegral.bigK(k2).multiply(4).divide(tScale);

        if (o1Scale.isZero()) {
            // special case where twoE * i3C = m2, then o2Scale is also zero
            // motion is exactly along one axis
            tRef = t0;
        } else {
            if (FastMath.abs(omega0Sorted.getX()).subtract(FastMath.abs(omega0Sorted.getY())).getReal() >= 0) {
                if (omega0Sorted.getX().getReal() >= 0) {
                    // omega is roughly towards +I
                    tRef = t0.subtract(jacobi.arcsn(omega0Sorted.getY().divide(o2Scale)).divide(tScale));
                } else {
                    // omega is roughly towards -I
                    tRef = t0.add(jacobi.arcsn(omega0Sorted.getY().divide(o2Scale)).divide(tScale).subtract(period.multiply(0.5)));
                }
            } else {
                if (omega0Sorted.getY().getReal() >= 0) {
                    // omega is roughly towards +J
                    tRef = t0.subtract(jacobi.arccn(omega0Sorted.getX().divide(o1Scale)).divide(tScale));
                } else {
                    // omega is roughly towards -J
                    tRef = t0.add(jacobi.arccn(omega0Sorted.getX().divide(o1Scale)).divide(tScale));
                }
            }
        }

        phiSlope           = FastMath.sqrt(m2).divide(i3C);
        phiQuadratureModel = computePhiQuadratureModel(t0);
        integOnePeriod     = phiQuadratureModel.getInterpolatedState(phiQuadratureModel.getFinalTime()).getPrimaryState()[0];

    }

    @SafeVarargs
    private static <T extends CalculusFieldElement<T>> T[] toArray(final T...elements) {
        return elements;
    }

    private static <T extends CalculusFieldElement<T>> T[] toArray(final Field<T> field, final int n, double d) {
        T[] array = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; ++i) {
            array[i] = field.getZero().newInstance(d);
        }
        return array;
    }

    private FieldDenseOutputModel<T> computePhiQuadratureModel(final T t0) {

        final T i1C = sortedInertia.getInertiaAxis1().getI();
        final T i2C = sortedInertia.getInertiaAxis2().getI();
        final T i3C = sortedInertia.getInertiaAxis3().getI();

        final T i32 = i3C.subtract(i2C);
        final T i31 = i3C.subtract(i1C);
        final T i21 = i2C.subtract(i1C);

        // coefficients for φ model
        final T b = phiSlope.multiply(i32).multiply(i31);
        final T c = i1C.multiply(i32);
        final T d = i3C.multiply(i21);

        // integrate the quadrature phi term on one period
        final DormandPrince853FieldIntegrator<T> integ = new DormandPrince853FieldIntegrator<>(t0.getField(),
                                                                                               1.0e-6 * period.getReal(),
                                                                                               1.0e-2 * period.getReal(),
                                                                                               phiSlope.getReal() * period.getReal() * 1.0e-13,
                                                                                               1.0e-13);
        final FieldDenseOutputModel<T> model = new FieldDenseOutputModel<>();
        integ.addStepHandler(model);

        integ.integrate(new FieldExpandableODE<T>(new FieldOrdinaryDifferentialEquation<T>() {

            /** {@inheritDoc} */
            @Override
            public int getDimension() {
                return 1;
            }

            /** {@inheritDoc} */
           @Override
            public T[] computeDerivatives(final T t, final T[] y) {
                final T sn = jacobi.valuesN(t.subtract(tRef).multiply(tScale)).sn();
                return toArray(b.divide(c.add(d.multiply(sn).multiply(sn))));
            }

        }), new FieldODEState<>(t0, toArray(t0.getField().getZero())), t0.add(period));

        return model;

    }

    public T[] computeTheoreticalState(T t) {

        final T t0            = getInitialTime();

        // angular velocity
        final FieldCopolarN<T> valuesN     = jacobi.valuesN(t.subtract(tRef).multiply(tScale));
        final FieldVector3D<T> omegaSorted = new FieldVector3D<>(o1Scale.multiply(valuesN.cn()),
                                                                 o2Scale.multiply(valuesN.sn()),
                                                                 o3Scale.multiply(valuesN.dn()));
        final FieldVector3D<T> omegaBody   = sortedToBody.applyTo(omegaSorted);

        // first Euler angles are directly linked to angular velocity
        final T   psi         = FastMath.atan2(sortedInertia.getInertiaAxis1().getI().multiply(omegaSorted.getX()),
                                               sortedInertia.getInertiaAxis2().getI().multiply(omegaSorted.getY()));
        final T   theta       = FastMath.acos(omegaSorted.getZ().divide(phiSlope));
        final T   phiLinear   = phiSlope.multiply(t.subtract(t0));

        // third Euler angle results from a quadrature
        final int nbPeriods   = (int) FastMath.floor(t.subtract(t0).divide(period)).getReal();
        final T tStartInteg   = t0.add(period.multiply(nbPeriods));
        final T integPartial  = phiQuadratureModel.getInterpolatedState(t.subtract(tStartInteg)).getPrimaryState()[0];
        final T phiQuadrature = integOnePeriod.multiply(nbPeriods).add(integPartial);
        final T phi           = phiLinear.add(phiQuadrature);

        // rotation between computation frame (aligned with momentum) and sorted computation frame
        // (it is simply the angles equations provided by Landau & Lifchitz)
        final FieldRotation<T> alignedToSorted = new FieldRotation<>(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                                                                     phi, theta, psi);

        // combine with offset rotation to get back from regular inertial frame to body frame
        FieldRotation<T> inertToBody = sortedToBody.applyTo(alignedToSorted.applyTo(inertToAligned));

        return toArray(omegaBody.getX(), omegaBody.getY(), omegaBody.getZ(),
                       inertToBody.getQ0(), inertToBody.getQ1(), inertToBody.getQ2(), inertToBody.getQ3());

    }

    public T[] doComputeDerivatives(T t, T[] y) {

        final  T[] yDot = MathArrays.buildArray(t.getField(), getDimension());

        // compute the derivatives using Euler equations
        final T[]   omega    = Arrays.copyOfRange(y, 0, 3);
        final T[]   minusOiO = FieldVector3D.crossProduct(new FieldVector3D<>(omega),
                                                          new FieldVector3D<>(inertiaTensor.operate(omega))).negate().toArray();
        final FieldVector<T> omegaDot = inertiaSolver.solve(MatrixUtils.createFieldVector(minusOiO));
        yDot[0] = omegaDot.getEntry(0);
        yDot[1] = omegaDot.getEntry(1);
        yDot[2] = omegaDot.getEntry(2);

        // compute the derivatives using Qdot = 0.5 * Omega_inertialframe * Q
        yDot[3] = y[0].multiply(y[4]).negate().subtract(y[1].multiply(y[5])).subtract(y[2].multiply(y[6])).multiply(0.5);
        yDot[4] = y[0].multiply(y[3]).         add(     y[2].multiply(y[5])).subtract(y[1].multiply(y[6])).multiply(0.5);
        yDot[5] = y[1].multiply(y[3]).         subtract(y[2].multiply(y[4])).add(     y[0].multiply(y[6])).multiply(0.5);
        yDot[6] = y[2].multiply(y[3]).         add(     y[1].multiply(y[4])).subtract(y[0].multiply(y[5])).multiply(0.5);

        return yDot;

    }

    /** Container for inertia of a 3D object.
     * <p>
     * Instances of this class are immutable
     * </p>
    */
   public static class Inertia<T extends CalculusFieldElement<T>> {

       /** Inertia along first axis. */
       private final InertiaAxis<T> iA1;

       /** Inertia along second axis. */
       private final InertiaAxis<T> iA2;

       /** Inertia along third axis. */
       private final InertiaAxis<T> iA3;

       /** Simple constructor from principal axes.
        * @param iA1 inertia along first axis
        * @param iA2 inertia along second axis
        * @param iA3 inertia along third axis
        */
       public Inertia(final InertiaAxis<T> iA1, final InertiaAxis<T> iA2, final InertiaAxis<T> iA3) {
           this.iA1 = iA1;
           this.iA2 = iA2;
           this.iA3 = iA3;
       }

       /** Swap axes 1 and 2.
        * @return inertia with swapped axes
        */
       public Inertia<T> swap12() {
           return new Inertia<>(iA2, iA1, iA3.negate());
       }

       /** Swap axes 1 and 3.
        * @return inertia with swapped axes
        */
       public Inertia<T> swap13() {
           return new Inertia<>(iA3, iA2.negate(), iA1);
       }

       /** Swap axes 2 and 3.
        * @return inertia with swapped axes
        */
       public Inertia<T> swap23() {
           return new Inertia<>(iA1.negate(), iA3, iA2);
       }

       /** Get inertia along first axis.
        * @return inertia along first axis
        */
       public InertiaAxis<T> getInertiaAxis1() {
           return iA1;
       }

       /** Get inertia along second axis.
        * @return inertia along second axis
        */
       public InertiaAxis<T> getInertiaAxis2() {
           return iA2;
       }

       /** Get inertia along third axis.
        * @return inertia along third axis
        */
       public InertiaAxis<T> getInertiaAxis3() {
           return iA3;
       }

   }

   /** Container for moment of inertia and associated inertia axis.
     * <p>
     * Instances of this class are immutable
     * </p>
     */
    public static class InertiaAxis<T extends CalculusFieldElement<T>> {

        /** Moment of inertia. */
        private final T i;

        /** Inertia axis. */
        private final FieldVector3D<T> a;

        /** Simple constructor to pair a moment of inertia with its associated axis.
         * @param i moment of inertia
         * @param a inertia axis
         */
        public InertiaAxis(final T i, final FieldVector3D<T> a) {
            this.i = i;
            this.a = a;
        }

        /** Reverse the inertia axis.
         * @return new container with reversed axis
         */
        public InertiaAxis<T> negate() {
            return new InertiaAxis<>(i, a.negate());
        }

        /** Get the moment of inertia.
         * @return moment of inertia
         */
        public T getI() {
            return i;
        }

        /** Get the inertia axis.
         * @return inertia axis
         */
        public FieldVector3D<T> getA() {
            return a;
        }

    }

}
