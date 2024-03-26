/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.geometry.euclidean.threed;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/** Enumerate representing a rotation order specification for Cardan or Euler angles.
 * <p>
 * Since Hipparchus 1.7 this class is an enumerate class.
 * </p>
 */
public enum RotationOrder {

    /** Set of Cardan angles.
     * this ordered set of rotations is around X, then around Y, then
     * around Z
     */
     XYZ(RotationStage.X, RotationStage.Y, RotationStage.Z),

    /** Set of Cardan angles.
     * this ordered set of rotations is around X, then around Z, then
     * around Y
     */
     XZY(RotationStage.X, RotationStage.Z, RotationStage.Y),

    /** Set of Cardan angles.
     * this ordered set of rotations is around Y, then around X, then
     * around Z
     */
     YXZ(RotationStage.Y, RotationStage.X, RotationStage.Z),

    /** Set of Cardan angles.
     * this ordered set of rotations is around Y, then around Z, then
     * around X
     */
     YZX(RotationStage.Y, RotationStage.Z, RotationStage.X),

    /** Set of Cardan angles.
     * this ordered set of rotations is around Z, then around X, then
     * around Y
     */
     ZXY(RotationStage.Z, RotationStage.X, RotationStage.Y),

    /** Set of Cardan angles.
     * this ordered set of rotations is around Z, then around Y, then
     * around X
     */
     ZYX(RotationStage.Z, RotationStage.Y, RotationStage.X),

    /** Set of Euler angles.
     * this ordered set of rotations is around X, then around Y, then
     * around X
     */
     XYX(RotationStage.X, RotationStage.Y, RotationStage.X),

    /** Set of Euler angles.
     * this ordered set of rotations is around X, then around Z, then
     * around X
     */
     XZX(RotationStage.X, RotationStage.Z, RotationStage.X),

    /** Set of Euler angles.
     * this ordered set of rotations is around Y, then around X, then
     * around Y
     */
     YXY(RotationStage.Y, RotationStage.X, RotationStage.Y),

    /** Set of Euler angles.
     * this ordered set of rotations is around Y, then around Z, then
     * around Y
     */
     YZY(RotationStage.Y, RotationStage.Z, RotationStage.Y),

    /** Set of Euler angles.
     * this ordered set of rotations is around Z, then around X, then
     * around Z
     */
     ZXZ(RotationStage.Z, RotationStage.X, RotationStage.Z),

    /** Set of Euler angles.
     * this ordered set of rotations is around Z, then around Y, then
     * around Z
     */
     ZYZ(RotationStage.Z, RotationStage.Y, RotationStage.Z);

    /** Switch to safe computation of asin/acos.
     * @since 3.1
     */
    private static final double SAFE_SWITCH = 0.999;

    /** Name of the rotations order. */
    private final String name;

    /** First stage. */
    private final RotationStage stage1;

    /** Second stage. */
    private final RotationStage stage2;

    /** Third stage. */
    private final RotationStage stage3;

    /** Missing stage (for Euler rotations). */
    private final RotationStage missing;

    /** Sign for direct order (i.e. X before Y, Y before Z, Z before X). */
    private final double sign;

    /**
     * /** Private constructor. This is a utility class that cannot be
     * instantiated by the user, so its only constructor is private.
     *
     * @param stage1 first stage
     * @param stage2 second stage
     * @param stage3 third stage
     */
    RotationOrder(final RotationStage stage1, final RotationStage stage2, final RotationStage stage3) {
        this.name   = stage1.name() + stage2.name() + stage3.name();
        this.stage1 = stage1;
        this.stage2 = stage2;
        this.stage3 = stage3;

        if (stage1 == stage3) {
            // Euler rotations
            if (stage1 != RotationStage.X && stage2 != RotationStage.X) {
                missing = RotationStage.X;
            } else if (stage1 != RotationStage.Y && stage2 != RotationStage.Y) {
                missing = RotationStage.Y;
            } else {
                missing = RotationStage.Z;
            }
        } else {
            // Cardan rotations
            missing = null;
        }

        // check if the first two rotations are in direct or indirect order
        final Vector3D a1 = stage1.getAxis();
        final Vector3D a2 = stage2.getAxis();
        final Vector3D a3 = missing == null ? stage3.getAxis() : missing.getAxis();
        sign = FastMath.copySign(1.0, Vector3D.dotProduct(a3, Vector3D.crossProduct(a1, a2)));

    }

    /** Get a string representation of the instance.
     * @return a string representation of the instance (in fact, its name)
     */
    @Override
    public String toString() {
        return name;
    }

    /** Get the axis of the first rotation.
     * @return axis of the first rotation
     */
    public Vector3D getA1() {
        return stage1.getAxis();
    }

    /** Get the axis of the second rotation.
     * @return axis of the second rotation
     */
    public Vector3D getA2() {
        return stage2.getAxis();
    }

    /** Get the axis of the third rotation.
     * @return axis of the third rotation
     */
    public Vector3D getA3() {
        return stage3.getAxis();
    }

    /** Get the rotation order corresponding to a string representation.
     * @param value name
     * @return a rotation order object
     * @since 1.7
     */
    public static RotationOrder getRotationOrder(final String value) {
        try {
            return RotationOrder.valueOf(value);
        } catch (IllegalArgumentException iae) {
            // Invalid value. An exception is thrown
            throw new MathIllegalStateException(iae,
                                                LocalizedGeometryFormats.INVALID_ROTATION_ORDER_NAME,
                                                value);
        }
    }

    /** Get the Cardan or Euler angles corresponding to the instance.
     * @param rotation rotation from which angles should be extracted
     * @param convention convention to use for the semantics of the angle
     * @return an array of three angles, in the order specified by the set
     * @since 3.1
     */
    public double[] getAngles(final Rotation rotation, final RotationConvention convention) {

        // Cardan/Euler angles rotations matrices have the following form, taking
        // RotationOrder.XYZ (φ about X, θ about Y, ψ about Z) and RotationConvention.FRAME_TRANSFORM
        // as an example:
        // [  cos θ cos ψ     cos φ sin ψ + sin φ sin θ cos ψ     sin φ sin ψ - cos φ sin θ cos ψ ]
        // [ -cos θ sin ψ     cos φ cos ψ - sin φ sin θ sin ψ     cos ψ sin φ + cos φ sin θ sin ψ ]
        // [  sin θ                       - sin φ cos θ                         cos φ cos θ       ]

        // One can see that there is a "simple" column (the first column in the example above) and a "simple"
        // row (the last row in the example above). In fact, the simple column always corresponds to the axis
        // of the first rotation in RotationConvention.FRAME_TRANSFORM convention (i.e. X axis here, hence
        // first column) and to the axis of the third rotation in RotationConvention.VECTOR_OPERATOR
        // convention. The "simple" row always corresponds to axis of the third rotation in
        // RotationConvention.FRAME_TRANSFORM convention (i.e. Z axis here, hence last row) and to the axis
        // of the first rotation in RotationConvention.VECTOR_OPERATOR convention.
        final Vector3D vCol = getColumnVector(rotation, convention);
        final Vector3D vRow = getRowVector(rotation, convention);

        // in the example above, we see that the components of the simple row
        // are [ sin θ, -sin φ cos θ, cos φ cos θ ], which means that if we select
        // θ strictly between -π/2 and +π/2 (i.e. cos θ > 0), then we can extract
        // φ as atan2(-Yrow, +Zrow), i.e. the coordinates along the axes of the
        // two final rotations in the Cardan sequence.

        // in the example above, we see that the components of the simple column
        // are [ cos θ cos ψ, -cos θ sin ψ, sin θ ], which means that if we select
        // θ strictly between -π/2 and +π/2 (i.e. cos θ > 0), then we can extract
        // ψ as atan2(-Ycol, +Xcol), i.e. the coordinates along the axes of the
        // two initial rotations in the Cardan sequence.

        // if we are exactly on the singularity (i.e. θ = ±π/2), then the matrix
        // degenerates to:
        // [  0    sin(ψ ± φ)    ∓ cos(ψ ± φ) ]
        // [  0    cos(ψ ± φ)    ± sin(ψ ± φ) ]
        // [ ±1       0              0        ]
        // so we can set θ=±π/2 copying the sign from ±1 term,
        // arbitrarily set φ=0
        // and extract ψ as atan2(XmidCol, YmidCol)

        // These results generalize to all sequences, with some adjustments for
        // Euler sequence where we need to replace one axis by the missing axis
        // and some sign adjustments depending on the rotation order being direct
        // (i.e. X before Y, Y before Z, Z before X) or indirect

        if (missing == null) {
            // this is a Cardan angles order

            if (convention == RotationConvention.FRAME_TRANSFORM) {
                final double s = stage2.getComponent(vRow) * -sign;
                final double c = stage3.getComponent(vRow);
                if (s * s + c * c > 1.0e-30) {
                    // regular case
                    return new double[] {
                        FastMath.atan2(s, c),
                        safeAsin(stage1.getComponent(vRow), stage2.getComponent(vRow), stage3.getComponent(vRow)) * sign,
                        FastMath.atan2(stage2.getComponent(vCol) * -sign, stage1.getComponent(vCol))
                    };
                } else {
                    // near singularity
                    final Vector3D midCol = rotation.applyTo(stage2.getAxis());
                    return new double[] {
                        0.0,
                        FastMath.copySign(MathUtils.SEMI_PI, stage1.getComponent(vRow) * sign),
                        FastMath.atan2(stage1.getComponent(midCol) * sign, stage2.getComponent(midCol))
                    };
                }
            } else {
                final double s = stage2.getComponent(vCol) * -sign;
                final double c = stage3.getComponent(vCol);
                if (s * s + c * c > 1.0e-30) {
                    // regular case
                    return new double[] {
                        FastMath.atan2(s, c),
                        safeAsin(stage3.getComponent(vRow), stage1.getComponent(vRow), stage2.getComponent(vRow)) * sign,
                        FastMath.atan2(stage2.getComponent(vRow) * -sign, stage1.getComponent(vRow))
                    };
                } else {
                    // near singularity
                    final Vector3D midRow = rotation.applyInverseTo(stage2.getAxis());
                    return new double[] {
                        0.0,
                        FastMath.copySign(MathUtils.SEMI_PI, stage3.getComponent(vRow) * sign),
                        FastMath.atan2(stage1.getComponent(midRow) * sign, stage2.getComponent(midRow))
                    };
                }
            }
        } else {
            // this is an Euler angles order
            if (convention == RotationConvention.FRAME_TRANSFORM) {
                final double s = stage2.getComponent(vRow);
                final double c = missing.getComponent(vRow) * -sign;
                if (s * s + c * c > 1.0e-30) {
                    // regular case
                    return new double[] {
                        FastMath.atan2(s, c),
                        safeAcos(stage1.getComponent(vRow), stage2.getComponent(vRow), missing.getComponent(vRow)),
                        FastMath.atan2(stage2.getComponent(vCol), missing.getComponent(vCol) * sign)
                    };
                } else {
                    // near singularity
                    final Vector3D midCol = rotation.applyTo(stage2.getAxis());
                    return new double[] {
                        FastMath.atan2(missing.getComponent(midCol) * -sign, stage2.getComponent(midCol)) *
                        FastMath.copySign(1, stage1.getComponent(vRow)),
                        stage1.getComponent(vRow) > 0 ? 0 : FastMath.PI,
                        0.0
                    };
                }
            } else {
                final double s = stage2.getComponent(vCol);
                final double c = missing.getComponent(vCol) * -sign;
                if (s * s + c * c > 1.0e-30) {
                    // regular case
                    return new double[] {
                        FastMath.atan2(s, c),
                        safeAcos(stage1.getComponent(vRow), stage2.getComponent(vRow), missing.getComponent(vRow)),
                        FastMath.atan2(stage2.getComponent(vRow), missing.getComponent(vRow) * sign)
                    };
                } else {
                    // near singularity
                    final Vector3D midRow = rotation.applyInverseTo(stage2.getAxis());
                    return new double[] {
                        FastMath.atan2(missing.getComponent(midRow) * -sign, stage2.getComponent(midRow)) *
                        FastMath.copySign(1, stage1.getComponent(vRow)),
                        stage1.getComponent(vRow) > 0 ? 0 : FastMath.PI,
                        0.0
                    };
                }
            }
        }
    }

    /** Get the Cardan or Euler angles corresponding to the instance.
     * @param <T> type of the field elements
     * @param rotation rotation from which angles should be extracted
     * @param convention convention to use for the semantics of the angle
     * @return an array of three angles, in the order specified by the set
     * @since 3.1
     */
    public <T extends CalculusFieldElement<T>> T[] getAngles(final FieldRotation<T> rotation, final RotationConvention convention) {

        // Cardan/Euler angles rotations matrices have the following form, taking
        // RotationOrder.XYZ (φ about X, θ about Y, ψ about Z) and RotationConvention.FRAME_TRANSFORM
        // as an example:
        // [  cos θ cos ψ     cos φ sin ψ + sin φ sin θ cos ψ     sin φ sin ψ - cos φ sin θ cos ψ ]
        // [ -cos θ sin ψ     cos φ cos ψ - sin φ sin θ sin ψ     cos ψ sin φ + cos φ sin θ sin ψ ]
        // [  sin θ                       - sin φ cos θ                         cos φ cos θ       ]

        // One can see that there is a "simple" column (the first column in the example above) and a "simple"
        // row (the last row in the example above). In fact, the simple column always corresponds to the axis
        // of the first rotation in RotationConvention.FRAME_TRANSFORM convention (i.e. X axis here, hence
        // first column) and to the axis of the third rotation in RotationConvention.VECTOR_OPERATOR
        // convention. The "simple" row always corresponds to axis of the third rotation in
        // RotationConvention.FRAME_TRANSFORM convention (i.e. Z axis here, hence last row) and to the axis
        // of the first rotation in RotationConvention.VECTOR_OPERATOR convention.
        final FieldVector3D<T> vCol = getColumnVector(rotation, convention);
        final FieldVector3D<T> vRow = getRowVector(rotation, convention);

        // in the example above, we see that the components of the simple row
        // are [ sin θ, -sin φ cos θ, cos φ cos θ ], which means that if we select
        // θ strictly between -π/2 and +π/2 (i.e. cos θ > 0), then we can extract
        // φ as atan2(-Yrow, +Zrow), i.e. the coordinates along the axes of the
        // two final rotations in the Cardan sequence.

        // in the example above, we see that the components of the simple column
        // are [ cos θ cos ψ, -cos θ sin ψ, sin θ ], which means that if we select
        // θ strictly between -π/2 and +π/2 (i.e. cos θ > 0), then we can extract
        // ψ as atan2(-Ycol, +Xcol), i.e. the coordinates along the axes of the
        // two initial rotations in the Cardan sequence.

        // if we are exactly on the singularity (i.e. θ = ±π/2), then the matrix
        // degenerates to:
        // [  0    sin(ψ ± φ)    ∓ cos(ψ ± φ) ]
        // [  0    cos(ψ ± φ)    ± sin(ψ ± φ) ]
        // [ ±1       0              0        ]
        // so we can set θ=±π/2 copying the sign from ±1 term,
        // arbitrarily set φ=0
        // and extract ψ as atan2(XmidCol, YmidCol)

        // These results generalize to all sequences, with some adjustments for
        // Euler sequence where we need to replace one axis by the missing axis
        // and some sign adjustments depending on the rotation order being direct
        // (i.e. X before Y, Y before Z, Z before X) or indirect

        if (missing == null) {
            // this is a Cardan angles order
            if (convention == RotationConvention.FRAME_TRANSFORM) {
                final T s = stage2.getComponent(vRow).multiply(-sign);
                final T c = stage3.getComponent(vRow);
                if (s.square().add(c.square()).getReal() > 1.0e-30) {
                    // regular case
                    return buildArray(FastMath.atan2(s, c),
                                      safeAsin(stage1.getComponent(vRow), stage2.getComponent(vRow), stage3.getComponent(vRow)).multiply(sign),
                                      FastMath.atan2(stage2.getComponent(vCol).multiply(-sign), stage1.getComponent(vCol)));
                } else {
                    // near singularity
                    final FieldVector3D<T> midCol = rotation.applyTo(stage2.getAxis());
                    return buildArray(s.getField().getZero(),
                                      FastMath.copySign(s.getPi().multiply(0.5), stage1.getComponent(vRow).multiply(sign)),
                                      FastMath.atan2(stage1.getComponent(midCol).multiply(sign), stage2.getComponent(midCol)));
                }
            } else {
                final T s = stage2.getComponent(vCol).multiply(-sign);
                final T c = stage3.getComponent(vCol);
                if (s.square().add(c.square()).getReal() > 1.0e-30) {
                    // regular case
                    return buildArray(FastMath.atan2(stage2.getComponent(vCol).multiply(-sign), stage3.getComponent(vCol)),
                                      safeAsin(stage3.getComponent(vRow), stage1.getComponent(vRow), stage2.getComponent(vRow)).multiply(sign),
                                      FastMath.atan2(stage2.getComponent(vRow).multiply(-sign), stage1.getComponent(vRow)));
                } else {
                    // near singularity
                    final FieldVector3D<T> midRow = rotation.applyInverseTo(stage2.getAxis());
                    return buildArray(s.getField().getZero(),
                                      FastMath.copySign(s.getPi().multiply(0.5), stage3.getComponent(vRow).multiply(sign)),
                                      FastMath.atan2(stage1.getComponent(midRow).multiply(sign), stage2.getComponent(midRow)));
                }
            }
        } else {
            // this is an Euler angles order
            if (convention == RotationConvention.FRAME_TRANSFORM) {
                final T s = stage2.getComponent(vRow);
                final T c = missing.getComponent(vRow).multiply(-sign);
                if (s.square().add(c.square()).getReal() > 1.0e-30) {
                    // regular case
                    return buildArray(FastMath.atan2(s, c),
                                      safeAcos(stage1.getComponent(vRow), stage2.getComponent(vRow), missing.getComponent(vRow)),
                                      FastMath.atan2(stage2.getComponent(vCol), missing.getComponent(vCol).multiply(sign)));
                } else {
                    // near singularity
                    final FieldVector3D<T> midCol = rotation.applyTo(stage2.getAxis());
                    return buildArray(FastMath.atan2(missing.getComponent(midCol).multiply(-sign), stage2.getComponent(midCol)).
                                      multiply(FastMath.copySign(1, stage1.getComponent(vRow).getReal())),
                                      stage1.getComponent(vRow).getReal() > 0 ? s.getField().getZero() : s.getPi(),
                                      s.getField().getZero());
                }
            } else {
                final T s = stage2.getComponent(vCol);
                final T c = missing.getComponent(vCol).multiply(-sign);
                if (s.square().add(c.square()).getReal() > 1.0e-30) {
                    // regular case
                    return buildArray(FastMath.atan2(s, c),
                                      safeAcos(stage1.getComponent(vRow), stage2.getComponent(vRow), missing.getComponent(vRow)),
                                      FastMath.atan2(stage2.getComponent(vRow), missing.getComponent(vRow).multiply(sign)));
                } else {
                    // near singularity
                    final FieldVector3D<T> midRow = rotation.applyInverseTo(stage2.getAxis());
                    return buildArray(FastMath.atan2(missing.getComponent(midRow).multiply(-sign), stage2.getComponent(midRow)).
                                      multiply(FastMath.copySign(1, stage1.getComponent(vRow).getReal())),
                                      stage1.getComponent(vRow).getReal() > 0 ? s.getField().getZero() : s.getPi(),
                                      s.getField().getZero());
                }
            }
        }
    }

    /** Get the simplest column vector from the rotation matrix.
     * @param rotation rotation
     * @param convention convention to use for the semantics of the angle
     * @return column vector
     * @since 3.1
     */
    private Vector3D getColumnVector(final Rotation rotation,
                                     final RotationConvention convention) {
        return rotation.applyTo(convention == RotationConvention.VECTOR_OPERATOR ?
                                stage3.getAxis() :
                                stage1.getAxis());
    }

    /** Get the simplest row vector from the rotation matrix.
     * @param rotation rotation
     * @param convention convention to use for the semantics of the angle
     * @return row vector
     * @since 3.1
     */
    private Vector3D getRowVector(final Rotation rotation,
                                  final RotationConvention convention) {
        return rotation.applyInverseTo(convention == RotationConvention.VECTOR_OPERATOR ?
                                       stage1.getAxis() :
                                       stage3.getAxis());
    }

    /** Get the simplest column vector from the rotation matrix.
     * @param <T> type of the field elements
     * @param rotation rotation
     * @param convention convention to use for the semantics of the angle
     * @return column vector
     * @since 3.1
     */
    private <T extends CalculusFieldElement<T>> FieldVector3D<T> getColumnVector(final FieldRotation<T> rotation,
                                                                                 final RotationConvention convention) {
        return rotation.applyTo(convention == RotationConvention.VECTOR_OPERATOR ?
                                stage3.getAxis() :
                                stage1.getAxis());
    }

    /** Get the simplest row vector from the rotation matrix.
     * @param <T> type of the field elements
     * @param rotation rotation
     * @param convention convention to use for the semantics of the angle
     * @return row vector
     * @since 3.1
     */
    private <T extends CalculusFieldElement<T>> FieldVector3D<T> getRowVector(final FieldRotation<T> rotation,
                                                                              final RotationConvention convention) {
        return rotation.applyInverseTo(convention == RotationConvention.VECTOR_OPERATOR ?
                                       stage1.getAxis() :
                                       stage3.getAxis());
    }

    /** Safe computation of acos(some vector coordinate) working around singularities.
     * @param cos  cosine coordinate
     * @param sin1 one of the sine coordinates
     * @param sin2 other sine coordinate
     * @return acos of the coordinate
     * @since 3.1
     */
    private static double safeAcos(final double cos, final double sin1, final double sin2) {
        if (cos < -SAFE_SWITCH) {
            return FastMath.PI - FastMath.asin(FastMath.sqrt(sin1 * sin1 + sin2 * sin2));
        } else if (cos > SAFE_SWITCH) {
            return FastMath.asin(FastMath.sqrt(sin1 * sin1 + sin2 * sin2));
        } else {
            return FastMath.acos(cos);
        }
    }

    /** Safe computation of asin(some vector coordinate) working around singularities.
     * @param sin sine coordinate
     * @param cos1 one of the cosine coordinates
     * @param cos2 other cosine coordinate
     * @return asin of the coordinate
     * @since 3.1
     */
    private static double safeAsin(final double sin, final double cos1, final double cos2) {
        if (sin < -SAFE_SWITCH) {
            return -FastMath.acos(FastMath.sqrt(cos1 * cos1 + cos2 * cos2));
        } else if (sin > SAFE_SWITCH) {
            return FastMath.acos(FastMath.sqrt(cos1 * cos1 + cos2 * cos2));
        } else {
            return FastMath.asin(sin);
        }
    }

    /** Safe computation of acos(some vector coordinate) working around singularities.
     * @param <T> type of the field elements
     * @param cos  cosine coordinate
     * @param sin1 one of the sine coordinates
     * @param sin2 other sine coordinate
     * @return acos of the coordinate
     * @since 3.1
     */
    private static <T extends CalculusFieldElement<T>> T safeAcos(final T cos,
                                                                  final T sin1,
                                                                  final T sin2) {
        if (cos.getReal() < -SAFE_SWITCH) {
            return FastMath.asin(FastMath.sqrt(sin1.square().add(sin2.square()))).subtract(sin1.getPi()).negate();
        } else if (cos.getReal() > SAFE_SWITCH) {
            return FastMath.asin(FastMath.sqrt(sin1.square().add(sin2.square())));
        } else {
            return FastMath.acos(cos);
        }
    }

    /** Safe computation of asin(some vector coordinate) working around singularities.
     * @param <T> type of the field elements
     * @param sin sine coordinate
     * @param cos1 one of the cosine coordinates
     * @param cos2 other cosine coordinate
     * @return asin of the coordinate
     * @since 3.1
     */
    private static <T extends CalculusFieldElement<T>> T safeAsin(final T sin,
                                                                  final T cos1,
                                                                  final T cos2) {
        if (sin.getReal() < -SAFE_SWITCH) {
            return FastMath.acos(FastMath.sqrt(cos1.square().add(cos2.square()))).negate();
        } else if (sin.getReal() > SAFE_SWITCH) {
            return FastMath.acos(FastMath.sqrt(cos1.square().add(cos2.square())));
        } else {
            return FastMath.asin(sin);
        }
    }

    /** Create a dimension 3 array.
     * @param <T> type of the field elements
     * @param a0 first array element
     * @param a1 second array element
     * @param a2 third array element
     * @return new array
     * @since 3.1
     */
    private static <T extends CalculusFieldElement<T>> T[] buildArray(final T a0, final T a1, final T a2) {
        final T[] array = MathArrays.buildArray(a0.getField(), 3);
        array[0] = a0;
        array[1] = a1;
        array[2] = a2;
        return array;
    }

}
