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

import java.util.HashMap;
import java.util.Map;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.LocalizedGeometryFormats;


/**
 * This class is a utility representing a rotation order specification
 * for Cardan or Euler angles specification.
 *
 * This class cannot be instanciated by the user. He can only use one
 * of the twelve predefined supported orders as an argument to either
 * the {@link Rotation#Rotation(RotationOrder, RotationConvention, double, double, double)}
 * constructor or the {@link Rotation#getAngles} method.
 *
 * Since Hipparchus 1.7 this class is an enumerate class.
 *
 */
public enum RotationOrder {

    /** Set of Cardan angles.
     * this ordered set of rotations is around X, then around Y, then
     * around Z
     */
     XYZ("XYZ", Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_K),

    /** Set of Cardan angles.
     * this ordered set of rotations is around X, then around Z, then
     * around Y
     */
     XZY("XZY", Vector3D.PLUS_I, Vector3D.PLUS_K, Vector3D.PLUS_J),

    /** Set of Cardan angles.
     * this ordered set of rotations is around Y, then around X, then
     * around Z
     */
     YXZ("YXZ", Vector3D.PLUS_J, Vector3D.PLUS_I, Vector3D.PLUS_K),

    /** Set of Cardan angles.
     * this ordered set of rotations is around Y, then around Z, then
     * around X
     */
     YZX("YZX", Vector3D.PLUS_J, Vector3D.PLUS_K, Vector3D.PLUS_I),

    /** Set of Cardan angles.
     * this ordered set of rotations is around Z, then around X, then
     * around Y
     */
     ZXY("ZXY", Vector3D.PLUS_K, Vector3D.PLUS_I, Vector3D.PLUS_J),

    /** Set of Cardan angles.
     * this ordered set of rotations is around Z, then around Y, then
     * around X
     */
     ZYX("ZYX", Vector3D.PLUS_K, Vector3D.PLUS_J, Vector3D.PLUS_I),

    /** Set of Euler angles.
     * this ordered set of rotations is around X, then around Y, then
     * around X
     */
     XYX("XYX", Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_I),

    /** Set of Euler angles.
     * this ordered set of rotations is around X, then around Z, then
     * around X
     */
     XZX("XZX", Vector3D.PLUS_I, Vector3D.PLUS_K, Vector3D.PLUS_I),

    /** Set of Euler angles.
     * this ordered set of rotations is around Y, then around X, then
     * around Y
     */
     YXY("YXY", Vector3D.PLUS_J, Vector3D.PLUS_I, Vector3D.PLUS_J),

    /** Set of Euler angles.
     * this ordered set of rotations is around Y, then around Z, then
     * around Y
     */
     YZY("YZY", Vector3D.PLUS_J, Vector3D.PLUS_K, Vector3D.PLUS_J),

    /** Set of Euler angles.
     * this ordered set of rotations is around Z, then around X, then
     * around Z
     */
     ZXZ("ZXZ", Vector3D.PLUS_K, Vector3D.PLUS_I, Vector3D.PLUS_K),

    /** Set of Euler angles.
     * this ordered set of rotations is around Z, then around Y, then
     * around Z
     */
     ZYZ("ZYZ", Vector3D.PLUS_K, Vector3D.PLUS_J, Vector3D.PLUS_K);

    /** Codes map. */
    private static final Map<String, RotationOrder> CODES_MAP = new HashMap<>();
    static {
        for (final RotationOrder type : values()) {
            CODES_MAP.put(type.toString(), type);
        }
    }

    /** Name of the rotations order. */
    private final String name;

    /** Axis of the first rotation. */
    private final Vector3D a1;

    /** Axis of the second rotation. */
    private final Vector3D a2;

    /** Axis of the third rotation. */
    private final Vector3D a3;

    /** Private constructor.
     * This is a utility class that cannot be instantiated by the user,
     * so its only constructor is private.
     * @param name name of the rotation order
     * @param a1 axis of the first rotation
     * @param a2 axis of the second rotation
     * @param a3 axis of the third rotation
     */
    RotationOrder(final String name, final Vector3D a1, final Vector3D a2, final Vector3D a3) {
        this.name = name;
        this.a1   = a1;
        this.a2   = a2;
        this.a3   = a3;
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
        return a1;
    }

    /** Get the axis of the second rotation.
     * @return axis of the second rotation
     */
    public Vector3D getA2() {
        return a2;
    }

    /** Get the axis of the second rotation.
     * @return axis of the second rotation
     */
    public Vector3D getA3() {
        return a3;
    }

    /**
     * Get the rotation order corresponding to a string representation.
     * @param value name
     * @return a rotation order object
     * @since 1.7
     */
    public static RotationOrder getRotationOrder(final String value) {
        final RotationOrder type = CODES_MAP.get(value);
        if (type == null) {
            // Invalid value. An exception is thrown
            throw new MathIllegalStateException(LocalizedGeometryFormats.INVALID_ROTATION_ORDER_NAME, value);
        }
        return type;
    }

}
