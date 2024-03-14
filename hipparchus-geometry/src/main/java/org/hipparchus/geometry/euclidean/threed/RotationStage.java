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

package org.hipparchus.geometry.euclidean.threed;

import org.hipparchus.CalculusFieldElement;

/** Enumerate for one stage of {@link RotationOrder}.
 * @since 3.1
 */
enum RotationStage {

    /** Rotation around X axis. */
    X {

        /** {@inheritDoc} */
        @Override
        public Vector3D getAxis() {
            return Vector3D.PLUS_I;
        }

        /** {@inheritDoc} */
        @Override
        public double getComponent(final Vector3D v) {
            return v.getX();
        }

        /** {@inheritDoc} */
        @Override
        public <T extends CalculusFieldElement<T>> T getComponent(final FieldVector3D<T> v) {
            return v.getX();
        }

    },

    /** Rotation around Y axis. */
    Y {

        /** {@inheritDoc} */
        @Override
        public Vector3D getAxis() {
            return Vector3D.PLUS_J;
        }

        /** {@inheritDoc} */
        @Override
        public double getComponent(final Vector3D v) {
            return v.getY();
        }

        /** {@inheritDoc} */
        @Override
        public <T extends CalculusFieldElement<T>> T getComponent(final FieldVector3D<T> v) {
            return v.getY();
        }

    },

    /** Rotation around Z axis. */
    Z {

        /** {@inheritDoc} */
        @Override
        public Vector3D getAxis() {
            return Vector3D.PLUS_K;
        }

        /** {@inheritDoc} */
        @Override
        public double getComponent(final Vector3D v) {
            return v.getZ();
        }

        /** {@inheritDoc} */
        @Override
        public <T extends CalculusFieldElement<T>> T getComponent(final FieldVector3D<T> v) {
            return v.getZ();
        }

    };

    /** Get the rotation axis.
     * @return rotation axis
     */
    public abstract Vector3D getAxis();

    /** Get vector component along axis.
     * @param v vector from which component should be retrieved
     * @return vector component along axis
     */
    public abstract double getComponent(Vector3D v);

    /** Get vector component along axis.
     * @param <T> type of the field elements
     * @param v vector from which component should be retrieved
     * @return vector component along axis
     */
    public abstract <T extends CalculusFieldElement<T>> T getComponent(FieldVector3D<T> v);

}
