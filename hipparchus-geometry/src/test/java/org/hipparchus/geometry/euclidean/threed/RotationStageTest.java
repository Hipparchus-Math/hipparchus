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
import org.hipparchus.Field;
import org.hipparchus.util.Binary64Field;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class RotationStageTest {

    @Test
    public void testAxis() {
        Assertions.assertEquals(0.0,
                            Vector3D.distance(Vector3D.PLUS_I, RotationStage.X.getAxis()),
                            1.0e-15);
        Assertions.assertEquals(0.0,
                            Vector3D.distance(Vector3D.PLUS_J, RotationStage.Y.getAxis()),
                            1.0e-15);
        Assertions.assertEquals(0.0,
                            Vector3D.distance(Vector3D.PLUS_K, RotationStage.Z.getAxis()),
                            1.0e-15);
    }

    @Test
    public void testComponent() {
        final Vector3D v = new Vector3D(1.0, 2.0, 3.0);
        Assertions.assertEquals( 1.0, RotationStage.X.getComponent(v), 1.0e-15);
        Assertions.assertEquals( 2.0, RotationStage.Y.getComponent(v), 1.0e-15);
        Assertions.assertEquals( 3.0, RotationStage.Z.getComponent(v), 1.0e-15);
    }

    @Test
    public void testFieldComponent() {
        doTestFieldComponent(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestFieldComponent(final Field<T> field) {
        final FieldVector3D<T> v = new FieldVector3D<>(field.getZero().newInstance(1.0),
                                                       field.getZero().newInstance(2.0),
                                                       field.getZero().newInstance(3.0));
        Assertions.assertEquals( 1.0, RotationStage.X.getComponent(v).getReal(), 1.0e-15);
        Assertions.assertEquals( 2.0, RotationStage.Y.getComponent(v).getReal(), 1.0e-15);
        Assertions.assertEquals( 3.0, RotationStage.Z.getComponent(v).getReal(), 1.0e-15);
    }

}
