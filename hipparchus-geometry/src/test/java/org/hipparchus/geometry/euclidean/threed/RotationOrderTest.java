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

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Field;

import static org.junit.jupiter.api.Assertions.assertEquals;


class RotationOrderTest {

    @Test
    void testName() {

    RotationOrder[] orders = {
      RotationOrder.XYZ, RotationOrder.XZY, RotationOrder.YXZ,
      RotationOrder.YZX, RotationOrder.ZXY, RotationOrder.ZYX,
      RotationOrder.XYX, RotationOrder.XZX, RotationOrder.YXY,
      RotationOrder.YZY, RotationOrder.ZXZ, RotationOrder.ZYZ
    };

    for (int i = 0; i < orders.length; ++i) {
      assertEquals(getFieldName(orders[i]), orders[i].toString());
    }

  }

    @Test
    void testIssue72() {
      for (RotationOrder order : RotationOrder.values()) {
          RotationOrder buildOrder = RotationOrder.getRotationOrder(order.toString());
          assertEquals(0.0, Vector3D.distance1(order.getA1(), buildOrder.getA1()), Double.MIN_VALUE);
          assertEquals(0.0, Vector3D.distance1(order.getA2(), buildOrder.getA2()), Double.MIN_VALUE);
          assertEquals(0.0, Vector3D.distance1(order.getA3(), buildOrder.getA3()), Double.MIN_VALUE);
      }
  }

    @Test
    void testIssue72InvalidName() {
      String wrongName = "BCE";
      try {
          RotationOrder.getRotationOrder(wrongName);
      } catch (MathIllegalStateException mise) {
          assertEquals(LocalizedGeometryFormats.INVALID_ROTATION_ORDER_NAME, mise.getSpecifier());
          assertEquals(wrongName, mise.getParts()[0]);
      }
  }


  private String getFieldName(RotationOrder order) {
    try {
      Field[] fields = RotationOrder.class.getFields();
      for (int i = 0; i < fields.length; ++i) {
        if (fields[i].get(null) == order) {
          return fields[i].getName();
        }
      }
    } catch (IllegalAccessException iae) {
      // ignored
    }
    return "unknown";
  }

}
