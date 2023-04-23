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

package org.hipparchus.migration.geometry.euclidean;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.migration.LocalizedMigrationFormats;

/** This class represents exceptions thrown while extractiong Cardan
 * or Euler angles from a rotation.

 * @deprecated as of 1.0, this exception is replaced by {@link MathIllegalStateException}
 */
@Deprecated
public class CardanEulerSingularityException
  extends MathIllegalStateException {

    /** Serializable version identifier */
    private static final long serialVersionUID = -1360952845582206770L;

    /**
     * Simple constructor.
     * build an exception with a default message.
     * @param isCardan if true, the rotation is related to Cardan angles,
     * if false it is related to EulerAngles
     */
    public CardanEulerSingularityException(boolean isCardan) {
        super(isCardan ? LocalizedMigrationFormats.CARDAN_ANGLES_SINGULARITY : LocalizedMigrationFormats.EULER_ANGLES_SINGULARITY);
    }

}
