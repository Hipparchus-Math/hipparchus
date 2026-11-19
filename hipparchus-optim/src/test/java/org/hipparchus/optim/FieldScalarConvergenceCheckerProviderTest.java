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
package org.hipparchus.optim;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.util.Binary64Field;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import static org.junit.jupiter.api.Assertions.assertEquals;

class FieldScalarConvergenceCheckerProviderTest {

    @ParameterizedTest
    @ValueSource(booleans = {true, false})
    void testAlwaysSame(final boolean value) {
        // GIVEN
        final Binary64Field field = Binary64Field.getInstance();
        // WHEN
        final FieldScalarConvergenceCheckerProvider checkerProvider = new FieldScalarConvergenceCheckerProvider() {
            @Override
            public <T extends CalculusFieldElement<T>> ConvergenceChecker<T> getChecker(Field<T> field) {
                return (iteration, previous, current) -> value;
            }
        };
        // THEN
        assertEquals(value, checkerProvider.getChecker(field).converged(0, field.getZero(), field.getOne()));
    }
}
