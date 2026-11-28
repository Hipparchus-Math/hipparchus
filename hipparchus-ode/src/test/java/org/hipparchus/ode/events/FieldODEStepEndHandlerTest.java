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
package org.hipparchus.ode.events;

import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.util.Binary64;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class FieldODEStepEndHandlerTest {

    @Test
    void testDefaultMethods() {
        FieldODEStepEndHandler<Binary64> handler = (state, forward) -> Action.CONTINUE;
        
        FieldODEStateAndDerivative<Binary64> s = new FieldODEStateAndDerivative<>(new Binary64(0), new Binary64[]{new Binary64(0)}, new Binary64[]{new Binary64(0)});
        
        // Test default init (does nothing)
        assertDoesNotThrow(() -> handler.init(s, new Binary64(10)));
        
        // Test default resetState (returns current state)
        FieldODEState<Binary64> reset = handler.resetState(s);
        assertSame(s, reset);
    }
}
