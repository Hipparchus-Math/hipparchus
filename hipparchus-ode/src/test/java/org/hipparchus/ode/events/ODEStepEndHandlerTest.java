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

import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class ODEStepEndHandlerTest {

    @Test
    void testDefaultMethods() {
        ODEStepEndHandler handler = (state, forward) -> Action.CONTINUE;
        
        ODEStateAndDerivative s = new ODEStateAndDerivative(0, new double[1], new double[1]);
        
        // Test default init (does nothing)
        assertDoesNotThrow(() -> handler.init(s, 10));
        
        // Test default resetState (returns current state)
        ODEState reset = handler.resetState(s);
        assertSame(s, reset);
    }
}
