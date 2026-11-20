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
