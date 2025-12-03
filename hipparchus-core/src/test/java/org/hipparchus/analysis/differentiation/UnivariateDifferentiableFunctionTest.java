package org.hipparchus.analysis.differentiation;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class UnivariateDifferentiableFunctionTest {

    @Test
    void testValue() {
        // GIVEN
        final UnivariateDifferentiableFunction function = new UnivariateDifferentiableFunction() {
            @Override
            public <T extends Derivative<T>> T value(T x) throws MathIllegalArgumentException {
                return x.square();
            }
        };
        // WHEN
        final double value = function.value(2.0);
        // THEN
        assertEquals(4.0, value);
    }
}
