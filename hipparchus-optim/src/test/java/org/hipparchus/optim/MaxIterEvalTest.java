package org.hipparchus.optim;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

class MaxIterEvalTest {

    @Test
    void testMaxIter() {
        MaxIter m = new MaxIter(100);
        assertEquals(100, m.getMaxIter());
        assertEquals(Integer.MAX_VALUE, MaxIter.unlimited().getMaxIter());
        assertThrows(MathIllegalArgumentException.class, () -> new MaxIter(0));
        assertThrows(MathIllegalArgumentException.class, () -> new MaxIter(-1));
    }

    @Test
    void testMaxEval() {
        MaxEval m = new MaxEval(100);
        assertEquals(100, m.getMaxEval());
        assertEquals(Integer.MAX_VALUE, MaxEval.unlimited().getMaxEval());
        assertThrows(MathIllegalArgumentException.class, () -> new MaxEval(0));
        assertThrows(MathIllegalArgumentException.class, () -> new MaxEval(-1));
    }
}
