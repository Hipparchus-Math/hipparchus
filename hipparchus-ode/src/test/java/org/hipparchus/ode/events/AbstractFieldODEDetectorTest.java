package org.hipparchus.ode.events;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class AbstractFieldODEDetectorTest {

    private static class MockFieldDetector<T extends CalculusFieldElement<T>> extends AbstractFieldODEDetector<MockFieldDetector<T>, T> {
        MockFieldDetector(FieldAdaptableInterval<T> maxCheck, int maxIter, BracketedRealFieldUnivariateSolver<T> solver, FieldODEEventHandler<T> handler) {
            super(maxCheck, maxIter, solver, handler);
        }

        @Override
        public T g(FieldODEStateAndDerivative<T> s) {
            return s.getTime().getField().getZero();
        }

        @Override
        protected MockFieldDetector<T> create(FieldAdaptableInterval<T> newMaxCheck, int newMaxIter, BracketedRealFieldUnivariateSolver<T> newSolver, FieldODEEventHandler<T> newHandler) {
            return new MockFieldDetector<>(newMaxCheck, newMaxIter, newSolver, newHandler);
        }
    }

    @Test
    void testBuilderMethods() {
        Binary64Field field = Binary64Field.getInstance();
        FieldAdaptableInterval<Binary64> maxCheck = FieldAdaptableInterval.of(100);
        int maxIter = 50;
        BracketedRealFieldUnivariateSolver<Binary64> solver = new FieldBracketingNthOrderBrentSolver<>(new Binary64(1e-10), new Binary64(1e-10), new Binary64(1e-10), 5);
        FieldODEEventHandler<Binary64> handler = (state, detector, increasing) -> Action.CONTINUE;

        MockFieldDetector<Binary64> detector = new MockFieldDetector<>(maxCheck, maxIter, solver, handler);

        assertEquals(100, detector.getMaxCheckInterval().currentInterval(null, true));
        assertEquals(maxIter, detector.getMaxIterationCount());
        assertSame(solver, detector.getSolver());
        assertSame(handler, detector.getHandler());

        MockFieldDetector<Binary64> d2 = detector.withMaxCheck(new Binary64(200));
        assertEquals(200, d2.getMaxCheckInterval().currentInterval(null, true));
        assertNotSame(detector, d2);

        MockFieldDetector<Binary64> d3 = detector.withMaxCheck(FieldAdaptableInterval.of(300));
        assertEquals(300, d3.getMaxCheckInterval().currentInterval(null, true));

        MockFieldDetector<Binary64> d4 = detector.withMaxIter(60);
        assertEquals(60, d4.getMaxIterationCount());

        MockFieldDetector<Binary64> d5 = detector.withThreshold(new Binary64(1.0e-12));
        assertTrue(d5.getSolver() instanceof FieldBracketingNthOrderBrentSolver);

        BracketedRealFieldUnivariateSolver<Binary64> newSolver = new FieldBracketingNthOrderBrentSolver<>(new Binary64(1e-10), new Binary64(1e-10), new Binary64(1e-10), 5);
        MockFieldDetector<Binary64> d6 = detector.withSolver(newSolver);
        assertSame(newSolver, d6.getSolver());

        FieldODEEventHandler<Binary64> newHandler = (state, d, inc) -> Action.STOP;
        MockFieldDetector<Binary64> d7 = detector.withHandler(newHandler);
        assertSame(newHandler, d7.getHandler());
    }

    @Test
    void testInitAndForward() {
        Binary64Field field = Binary64Field.getInstance();
        FieldODEEventHandler<Binary64> handler = (state, detector, increasing) -> Action.CONTINUE;
        MockFieldDetector<Binary64> detector = new MockFieldDetector<>(FieldAdaptableInterval.of(100), 50, 
                new FieldBracketingNthOrderBrentSolver<>(new Binary64(1e-10), new Binary64(1e-10), new Binary64(1e-10), 5), handler);
        FieldODEStateAndDerivative<Binary64> s0 = new FieldODEStateAndDerivative<>(new Binary64(0), new Binary64[]{new Binary64(0)}, new Binary64[]{new Binary64(0)});
        
        detector.init(s0, new Binary64(10));
        assertTrue(detector.isForward());

        detector.init(s0, new Binary64(-10));
        assertFalse(detector.isForward());
    }
}
