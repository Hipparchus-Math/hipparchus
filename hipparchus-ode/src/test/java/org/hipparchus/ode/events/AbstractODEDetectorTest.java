package org.hipparchus.ode.events;

import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class AbstractODEDetectorTest {

    private static class MockDetector extends AbstractODEDetector<MockDetector> {
        MockDetector(AdaptableInterval maxCheck, int maxIter, BracketedUnivariateSolver<UnivariateFunction> solver, ODEEventHandler handler) {
            super(maxCheck, maxIter, solver, handler);
        }

        @Override
        public double g(ODEStateAndDerivative s) {
            return 0;
        }

        @Override
        protected MockDetector create(AdaptableInterval newMaxCheck, int newMaxIter, BracketedUnivariateSolver<UnivariateFunction> newSolver, ODEEventHandler newHandler) {
            return new MockDetector(newMaxCheck, newMaxIter, newSolver, newHandler);
        }
    }

    @Test
    void testBuilderMethods() {
        AdaptableInterval maxCheck = AdaptableInterval.of(100);
        int maxIter = 50;
        BracketedUnivariateSolver<UnivariateFunction> solver = new BracketingNthOrderBrentSolver();
        ODEEventHandler handler = (state, detector, increasing) -> Action.CONTINUE;

        MockDetector detector = new MockDetector(maxCheck, maxIter, solver, handler);

        assertEquals(100, detector.getMaxCheckInterval().currentInterval(null, true));
        assertEquals(maxIter, detector.getMaxIterationCount());
        assertSame(solver, detector.getSolver());
        assertSame(handler, detector.getHandler());

        MockDetector d2 = detector.withMaxCheck(200);
        assertEquals(200, d2.getMaxCheckInterval().currentInterval(null, true));
        assertNotSame(detector, d2);

        MockDetector d3 = detector.withMaxCheck(AdaptableInterval.of(300));
        assertEquals(300, d3.getMaxCheckInterval().currentInterval(null, true));

        MockDetector d4 = detector.withMaxIter(60);
        assertEquals(60, d4.getMaxIterationCount());

        MockDetector d5 = detector.withThreshold(1.0e-10);
        assertTrue(d5.getSolver() instanceof BracketingNthOrderBrentSolver);

        BracketedUnivariateSolver<UnivariateFunction> newSolver = new BracketingNthOrderBrentSolver();
        MockDetector d6 = detector.withSolver(newSolver);
        assertSame(newSolver, d6.getSolver());

        ODEEventHandler newHandler = (state, d, inc) -> Action.STOP;
        MockDetector d7 = detector.withHandler(newHandler);
        assertSame(newHandler, d7.getHandler());
    }

    @Test
    void testInitAndForward() {
        ODEEventHandler handler = (state, detector, increasing) -> Action.CONTINUE;
        MockDetector detector = new MockDetector(AdaptableInterval.of(100), 50, new BracketingNthOrderBrentSolver(), handler);
        ODEStateAndDerivative s0 = new ODEStateAndDerivative(0, new double[1], new double[1]);
        
        detector.init(s0, 10);
        assertTrue(detector.isForward());

        detector.init(s0, -10);
        assertFalse(detector.isForward());
    }
}
