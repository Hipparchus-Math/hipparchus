package org.hipparchus.stat.descriptive.vector;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.stat.descriptive.moment.Mean;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class VectorialCoverageTest {

    @Test
    public void testVectorialCovariance() {
        VectorialCovariance cov = new VectorialCovariance(2, true);
        assertEquals(0, cov.getN());
        assertNotNull(cov.getResult());
        assertEquals(0.0, cov.getResult().getEntry(0, 0));

        cov.increment(new double[] { 1.0, 2.0 });
        assertEquals(1, cov.getN());
        assertEquals(0.0, cov.getResult().getEntry(0, 0)); // n=1, result is zero matrix

        cov.increment(new double[] { 2.0, 3.0 });
        assertEquals(2, cov.getN());
        assertTrue(cov.getResult().getEntry(0, 0) > 0);

        VectorialCovariance cov2 = new VectorialCovariance(2, true);
        cov2.increment(new double[] { 1.0, 2.0 });
        cov2.increment(new double[] { 2.0, 3.0 });

        assertEquals(cov, cov2);
        assertEquals(cov.hashCode(), cov2.hashCode());
        assertNotEquals(cov, new Object());
        
        cov.clear();
        assertEquals(0, cov.getN());
        assertEquals(0.0, cov.getResult().getEntry(0, 0));
        assertNotEquals(cov, cov2);

        VectorialCovariance biased = new VectorialCovariance(2, false);
        biased.increment(new double[] { 1.0, 2.0 });
        biased.increment(new double[] { 2.0, 3.0 });
        assertNotEquals(cov2, biased);
    }

    @Test
    public void testVectorialStorelessStatistic() {
        VectorialStorelessStatistic stat = new VectorialStorelessStatistic(2, new Mean());
        assertEquals(2, stat.getDimension());
        assertEquals(0, stat.getN());

        stat.increment(new double[] { 1.0, 2.0 });
        stat.increment(new double[] { 3.0, 4.0 });
        assertEquals(2, stat.getN());
        assertArrayEquals(new double[] { 2.0, 3.0 }, stat.getResult());

        VectorialStorelessStatistic stat2 = new VectorialStorelessStatistic(2, new Mean());
        stat2.increment(new double[] { 1.0, 2.0 });
        stat2.increment(new double[] { 3.0, 4.0 });

        assertEquals(stat, stat2);
        assertEquals(stat.hashCode(), stat2.hashCode());
        assertNotEquals(stat, new Object());

        stat.clear();
        assertEquals(0, stat.getN());
        assertNotEquals(stat, stat2);
        
        assertThrows(MathIllegalArgumentException.class, () -> new VectorialStorelessStatistic(0, new Mean()));
    }
}
