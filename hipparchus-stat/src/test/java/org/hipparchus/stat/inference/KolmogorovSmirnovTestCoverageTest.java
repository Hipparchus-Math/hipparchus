package org.hipparchus.stat.inference;

import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.FieldMatrix;
import org.hipparchus.fraction.BigFraction;
import org.junit.jupiter.api.Test;
import java.lang.reflect.Method;

import static org.junit.jupiter.api.Assertions.*;

public class KolmogorovSmirnovTestCoverageTest {

    @Test
    public void testCdfExact() {
        KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
        int n = 5;
        double d = 0.4;
        // This will call cdf(d, n, true) -> exactK(d, n) -> createExactH(d, n)
        double p = ksTest.cdfExact(d, n);
        assertTrue(p > 0 && p < 1);
        
        // n=5, d=0.4. ninv = 0.2. ninvhalf = 0.1.
        // d > ninv.
        
        // Try very small d to cover d <= ninvhalf
        assertEquals(0.0, ksTest.cdf(0.05, 5, true));
        
        // Try d between ninvhalf and ninv to cover the second branch
        // n=5, ninv=0.2, ninvhalf=0.1. d=0.15
        double pSmall = ksTest.cdf(0.15, 5, true);
        assertTrue(pSmall > 0);
        
        // Try d >= 1
        assertEquals(1.0, ksTest.cdf(1.1, 5, true));
    }
    
    @Test
    public void testCdfLargeD() {
        KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
        int n = 5;
        // n=5, ninv=0.2. 1-ninv = 0.8. 0.8 <= d < 1
        double d = 0.9;
        double p = ksTest.cdf(d, n, false);
        assertEquals(1 - 2 * Math.pow(1 - d, n), p, 1e-10);
    }

    @Test
    public void testPelzGoodEdgeCases() {
        KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
        // pelzGood is called for n > 140
        int n = 200;
        double d = 0.1;
        double p = ksTest.pelzGood(d, n);
        assertTrue(p > 0 && p < 1);
    }
    
    @Test
    public void testCheckArray() throws Exception {
        KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
        Method checkArray = KolmogorovSmirnovTest.class.getDeclaredMethod("checkArray", double[].class);
        checkArray.setAccessible(true);
        
        // Null array
        try {
            checkArray.invoke(ksTest, (Object) null);
            fail("Expected NullArgumentException");
        } catch (Exception e) {
            // Success
        }
        
        // Short array
        try {
            checkArray.invoke(ksTest, (Object) new double[] { 1.0 });
            fail("Expected MathIllegalArgumentException");
        } catch (Exception e) {
            // Success
        }
    }

    @Test
    public void testOneSampleTestWithAlpha() {
        KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
        RealDistribution dist = new NormalDistribution(0, 1);
        double[] data = { 0.1, 0.2, -0.1, 0.3, 0.5 };
        boolean rejected = ksTest.kolmogorovSmirnovTest(dist, data, 0.05);
        assertFalse(rejected); // Should not be rejected at 0.05 alpha
    }

    @Test
    public void testTwoSampleStrict() {
        KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
        // Use overlapping samples to have a non-extreme p-value where strictness might matter
        double[] x = { 0.1, 0.3, 0.5, 0.7, 0.9 };
        double[] y = { 0.2, 0.3, 0.4, 0.6, 0.8 };
        double pStrict = ksTest.kolmogorovSmirnovTest(x, y, true);
        double pNonStrict = ksTest.kolmogorovSmirnovTest(x, y, false);
        // They should be different because of the tie at 0.3 or just because of how exactP handles equality
        assertNotEquals(pStrict, pNonStrict);
    }
}
