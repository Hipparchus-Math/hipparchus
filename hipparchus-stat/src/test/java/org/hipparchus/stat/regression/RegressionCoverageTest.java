package org.hipparchus.stat.regression;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class RegressionCoverageTest {

    @Test
    public void testRegressionResultsCoverage() {
        double[] params = new double[] { 1.0, 2.0 };
        // Symmetric compressed variance-covariance matrix: [var(p0), cov(p0, p1), var(p1)]
        double[][] varcov = new double[][] { { 0.1, 0.05 }, { 0.05, 0.2 } };
        long nobs = 100;
        int rank = 2;
        double sumy = 50.0;
        double sumysq = 200.0;
        double sse = 20.0;
        boolean containsConstant = true;
        boolean copyData = true;

        RegressionResults results = new RegressionResults(params, varcov, false, nobs, rank, sumy, sumysq, sse, containsConstant, copyData);

        assertEquals(1.0, results.getParameterEstimate(0));
        assertEquals(2.0, results.getParameterEstimate(1));
        assertArrayEquals(params, results.getParameterEstimates());
        assertEquals(Math.sqrt(0.1), results.getStdErrorOfEstimate(0), 1e-10);
        assertEquals(Math.sqrt(0.2), results.getStdErrorOfEstimate(1), 1e-10);
        assertNotNull(results.getStdErrorOfEstimates());
        assertEquals(0.05, results.getCovarianceOfParameters(0, 1));
        assertEquals(2, results.getNumberOfParameters());
        assertEquals(nobs, results.getN());
        assertEquals(sumysq - (sumy * sumy / nobs), results.getTotalSumSquares(), 1e-10);
        assertEquals(results.getTotalSumSquares() - sse, results.getRegressionSumSquares(), 1e-10);
        assertEquals(sse, results.getErrorSumSquares());
        assertEquals(sse / (nobs - rank), results.getMeanSquareError(), 1e-10);
        assertTrue(results.getRSquared() > 0);
        assertTrue(results.getAdjustedRSquared() > 0);
        assertTrue(results.hasIntercept());
    }

    @Test
    public void testRegressionResultsSymmetric() {
        double[] params = new double[] { 1.0, 2.0 };
        // Symmetric compressed: var(p0), cov(p0, p1), var(p1)
        double[][] varcov = new double[][] { { 0.1, 0.05, 0.2 } };
        RegressionResults results = new RegressionResults(params, varcov, true, 100, 2, 50.0, 200.0, 20.0, true, true);
        assertEquals(0.05, results.getCovarianceOfParameters(0, 1));
        assertEquals(0.2, results.getCovarianceOfParameters(1, 1));
    }

    @Test
    public void testRegressionResultsNull() {
        // Constructor with empty arrays
        RegressionResults results = new RegressionResults(new double[0], new double[0][0], true, 0, 0, 0, 0, 0, false, false);
        assertEquals(0, results.getNumberOfParameters());
        assertArrayEquals(new double[0], results.getParameterEstimates());
        // getParameterEstimate(0) will fail because range is [0, -1]
        assertThrows(MathIllegalArgumentException.class, () -> results.getParameterEstimate(0));
        assertThrows(MathIllegalArgumentException.class, () -> results.getStdErrorOfEstimate(0));
        assertArrayEquals(new double[0], results.getStdErrorOfEstimates());
        assertThrows(MathIllegalArgumentException.class, () -> results.getCovarianceOfParameters(0, 1));
    }

    @Test
    public void testMillerUpdatingRegressionAdditional() {
        MillerUpdatingRegression miller = new MillerUpdatingRegression(3, true);
        // orderOfRegressors includes the constant if included
        assertEquals(4, miller.getOrderOfRegressors().length);
        assertTrue(miller.hasIntercept());
        
        miller.addObservation(new double[] { 1.0, 2.0, 3.0 }, 10.0);
        miller.addObservation(new double[] { 2.0, 1.0, 4.0 }, 12.0);
        miller.addObservation(new double[] { 3.0, 3.0, 5.0 }, 15.0);
        miller.addObservation(new double[] { 4.0, 2.0, 6.0 }, 18.0);
        miller.addObservation(new double[] { 5.0, 3.0, 7.0 }, 20.0);
        
        RegressionResults results = miller.regress();
        assertNotNull(results);
        assertEquals(4, results.getNumberOfParameters()); // 3 variables + intercept
        
        double[] pcorr = miller.getPartialCorrelations(0);
        assertNotNull(pcorr);
        
        double diag = miller.getDiagonalOfHatMatrix(new double[] { 1.0, 2.0, 3.0 });
        assertTrue(diag >= 0);
    }
}
