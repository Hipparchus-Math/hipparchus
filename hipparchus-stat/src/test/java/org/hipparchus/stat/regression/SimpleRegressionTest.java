/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.stat.regression;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.ISAACRandom;
import org.hipparchus.stat.LocalizedStatFormats;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Random;


/**
 * Test cases for the TestStatistic class.
 *
 */

public final class SimpleRegressionTest {

    /*
     * NIST "Norris" refernce data set from
     * http://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Norris.dat
     * Strangely, order is {y,x}
     */
    private double[][] data = { { 0.1, 0.2 }, {338.8, 337.4 }, {118.1, 118.2 },
            {888.0, 884.6 }, {9.2, 10.1 }, {228.1, 226.5 }, {668.5, 666.3 }, {998.5, 996.3 },
            {449.1, 448.6 }, {778.9, 777.0 }, {559.2, 558.2 }, {0.3, 0.4 }, {0.1, 0.6 }, {778.1, 775.5 },
            {668.8, 666.9 }, {339.3, 338.0 }, {448.9, 447.5 }, {10.8, 11.6 }, {557.7, 556.0 },
            {228.3, 228.1 }, {998.0, 995.8 }, {888.8, 887.6 }, {119.6, 120.2 }, {0.3, 0.3 },
            {0.6, 0.3 }, {557.6, 556.8 }, {339.3, 339.1 }, {888.0, 887.2 }, {998.5, 999.0 },
            {778.9, 779.0 }, {10.2, 11.1 }, {117.6, 118.3 }, {228.9, 229.2 }, {668.4, 669.1 },
            {449.2, 448.9 }, {0.2, 0.5 }
    };

    /*
     * Correlation example from
     * http://www.xycoon.com/correlation.htm
     */
    private double[][] corrData = { { 101.0, 99.2 }, {100.1, 99.0 }, {100.0, 100.0 },
            {90.6, 111.6 }, {86.5, 122.2 }, {89.7, 117.6 }, {90.6, 121.1 }, {82.8, 136.0 },
            {70.1, 154.2 }, {65.4, 153.6 }, {61.3, 158.5 }, {62.5, 140.6 }, {63.6, 136.2 },
            {52.6, 168.0 }, {59.7, 154.3 }, {59.5, 149.0 }, {61.3, 165.5 }
    };

    /*
     * From Moore and Mcabe, "Introduction to the Practice of Statistics"
     * Example 10.3
     */
    private double[][] infData = { { 15.6, 5.2 }, {26.8, 6.1 }, {37.8, 8.7 }, {36.4, 8.5 },
            {35.5, 8.8 }, {18.6, 4.9 }, {15.3, 4.5 }, {7.9, 2.5 }, {0.0, 1.1 }
    };

    /*
     * Points to remove in the remove tests
     */
    private double[][] removeSingle = {infData[1]};
    private double[][] removeMultiple = { infData[1], infData[2] };
    private double removeX = infData[0][0];
    private double removeY = infData[0][1];


    /*
     * Data with bad linear fit
     */
    private double[][] infData2 = { { 1, 1 }, {2, 0 }, {3, 5 }, {4, 2 },
            {5, -1 }, {6, 12 }
    };


    /*
     * Data from NIST NOINT1
     */
    private double[][] noint1 = {
        {130.0,60.0},
        {131.0,61.0},
        {132.0,62.0},
        {133.0,63.0},
        {134.0,64.0},
        {135.0,65.0},
        {136.0,66.0},
        {137.0,67.0},
        {138.0,68.0},
        {139.0,69.0},
        {140.0,70.0}
    };

    /*
     * Data from NIST NOINT2
     *
     */
    private double[][] noint2 = {
        {3.0,4},
        {4,5},
        {4,6}
    };


    /**
     * Test that the SimpleRegression objects generated from combining two
     * SimpleRegression objects created from subsets of data are identical to
     * SimpleRegression objects created from the combined data.
     */
    @Test
    public void testAppend() {
        check(false);
        check(true);
    }

    /**
     * Checks that adding data to a single model gives the same result
     * as adding "parts" of the dataset to smaller models and using append
     * to aggregate the smaller models.
     *
     * @param includeIntercept
     */
    private void check(boolean includeIntercept) {
        final int sets = 2;
        final ISAACRandom rand = new ISAACRandom(10L);// Seed can be changed
        final SimpleRegression whole = new SimpleRegression(includeIntercept);// regression of the whole set
        final SimpleRegression parts = new SimpleRegression(includeIntercept);// regression with parts.

        for (int s = 0; s < sets; s++) {// loop through each subset of data.
            final double coef = rand.nextDouble();
            final SimpleRegression sub = new SimpleRegression(includeIntercept);// sub regression
            for (int i = 0; i < 5; i++) { // loop through individual samlpes.
                final double x = rand.nextDouble();
                final double y = x * coef + rand.nextDouble();// some noise
                sub.addData(x, y);
                whole.addData(x, y);
            }
            parts.append(sub);
            Assertions.assertTrue(equals(parts, whole, 1E-6));
        }
    }

    /**
     * Returns true iff the statistics reported by model1 are all within tol of
     * those reported by model2.
     *
     * @param model1 first model
     * @param model2 second model
     * @param tol tolerance
     * @return true if the two models report the same regression stats
     */
    private boolean equals(SimpleRegression model1, SimpleRegression model2, double tol) {
        if (model1.getN() != model2.getN()) {
            return false;
        }
        if (FastMath.abs(model1.getIntercept() - model2.getIntercept()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getInterceptStdErr() - model2.getInterceptStdErr()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getMeanSquareError() - model2.getMeanSquareError()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getR() - model2.getR()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getRegressionSumSquares() - model2.getRegressionSumSquares()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getRSquare() - model2.getRSquare()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getSignificance() - model2.getSignificance()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getSlope() - model2.getSlope()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getSlopeConfidenceInterval() - model2.getSlopeConfidenceInterval()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getSlopeStdErr() - model2.getSlopeStdErr()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getSumOfCrossProducts() - model2.getSumOfCrossProducts()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getSumSquaredErrors() - model2.getSumSquaredErrors()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getTotalSumSquares() - model2.getTotalSumSquares()) > tol) {
            return false;
        }
        if (FastMath.abs(model1.getXSumSquares() - model2.getXSumSquares()) > tol) {
            return false;
        }
        return true;
    }

    @Test
    public void testRegressIfaceMethod(){
        final SimpleRegression regression = new SimpleRegression(true);
        final UpdatingMultipleLinearRegression iface = regression;
        final SimpleRegression regressionNoint = new SimpleRegression( false );
        final SimpleRegression regressionIntOnly= new SimpleRegression( false );
        for (int i = 0; i < data.length; i++) {
            iface.addObservation( new double[]{data[i][1]}, data[i][0]);
            regressionNoint.addData(data[i][1], data[i][0]);
            regressionIntOnly.addData(1.0, data[i][0]);
        }

        //should not be null
        final RegressionResults fullReg = iface.regress( );
        Assertions.assertNotNull(fullReg);
        Assertions.assertEquals(regression.getIntercept(), fullReg.getParameterEstimate(0), 1.0e-16, "intercept");
        Assertions.assertEquals(regression.getInterceptStdErr(), fullReg.getStdErrorOfEstimate(0),1.0E-16,"intercept std err");
        Assertions.assertEquals(regression.getSlope(), fullReg.getParameterEstimate(1), 1.0e-16, "slope");
        Assertions.assertEquals(regression.getSlopeStdErr(), fullReg.getStdErrorOfEstimate(1),1.0E-16,"slope std err");
        Assertions.assertEquals(regression.getN(), fullReg.getN(), "number of observations");
        Assertions.assertEquals(regression.getRSquare(), fullReg.getRSquared(), 1.0E-16, "r-square");
        Assertions.assertEquals(regression.getRegressionSumSquares(), fullReg.getRegressionSumSquares() ,1.0E-16,"SSR");
        Assertions.assertEquals(regression.getMeanSquareError(), fullReg.getMeanSquareError() ,1.0E-16,"MSE");
        Assertions.assertEquals(regression.getSumSquaredErrors(), fullReg.getErrorSumSquares() ,1.0E-16,"SSE");


        final RegressionResults noInt   = iface.regress( new int[]{1} );
        Assertions.assertNotNull(noInt);
        Assertions.assertEquals(regressionNoint.getSlope(), noInt.getParameterEstimate(0), 1.0e-12, "slope");
        Assertions.assertEquals(regressionNoint.getSlopeStdErr(), noInt.getStdErrorOfEstimate(0),1.0E-16,"slope std err");
        Assertions.assertEquals(regressionNoint.getN(), noInt.getN(), "number of observations");
        Assertions.assertEquals(regressionNoint.getRSquare(), noInt.getRSquared(), 1.0E-16, "r-square");
        Assertions.assertEquals(regressionNoint.getRegressionSumSquares(), noInt.getRegressionSumSquares() ,1.0E-8,"SSR");
        Assertions.assertEquals(regressionNoint.getMeanSquareError(), noInt.getMeanSquareError() ,1.0E-16,"MSE");
        Assertions.assertEquals(regressionNoint.getSumSquaredErrors(), noInt.getErrorSumSquares() ,1.0E-16,"SSE");

        final RegressionResults onlyInt = iface.regress( new int[]{0} );
        Assertions.assertNotNull(onlyInt);
        Assertions.assertEquals(regressionIntOnly.getSlope(), onlyInt.getParameterEstimate(0), 1.0e-12, "slope");
        Assertions.assertEquals(regressionIntOnly.getSlopeStdErr(), onlyInt.getStdErrorOfEstimate(0),1.0E-12,"slope std err");
        Assertions.assertEquals(regressionIntOnly.getN(), onlyInt.getN(), "number of observations");
        Assertions.assertEquals(regressionIntOnly.getRSquare(), onlyInt.getRSquared(), 1.0E-14, "r-square");
        Assertions.assertEquals(regressionIntOnly.getSumSquaredErrors(), onlyInt.getErrorSumSquares() ,1.0E-8,"SSE");
        Assertions.assertEquals(regressionIntOnly.getRegressionSumSquares(), onlyInt.getRegressionSumSquares() ,1.0E-8,"SSR");
        Assertions.assertEquals(regressionIntOnly.getMeanSquareError(), onlyInt.getMeanSquareError() ,1.0E-8,"MSE");

    }

    /**
     * Verify that regress generates exceptions as advertised for bad model specifications.
     */
    @Test
    public void testRegressExceptions() {
        // No intercept
        final SimpleRegression noIntRegression = new SimpleRegression(false);
        noIntRegression.addData(noint2[0][1], noint2[0][0]);
        noIntRegression.addData(noint2[1][1], noint2[1][0]);
        noIntRegression.addData(noint2[2][1], noint2[2][0]);
        try { // null array
            noIntRegression.regress(null);
            Assertions.fail("Expecting MathIllegalArgumentException for null array");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try { // empty array
            noIntRegression.regress(new int[] {});
            Assertions.fail("Expecting MathIllegalArgumentException for empty array");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try { // more than 1 regressor
            noIntRegression.regress(new int[] {0, 1});
            Assertions.fail("Expecting MathIllegalArgumentException - too many regressors");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try { // invalid regressor
            noIntRegression.regress(new int[] {1});
            Assertions.fail("Expecting MathIllegalArgumentException - invalid regression");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }

        // With intercept
        final SimpleRegression regression = new SimpleRegression(true);
        regression.addData(noint2[0][1], noint2[0][0]);
        regression.addData(noint2[1][1], noint2[1][0]);
        regression.addData(noint2[2][1], noint2[2][0]);
        try { // null array
            regression.regress(null);
            Assertions.fail("Expecting MathIllegalArgumentException for null array");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try { // empty array
            regression.regress(new int[] {});
            Assertions.fail("Expecting MathIllegalArgumentException for empty array");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try { // more than 2 regressors
            regression.regress(new int[] {0, 1, 2});
            Assertions.fail("Expecting MathIllegalArgumentException - too many regressors");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try { // wrong order
            regression.regress(new int[] {1,0});
            Assertions.fail("Expecting MathIllegalArgumentException - invalid regression");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try { // out of range
            regression.regress(new int[] {3,4});
            Assertions.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try { // out of range
            regression.regress(new int[] {0,2});
            Assertions.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try { // out of range
            regression.regress(new int[] {2});
            Assertions.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
    }

    @Test
    public void testNoInterceot_noint2(){
         SimpleRegression regression = new SimpleRegression(false);
         regression.addData(noint2[0][1], noint2[0][0]);
         regression.addData(noint2[1][1], noint2[1][0]);
         regression.addData(noint2[2][1], noint2[2][0]);
         Assertions.assertEquals(0, regression.getIntercept(), 0, "intercept");
         Assertions.assertEquals(0.727272727272727,
                 regression.getSlope(), 10E-12, "slope");
         Assertions.assertEquals(0.420827318078432E-01,
                regression.getSlopeStdErr(),10E-12,"slope std err");
        Assertions.assertEquals(3, regression.getN(), "number of observations");
        Assertions.assertEquals(0.993348115299335,
            regression.getRSquare(), 10E-12, "r-square");
        Assertions.assertEquals(40.7272727272727,
            regression.getRegressionSumSquares(), 10E-9, "SSR");
        Assertions.assertEquals(0.136363636363636,
            regression.getMeanSquareError(), 10E-10, "MSE");
        Assertions.assertEquals(0.272727272727273,
            regression.getSumSquaredErrors(),10E-9,"SSE");
    }

    @Test
    public void testNoIntercept_noint1(){
        SimpleRegression regression = new SimpleRegression(false);
        for (int i = 0; i < noint1.length; i++) {
            regression.addData(noint1[i][1], noint1[i][0]);
        }
        Assertions.assertEquals(0, regression.getIntercept(), 0, "intercept");
        Assertions.assertEquals(2.07438016528926, regression.getSlope(), 10E-12, "slope");
        Assertions.assertEquals(0.165289256198347E-01,
                regression.getSlopeStdErr(),10E-12,"slope std err");
        Assertions.assertEquals(11, regression.getN(), "number of observations");
        Assertions.assertEquals(0.999365492298663,
            regression.getRSquare(), 10E-12, "r-square");
        Assertions.assertEquals(200457.727272727,
            regression.getRegressionSumSquares(), 10E-9, "SSR");
        Assertions.assertEquals(12.7272727272727,
            regression.getMeanSquareError(), 10E-10, "MSE");
        Assertions.assertEquals(127.272727272727,
            regression.getSumSquaredErrors(),10E-9,"SSE");

    }

    @Test
    public void testNorris() {
        SimpleRegression regression = new SimpleRegression();
        for (int i = 0; i < data.length; i++) {
            regression.addData(data[i][1], data[i][0]);
        }
        // Tests against certified values from
        // http://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Norris.dat
        Assertions.assertEquals(1.00211681802045, regression.getSlope(), 10E-12, "slope");
        Assertions.assertEquals(0.429796848199937E-03,
                regression.getSlopeStdErr(),10E-12,"slope std err");
        Assertions.assertEquals(36, regression.getN(), "number of observations");
        Assertions.assertEquals( -0.262323073774029,
            regression.getIntercept(),10E-12,"intercept");
        Assertions.assertEquals(0.232818234301152,
            regression.getInterceptStdErr(),10E-12,"std err intercept");
        Assertions.assertEquals(0.999993745883712,
            regression.getRSquare(), 10E-12, "r-square");
        Assertions.assertEquals(4255954.13232369,
            regression.getRegressionSumSquares(), 10E-9, "SSR");
        Assertions.assertEquals(0.782864662630069,
            regression.getMeanSquareError(), 10E-10, "MSE");
        Assertions.assertEquals(26.6173985294224,
            regression.getSumSquaredErrors(),10E-9,"SSE");
        // ------------  End certified data tests

        Assertions.assertEquals( -0.262323073774029,
            regression.predict(0), 10E-12, "predict(0)");
        Assertions.assertEquals(1.00211681802045 - 0.262323073774029,
            regression.predict(1), 10E-12, "predict(1)");
    }

    @Test
    public void testCorr() {
        SimpleRegression regression = new SimpleRegression();
        regression.addData(corrData);
        Assertions.assertEquals(17, regression.getN(), "number of observations");
        Assertions.assertEquals(.896123, regression.getRSquare(), 10E-6, "r-square");
        Assertions.assertEquals(-0.94663767742, regression.getR(), 1E-10, "r");
    }

    @Test
    public void testNaNs() {
        SimpleRegression regression = new SimpleRegression();
        Assertions.assertTrue(Double.isNaN(regression.getIntercept()), "intercept not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getSlope()), "slope not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getSlopeStdErr()), "slope std err not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getInterceptStdErr()), "intercept std err not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getMeanSquareError()), "MSE not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getR()), "e not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getRSquare()), "r-square not NaN");
        Assertions.assertTrue( Double.isNaN(regression.getRegressionSumSquares()), "RSS not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getSumSquaredErrors()),"SSE not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getTotalSumSquares()), "SSTO not NaN");
        Assertions.assertTrue(Double.isNaN(regression.predict(0)), "predict not NaN");

        regression.addData(1, 2);
        regression.addData(1, 3);

        // No x variation, so these should still blow...
        Assertions.assertTrue(Double.isNaN(regression.getIntercept()), "intercept not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getSlope()), "slope not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getSlopeStdErr()), "slope std err not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getInterceptStdErr()), "intercept std err not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getMeanSquareError()), "MSE not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getR()), "e not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getRSquare()), "r-square not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getRegressionSumSquares()), "RSS not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getSumSquaredErrors()), "SSE not NaN");
        Assertions.assertTrue(Double.isNaN(regression.predict(0)), "predict not NaN");

        // but SSTO should be OK
        Assertions.assertFalse(Double.isNaN(regression.getTotalSumSquares()), "SSTO NaN");

        regression = new SimpleRegression();

        regression.addData(1, 2);
        regression.addData(3, 3);

        // All should be OK except MSE, s(b0), s(b1) which need one more df
        Assertions.assertFalse(Double.isNaN(regression.getIntercept()), "interceptNaN");
        Assertions.assertFalse(Double.isNaN(regression.getSlope()), "slope NaN");
        Assertions.assertTrue(Double.isNaN(regression.getSlopeStdErr()), "slope std err not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getInterceptStdErr()), "intercept std err not NaN");
        Assertions.assertTrue(Double.isNaN(regression.getMeanSquareError()), "MSE not NaN");
        Assertions.assertFalse(Double.isNaN(regression.getR()), "r NaN");
        Assertions.assertFalse(Double.isNaN(regression.getRSquare()), "r-square NaN");
        Assertions.assertFalse(Double.isNaN(regression.getRegressionSumSquares()), "RSS NaN");
        Assertions.assertFalse(Double.isNaN(regression.getSumSquaredErrors()), "SSE NaN");
        Assertions.assertFalse(Double.isNaN(regression.getTotalSumSquares()), "SSTO NaN");
        Assertions.assertFalse(Double.isNaN(regression.predict(0)), "predict NaN");

        regression.addData(1, 4);

        // MSE, MSE, s(b0), s(b1) should all be OK now
        Assertions.assertFalse(Double.isNaN(regression.getMeanSquareError()), "MSE NaN");
        Assertions.assertFalse(Double.isNaN(regression.getSlopeStdErr()), "slope std err NaN");
        Assertions.assertFalse(Double.isNaN(regression.getInterceptStdErr()), "intercept std err NaN");
    }

    @Test
    public void testClear() {
        SimpleRegression regression = new SimpleRegression();
        regression.addData(corrData);
        Assertions.assertEquals(17, regression.getN(), "number of observations");
        regression.clear();
        Assertions.assertEquals(0, regression.getN(), "number of observations");
        regression.addData(corrData);
        Assertions.assertEquals(.896123, regression.getRSquare(), 10E-6, "r-square");
        regression.addData(data);
        Assertions.assertEquals(53, regression.getN(), "number of observations");
    }

    @Test
    public void testInference() {
        //----------  verified against R, version 1.8.1 -----
        // infData
        SimpleRegression regression = new SimpleRegression();
        regression.addData(infData);
        Assertions.assertEquals(0.011448491,
                regression.getSlopeStdErr(), 1E-10, "slope std err");
        Assertions.assertEquals(0.286036932,
                regression.getInterceptStdErr(),1E-8,"std err intercept");
        Assertions.assertEquals(4.596e-07,
                regression.getSignificance(),1E-8,"significance");
        Assertions.assertEquals(0.0270713794287,
                regression.getSlopeConfidenceInterval(),1E-8,"slope conf interval half-width");
        // infData2
        regression = new SimpleRegression();
        regression.addData(infData2);
        Assertions.assertEquals(1.07260253,
                regression.getSlopeStdErr(), 1E-8, "slope std err");
        Assertions.assertEquals(4.17718672,
                regression.getInterceptStdErr(),1E-8,"std err intercept");
        Assertions.assertEquals(0.261829133982,
                regression.getSignificance(),1E-11,"significance");
        Assertions.assertEquals(2.97802204827,
                regression.getSlopeConfidenceInterval(),1E-8,"slope conf interval half-width");
        //------------- End R-verified tests -------------------------------

        //FIXME: get a real example to test against with alpha = .01
        Assertions.assertTrue(regression.getSlopeConfidenceInterval() < regression.getSlopeConfidenceInterval(0.01),
                "tighter means wider");

        try {
            regression.getSlopeConfidenceInterval(1);
            Assertions.fail("expecting MathIllegalArgumentException for alpha = 1");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }

    }

    @Test
    public void testPerfect() {
        SimpleRegression regression = new SimpleRegression();
        int n = 100;
        for (int i = 0; i < n; i++) {
            regression.addData(((double) i) / (n - 1), i);
        }
        Assertions.assertEquals(0.0, regression.getSignificance(), 1.0e-5);
        Assertions.assertTrue(regression.getSlope() > 0.0);
        Assertions.assertTrue(regression.getSumSquaredErrors() >= 0.0);
    }

    @Test
    public void testPerfect2() {
        SimpleRegression regression = new SimpleRegression();
        regression.addData(0, 0);
        regression.addData(1, 1);
        regression.addData(2, 2);
        Assertions.assertEquals(0.0, regression.getSlopeStdErr(), 0.0);
        Assertions.assertEquals(0.0, regression.getSignificance(), Double.MIN_VALUE);
        Assertions.assertEquals(1, regression.getRSquare(), Double.MIN_VALUE);
    }

    @Test
    public void testPerfectNegative() {
        SimpleRegression regression = new SimpleRegression();
        int n = 100;
        for (int i = 0; i < n; i++) {
            regression.addData(- ((double) i) / (n - 1), i);
        }

        Assertions.assertEquals(0.0, regression.getSignificance(), 1.0e-5);
        Assertions.assertTrue(regression.getSlope() < 0.0);
    }

    @Test
    public void testRandom() {
        SimpleRegression regression = new SimpleRegression();
        Random random = new Random(1);
        int n = 100;
        for (int i = 0; i < n; i++) {
            regression.addData(((double) i) / (n - 1), random.nextDouble());
        }

        Assertions.assertTrue( 0.0 < regression.getSignificance()
                    && regression.getSignificance() < 1.0);
    }


    // Jira MATH-85 = Bugzilla 39432
    @Test
    public void testSSENonNegative() {
        double[] y = { 8915.102, 8919.302, 8923.502 };
        double[] x = { 1.107178495E2, 1.107264895E2, 1.107351295E2 };
        SimpleRegression reg = new SimpleRegression();
        for (int i = 0; i < x.length; i++) {
            reg.addData(x[i], y[i]);
        }
        Assertions.assertTrue(reg.getSumSquaredErrors() >= 0.0);
    }

    // Test remove X,Y (single observation)
    @Test
    public void testRemoveXY() {
        // Create regression with inference data then remove to test
        SimpleRegression regression = new SimpleRegression();
        Assertions.assertTrue(regression.hasIntercept());
        regression.addData(infData);
        regression.removeData(removeX, removeY);
        regression.addData(removeX, removeY);
        // Use the inference assertions to make sure that everything worked
        Assertions.assertEquals(0.011448491,
                regression.getSlopeStdErr(), 1E-10, "slope std err");
        Assertions.assertEquals(0.286036932,
                regression.getInterceptStdErr(),1E-8,"std err intercept");
        Assertions.assertEquals(4.596e-07,
                regression.getSignificance(),1E-8,"significance");
        Assertions.assertEquals(0.0270713794287,
                regression.getSlopeConfidenceInterval(),1E-8,"slope conf interval half-width");
     }

    // Test remove single observation in array
    @Test
    public void testRemoveSingle() {
        // Create regression with inference data then remove to test
        SimpleRegression regression = new SimpleRegression();
        Assertions.assertTrue(regression.hasIntercept());
        regression.addData(infData);
        regression.removeData(removeSingle);
        regression.addData(removeSingle);
        // Use the inference assertions to make sure that everything worked
        Assertions.assertEquals(0.011448491,
                regression.getSlopeStdErr(), 1E-10, "slope std err");
        Assertions.assertEquals(0.286036932,
                regression.getInterceptStdErr(),1E-8,"std err intercept");
        Assertions.assertEquals(4.596e-07,
                regression.getSignificance(),1E-8,"significance");
        Assertions.assertEquals(0.0270713794287,
                regression.getSlopeConfidenceInterval(),1E-8,"slope conf interval half-width");
     }

    // Test remove multiple observations
    @Test
    public void testRemoveMultiple() {
        // Create regression with inference data then remove to test
        SimpleRegression regression = new SimpleRegression();
        Assertions.assertTrue(regression.hasIntercept());
        regression.addData(infData);
        regression.removeData(removeMultiple);
        regression.addData(removeMultiple);
        // Use the inference assertions to make sure that everything worked
        Assertions.assertEquals(0.011448491,
                regression.getSlopeStdErr(), 1E-10, "slope std err");
        Assertions.assertEquals(0.286036932,
                regression.getInterceptStdErr(),1E-8,"std err intercept");
        Assertions.assertEquals(4.596e-07,
                regression.getSignificance(),1E-8,"significance");
        Assertions.assertEquals(0.0270713794287,
                regression.getSlopeConfidenceInterval(),1E-8,"slope conf interval half-width");
     }

    // Test remove multiple observations
    @Test
    public void testRemoveMultipleNoIntercept() {
        // Create regression with inference data then remove to test
        SimpleRegression regression = new SimpleRegression(false);
        Assertions.assertFalse(regression.hasIntercept());
        Assertions.assertEquals(0.0, regression.getIntercept(), 1.0e-15);
        regression.addData(infData);
        Assertions.assertEquals(0.30593, regression.predict(1.25), 1.0e-5);
        regression.removeData(removeMultiple);
        regression.addData(removeMultiple);
        // Use the inference assertions to make sure that everything worked
        Assertions.assertEquals(0.0103629732,
                regression.getSlopeStdErr(), 1E-10, "slope std err");
        Assertions.assertTrue(Double.isNaN(regression.getInterceptStdErr()), "std err intercept");
        Assertions.assertEquals(6.199e-08,
                regression.getSignificance(),1E-10,"significance");
        Assertions.assertEquals(0.02450454,
                regression.getSlopeConfidenceInterval(),1E-8,"slope conf interval half-width");
     }

    // Remove observation when empty
    @Test
    public void testRemoveObsFromEmpty() {
        SimpleRegression regression = new SimpleRegression();
        regression.removeData(removeX, removeY);
        Assertions.assertEquals(0, regression.getN());
    }

    // Remove single observation to empty
    @Test
    public void testRemoveObsFromSingle() {
        SimpleRegression regression = new SimpleRegression();
        regression.addData(removeX, removeY);
        regression.removeData(removeX, removeY);
        Assertions.assertEquals(0, regression.getN());
    }

    // Remove multiple observations to empty
    @Test
    public void testRemoveMultipleToEmpty() {
        SimpleRegression regression = new SimpleRegression();
        regression.addData(removeMultiple);
        regression.removeData(removeMultiple);
        Assertions.assertEquals(0, regression.getN());
    }

    // Remove multiple observations past empty (i.e. size of array > n)
    @Test
    public void testRemoveMultiplePastEmpty() {
        SimpleRegression regression = new SimpleRegression();
        regression.addData(removeX, removeY);
        regression.removeData(removeMultiple);
        Assertions.assertEquals(0, regression.getN());
    }

    @Test
    public void testWrongDimensions() {
        try {
            new SimpleRegression().addData(new double[1][1]);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedStatFormats.INVALID_REGRESSION_OBSERVATION, miae.getSpecifier());
        }
        try {
            new SimpleRegression().addObservation(null, 0.0);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedStatFormats.INVALID_REGRESSION_OBSERVATION, miae.getSpecifier());
        }
        try {
            new SimpleRegression().addObservation(new double[0], 0.0);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedStatFormats.INVALID_REGRESSION_OBSERVATION, miae.getSpecifier());
        }
        try {
            new SimpleRegression().addObservations(new double[][] { null, null }, new double[2]);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedStatFormats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS, miae.getSpecifier());
        }
        try {
            new SimpleRegression().addObservations(new double[][] { new double[0], new double[0] }, new double[2]);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedStatFormats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS, miae.getSpecifier());
        }
    }

    @Test
    public void testFewPoints() {
        SimpleRegression sr = new SimpleRegression();
        sr.addObservations(new double[][] { new double[] { 1.0, 1.5 }}, new double[] { 1.0 });
        Assertions.assertEquals(1, sr.getN());
        Assertions.assertTrue(Double.isNaN(sr.getXSumSquares()));
        sr.addObservations(new double[][] { new double[] { 1.0, 1.5 }}, new double[] { 1.0 });
        Assertions.assertEquals(2, sr.getN());
        Assertions.assertFalse(Double.isNaN(sr.getXSumSquares()));
        Assertions.assertTrue(Double.isNaN(sr.getSlopeConfidenceInterval()));
        Assertions.assertTrue(Double.isNaN(sr.getSignificance()));
        try {
            sr.regress();
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedStatFormats.NOT_ENOUGH_DATA_REGRESSION, miae.getSpecifier());
        }
        sr.addObservations(new double[][] { new double[] { 1.0, 1.5 }}, new double[] { 1.0 });
        RegressionResults results = sr.regress();
        Assertions.assertTrue(Double.isNaN(results.getParameterEstimate(1)));
        results = sr.regress(new int[] { 1 });
        Assertions.assertEquals(1.0, results.getParameterEstimate(0), 1.0e-15);
        sr.addObservations(new double[][] { new double[] { 2.0, 2.5 }}, new double[] { 2.0 });
        results = sr.regress();
        Assertions.assertFalse(Double.isNaN(results.getParameterEstimate(1)));
        sr.addObservations(new double[][] { new double[] { Double.NaN, Double.NaN }}, new double[] { Double.NaN });
        results = sr.regress(new int[] { 1 });
        Assertions.assertTrue(Double.isNaN(results.getParameterEstimate(0)));

    }

    @Test
    public void testFewPointsWithoutIntercept() {
        SimpleRegression sr = new SimpleRegression(false);
        sr.addObservations(new double[][] { new double[] { 1.0, 1.5 }}, new double[] { 1.0 });
        Assertions.assertEquals(1, sr.getN());
        Assertions.assertTrue(Double.isNaN(sr.getXSumSquares()));
        try {
            sr.regress();
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedStatFormats.NOT_ENOUGH_DATA_REGRESSION, miae.getSpecifier());
        }
        sr.addObservations(new double[][] { new double[] { 1.0, 1.5 }}, new double[] { 1.0 });
        RegressionResults results = sr.regress();
        Assertions.assertFalse(Double.isNaN(results.getParameterEstimate(0)));
        sr.addObservations(new double[][] { new double[] { Double.NaN, 1.0 }}, new double[] { 2.0 });
        results = sr.regress();
        Assertions.assertTrue(Double.isNaN(results.getParameterEstimate(0)));
    }

}
