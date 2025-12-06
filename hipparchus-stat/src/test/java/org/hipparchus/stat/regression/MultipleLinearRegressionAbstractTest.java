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
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;


public abstract class MultipleLinearRegressionAbstractTest {

    protected AbstractMultipleLinearRegression regression;

    @BeforeEach
    public void setUp(){
        regression = createRegression();
    }

    protected abstract AbstractMultipleLinearRegression createRegression();

    protected abstract int getNumberOfRegressors();

    protected abstract int getSampleSize();

    @Test
    public void canEstimateRegressionParameters(){
        double[] beta = regression.estimateRegressionParameters();
        assertEquals(getNumberOfRegressors(), beta.length);
    }

    @Test
    public void canEstimateResiduals(){
        double[] e = regression.estimateResiduals();
        assertEquals(getSampleSize(), e.length);
    }

    @Test
    public void canEstimateRegressionParametersVariance(){
        double[][] variance = regression.estimateRegressionParametersVariance();
        assertEquals(getNumberOfRegressors(), variance.length);
    }

    @Test
    public void canEstimateRegressandVariance(){
        if (getSampleSize() > getNumberOfRegressors()) {
            double variance = regression.estimateRegressandVariance();
            assertTrue(variance > 0.0);
        }
    }

    /**
     * Verifies that newSampleData methods consistently insert unitary columns
     * in design matrix.  Confirms the fix for MATH-411.
     */
    @Test
    public void testNewSample() {
        double[] design = new double[] {
          1, 19, 22, 33,
          2, 20, 30, 40,
          3, 25, 35, 45,
          4, 27, 37, 47
        };
        double[] y = new double[] {1, 2, 3, 4};
        double[][] x = new double[][] {
          {19, 22, 33},
          {20, 30, 40},
          {25, 35, 45},
          {27, 37, 47}
        };
        AbstractMultipleLinearRegression regression = createRegression();
        regression.newSampleData(design, 4, 3);
        RealMatrix flatX = regression.getX().copy();
        RealVector flatY = regression.getY().copy();
        regression.newXSampleData(x);
        regression.newYSampleData(y);
        assertEquals(flatX, regression.getX());
        assertEquals(flatY, regression.getY());

        // No intercept
        regression.setNoIntercept(true);
        regression.newSampleData(design, 4, 3);
        flatX = regression.getX().copy();
        flatY = regression.getY().copy();
        regression.newXSampleData(x);
        regression.newYSampleData(y);
        assertEquals(flatX, regression.getX());
        assertEquals(flatY, regression.getY());
    }

    @Test
    public void testNewSampleNullData() {
        assertThrows(NullArgumentException.class, () -> {
            double[] data = null;
            createRegression().newSampleData(data, 2, 3);
        });
    }

    @Test
    public void testNewSampleInvalidData() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            double[] data = new double[]{1, 2, 3, 4};
            createRegression().newSampleData(data, 2, 3);
        });
    }

    @Test
    public void testNewSampleInsufficientData() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            double[] data = new double[]{1, 2, 3, 4};
            createRegression().newSampleData(data, 1, 3);
        });
    }

    @Test
    public void testXSampleDataNull() {
        assertThrows(NullArgumentException.class, () -> createRegression().newXSampleData(null));
    }

    @Test
    public void testYSampleDataNull() {
        assertThrows(NullArgumentException.class, () -> createRegression().newYSampleData(null));
    }

}
