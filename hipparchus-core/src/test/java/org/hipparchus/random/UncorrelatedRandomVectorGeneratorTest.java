//Licensed to the Apache Software Foundation (ASF) under one
//or more contributor license agreements.  See the NOTICE file
//distributed with this work for additional information
//regarding copyright ownership.  The ASF licenses this file
//to you under the Apache License, Version 2.0 (the
//"License"); you may not use this file except in compliance
//with the License.  You may obtain a copy of the License at

//http://www.apache.org/licenses/LICENSE-2.0

//Unless required by applicable law or agreed to in writing,
//software distributed under the License is distributed on an
//"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
//KIND, either express or implied.  See the License for the
//specific language governing permissions and limitations
//under the License.

package org.hipparchus.random;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;

public class UncorrelatedRandomVectorGeneratorTest {
    private double[] mean;
    private double[] standardDeviation;
    private UncorrelatedRandomVectorGenerator generator;

    public UncorrelatedRandomVectorGeneratorTest() {
        mean              = new double[] {0.0, 1.0, -3.0, 2.3};
        standardDeviation = new double[] {1.0, 2.0, 10.0, 0.1};
        RandomGenerator rg = new JDKRandomGenerator();
        rg.setSeed(17399225432l);
        generator =
            new UncorrelatedRandomVectorGenerator(mean, standardDeviation,
                    new GaussianRandomGenerator(rg));
    }

    @Test
    public void testMeanAndCorrelation() {
        final int n = generator.nextVector().length;
        final double[] estimatedMean = new double[generator.nextVector().length];
        final RealMatrix matrix = new Array2DRowRealMatrix(10000, n);
        for (int i = 0; i < 10000; ++i) {
            double[] v = generator.nextVector();
            matrix.setRow(i, v);
        }

        for (int i = 0; i < n; i++) {
            estimatedMean[i] = UnitTestUtils.mean(matrix.getColumn(i));
        }

        // double[] estimatedMean = meanStat.getResult();
        double scale;
        RealMatrix estimatedCorrelation = UnitTestUtils.covarianceMatrix(matrix);
        //RealMatrix estimatedCorrelation = covStat.getResult();
        for (int i = 0; i < estimatedMean.length; ++i) {
            Assert.assertEquals(mean[i], estimatedMean[i], 0.07);
            for (int j = 0; j < i; ++j) {
                scale = standardDeviation[i] * standardDeviation[j];
                Assert.assertEquals(0, estimatedCorrelation.getEntry(i, j) / scale, 0.03);
            }
            scale = standardDeviation[i] * standardDeviation[i];
            Assert.assertEquals(1, estimatedCorrelation.getEntry(i, i) / scale, 0.025);
        }
    }
}
