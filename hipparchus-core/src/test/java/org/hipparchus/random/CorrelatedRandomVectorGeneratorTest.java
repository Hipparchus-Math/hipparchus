//Licensed to the Apache Software Foundation (ASF) under one
//or more contributor license agreements.  See the NOTICE file
//distributed with this work for additional information
//regarding copyright ownership.  The ASF licenses this file
//to you under the Apache License, Version 2.0 (the
//"License"); you may not use this file except in compliance
//with the License.  You may obtain a copy of the License at

//https://www.apache.org/licenses/LICENSE-2.0

//Unless required by applicable law or agreed to in writing,
//software distributed under the License is distributed on an
//"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
//KIND, either express or implied.  See the License for the
//specific language governing permissions and limitations
//under the License.

package org.hipparchus.random;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

public class CorrelatedRandomVectorGeneratorTest {
    private double[] mean;
    private RealMatrix covariance;
    private CorrelatedRandomVectorGenerator generator;

    public CorrelatedRandomVectorGeneratorTest() {
        mean = new double[] { 0.0, 1.0, -3.0, 2.3 };

        RealMatrix b = MatrixUtils.createRealMatrix(4, 3);
        int counter = 0;
        for (int i = 0; i < b.getRowDimension(); ++i) {
            for (int j = 0; j < b.getColumnDimension(); ++j) {
                b.setEntry(i, j, 1.0 + 0.1 * ++counter);
            }
        }
        RealMatrix bbt = b.multiplyTransposed(b);
        covariance = MatrixUtils.createRealMatrix(mean.length, mean.length);
        for (int i = 0; i < covariance.getRowDimension(); ++i) {
            covariance.setEntry(i, i, bbt.getEntry(i, i));
            for (int j = 0; j < covariance.getColumnDimension(); ++j) {
                double s = bbt.getEntry(i, j);
                covariance.setEntry(i, j, s);
                covariance.setEntry(j, i, s);
            }
        }

        RandomGenerator rg = new JDKRandomGenerator();
        rg.setSeed(17399225432l);
        GaussianRandomGenerator rawGenerator = new GaussianRandomGenerator(rg);
        generator = new CorrelatedRandomVectorGenerator(mean,
                                                        covariance,
                                                        1.0e-12 * covariance.getNorm1(),
                                                        rawGenerator);
    }

    @Test
    public void testRank() {
        Assertions.assertEquals(2, generator.getRank());
    }

    @Test
    public void testMath226() {
        double[] mean = { 1, 1, 10, 1 };
        double[][] cov = {
                { 1, 3, 2, 6 },
                { 3, 13, 16, 2 },
                { 2, 16, 38, -1 },
                { 6, 2, -1, 197 }
        };
        RealMatrix covRM = MatrixUtils.createRealMatrix(cov);
        JDKRandomGenerator jg = new JDKRandomGenerator();
        jg.setSeed(5322145245211l);
        NormalizedRandomGenerator rg = new GaussianRandomGenerator(jg);
        CorrelatedRandomVectorGenerator sg =
            new CorrelatedRandomVectorGenerator(mean, covRM, 0.00001, rg);

        double[] min = new double[mean.length];
        Arrays.fill(min, Double.POSITIVE_INFINITY);
        double[] max = new double[mean.length];
        Arrays.fill(max, Double.NEGATIVE_INFINITY);
        for (int i = 0; i < 10; i++) {
            double[] generated = sg.nextVector();
            for (int j = 0; j < generated.length; ++j) {
                min[j] = FastMath.min(min[j], generated[j]);
                max[j] = FastMath.max(max[j], generated[j]);
            }
        }
        for (int j = 0; j < min.length; ++j) {
            Assertions.assertTrue(max[j] - min[j] > 2.0);
        }

    }

    @Test
    public void testRootMatrix() {
        RealMatrix b = generator.getRootMatrix();
        RealMatrix bbt = b.multiplyTransposed(b);
        for (int i = 0; i < covariance.getRowDimension(); ++i) {
            for (int j = 0; j < covariance.getColumnDimension(); ++j) {
                Assertions.assertEquals(covariance.getEntry(i, j), bbt.getEntry(i, j), 1.0e-12);
            }
        }
    }

    @Test
    public void testMeanAndCovariance() {

        final double[] meanStat = new double[mean.length];
        final RealMatrix matrix = new Array2DRowRealMatrix(5000, mean.length);
        for (int i = 0; i < 5000; ++i) {
            double[] v = generator.nextVector();
            matrix.setRow(i, v);
        }

        for (int i = 0; i < mean.length; i++) {
            meanStat[i] = UnitTestUtils.mean(matrix.getColumn(i));
        }

        RealMatrix estimatedCovariance = UnitTestUtils.covarianceMatrix(matrix);
        for (int i = 0; i < meanStat.length; ++i) {
            Assertions.assertEquals(mean[i], meanStat[i], 0.07);
            for (int j = 0; j <= i; ++j) {
                Assertions.assertEquals(covariance.getEntry(i, j),
                                    estimatedCovariance.getEntry(i, j),
                                    0.1 * (1.0 + FastMath.abs(mean[i])) * (1.0 + FastMath.abs(mean[j])));
            }
        }

    }

    @Test
    public void testSampleWithZeroCovariance() {
        final double[][] covMatrix1 = new double[][]{
                {0.013445532, 0.010394690, 0.009881156, 0.010499559},
                {0.010394690, 0.023006616, 0.008196856, 0.010732709},
                {0.009881156, 0.008196856, 0.019023866, 0.009210099},
                {0.010499559, 0.010732709, 0.009210099, 0.019107243}
        };

        final double[][] covMatrix2 = new double[][]{
                {0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 0.013445532, 0.010394690, 0.009881156, 0.010499559},
                {0.0, 0.010394690, 0.023006616, 0.008196856, 0.010732709},
                {0.0, 0.009881156, 0.008196856, 0.019023866, 0.009210099},
                {0.0, 0.010499559, 0.010732709, 0.009210099, 0.019107243}
        };

        final double[][] covMatrix3 = new double[][]{
                {0.013445532, 0.010394690, 0.0, 0.009881156, 0.010499559},
                {0.010394690, 0.023006616, 0.0, 0.008196856, 0.010732709},
                {0.0, 0.0, 0.0, 0.0, 0.0},
                {0.009881156, 0.008196856, 0.0, 0.019023866, 0.009210099},
                {0.010499559, 0.010732709, 0.0, 0.009210099, 0.019107243}
        };

        testSampler(covMatrix1, 10000, 0.001);
        testSampler(covMatrix2, 10000, 0.001);
        testSampler(covMatrix3, 10000, 0.001);

    }

    private CorrelatedRandomVectorGenerator createSampler(double[][] cov) {
        RealMatrix matrix = new Array2DRowRealMatrix(cov);
        double small = 10e-12 * matrix.getNorm1();
        return new CorrelatedRandomVectorGenerator(
                new double[cov.length],
                matrix,
                small,
                new GaussianRandomGenerator(new Well1024a(0x366a26b94e520f41l)));
    }

    private void testSampler(final double[][] covMatrix, int samples, double epsilon) {
        CorrelatedRandomVectorGenerator sampler = createSampler(covMatrix);
        RealMatrix matrix = new Array2DRowRealMatrix(samples, covMatrix.length);

        for (int i = 0; i < samples; ++i) {
            matrix.setRow(i, sampler.nextVector());
        }

        RealMatrix sampleCov = UnitTestUtils.covarianceMatrix(matrix);
        for (int r = 0; r < covMatrix.length; ++r) {
            UnitTestUtils.assertEquals(covMatrix[r], sampleCov.getColumn(r), epsilon);
        }

    }

}
