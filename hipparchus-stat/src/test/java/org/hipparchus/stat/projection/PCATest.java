/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.stat.projection;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.stat.LocalizedStatFormats;
import org.junit.Assert;
import org.junit.Test;

public class PCATest {

    // example from:
    // https://towardsdatascience.com/the-mathematics-behind-principal-component-analysis-fff2d7f4b643
    // https://stattrek.com/matrix-algebra/covariance-matrix
    private static final double[][] SCORES = {
            {90, 60, 90},
            {90, 90, 30},
            {60, 60, 60},
            {60, 60, 90},
            {30, 30, 30},
    };

    public static final double[] EXPECTED_MEAN = {66, 60, 60};
    public static final double[] EXPECTED_VARIANCE = {1137.587441, 786.387983, 56.024575};

    private static final double[][] EXPECTED_COMPONENTS = {
            {  0.6558023,  0.385999 },
            {  0.4291978,  0.516366 },
            {  0.6210577, -0.7644414 }
    };

    /**
     * This is the expected value (give or take sign) when centering (covariance)
     * but no scaling is applied. In general, components with the same values but
     * differing sign are
     * <a href="https://stats.stackexchange.com/questions/30348/is-it-acceptable-to-reverse-a-sign-of-a-principal-component-score">equivalent</a>.
     *
     * The result has been cross-checked with
     * <a href="https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html">sklearn.decomposition.PCA</a>
     * <a href="https://javadoc.io/static/nz.ac.waikato.cms.weka/weka-dev/3.9.4/weka/attributeSelection/PrincipalComponents.html>weka.attributeSelection.PrincipalComponents</a>
     * (with the <code>centerData</code> option set to <code>true</code>)
     * <a href="https://github.com/datumbox/datumbox-framework/blob/develop/datumbox-framework-core/src/main/java/com/datumbox/framework/core/machinelearning/featureselection/PCA.java">com.datumbox.framework.core.machinelearning.featureselection.PCA</a>
     * (but for whatever reason datumbox does the transform on the unnormalized original data rather than normalized data - normalizing manually achieves the result below)
     * <a href="https://www.javadoc.io/doc/com.github.haifengl/smile-core/latest/smile/feature/extraction/PCA.html">smile.feature.extraction.PCA</a>
     * (using the PCA.fit method)
     * <a href="https://au.mathworks.com/help/stats/pca.html">pca from matlab</a>
     */
    private static final double[][] EXPECTED_COV = {
            {  34.3709848, -13.6692708 },
            {   9.9834573,  47.6882055 },
            {  -3.9348135,  -2.3159927 },
            {  14.6969171, -25.2492347 },
            { -55.1165457,  -6.4537072 },
    };

    /**
     * This is the expected value give or take sign when centering and scaling (correlation).
     *
     * The result has been cross-checked with
     * <a href="https://javadoc.io/static/nz.ac.waikato.cms.weka/weka-dev/3.9.4/weka/attributeSelection/PrincipalComponents.html></a>
     * <a href="https://au.mathworks.com/help/stats/pca.html">pca from matlab</a>
     * (using the 'VariableWeights','variance' option)
     */
    private static final double[][] EXPECTED_COR = {
            {  0.9118256, -0.942809  },
            {  1.3832302,  1.4142136 },
            { -0.1690309,  0.0       },
            {  0.0666714, -0.942809  },
            { -2.1926964,  0.4714045 },
    };

    /**
     * This is the expected value give or take sign when centering and scaling (correlation)
     * is applied with no bias adjustment. In general, bias adjustment is more accurate but
     * alters most PCA machine learning models by an insignificant amount so is often not accounted for.
     *
     * The result has been cross-checked with
     * <a href="https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html">sklearn.decomposition.PCA</a>
     * (but with prior preprocessing using <a href="https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html">sklearn.preprocessing.StandardScaler</a>)
     * <a href="https://www.javadoc.io/doc/com.github.haifengl/smile-core/latest/smile/feature/extraction/PCA.html">smile.feature.extraction.PCA</a>
     * (using the PCA.cor method)
     */
    public static final double[][] EXPECTED_COR_NO_BIAS = {
            {  1.0194521, -1.0540926 },
            {  1.5464984,  1.5811388 },
            { -0.1889822,  0.0       },
            {  0.0745409, -1.0540926 },
            { -2.451509,   0.5270463 },
    };
    public static final double DELTA = 0.000001;

    @Test
    public void defaultSettings() {
        PCA pca = new PCA(2);
        Assert.assertEquals(2, pca.getNumComponents());
        Assert.assertFalse(pca.isScale());
        Assert.assertTrue(pca.isBiasCorrection());
    }

    @Test
    public void covariance() {
        PCA pca = new PCA(2);
        double[][] actual = pca.fitAndTransform(SCORES);
        Assert.assertArrayEquals(EXPECTED_MEAN, pca.getCenter(), DELTA);
        Assert.assertArrayEquals(EXPECTED_VARIANCE, pca.getVariance(), DELTA);
        assertExpected(EXPECTED_COMPONENTS, pca.getComponents());
        assertExpected(EXPECTED_COV, actual);

        // calling fit and transform individually should be the same as combo method
        pca = new PCA(2);
        actual = pca.fit(SCORES).transform(SCORES);
        Assert.assertArrayEquals(EXPECTED_MEAN, pca.getCenter(), DELTA);
        assertExpected(EXPECTED_COV, actual);
    }

    @Test
    public void correlation() {
        PCA pca = new PCA(2, true, true);
        double[][] actual = pca.fitAndTransform(SCORES);
        assertExpected(EXPECTED_COR, actual);
    }

    @Test
    public void correlationNoBias() {
        PCA pca = new PCA(2, true, false);
        double[][] actual = pca.fitAndTransform(SCORES);
        assertExpected(EXPECTED_COR_NO_BIAS, actual);
    }

    @Test
    public void transformWithoutFit() {
        PCA pca = new PCA(2);
        try {
            pca.transform(SCORES);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedStatFormats.ILLEGAL_STATE_PCA, mise.getSpecifier());
            Assert.assertEquals("transform", mise.getParts()[0]);
        }
    }

    private static void assertExpected(double[][] expected, double[][] actual) {
        for (int i = 0; i < expected.length; i++) {
            double[] e = expected[i];
            double[] t = actual[i];
            Assert.assertArrayEquals(e, t, DELTA);
        }
    }
}
