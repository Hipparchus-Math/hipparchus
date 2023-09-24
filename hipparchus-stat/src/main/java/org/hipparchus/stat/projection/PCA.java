/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.stat.projection;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.linear.EigenDecompositionSymmetric;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.stat.LocalizedStatFormats;
import org.hipparchus.stat.StatUtils;
import org.hipparchus.stat.correlation.Covariance;
import org.hipparchus.stat.descriptive.moment.StandardDeviation;

/**
 * Principal component analysis (PCA) is a statistical technique for reducing the dimensionality of a dataset.
 * <a href="https://en.wikipedia.org/wiki/Principal_component_analysis">PCA</a> can be thought of as a
 * projection or scaling of the data to reduce the number of dimensions but done in a way
 * that preserves as much information as possible.
 * @since 3.0
 */
public class PCA {
    /**
     * The number of components (reduced dimensions) for this projection.
     */
    private final int numC;

    /**
     * Whether to scale (standardize) the input data as well as center (normalize).
     */
    private final boolean scale;

    /**
     * Whether to correct for bias when standardizing. Ignored when only centering.
     */
    private final boolean biasCorrection;

    /**
     * The by column (feature) averages (means) from the fitted data.
     */
    private double[] center;

    /**
     * The by column (feature) standard deviates from the fitted data.
     */
    private double[] std;

    /**
     * The eigenValues (variance) of our projection model.
     */
    private double[] eigenValues;

    /**
     * The eigenVectors (components) of our projection model.
     */
    private RealMatrix principalComponents;

    /**
     * Utility class when scaling.
     */
    private final StandardDeviation sd;

    /**
     * Create a PCA with the ability to adjust scaling parameters.
     *
     * @param numC the number of components
     * @param scale whether to also scale (correlation) rather than just center (covariance)
     * @param biasCorrection whether to adjust for bias when scaling
     */
    public PCA(int numC, boolean scale, boolean biasCorrection) {
        this.numC = numC;
        this.scale = scale;
        this.biasCorrection = biasCorrection;
        sd = scale ? new StandardDeviation(biasCorrection) : null;
    }

    /**
     * A default PCA will center but not scale.
     *
     * @param numC the number of components
     */
    public PCA(int numC) {
        this(numC, false, true);
    }

    /** GEt number of components.
     * @return the number of components
     */
    public int getNumComponents() {
        return numC;
    }

    /** Check whether scaling (correlation) or no scaling (covariance) is used.
     * @return whether scaling (correlation) or no scaling (covariance) is used
     */
    public boolean isScale() {
        return scale;
    }

    /** Check whether scaling (correlation), if in use, adjusts for bias.
     * @return whether scaling (correlation), if in use, adjusts for bias
     */
    public boolean isBiasCorrection() {
        return biasCorrection;
    }

    /** Get principal component variances.
     * @return the principal component variances, ordered from largest to smallest, which are the eigenvalues of the covariance or correlation matrix of the fitted data
     */
    public double[] getVariance() {
        validateState("getVariance");
        return eigenValues.clone();
    }

    /** Get by column center (or mean) of the fitted data.
     * @return the by column center (or mean) of the fitted data
     */
    public double[] getCenter() {
        validateState("getCenter");
        return center.clone();
    }

    /**
     * Returns the principal components of our projection model.
     * These are the eigenvectors of our covariance/correlation matrix.
     *
     * @return the principal components
     */
    public double[][] getComponents() {
        validateState("getComponents");
        return principalComponents.getData();
    }

    /**
     * Fit our model to the data and then transform it to the reduced dimensions.
     *
     * @param data the input data
     * @return the fitted data
     */
    public double[][] fitAndTransform(double[][] data) {
        center = null;
        RealMatrix normalizedM = getNormalizedMatrix(data);
        calculatePrincipalComponents(normalizedM);
        return normalizedM.multiply(principalComponents).getData();
    }

    /**
     * Transform the supplied data using our projection model.
     *
     * @param data the input data
     * @return the fitted data
     */
    public double[][] transform(double[][] data) {
        validateState("transform");
        RealMatrix normalizedM = getNormalizedMatrix(data);
        return normalizedM.multiply(principalComponents).getData();
    }

    /**
     * Fit our model to the data, ready for subsequence transforms.
     *
     * @param data the input data
     * @return this
     */
    public PCA fit(double[][] data) {
        center = null;
        RealMatrix normalized = getNormalizedMatrix(data);
        calculatePrincipalComponents(normalized);
        return this;
    }

    /** Check if the state allows an operation to be performed.
     * @param from name of the operation
     * @exception MathIllegalStateException if the state does not allows operation
     */
    private void validateState(String from) {
        if (center == null) {
            throw new MathIllegalStateException(LocalizedStatFormats.ILLEGAL_STATE_PCA, from);
        }

    }

    /** Compute eigenvalues and principal components.
     * <p>
     * The results are stored in the instance itself
     * <p>
     * @param normalizedM normalized matrix
     */
    private void calculatePrincipalComponents(RealMatrix normalizedM) {
        RealMatrix covarianceM = new Covariance(normalizedM).getCovarianceMatrix();
        EigenDecompositionSymmetric decomposition = new EigenDecompositionSymmetric(covarianceM);
        eigenValues = decomposition.getEigenvalues();
        principalComponents = MatrixUtils.createRealMatrix(eigenValues.length, numC);
        for (int c = 0; c < numC; c++) {
            for (int f = 0; f < eigenValues.length; f++) {
                principalComponents.setEntry(f, c, decomposition.getEigenvector(c).getEntry(f));
            }
        }
    }

    /**
     * This will either normalize (center) or
     * standardize (center plus scale) the input data.
     *
     * @param input the input data
     * @return the normalized (or standardized) matrix
     */
    private RealMatrix getNormalizedMatrix(double[][] input) {
        int numS = input.length;
        int numF = input[0].length;
        boolean calculating = center == null;
        if (calculating) {
            center = new double[numF];
            if (scale) {
                std = new double[numF];
            }
        }

        double[][] normalized = new double[numS][numF];
        for (int f = 0; f < numF; f++) {
            if (calculating) {
                calculateNormalizeParameters(input, numS, f);
            }
            for (int s = 0; s < numS; s++) {
                normalized[s][f] = input[s][f] - center[f];
            }
            if (scale) {
                for (int s = 0; s < numS; s++) {
                    normalized[s][f] /= std[f];
                }
            }
        }

        return MatrixUtils.createRealMatrix(normalized);
    }

    /** compute normalized parameters.
     * @param input the input data
     * @param numS number of data rows
     * @param f index of the component
     */
    private void calculateNormalizeParameters(double[][] input, int numS, int f) {
        double[] column = new double[numS];
        for (int s = 0; s < numS; s++) {
            column[s] = input[s][f];
        }
        center[f] = StatUtils.mean(column);
        if (scale) {
            std[f] = sd.evaluate(column, center[f]);
        }
    }
}
