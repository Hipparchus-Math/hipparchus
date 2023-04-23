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
package org.hipparchus.migration;

import java.util.Locale;
import java.util.MissingResourceException;
import java.util.ResourceBundle;

import org.hipparchus.exception.Localizable;
import org.hipparchus.exception.UTF8Control;

/**
 * Enumeration for localized messages formats used in exceptions messages.
 * <p>
 * The constants in this enumeration represent the available
 * formats as localized strings. These formats are intended to be
 * localized using simple properties files, using the constant
 * name as the key and the property value as the message format.
 * The source English format is provided in the constants themselves
 * to serve both as a reminder for developers to understand the parameters
 * needed by each format, as a basis for translators to create
 * localized properties files, and as a default format if some
 * translation is missing.
 * </p>
 * @deprecated these formats are not used at all, they are intended only
 * as a migration help from Apache Commons Math
 */
@Deprecated
public enum LocalizedMigrationFormats implements Localizable {

    /** ARGUMENT_OUTSIDE_DOMAIN. */
    ARGUMENT_OUTSIDE_DOMAIN("Argument {0} outside domain [{1} ; {2}]"),

    /** ASSYMETRIC_EIGEN_NOT_SUPPORTED. */
    ASSYMETRIC_EIGEN_NOT_SUPPORTED("eigen decomposition of assymetric matrices not supported yet"),

    /** BOBYQA_BOUND_DIFFERENCE_CONDITION. */
    BOBYQA_BOUND_DIFFERENCE_CONDITION("the difference between the upper and lower bound must be larger than twice the initial trust region radius ({0})"),

    /** CANNOT_CLEAR_STATISTIC_CONSTRUCTED_FROM_EXTERNAL_MOMENTS. */
    CANNOT_CLEAR_STATISTIC_CONSTRUCTED_FROM_EXTERNAL_MOMENTS("statistics constructed from external moments cannot be cleared"),

    /** CANNOT_FORMAT_INSTANCE_AS_3D_VECTOR. */
    CANNOT_FORMAT_INSTANCE_AS_3D_VECTOR("cannot format a {0} instance as a 3D vector"),

    /** CANNOT_FORMAT_INSTANCE_AS_REAL_VECTOR. */
    CANNOT_FORMAT_INSTANCE_AS_REAL_VECTOR("cannot format a {0} instance as a real vector"),

    /** CANNOT_INCREMENT_STATISTIC_CONSTRUCTED_FROM_EXTERNAL_MOMENTS. */
    CANNOT_INCREMENT_STATISTIC_CONSTRUCTED_FROM_EXTERNAL_MOMENTS("statistics constructed from external moments cannot be incremented"),

    /** CANNOT_RETRIEVE_AT_NEGATIVE_INDEX. */
    CANNOT_RETRIEVE_AT_NEGATIVE_INDEX("elements cannot be retrieved from a negative array index {0}"),

    /** CANNOT_SET_AT_NEGATIVE_INDEX. */
    CANNOT_SET_AT_NEGATIVE_INDEX("cannot set an element at a negative index {0}"),

    /** CANNOT_TRANSFORM_TO_DOUBLE. */
    CANNOT_TRANSFORM_TO_DOUBLE("Conversion Exception in Transformation: {0}"),

    /** CARDAN_ANGLES_SINGULARITY. */
    CARDAN_ANGLES_SINGULARITY("Cardan angles singularity"),

    /** CLASS_DOESNT_IMPLEMENT_COMPARABLE. */
    CLASS_DOESNT_IMPLEMENT_COMPARABLE("class ({0}) does not implement Comparable"),

    /** COLUMN_INDEX_OUT_OF_RANGE. */
    COLUMN_INDEX_OUT_OF_RANGE("column index {0} out of allowed range [{1}, {2}]"),

    /** DIMENSIONS_MISMATCH_SIMPLE. */
    DIMENSIONS_MISMATCH_SIMPLE("{0} != {1}"),

    /** EMPTY_STRING_FOR_IMAGINARY_CHARACTER. */
    EMPTY_STRING_FOR_IMAGINARY_CHARACTER("empty string for imaginary character"),

    /** EULER_ANGLES_SINGULARITY. */
    EULER_ANGLES_SINGULARITY("Euler angles singularity"),

    /** FUNCTION_NOT_DIFFERENTIABLE. */
    FUNCTION_NOT_DIFFERENTIABLE("function is not differentiable"),

    /** FUNCTION_NOT_POLYNOMIAL. */
    FUNCTION_NOT_POLYNOMIAL("function is not polynomial"),

    /** INDEX_OUT_OF_RANGE. */
    INDEX_OUT_OF_RANGE("index {0} out of allowed range [{1}, {2}]"),

    /** INVALID_BRACKETING_PARAMETERS. */
    INVALID_BRACKETING_PARAMETERS("invalid bracketing parameters:  lower bound={0},  initial={1}, upper bound={2}"),

    /** INVALID_INTERVAL_INITIAL_VALUE_PARAMETERS. */
    INVALID_INTERVAL_INITIAL_VALUE_PARAMETERS("invalid interval, initial value parameters:  lower={0}, initial={1}, upper={2}"),

    /** INVALID_ITERATIONS_LIMITS. */
    INVALID_ITERATIONS_LIMITS("invalid iteration limits: min={0}, max={1}"),

    /** INVALID_REGRESSION_ARRAY. */
    INVALID_REGRESSION_ARRAY("input data array length = {0} does not match the number of observations = {1} and the number of regressors = {2}"),

    /** ITERATOR_EXHAUSTED. */
    ITERATOR_EXHAUSTED("iterator exhausted"),

    /** LOESS_EXPECTS_AT_LEAST_ONE_POINT. */
    LOESS_EXPECTS_AT_LEAST_ONE_POINT("Loess expects at least 1 point"),

    /** MAP_MODIFIED_WHILE_ITERATING. */
    MAP_MODIFIED_WHILE_ITERATING("map has been modified while iterating"),

    /** MAX_ITERATIONS_EXCEEDED. */
    MAX_ITERATIONS_EXCEEDED("maximal number of iterations ({0}) exceeded"),

    /** MISMATCHED_LOESS_ABSCISSA_ORDINATE_ARRAYS. */
    MISMATCHED_LOESS_ABSCISSA_ORDINATE_ARRAYS("Loess expects the abscissa and ordinate arrays to be of the same size, but got {0} abscissae and {1} ordinatae"),

    /** NEGATIVE_BRIGHTNESS_EXPONENT. */
    NEGATIVE_BRIGHTNESS_EXPONENT("brightness exponent should be positive or null, but got {0}"),

    /** NEGATIVE_ELEMENT_AT_2D_INDEX. */
    NEGATIVE_ELEMENT_AT_2D_INDEX("element ({0}, {1}) is negative: {2}"),

    /** NEGATIVE_NUMBER_OF_SUCCESSES. */
    NEGATIVE_NUMBER_OF_SUCCESSES("number of successes must be non-negative ({0})"),

    /** NEGATIVE_NUMBER_OF_TRIALS. */
    NEGATIVE_NUMBER_OF_TRIALS("number of trials must be non-negative ({0})"),

    /** NON_POSITIVE_DEFINITE_MATRIX. */
    NON_POSITIVE_DEFINITE_MATRIX("not positive definite matrix: diagonal element at ({1},{1}) is smaller than {2} ({0})"),

    /** NON_POSITIVE_MICROSPHERE_ELEMENTS. */
    NON_POSITIVE_MICROSPHERE_ELEMENTS("number of microsphere elements must be positive, but got {0}"),

    /** NON_POSITIVE_POLYNOMIAL_DEGREE. */
    NON_POSITIVE_POLYNOMIAL_DEGREE("polynomial degree must be positive: degree={0}"),

    /** NON_REAL_FINITE_ABSCISSA. */
    NON_REAL_FINITE_ABSCISSA("all abscissae must be finite real numbers, but {0}-th is {1}"),

    /** NON_REAL_FINITE_ORDINATE. */
    NON_REAL_FINITE_ORDINATE("all ordinatae must be finite real numbers, but {0}-th is {1}"),

    /** NON_REAL_FINITE_WEIGHT. */
    NON_REAL_FINITE_WEIGHT("all weights must be finite real numbers, but {0}-th is {1}"),

    /** NOT_ADDITION_COMPATIBLE_MATRICES. */
    NOT_ADDITION_COMPATIBLE_MATRICES("{0}x{1} and {2}x{3} matrices are not addition compatible"),

    /** NOT_DECREASING_NUMBER_OF_POINTS. */
    NOT_DECREASING_NUMBER_OF_POINTS("points {0} and {1} are not decreasing ({2} < {3})"),

    /** NOT_INCREASING_NUMBER_OF_POINTS. */
    NOT_INCREASING_NUMBER_OF_POINTS("points {0} and {1} are not increasing ({2} > {3})"),

    /** NOT_MULTIPLICATION_COMPATIBLE_MATRICES. */
    NOT_MULTIPLICATION_COMPATIBLE_MATRICES("{0}x{1} and {2}x{3} matrices are not multiplication compatible"),

    /** NOT_POSITIVE_DEGREES_OF_FREEDOM. */
    NOT_POSITIVE_DEGREES_OF_FREEDOM("degrees of freedom must be positive ({0})"),

    /** NOT_POSITIVE_ELEMENT_AT_INDEX. */
    NOT_POSITIVE_ELEMENT_AT_INDEX("element {0} is not positive: {1}"),

    /** NOT_POSITIVE_LENGTH. */
    NOT_POSITIVE_LENGTH("length must be positive ({0})"),

    /** NOT_POSITIVE_MEAN. */
    NOT_POSITIVE_MEAN("mean must be positive ({0})"),

    /** NOT_POSITIVE_PERMUTATION. */
    NOT_POSITIVE_PERMUTATION("permutation k ({0}) must be positive"),

    /** NOT_POSITIVE_POISSON_MEAN. */
    NOT_POSITIVE_POISSON_MEAN("the Poisson mean must be positive ({0})"),

    /** NOT_POSITIVE_POPULATION_SIZE. */
    NOT_POSITIVE_POPULATION_SIZE("population size must be positive ({0})"),

    /** NOT_POSITIVE_ROW_DIMENSION. */
    NOT_POSITIVE_ROW_DIMENSION("invalid row dimension: {0} (must be positive)"),

    /** NOT_POSITIVE_SAMPLE_SIZE. */
    NOT_POSITIVE_SAMPLE_SIZE("sample size must be positive ({0})"),

    /** NOT_POSITIVE_SHAPE. */
    NOT_POSITIVE_SHAPE("shape must be positive ({0})"),

    /** NOT_POSITIVE_STANDARD_DEVIATION. */
    NOT_POSITIVE_STANDARD_DEVIATION("standard deviation must be positive ({0})"),

    /** NOT_POSITIVE_UPPER_BOUND. */
    NOT_POSITIVE_UPPER_BOUND("upper bound must be positive ({0})"),

    /** NOT_STRICTLY_DECREASING_NUMBER_OF_POINTS. */
    NOT_STRICTLY_DECREASING_NUMBER_OF_POINTS("points {0} and {1} are not strictly decreasing ({2} <= {3})"),

    /** NOT_STRICTLY_INCREASING_KNOT_VALUES. */
    NOT_STRICTLY_INCREASING_KNOT_VALUES("knot values must be strictly increasing"),

    /** NOT_STRICTLY_INCREASING_NUMBER_OF_POINTS. */
    NOT_STRICTLY_INCREASING_NUMBER_OF_POINTS("points {0} and {1} are not strictly increasing ({2} >= {3})"),

    /** NOT_SUBTRACTION_COMPATIBLE_MATRICES. */
    NOT_SUBTRACTION_COMPATIBLE_MATRICES("{0}x{1} and {2}x{3} matrices are not subtraction compatible"),

    /** NOT_SYMMETRIC_MATRIX. */
    NOT_SYMMETRIC_MATRIX("not symmetric matrix"),

    /** NO_BIN_SELECTED. */
    NO_BIN_SELECTED("no bin selected"),

    /** NO_DEGREES_OF_FREEDOM. */
    NO_DEGREES_OF_FREEDOM("no degrees of freedom ({0} measurements, {1} parameters)"),

    /** NO_DENSITY_FOR_THIS_DISTRIBUTION. */
    NO_DENSITY_FOR_THIS_DISTRIBUTION("This distribution does not have a density function implemented"),

    /** NO_RESULT_AVAILABLE. */
    NO_RESULT_AVAILABLE("no result available"),

    /** NO_SUCH_MATRIX_ENTRY. */
    NO_SUCH_MATRIX_ENTRY("no entry at indices ({0}, {1}) in a {2}x{3} matrix"),

    /** N_POINTS_GAUSS_LEGENDRE_INTEGRATOR_NOT_SUPPORTED. */
    N_POINTS_GAUSS_LEGENDRE_INTEGRATOR_NOT_SUPPORTED("{0} points Legendre-Gauss integrator not supported, number of points must be in the {1}-{2} range"),

    /** OBJECT_TRANSFORMATION. */
    OBJECT_TRANSFORMATION("conversion exception in transformation"),

    /** OBSERVED_COUNTS_ALL_ZERO. */
    OBSERVED_COUNTS_ALL_ZERO("observed counts are all 0 in observed array {0}"),

    /** OUT_OF_ORDER_ABSCISSA_ARRAY. */
    OUT_OF_ORDER_ABSCISSA_ARRAY("the abscissae array must be sorted in a strictly increasing order, but the {0}-th element is {1} whereas {2}-th is {3}"),

    /** OUT_OF_RANGE_RIGHT. */
    OUT_OF_RANGE_RIGHT("{0} out of [{1}, {2}) range"),

    /** POLYNOMIAL_INTERPOLANTS_MISMATCH_SEGMENTS. */
    POLYNOMIAL_INTERPOLANTS_MISMATCH_SEGMENTS("number of polynomial interpolants must match the number of segments ({0} != {1} - 1)"),

    /** POWER_NEGATIVE_PARAMETERS. */
    POWER_NEGATIVE_PARAMETERS("cannot raise an integral value to a negative power ({0}^{1})"),

    /** ROW_INDEX_OUT_OF_RANGE. */
    ROW_INDEX_OUT_OF_RANGE("row index {0} out of allowed range [{1}, {2}]"),

    /** SAME_SIGN_AT_ENDPOINTS. */
    SAME_SIGN_AT_ENDPOINTS("function values at endpoints do not have different signs, endpoints: [{0}, {1}], values: [{2}, {3}]"),

    /** UNABLE_TO_BRACKET_OPTIMUM_IN_LINE_SEARCH. */
    UNABLE_TO_BRACKET_OPTIMUM_IN_LINE_SEARCH("unable to bracket optimum in line search"),

    /** UNABLE_TO_COMPUTE_COVARIANCE_SINGULAR_PROBLEM. */
    UNABLE_TO_COMPUTE_COVARIANCE_SINGULAR_PROBLEM("unable to compute covariances: singular problem"),

    /** UNABLE_TO_FIRST_GUESS_HARMONIC_COEFFICIENTS. */
    UNABLE_TO_FIRST_GUESS_HARMONIC_COEFFICIENTS("unable to first guess the harmonic coefficients"),

    /** UNPARSEABLE_3D_VECTOR. */
    UNPARSEABLE_3D_VECTOR("unparseable 3D vector: \"{0}\""),

    /** UNPARSEABLE_COMPLEX_NUMBER. */
    UNPARSEABLE_COMPLEX_NUMBER("unparseable complex number: \"{0}\""),

    /** UNPARSEABLE_REAL_VECTOR. */
    UNPARSEABLE_REAL_VECTOR("unparseable real vector: \"{0}\""),

    /** UNSUPPORTED_EXPANSION_MODE. */
    UNSUPPORTED_EXPANSION_MODE("unsupported expansion mode {0}, supported modes are {1} ({2}) and {3} ({4})"),

    /** VECTOR_LENGTH_MISMATCH. */
    VECTOR_LENGTH_MISMATCH("vector length mismatch: got {0} but expected {1}"),

    /** WRONG_BLOCK_LENGTH. */
    WRONG_BLOCK_LENGTH("wrong array shape (block length = {0}, expected {1})");

    /** Source English format. */
    private final String sourceFormat;

    /** Simple constructor.
     * @param sourceFormat source English format to use when no
     * localized version is available
     */
    LocalizedMigrationFormats(final String sourceFormat) {
        this.sourceFormat = sourceFormat;
    }

    /** {@inheritDoc} */
    @Override
    public String getSourceString() {
        return sourceFormat;
    }

    /** {@inheritDoc} */
    @Override
    public String getLocalizedString(final Locale locale) {
        try {
            final String path = LocalizedMigrationFormats.class.getName().replaceAll("\\.", "/");
            ResourceBundle bundle =
                    ResourceBundle.getBundle("assets/" + path, locale, new UTF8Control());
            if (bundle.getLocale().getLanguage().equals(locale.getLanguage())) {
                final String translated = bundle.getString(name());
                if ((translated != null) &&
                    (translated.length() > 0) &&
                    (!translated.toLowerCase(locale).contains("missing translation"))) {
                    // the value of the resource is the translated format
                    return translated;
                }
            }

        } catch (MissingResourceException mre) { // NOPMD
            // do nothing here
        }

        // either the locale is not supported or the resource is unknown
        // don't translate and fall back to using the source format
        return sourceFormat;

    }

}
