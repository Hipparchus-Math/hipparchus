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
package org.hipparchus.geometry;

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
 */
public enum LocalizedGeometryFormats implements Localizable {

    /** CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR. */
    CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR("cannot normalize a zero norm vector"),

    /** CLOSE_VERTICES. */
    CLOSE_VERTICES("too close vertices near point ({0}, {1}, {2})"),

    /** CLOSEST_ORTHOGONAL_MATRIX_HAS_NEGATIVE_DETERMINANT. */
    CLOSEST_ORTHOGONAL_MATRIX_HAS_NEGATIVE_DETERMINANT("the closest orthogonal matrix has a negative determinant {0}"),

    /** CROSSING_BOUNDARY_LOOPS. */
    CROSSING_BOUNDARY_LOOPS("some outline boundary loops cross each other"),

    /** EDGE_CONNECTED_TO_ONE_FACET. */
    EDGE_CONNECTED_TO_ONE_FACET("edge joining points ({0}, {1}, {2}) and ({3}, {4}, {5}) is connected to one facet only"),

    /** FACET_ORIENTATION_MISMATCH. */
    FACET_ORIENTATION_MISMATCH("facets orientation mismatch around edge joining points ({0}, {1}, {2}) and ({3}, {4}, {5})"),

    /** INCONSISTENT_STATE_AT_2_PI_WRAPPING. */
    INCONSISTENT_STATE_AT_2_PI_WRAPPING("inconsistent state at 2\u03c0 wrapping"),

    /** NON_INVERTIBLE_TRANSFORM. */
    NON_INVERTIBLE_TRANSFORM("non-invertible affine transform collapses some lines into single points"),

    /** NOT_CONVEX. */
    NOT_CONVEX("vertices do not form a convex hull in CCW winding"),

    /** NOT_CONVEX_HYPERPLANES. */
    NOT_CONVEX_HYPERPLANES("hyperplanes do not define a convex region"),

    /** NOT_SUPPORTED_IN_DIMENSION_N. */
    NOT_SUPPORTED_IN_DIMENSION_N("method not supported in dimension {0}"),

    /** OUTLINE_BOUNDARY_LOOP_OPEN. */
    OUTLINE_BOUNDARY_LOOP_OPEN("an outline boundary loop is open"),

    /** FACET_WITH_SEVERAL_BOUNDARY_LOOPS. */
    FACET_WITH_SEVERAL_BOUNDARY_LOOPS("a facet has several boundary loops"),

    /** OUT_OF_PLANE. */
    OUT_OF_PLANE("point ({0}, {1}, {2}) is out of plane"),

    /** ROTATION_MATRIX_DIMENSIONS. */
    ROTATION_MATRIX_DIMENSIONS("a {0}x{1} matrix cannot be a rotation matrix"),

    /** UNABLE_TO_ORTHOGONOLIZE_MATRIX. */
    UNABLE_TO_ORTHOGONOLIZE_MATRIX("unable to orthogonalize matrix in {0} iterations"),

    /** ZERO_NORM_FOR_ROTATION_AXIS. */
    ZERO_NORM_FOR_ROTATION_AXIS("zero norm for rotation axis"),

    /** ZERO_NORM_FOR_ROTATION_DEFINING_VECTOR. */
    ZERO_NORM_FOR_ROTATION_DEFINING_VECTOR("zero norm for rotation defining vector"),

    /** TOO_SMALL_TOLERANCE. */
    TOO_SMALL_TOLERANCE("tolerance {0,number,0.00000E00} is not computationally feasible, it is smaller than {1} ({2,number,0.00000E00})"),

    /** INVALID_ROTATION_ORDER_NAME. */
    INVALID_ROTATION_ORDER_NAME("the value {0} does not correspond to a rotation order");

    /** Source English format. */
    private final String sourceFormat;

    /** Simple constructor.
     * @param sourceFormat source English format to use when no
     * localized version is available
     */
    LocalizedGeometryFormats(final String sourceFormat) {
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
            final String path = LocalizedGeometryFormats.class.getName().replaceAll("\\.", "/");
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
