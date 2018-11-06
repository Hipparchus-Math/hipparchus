/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.complex;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Comparator for Complex Numbers.
 *
 */
public class ComplexComparator implements Comparator<Complex>, Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20171113L;

    /** Compare two complex numbers, using real ordering as the primary sort order and
     * imaginary ordering as the secondary sort order.
     * @param o1 first complex number
     * @param o2 second complex number
     * @return a negative value if o1 real part is less than o2 real part
     * or if real parts are equal and o1 imaginary part is less than o2 imaginary part
     */
    @Override
    public int compare(Complex o1, Complex o2) {
        if (o1 == null) {
            return o2 == null ? 0 : -1;
        } else if (o2 == null) {
            return 1;
        }

        final int cR = Double.compare(o1.getReal(), o2.getReal());
        if (cR == 0) {
            return Double.compare(o1.getImaginary(),o2.getImaginary());
        } else {
            return cR;
        }
    }

}
