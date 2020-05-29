/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus;

import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Interface representing a <a href="http://mathworld.wolfram.com/RealNumber.html">real</a>
 * <a href="http://mathworld.wolfram.com/Field.html">field</a>.
 * @param <T> the type of the field elements
 * @see FieldElement
 */
public interface RealFieldElement<T> extends CalculusFieldElement<T> {

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     */
    T remainder(double a);

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     * @exception MathIllegalArgumentException if number of free parameters or orders are inconsistent
     */
    T remainder(T a)
        throws MathIllegalArgumentException;

    /** absolute value.
     * @return abs(this)
     */
    T abs();

}
