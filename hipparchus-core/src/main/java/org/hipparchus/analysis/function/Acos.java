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

package org.hipparchus.analysis.function;

import org.hipparchus.analysis.differentiation.Derivative;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.ExtendedUnivariateDifferentiableFunction;
import org.hipparchus.util.FastMath;

/**
 * Arc-cosine function.
 *
 */
public class Acos implements ExtendedUnivariateDifferentiableFunction {
    /** {@inheritDoc} */
    @Override
    public double value(double x) {
        return FastMath.acos(x);
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure value(final DerivativeStructure t) {
        return t.acos();
    }

    /** {@inheritDoc} */
    @Override
    public <T extends Derivative<T>> T value(T x) {
        return x.acos();
    }

}
