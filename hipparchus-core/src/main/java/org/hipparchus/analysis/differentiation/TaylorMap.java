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
package org.hipparchus.analysis.differentiation;

/** Container for a Taylor map.
 * <p>
 * A dimension n Taylor map is a set of n {@link DerivativeStructure}
 * \((f_1, f_2, \ldots, f_n)\) depending on n parameters \((p_1, p_2, \ldots, p_n)\).
 * </p>
 * @since 2.2
 */
public class TaylorMap {

    /** Mapping functions. */
    private final DerivativeStructure[] functions;

    /** Simple constructor.
     * @param functions functions composing the map
     */
    public TaylorMap(final DerivativeStructure[] functions) {
        this.functions = functions.clone();
    }

    /** Get the dimension of the map.
     * @return dimension of the map
     */
    public int getDimension() {
        return functions.length;
    }

    /** Evaluate Taylor expansion of the map at some offset.
     * @param deltaP parameters offsets \((\Delta p_1, \Delta p_2, \ldots, \Delta p_n)\)
     * @return value of the Taylor expansion at \((p_1 + \Delta p_1, p_2 + \Delta p_2, \ldots, p_n + \Delta p_n)\)
     */
    public double[] value(final double[] deltaP) {
        final double[] value = new double[functions.length];
        for (int i = 0; i < functions.length; ++i) {
            value[i] = functions[i].taylor(deltaP);
        }
        return value;
    }

    /** Invert the instance.
     * <p>
     * Consider {@link #value(double[]) Taylor expansion} of the map with
     * small parameters offsets \((\Delta p_1, \Delta p_2, \ldots, \Delta p_n)\)
     * which leads to evaluation offsets \((f_1 + df_1, f_2 + df_2, \ldots, f_n + df_n)\).
     * The map inversion defines a Taylor map that computes \((\Delta p_1,
     * \Delta p_2, \ldots, \Delta p_n)\) from \((df_1, df_2, \ldots, df_n)\).
     * @return inverted map
     */
    public TaylorMap invert() {

        final double[][] map = new double[functions.length][];
        for (int i = 0; i < map.length; ++i) {
            map[i] = functions[i].getAllDerivatives();
        }

        final DSFactory  factory   = functions[0].getFactory();
        final double[][] invertedF = factory.getCompiler().invertMap(map);

        final DerivativeStructure[] inverted = new DerivativeStructure[functions.length];
        for (int i = 0; i < map.length; ++i) {
            inverted[i] = factory.build(invertedF[i]);
        }

        return new TaylorMap(inverted);

    }

}
