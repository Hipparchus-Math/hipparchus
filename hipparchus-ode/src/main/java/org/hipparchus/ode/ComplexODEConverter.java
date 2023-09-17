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

package org.hipparchus.ode;

import org.hipparchus.complex.Complex;

/** This class converts {@link ComplexOrdinaryDifferentialEquation complex Ordinary
 * Differential Equations} into {@link OrdinaryDifferentialEquation real ones}.
 *
 * <p>This class is a wrapper around a {@link ComplexOrdinaryDifferentialEquation} which
 * allow to use a {@link ODEIntegrator} to integrate it.</p>
 *
 * <p>The transformation is done by changing the n dimension state
 * vector to a 2n dimension vector, where the even components are
 * real parts and odd components are imaginary parts.</p>
 *
 * <p>One should be aware that the data is duplicated during the
 * transformation process and that for each call to {@link
 * OrdinaryDifferentialEquation#computeDerivatives(double, double[])
 * computeDerivatives}, this wrapper does copy 4n scalars : 2n before
 * the call to {@link
 * OrdinaryDifferentialEquation#computeDerivatives(double, double[])
 * computeDerivatives} in order to dispatch the y state vector,
 * and 2n after the call to gather zDot. Since the underlying problem
 * by itself perhaps also needs to copy data and dispatch the arrays
 * into domain objects, this has an impact on both memory and CPU usage.
 * The only way to avoid this duplication is to perform the transformation
 * at the problem level, i.e. to implement the problem as a first order one
 * and then avoid using this class.</p>
 *
 * <p>
 * The proper way to use the converter is as follows:
 * </p>
 * <pre>
 *   ODEIntegrator                       integrator       = ...build some integrator...;
 *   ComplexOrdinaryDifferentialEquation complexEquations = ...set up the complex problem...;
 *   ComplexODEState                     initialState     = ...set up initial state...;
 *   ComplexODEConverter                 converter        = new ComplexODEConverter();
 *   ComplexODEStateAndDerivative        finalstate       =
 *      converter.convertStateAndDerivative(integrator.integrate(converter.convertEquations(complexEquations),
 *                                                               converter.convertState(initialState),
 *                                                               t);
 * </pre>
 * <p>
 * If there are {@link ComplexSecondaryODE complex secondary equations}, they must be converted
 * too and both the converted primary equations and converted secondary equations must be
 * combined together using {@link ExpandableODE ExpandableODE} as usual for regular real equations.
 * </p>
 *
 * @see ComplexOrdinaryDifferentialEquation
 * @see OrdinaryDifferentialEquation
 * @since 1.4
 */

public class ComplexODEConverter {

    /** Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 3.0
     */
    public ComplexODEConverter() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** Convert an equations set.
     * @param equations equations to convert
     * @return converted equations
     */
    public OrdinaryDifferentialEquation convertEquations(final ComplexOrdinaryDifferentialEquation equations) {
        return new OrdinaryDifferentialEquation() {

            /** {@inheritDoc}
             * <p>The dimension of the real problem is twice the
             * dimension of the underlying complex problem.</p>
             * @return dimension of the problem
             */
            @Override
            public int getDimension() {
                return 2 * equations.getDimension();
            }

            /** {@inheritDoc} */
            @Override
            public void init(final double t0, final double[] y0, final double finalTime) {
                equations.init(t0, convert(y0), finalTime);
            }

            /** {@inheritDoc} */
            @Override
            public double[] computeDerivatives(final double t, final double[] y) {
                return convert(equations.computeDerivatives(t, convert(y)));
            }

        };
    }

    /** Convert a secondary equations set.
     * @param equations equations to convert
     * @return converted equations
     */
    public SecondaryODE convertSecondaryEquations(final ComplexSecondaryODE equations) {
        return new SecondaryODE() {

            /** {@inheritDoc}
             * <p>The dimension of the real problem is twice the
             * dimension of the underlying complex problem.</p>
             * @return dimension of the problem
             */
            @Override
            public int getDimension() {
                return 2 * equations.getDimension();
            }

            /** {@inheritDoc} */
            @Override
            public void init(final double t0, final double[] primary0, final double[] secondary0, final double finalTime) {
                equations.init(t0, convert(primary0), convert(secondary0), finalTime);
            }

            /** {@inheritDoc} */
            @Override
            public double[] computeDerivatives(final double t, final double[] primary, final double[] primaryDot,
                                               final double[] secondary) {
                return convert(equations.computeDerivatives(t, convert(primary), convert(primaryDot), convert(secondary)));
            }

        };
    }

    /** Convert a complex state (typically the initial state).
     * @param state state to convert
     * @return converted state
     */
    public ODEState convertState(final ComplexODEState state) {
        final double[][] secondary = new double[state.getNumberOfSecondaryStates()][];
        for (int index = 0; index < secondary.length; ++index) {
            secondary[index] = convert(state.getSecondaryState(index + 1));
        }
        return new ODEState(state.getTime(),
                            convert(state.getPrimaryState()),
                            secondary);
    }

    /** Convert a real state and derivatives (typically the final state or some intermediate state for
     * step handling or event handling).
     * @param state state to convert
     * @return converted state
     */
    public ComplexODEStateAndDerivative convertState(final ODEStateAndDerivative state) {
        final Complex[][] secondary           = new Complex[state.getNumberOfSecondaryStates()][];
        final Complex[][] secondaryDerivative = new Complex[state.getNumberOfSecondaryStates()][];
        for (int index = 0; index < secondary.length; ++index) {
            secondary[index]           = convert(state.getSecondaryState(index + 1));
            secondaryDerivative[index] = convert(state.getSecondaryDerivative(index + 1));
        }
        return new ComplexODEStateAndDerivative(state.getTime(),
                                                convert(state.getPrimaryState()),
                                                convert(state.getPrimaryDerivative()),
                                                secondary, secondaryDerivative);
    }

    /** Convert a real array into a complex array.
     * @param a array to convert
     * @return converted array
     */
    private Complex[] convert(final double[] a) {
        final Complex[] converted = new Complex[a.length / 2];
        for (int i = 0; i < converted.length; ++i) {
            converted[i] = new Complex(a[2 * i], a[2 * i + 1]);
        }
        return converted;
    }

    /** Convert a complex array into a real array.
     * @param a array to convert
     * @return converted array
     */
    private double[] convert(final Complex[] a) {
        final double[] converted = new double[a.length * 2];
        for (int i = 0; i < a.length; ++i) {
            converted[2 * i]     = a[i].getReal();
            converted[2 * i + 1] = a[i].getImaginary();
        }
        return converted;
    }

}
