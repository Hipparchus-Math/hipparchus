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
package org.hipparchus.migration.ode;

import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.NamedParameterJacobianProvider;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.ParameterConfiguration;
import org.hipparchus.ode.ParametersController;
import org.hipparchus.ode.SecondaryODE;

/**
 * This class defines a set of {@link SecondaryODE secondary equations} to
 * compute the Jacobian matrices with respect to the initial state vector and, if
 * any, to some parameters of the primary ODE set.
 * <p>
 * It is intended to be packed into an {@link ExpandableODE}
 * in conjunction with a primary set of ODE, which may be:</p>
 * <ul>
 * <li>a {@link FirstOrderDifferentialEquations}</li>
 * <li>a {@link MainStateJacobianProvider}</li>
 * </ul>
 * <p>In order to compute Jacobian matrices with respect to some parameters of the
 * primary ODE set, the following parameter Jacobian providers may be set:</p>
 * <ul>
 * <li>a {@link ParametersController}</li>
 * </ul>
 *
 * @see ExpandableODE
 * @see FirstOrderDifferentialEquations
 * @see MainStateJacobianProvider
 * @see NamedParameterJacobianProvider
 * @see ParametersController
 * @deprecated as of 1.0, replaced with {@link org.hipparchus.ode.VariationalEquation}
 */
@Deprecated
public class JacobianMatrices {

    /** Expandable first order differential equation. */
    private ExpandableODE efode;

    /** Index of the instance in the expandable set. */
    private int index;

    /** FODE with exact primary Jacobian computation skill. */
    private MainStateJacobianProvider jode;

    /** FODE without exact parameter Jacobian computation skill. */
    private ParametersController parametersController;

    /** Primary state vector dimension. */
    private int stateDim;

    /** Selected parameters for parameter Jacobian computation. */
    private MutableParameterConfiguration[] selectedParameters;

    /** FODE with exact parameter Jacobian computation skill. */
    private List<NamedParameterJacobianProvider> jacobianProviders;

    /** Parameters dimension. */
    private int paramDim;

    /** Boolean for selected parameters consistency. */
    private boolean dirtyParameter;

    /** State and parameters Jacobian matrices in a row. */
    private double[] matricesData;

    /** Simple constructor for a secondary equations set computing Jacobian matrices.
     * <p>
     * Parameters must belong to the supported ones given by {@link
     * org.hipparchus.ode.Parameterizable#getParametersNames()}, so the primary set of differential
     * equations must be {@link org.hipparchus.ode.Parameterizable}.
     * </p>
     * <p>Note that each selection clears the previous selected parameters.</p>
     *
     * @param fode the primary first order differential equations set to extend
     * @param hY step used for finite difference computation with respect to state vector
     * @param parameters parameters to consider for Jacobian matrices processing
     * (may be null if parameters Jacobians is not desired)
     * @exception MathIllegalArgumentException if there is a dimension mismatch between
     * the steps array {@code hY} and the equation dimension
     */
    public JacobianMatrices(final OrdinaryDifferentialEquation fode, final double[] hY,
                            final String... parameters)
        throws MathIllegalArgumentException {
        this(new MainStateJacobianWrapper(fode, hY), parameters);
    }

    /** Simple constructor for a secondary equations set computing Jacobian matrices.
     * <p>
     * Parameters must belong to the supported ones given by {@link
     * org.hipparchus.ode.Parameterizable#getParametersNames()}, so the primary set of differential
     * equations must be {@link org.hipparchus.ode.Parameterizable}.
     * </p>
     * <p>Note that each selection clears the previous selected parameters.</p>
     *
     * @param jode the primary first order differential equations set to extend
     * @param parameters parameters to consider for Jacobian matrices processing
     * (may be null if parameters Jacobians is not desired)
     */
    public JacobianMatrices(final MainStateJacobianProvider jode,
                            final String... parameters) {

        this.efode = null;
        this.index = -1;

        this.jode = jode;
        this.parametersController = null;

        this.stateDim = jode.getDimension();

        if (parameters == null) {
            selectedParameters = null;
            paramDim = 0;
        } else {
            this.selectedParameters = new MutableParameterConfiguration[parameters.length];
            for (int i = 0; i < parameters.length; ++i) {
                selectedParameters[i] = new MutableParameterConfiguration(parameters[i], Double.NaN);
            }
            paramDim = parameters.length;
        }
        this.dirtyParameter = false;

        this.jacobianProviders = new ArrayList<>();

        // set the default initial state Jacobian to the identity
        // and the default initial parameters Jacobian to the null matrix
        matricesData = new double[(stateDim + paramDim) * stateDim];
        for (int i = 0; i < stateDim; ++i) {
            matricesData[i * (stateDim + 1)] = 1.0;
        }

    }

    /** Register the variational equations for the Jacobians matrices to the expandable set.
     * <p>
     * This method must be called <em>before {@link #setUpInitialState(ODEState)}</em>
     * </p>
     * @param expandable expandable set into which variational equations should be registered
     * @throws MathIllegalArgumentException if the dimension of the partial state does not
     * match the selected equations set dimension
     * @exception MismatchedEquations if the primary set of the expandable set does
     * not match the one used to build the instance
     * @see ExpandableODE#addSecondaryEquations(SecondaryODE)
     * @see #setUpInitialState(ODEState)
     */
    public void registerVariationalEquations(final ExpandableODE expandable)
        throws MathIllegalArgumentException, MismatchedEquations {

        // safety checks
        final OrdinaryDifferentialEquation ode = (jode instanceof MainStateJacobianWrapper) ?
                                                 ((MainStateJacobianWrapper) jode).ode :
                                                 jode;
        if (expandable.getPrimary() != ode) {
            throw new MismatchedEquations();
        }

        efode = expandable;
        index = efode.addSecondaryEquations(new JacobiansSecondaryODE());

    }

    /** Set up initial state.
     * <p>
     * This method inserts the initial Jacobian matrices data into
     * an {@link ODEState ODE state} by overriding the additional
     * state components corresponding to the instance. It must be
     * called prior to integrate the equations.
     * </p>
     * <p>
     * This method must be called <em>after</em> {@link
     * #registerVariationalEquations(ExpandableODE)},
     * {@link #setInitialMainStateJacobian(double[][])} and
     * {@link #setInitialParameterJacobian(String, double[])}.
     * </p>
     * @param initialState initial state, without the initial Jacobians
     * matrices
     * @return a new instance of initial state, with the initial Jacobians
     * matrices properly initialized
     */
    public ODEState setUpInitialState(final ODEState initialState) { // NOPMD - PMD false positive

        // insert the matrices data into secondary states
        final double[][] secondary = new double[efode.getMapper().getNumberOfEquations() - 1][];
        for (int i = 0; i < initialState.getNumberOfSecondaryStates(); ++i) {
            if (i + 1 != index) {
                secondary[i] = initialState.getSecondaryState(i + 1);
            }
        }
        secondary[index - 1] = matricesData;

        // create an updated initial state
        return new ODEState(initialState.getTime(), initialState.getPrimaryState(), secondary);

    }

    /** Add a parameter Jacobian provider.
     * @param provider the parameter Jacobian provider to compute exactly the parameter Jacobian matrix
     */
    public void addParameterJacobianProvider(final NamedParameterJacobianProvider provider) {
        jacobianProviders.add(provider);
    }

    /** Set a parameter Jacobian provider.
     * @param pc the controller to compute the parameter Jacobian matrix using finite differences
     * @deprecated as of 1.0, replaced with {@link #setParametersController(ParametersController)}
     */
    @Deprecated
    public void setParameterizedODE(final ParametersController pc) {
        setParametersController(pc);
    }

    /** Set a parameter Jacobian provider.
     * @param parametersController the controller to compute the parameter Jacobian matrix using finite differences
     */
    public void setParametersController(final ParametersController parametersController) {
        this.parametersController = parametersController;
        dirtyParameter = true;
    }

    /** Set the step associated to a parameter in order to compute by finite
     *  difference the Jacobian matrix.
     * <p>
     * Needed if and only if the primary ODE set is a {@link ParametersController}.
     * </p>
     * <p>
     * Given a non zero parameter value pval for the parameter, a reasonable value
     * for such a step is {@code pval * FastMath.sqrt(Precision.EPSILON)}.
     * </p>
     * <p>
     * A zero value for such a step doesn't enable to compute the parameter Jacobian matrix.
     * </p>
     * @param parameter parameter to consider for Jacobian processing
     * @param hP step for Jacobian finite difference computation w.r.t. the specified parameter
     * @see ParametersController
     * @exception MathIllegalArgumentException if the parameter is not supported
     */
    public void setParameterStep(final String parameter, final double hP)
        throws MathIllegalArgumentException {

        for (MutableParameterConfiguration param: selectedParameters) {
            if (parameter.equals(param.getParameterName())) {
                param.setHP(hP);
                dirtyParameter = true;
                return;
            }
        }

        throw new MathIllegalArgumentException(LocalizedODEFormats.UNKNOWN_PARAMETER, parameter);

    }

    /** Set the initial value of the Jacobian matrix with respect to state.
     * <p>
     * If this method is not called, the initial value of the Jacobian
     * matrix with respect to state is set to identity.
     * </p>
     * <p>
     * This method must be called <em>before {@link #setUpInitialState(ODEState)}</em>
     * </p>
     * @param dYdY0 initial Jacobian matrix w.r.t. state
     * @exception MathIllegalArgumentException if matrix dimensions are incorrect
     */
    public void setInitialMainStateJacobian(final double[][] dYdY0)
        throws MathIllegalArgumentException {

        // Check dimensions
        checkDimension(stateDim, dYdY0);
        checkDimension(stateDim, dYdY0[0]);

        // store the matrix in row major order as a single dimension array
        int i = 0;
        for (final double[] row : dYdY0) {
            System.arraycopy(row, 0, matricesData, i, stateDim);
            i += stateDim;
        }

    }

    /** Set the initial value of a column of the Jacobian matrix with respect to one parameter.
     * <p>
     * If this method is not called for some parameter, the initial value of
     * the column of the Jacobian matrix with respect to this parameter is set to zero.
     * </p>
     * <p>
     * This method must be called <em>before {@link #setUpInitialState(ODEState)}</em>
     * </p>
     * @param pName parameter name
     * @param dYdP initial Jacobian column vector with respect to the parameter
     * @exception MathIllegalArgumentException if a parameter is not supported
     * @throws MathIllegalArgumentException if the column vector does not match state dimension
     */
    public void setInitialParameterJacobian(final String pName, final double[] dYdP)
        throws MathIllegalArgumentException {

        // Check dimensions
        checkDimension(stateDim, dYdP);

        // store the column in a global single dimension array
        int i = stateDim * stateDim;
        for (MutableParameterConfiguration param: selectedParameters) {
            if (pName.equals(param.getParameterName())) {
                System.arraycopy(dYdP, 0, matricesData, i, stateDim);
                return;
            }
            i += stateDim;
        }

        throw new MathIllegalArgumentException(LocalizedODEFormats.UNKNOWN_PARAMETER, pName);

    }

    /** Extract the Jacobian matrix with respect to state.
     * @param state state from which to extract Jacobian matrix
     * @return Jacobian matrix dY/dY0 with respect to state.
     */
    public double[][] extractMainSetJacobian(final ODEState state) {

        // get current state for this set of equations from the expandable fode
        final double[] p = state.getSecondaryState(index);

        final double[][] dYdY0 = new double[stateDim][stateDim];
        int j = 0;
        for (int i = 0; i < stateDim; i++) {
            System.arraycopy(p, j, dYdY0[i], 0, stateDim);
            j += stateDim;
        }

        return dYdY0;

    }

    /** Extract the Jacobian matrix with respect to one parameter.
     * @param state state from which to extract Jacobian matrix
     * @param pName name of the parameter for the computed Jacobian matrix
     * @return Jacobian matrix dY/dP with respect to the named parameter
     */
    public double[] extractParameterJacobian(final ODEState state, final String pName) {

        // get current state for this set of equations from the expandable fode
        final double[] p = state.getSecondaryState(index);

        final double[] dYdP = new double[stateDim];
        int i = stateDim * stateDim;
        for (MutableParameterConfiguration param: selectedParameters) {
            if (param.getParameterName().equals(pName)) {
                System.arraycopy(p, i, dYdP, 0, stateDim);
                break;
            }
            i += stateDim;
        }

        return dYdP;

    }

    /** Check array dimensions.
     * @param expected expected dimension
     * @param array (may be null if expected is 0)
     * @throws MathIllegalArgumentException if the array dimension does not match the expected one
     */
    private void checkDimension(final int expected, final Object array)
        throws MathIllegalArgumentException {
        int arrayDimension = (array == null) ? 0 : Array.getLength(array);
        if (arrayDimension != expected) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   arrayDimension, expected);
        }
    }

    /** Local implementation of secondary equations.
     * <p>
     * This class is an inner class to ensure proper scheduling of calls
     * by forcing the use of {@link JacobianMatrices#registerVariationalEquations(ExpandableODE)}.
     * </p>
     */
    private class JacobiansSecondaryODE implements SecondaryODE {

        /** {@inheritDoc} */
        @Override
        public int getDimension() {
            return stateDim * (stateDim + paramDim);
        }

        /** {@inheritDoc} */
        @Override
        public double[] computeDerivatives(final double t, final double[] y, final double[] yDot,
                                           final double[] z)
            throws MathIllegalArgumentException, MathIllegalStateException {

            try {

                // Lazy initialization
                Constructor<ParameterConfiguration> configCtr =
                                ParameterConfiguration.class.getDeclaredConstructor(String.class, Double.TYPE);
                configCtr.setAccessible(true); // NOPMD
                @SuppressWarnings("unchecked")
                Constructor<NamedParameterJacobianProvider> providerCtr =
                (Constructor<NamedParameterJacobianProvider>)
                Class.forName("org.hipparchus.ode.ParameterJacobianWrapper").getDeclaredConstructor(OrdinaryDifferentialEquation.class,
                                                                                                    double[].class,
                                                                                                    ParametersController.class,
                                                                                                    ParameterConfiguration[].class);
                providerCtr.setAccessible(true); // NOPMD
                if (dirtyParameter && (paramDim != 0)) {
                    ParameterConfiguration [] immutable = new ParameterConfiguration[selectedParameters.length];
                    for (int i = 0; i < selectedParameters.length; ++i) {
                        immutable[i] = configCtr.newInstance(selectedParameters[i].getParameterName(),
                                                             selectedParameters[i].getHP());
                    }
                    jacobianProviders.add(providerCtr.newInstance(jode, new double[jode.getDimension()],
                                                                  parametersController, immutable));
                    dirtyParameter = false;
                }

            } catch (InstantiationException | IllegalAccessException | IllegalArgumentException |
                     InvocationTargetException | NoSuchMethodException | SecurityException | ClassNotFoundException e) {
                throw new MathIllegalStateException(e, LocalizedCoreFormats.SIMPLE_MESSAGE, e.getLocalizedMessage());
            }

            // variational equations:
            // from d[dy/dt]/dy0 and d[dy/dt]/dp to d[dy/dy0]/dt and d[dy/dp]/dt

            // compute Jacobian matrix with respect to primary state
            double[][] dFdY = jode.computeMainStateJacobian(t, y, yDot);

            // Dispatch Jacobian matrix in the compound secondary state vector
            final double[] zDot = new double[z.length];
            for (int i = 0; i < stateDim; ++i) {
                final double[] dFdYi = dFdY[i];
                for (int j = 0; j < stateDim; ++j) {
                    double s = 0;
                    final int startIndex = j;
                    int zIndex = startIndex;
                    for (int l = 0; l < stateDim; ++l) {
                        s += dFdYi[l] * z[zIndex];
                        zIndex += stateDim;
                    }
                    zDot[startIndex + i * stateDim] = s;
                }
            }

            if (paramDim != 0) {
                // compute Jacobian matrices with respect to parameters
                int startIndex = stateDim * stateDim;
                for (MutableParameterConfiguration param: selectedParameters) {
                    boolean found = false;
                    for (int k = 0 ; (!found) && (k < jacobianProviders.size()); ++k) {
                        final NamedParameterJacobianProvider provider = jacobianProviders.get(k);
                        if (provider.isSupported(param.getParameterName())) {
                            final double[] dFdP =
                                            provider.computeParameterJacobian(t, y, yDot,
                                                                              param.getParameterName());
                            for (int i = 0; i < stateDim; ++i) {
                                final double[] dFdYi = dFdY[i];
                                int zIndex = startIndex;
                                double s = dFdP[i];
                                for (int l = 0; l < stateDim; ++l) {
                                    s += dFdYi[l] * z[zIndex];
                                    zIndex++;
                                }
                                zDot[startIndex + i] = s;
                            }
                            found = true;
                        }
                    }
                    if (! found) {
                        Arrays.fill(zDot, startIndex, startIndex + stateDim, 0.0);
                    }
                    startIndex += stateDim;
                }
            }

            return zDot;

        }
    }

    /** Wrapper class to compute jacobian matrices by finite differences for ODE
     *  which do not compute them by themselves.
     */
    private static class MainStateJacobianWrapper implements MainStateJacobianProvider {

        /** Raw ODE without jacobians computation skill to be wrapped into a MainStateJacobianProvider. */
        private final OrdinaryDifferentialEquation ode;

        /** Steps for finite difference computation of the jacobian df/dy w.r.t. state. */
        private final double[] hY;

        /** Wrap a {@link FirstOrderDifferentialEquations} into a {@link MainStateJacobianProvider}.
         * @param ode original ODE problem, without jacobians computation skill
         * @param hY step sizes to compute the jacobian df/dy
         * @exception MathIllegalArgumentException if there is a dimension mismatch between
         * the steps array {@code hY} and the equation dimension
         */
        MainStateJacobianWrapper(final OrdinaryDifferentialEquation ode,
                                 final double[] hY)
            throws MathIllegalArgumentException {
            this.ode = ode;
            this.hY = hY.clone();
            if (hY.length != ode.getDimension()) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                       ode.getDimension(), hY.length);
            }
        }

        /** {@inheritDoc} */
        @Override
        public int getDimension() {
            return ode.getDimension();
        }

        /** {@inheritDoc} */
        @Override
        public double[] computeDerivatives(double t, double[] y)
            throws MathIllegalArgumentException, MathIllegalStateException {
            return ode.computeDerivatives(t, y);
        }

        /** {@inheritDoc} */
        @Override
        public double[][] computeMainStateJacobian(double t, double[] y, double[] yDot)
            throws MathIllegalArgumentException, MathIllegalStateException {

            final int n = ode.getDimension();
            final double[][] dFdY = new double[n][n];
            for (int j = 0; j < n; ++j) {
                final double savedYj = y[j];
                y[j] += hY[j];
                final double[] tmpDot = ode.computeDerivatives(t, y);
                for (int i = 0; i < n; ++i) {
                    dFdY[i][j] = (tmpDot[i] - yDot[i]) / hY[j];
                }
                y[j] = savedYj;
            }
            return dFdY;
        }

    }

    /**
     * Special exception for equations mismatch.
     */
    public static class MismatchedEquations extends MathIllegalArgumentException {

        /** Serializable UID. */
        private static final long serialVersionUID = 20120902L;

        /** Simple constructor. */
        public MismatchedEquations() {
            super(LocalizedODEFormats.UNMATCHED_ODE_IN_EXPANDED_SET);
        }

    }

    /** Selected parameter for parameter Jacobian computation. */
    private static class MutableParameterConfiguration {

        /** Parameter name. */
        private String parameterName;

        /** Parameter step for finite difference computation. */
        private double hP;

        /** Parameter name and step pair constructor.
         * @param parameterName parameter name
         * @param hP parameter step
         */
        MutableParameterConfiguration(final String parameterName, final double hP) {
            this.parameterName = parameterName;
            this.hP = hP;
        }

        /** Get parameter name.
         * @return parameterName parameter name
         */
        public String getParameterName() {
            return parameterName;
        }

        /** Get parameter step.
         * @return hP parameter step
         */
        public double getHP() {
            return hP;
        }

        /** Set parameter step.
         * @param hParam parameter step
         */
        public void setHP(final double hParam) {
            this.hP = hParam;
        }

    }

}

