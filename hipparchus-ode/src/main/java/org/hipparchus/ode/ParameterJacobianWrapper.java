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
package org.hipparchus.ode;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;

/** Wrapper class to compute Jacobian matrices by finite differences for ODE
 *  which do not compute them by themselves.
 *
 */
class ParameterJacobianWrapper implements NamedParameterJacobianProvider {

    /** Main ODE set. */
    private final OrdinaryDifferentialEquation fode;

    /** Raw ODE without Jacobian computation skill to be wrapped into a ParameterJacobianProvider. */
    private final ParametersController pode;

    /** Steps for finite difference computation of the Jacobian df/dp w.r.t. parameters. */
    private final Map<String, Double> hParam;

    /** Wrap a {@link ParametersController} into a {@link NamedParameterJacobianProvider}.
     * @param fode main first order differential equations set
     * @param pode secondary problem, without parameter Jacobian computation skill
     * @param paramsAndSteps parameters and steps to compute the Jacobians df/dp
     * @see JacobianMatrices#setParameterStep(String, double)
     */
    ParameterJacobianWrapper(final OrdinaryDifferentialEquation fode,
                             final ParametersController pode,
                             final ParameterConfiguration[] paramsAndSteps) {
        this.fode = fode;
        this.pode = pode;
        this.hParam = new HashMap<String, Double>();

        // set up parameters for jacobian computation
        for (final ParameterConfiguration param : paramsAndSteps) {
            final String name = param.getParameterName();
            if (pode.isSupported(name)) {
                hParam.put(name, param.getHP());
            }
        }
    }

    /** {@inheritDoc} */
    @Override
    public Collection<String> getParametersNames() {
        return pode.getParametersNames();
    }

    /** {@inheritDoc} */
    @Override
    public boolean isSupported(String name) {
        return pode.isSupported(name);
    }

    /** {@inheritDoc} */
    @Override
    public double[] computeParameterJacobian(final double t, final double[] y,
                                             final double[] yDot, final String paramName)
        throws MathIllegalArgumentException, MathIllegalStateException {

        final int n = fode.getDimension();
        final double[] dFdP = new double[n];
        if (pode.isSupported(paramName)) {

            // compute the jacobian df/dp w.r.t. parameter
            final double p  = pode.getParameter(paramName);
            final double hP = hParam.get(paramName);
            pode.setParameter(paramName, p + hP);
            final double[] tmpDot = fode.computeDerivatives(t, y);
            for (int i = 0; i < n; ++i) {
                dFdP[i] = (tmpDot[i] - yDot[i]) / hP;
            }
            pode.setParameter(paramName, p);
        } else {
            Arrays.fill(dFdP, 0, n, 0.0);
        }

        return dFdP;

    }

}
