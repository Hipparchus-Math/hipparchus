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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.AbstractIntegrator;
import org.hipparchus.ode.AbstractParameterizable;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.NamedParameterJacobianProvider;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.ParametersController;
import org.hipparchus.ode.nonstiff.DormandPrince54Integrator;
import org.hipparchus.stat.descriptive.StreamingStatistics;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

@Deprecated
public class JacobianMatricesTest {

    @Test
    public void testLowAccuracyExternalDifferentiation()
        throws MathIllegalArgumentException, MathIllegalStateException {
        // this test does not really test JacobianMatrices,
        // it only shows that WITHOUT this class, attempting to recover
        // the jacobians from external differentiation on simple integration
        // results with low accuracy gives very poor results. In fact,
        // the curves dy/dp = g(b) when b varies from 2.88 to 3.08 are
        // essentially noise.
        // This test is taken from Hairer, Norsett and Wanner book
        // Solving Ordinary Differential Equations I (Nonstiff problems),
        // the curves dy/dp = g(b) are in figure 6.5
        ODEIntegrator integ =
            new DormandPrince54Integrator(1.0e-8, 100.0, new double[] { 1.0e-4, 1.0e-4 }, new double[] { 1.0e-4, 1.0e-4 });
        double hP = 1.0e-12;
        StreamingStatistics residualsP0 = new StreamingStatistics();
        StreamingStatistics residualsP1 = new StreamingStatistics();
        for (double b = 2.88; b < 3.08; b += 0.001) {
            Brusselator brusselator = new Brusselator(b);
            double[] y = { 1.3, b };
            y = integ.integrate(brusselator, new ODEState(0, y), 20.0).getPrimaryState();
            double[] yP = { 1.3, b + hP };
            yP = integ.integrate(brusselator, new ODEState(0, yP), 20.0).getPrimaryState();
            residualsP0.addValue((yP[0] - y[0]) / hP - brusselator.dYdP0());
            residualsP1.addValue((yP[1] - y[1]) / hP - brusselator.dYdP1());
        }
        Assert.assertTrue((residualsP0.getMax() - residualsP0.getMin()) > 500);
        Assert.assertTrue(residualsP0.getStandardDeviation() > 30);
        Assert.assertTrue((residualsP1.getMax() - residualsP1.getMin()) > 700);
        Assert.assertTrue(residualsP1.getStandardDeviation() > 40);
    }

    @Test
    public void testHighAccuracyExternalDifferentiation()
        throws MathIllegalArgumentException, MathIllegalStateException {
        ODEIntegrator integ =
            new DormandPrince54Integrator(1.0e-8, 100.0, new double[] { 1.0e-10, 1.0e-10 }, new double[] { 1.0e-10, 1.0e-10 });
        double hP = 1.0e-12;
        StreamingStatistics residualsP0 = new StreamingStatistics();
        StreamingStatistics residualsP1 = new StreamingStatistics();
        for (double b = 2.88; b < 3.08; b += 0.001) {
            ParamBrusselator brusselator = new ParamBrusselator(b);
            double[] y = { 1.3, b };
            y = integ.integrate(brusselator, new ODEState(0, y), 20.0).getPrimaryState();
            double[] yP = { 1.3, b + hP };
            brusselator.setParameter("b", b + hP);
            yP = integ.integrate(brusselator, new ODEState(0, yP), 20.0).getPrimaryState();
            residualsP0.addValue((yP[0] - y[0]) / hP - brusselator.dYdP0());
            residualsP1.addValue((yP[1] - y[1]) / hP - brusselator.dYdP1());
        }
        Assert.assertTrue((residualsP0.getMax() - residualsP0.getMin()) > 0.02);
        Assert.assertTrue((residualsP0.getMax() - residualsP0.getMin()) < 0.03);
        Assert.assertTrue(residualsP0.getStandardDeviation() > 0.003);
        Assert.assertTrue(residualsP0.getStandardDeviation() < 0.004);
        Assert.assertTrue((residualsP1.getMax() - residualsP1.getMin()) > 0.04);
        Assert.assertTrue((residualsP1.getMax() - residualsP1.getMin()) < 0.05);
        Assert.assertTrue(residualsP1.getStandardDeviation() > 0.007);
        Assert.assertTrue(residualsP1.getStandardDeviation() < 0.008);
    }

    @Test
    public void testWrongParameterName() {
        final String name = "an-unknown-parameter";
        try {
            ParamBrusselator brusselator = new ParamBrusselator(2.9);
            brusselator.setParameter(name, 3.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException upe) {
            Assert.assertEquals(LocalizedODEFormats.UNKNOWN_PARAMETER, upe.getSpecifier());
            Assert.assertEquals(name, upe.getParts()[0]);
        }
    }

    @Test
    public void testInternalDifferentiation()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        AbstractIntegrator integ =
                        new DormandPrince54Integrator(1.0e-8, 100.0, new double[] { 1.0e-4, 1.0e-4 }, new double[] { 1.0e-4, 1.0e-4 });
        double hP = 1.0e-12;
        double hY = 1.0e-12;
        StreamingStatistics residualsP0 = new StreamingStatistics();
        StreamingStatistics residualsP1 = new StreamingStatistics();
        for (double b = 2.88; b < 3.08; b += 0.001) {
                ParamBrusselator brusselator = new ParamBrusselator(b);
                brusselator.setParameter(ParamBrusselator.B, b);
            double[] z = { 1.3, b };

            JacobianMatrices jacob = new JacobianMatrices(brusselator, new double[] { hY, hY }, ParamBrusselator.B);
            jacob.setParametersController(brusselator);
            jacob.setParameterStep(ParamBrusselator.B, hP);
            jacob.setInitialParameterJacobian(ParamBrusselator.B, new double[] { 0.0, 1.0 });

            ExpandableODE efode = new ExpandableODE(brusselator);
            jacob.registerVariationalEquations(efode);

            integ.setMaxEvaluations(5000);
            final ODEState initialState = jacob.setUpInitialState(new ODEState(0, z));
            final ODEStateAndDerivative finalState = integ.integrate(efode, initialState, 20.0);
            final double[]   dZdP  = jacob.extractParameterJacobian(finalState, ParamBrusselator.B);
//            Assert.assertEquals(5000, integ.getMaxEvaluations());
//            Assert.assertTrue(integ.getEvaluations() > 1500);
//            Assert.assertTrue(integ.getEvaluations() < 2100);
//            Assert.assertEquals(4 * integ.getEvaluations(), integ.getEvaluations());
            residualsP0.addValue(dZdP[0] - brusselator.dYdP0());
            residualsP1.addValue(dZdP[1] - brusselator.dYdP1());
        }
        Assert.assertTrue((residualsP0.getMax() - residualsP0.getMin()) < 0.02);
        Assert.assertTrue(residualsP0.getStandardDeviation() < 0.003);
        Assert.assertTrue((residualsP1.getMax() - residualsP1.getMin()) < 0.05);
        Assert.assertTrue(residualsP1.getStandardDeviation() < 0.01);
    }

    @Test
    public void testAnalyticalDifferentiation()
        throws MathIllegalArgumentException, MathIllegalStateException {
        AbstractIntegrator integ =
            new DormandPrince54Integrator(1.0e-8, 100.0, new double[] { 1.0e-4, 1.0e-4 }, new double[] { 1.0e-4, 1.0e-4 });
        StreamingStatistics residualsP0 = new StreamingStatistics();
        StreamingStatistics residualsP1 = new StreamingStatistics();
        for (double b = 2.88; b < 3.08; b += 0.001) {
            Brusselator brusselator = new Brusselator(b);
            double[] z = { 1.3, b };

            JacobianMatrices jacob = new JacobianMatrices(brusselator, Brusselator.B);
            jacob.addParameterJacobianProvider(brusselator);
            jacob.setInitialParameterJacobian(Brusselator.B, new double[] { 0.0, 1.0 });

            ExpandableODE efode = new ExpandableODE(brusselator);
            jacob.registerVariationalEquations(efode);

            integ.setMaxEvaluations(5000);
            final ODEState initialState = jacob.setUpInitialState(new ODEState(0, z));
            final ODEStateAndDerivative finalState = integ.integrate(efode, initialState, 20.0);
            final double[] dZdP = jacob.extractParameterJacobian(finalState, Brusselator.B);
//            Assert.assertEquals(5000, integ.getMaxEvaluations());
//            Assert.assertTrue(integ.getEvaluations() > 350);
//            Assert.assertTrue(integ.getEvaluations() < 510);
            residualsP0.addValue(dZdP[0] - brusselator.dYdP0());
            residualsP1.addValue(dZdP[1] - brusselator.dYdP1());
        }
        Assert.assertTrue((residualsP0.getMax() - residualsP0.getMin()) < 0.014);
        Assert.assertTrue(residualsP0.getStandardDeviation() < 0.003);
        Assert.assertTrue((residualsP1.getMax() - residualsP1.getMin()) < 0.05);
        Assert.assertTrue(residualsP1.getStandardDeviation() < 0.01);
    }

    @Test
    public void testFinalResult()
        throws MathIllegalArgumentException, MathIllegalStateException {

        AbstractIntegrator integ =
            new DormandPrince54Integrator(1.0e-8, 100.0, new double[] { 1.0e-10, 1.0e-10 }, new double[] { 1.0e-10, 1.0e-10 });
        double[] y = new double[] { 0.0, 1.0 };
        Circle circle = new Circle(y, 1.0, 1.0, 0.1);

        JacobianMatrices jacob = new JacobianMatrices(circle, Circle.CX, Circle.CY, Circle.OMEGA);
        jacob.addParameterJacobianProvider(circle);
        jacob.setInitialMainStateJacobian(circle.exactDyDy0(0));
        jacob.setInitialParameterJacobian(Circle.CX, circle.exactDyDcx(0));
        jacob.setInitialParameterJacobian(Circle.CY, circle.exactDyDcy(0));
        jacob.setInitialParameterJacobian(Circle.OMEGA, circle.exactDyDom(0));

        ExpandableODE efode = new ExpandableODE(circle);
        jacob.registerVariationalEquations(efode);

        integ.setMaxEvaluations(5000);

        double t = 18 * FastMath.PI;
        final ODEState initialState = jacob.setUpInitialState(new ODEState(0, y));
        final ODEStateAndDerivative finalState = integ.integrate(efode, initialState, t);
        y = finalState.getPrimaryState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(circle.exactY(t)[i], y[i], 1.0e-9);
        }

        double[][] dydy0 = jacob.extractMainSetJacobian(finalState);
        for (int i = 0; i < dydy0.length; ++i) {
            for (int j = 0; j < dydy0[i].length; ++j) {
                Assert.assertEquals(circle.exactDyDy0(t)[i][j], dydy0[i][j], 1.0e-9);
            }
        }
        double[] dydcx = jacob.extractParameterJacobian(finalState, Circle.CX);
        for (int i = 0; i < dydcx.length; ++i) {
            Assert.assertEquals(circle.exactDyDcx(t)[i], dydcx[i], 1.0e-7);
        }
        double[] dydcy = jacob.extractParameterJacobian(finalState, Circle.CY);
        for (int i = 0; i < dydcy.length; ++i) {
            Assert.assertEquals(circle.exactDyDcy(t)[i], dydcy[i], 1.0e-7);
        }
        double[] dydom = jacob.extractParameterJacobian(finalState, Circle.OMEGA);
        for (int i = 0; i < dydom.length; ++i) {
            Assert.assertEquals(circle.exactDyDom(t)[i], dydom[i], 1.0e-7);
        }
    }

    @Test
    public void testParameterizable()
        throws MathIllegalArgumentException, MathIllegalStateException {

        AbstractIntegrator integ =
            new DormandPrince54Integrator(1.0e-8, 100.0, new double[] { 1.0e-10, 1.0e-10 }, new double[] { 1.0e-10, 1.0e-10 });
        double[] y = new double[] { 0.0, 1.0 };
        ParameterizedCircle pcircle = new ParameterizedCircle(y, 1.0, 1.0, 0.1);

        double hP = 1.0e-12;
        double hY = 1.0e-12;

        JacobianMatrices jacob = new JacobianMatrices(pcircle, new double[] { hY, hY },
                                                      ParameterizedCircle.CX, ParameterizedCircle.CY,
                                                      ParameterizedCircle.OMEGA);
        jacob.setParametersController(pcircle);
        jacob.setParameterStep(ParameterizedCircle.CX,    hP);
        jacob.setParameterStep(ParameterizedCircle.CY,    hP);
        jacob.setParameterStep(ParameterizedCircle.OMEGA, hP);
        jacob.setInitialMainStateJacobian(pcircle.exactDyDy0(0));
        jacob.setInitialParameterJacobian(ParameterizedCircle.CX, pcircle.exactDyDcx(0));
        jacob.setInitialParameterJacobian(ParameterizedCircle.CY, pcircle.exactDyDcy(0));
        jacob.setInitialParameterJacobian(ParameterizedCircle.OMEGA, pcircle.exactDyDom(0));

        ExpandableODE efode = new ExpandableODE(pcircle);
        jacob.registerVariationalEquations(efode);

        integ.setMaxEvaluations(50000);

        double t = 18 * FastMath.PI;
        final ODEState initialState = jacob.setUpInitialState(new ODEState(0, y));
        final ODEStateAndDerivative finalState = integ.integrate(efode, initialState, t);
        y = finalState.getPrimaryState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(pcircle.exactY(t)[i], y[i], 1.0e-9);
        }

        double[][] dydy0 = jacob.extractMainSetJacobian(finalState);
        for (int i = 0; i < dydy0.length; ++i) {
            for (int j = 0; j < dydy0[i].length; ++j) {
                Assert.assertEquals(pcircle.exactDyDy0(t)[i][j], dydy0[i][j], 5.0e-4);
            }
        }

        double[] dydp0 = jacob.extractParameterJacobian(finalState, ParameterizedCircle.CX);
        for (int i = 0; i < dydp0.length; ++i) {
            Assert.assertEquals(pcircle.exactDyDcx(t)[i], dydp0[i], 5.0e-4);
        }

        double[] dydp1 =  jacob.extractParameterJacobian(finalState, ParameterizedCircle.OMEGA);
        for (int i = 0; i < dydp1.length; ++i) {
            Assert.assertEquals(pcircle.exactDyDom(t)[i], dydp1[i], 1.0e-2);
        }
    }

    private static class Brusselator extends AbstractParameterizable
        implements MainStateJacobianProvider, NamedParameterJacobianProvider {

        public static final String B = "b";

        private double b;

        public Brusselator(double b) {
            super(B);
            this.b = b;
        }

        @Override
        public int getDimension() {
            return 2;
        }

        @Override
        public double[] computeDerivatives(double t, double[] y) {
            double prod = y[0] * y[0] * y[1];
            return new double[] {
                1 + prod - (b + 1) * y[0],
                b * y[0] - prod
            };
        }

        @Override
        public double[][] computeMainStateJacobian(double t, double[] y, double[] yDot) {
            double p = 2 * y[0] * y[1];
            double y02 = y[0] * y[0];
            return new double[][] {
                { p - (1 + b),  y02 },
                { b - p,  -y02}
            };
        }

        @Override
        public double[] computeParameterJacobian(double t, double[] y, double[] yDot,
                                             String paramName) {
            if (isSupported(paramName)) {
                return new double[] { -y[0], y[0] };
            } else {
                return new double[] { 0, 0 };
            }
        }

        public double dYdP0() {
            return -1088.232716447743 + (1050.775747149553 + (-339.012934631828 + 36.52917025056327 * b) * b) * b;
        }

        public double dYdP1() {
            return 1502.824469929139 + (-1438.6974831849952 + (460.959476642384 - 49.43847385647082 * b) * b) * b;
        }

    }

    private static class ParamBrusselator extends AbstractParameterizable
        implements OrdinaryDifferentialEquation, ParametersController {

        public static final String B = "b";

        private double b;

        public ParamBrusselator(double b) {
            super(B);
            this.b = b;
        }

        @Override
        public int getDimension() {
            return 2;
        }

        /** {@inheritDoc} */
        @Override
        public double getParameter(final String name)
            throws MathIllegalArgumentException {
            complainIfNotSupported(name);
            return b;
        }

        /** {@inheritDoc} */
        @Override
        public void setParameter(final String name, final double value)
            throws MathIllegalArgumentException {
            complainIfNotSupported(name);
            b = value;
        }

        @Override
        public double[] computeDerivatives(double t, double[] y) {
            double prod = y[0] * y[0] * y[1];
            return new double[] {
                1 + prod - (b + 1) * y[0],
                b * y[0] - prod
            };
        }

        public double dYdP0() {
            return -1088.232716447743 + (1050.775747149553 + (-339.012934631828 + 36.52917025056327 * b) * b) * b;
        }

        public double dYdP1() {
            return 1502.824469929139 + (-1438.6974831849952 + (460.959476642384 - 49.43847385647082 * b) * b) * b;
        }

    }

    /** ODE representing a point moving on a circle with provided center and angular rate. */
    private static class Circle extends AbstractParameterizable
        implements MainStateJacobianProvider, NamedParameterJacobianProvider {

        public static final String CX = "cx";
        public static final String CY = "cy";
        public static final String OMEGA = "omega";

        private final double[] y0;
        private double cx;
        private double cy;
        private double omega;

        public Circle(double[] y0, double cx, double cy, double omega) {
            super(CX,CY,OMEGA);
            this.y0    = y0.clone();
            this.cx    = cx;
            this.cy    = cy;
            this.omega = omega;
        }

        @Override
        public int getDimension() {
            return 2;
        }

        @Override
        public double[] computeDerivatives(double t, double[] y) {
            return new double[] {
                omega * (cy - y[1]),
                omega * (y[0] - cx)
            };
        }

        @Override
        public double[][] computeMainStateJacobian(double t, double[] y, double[] yDot) {
            return new double[][] {
                { 0, -omega },
                { omega, 0 }
            };
        }

        @Override
        public double[] computeParameterJacobian(double t, double[] y, double[] yDot, String paramName)
            throws MathIllegalArgumentException {
            complainIfNotSupported(paramName);
            if (paramName.equals(CX)) {
                return new double[] { 0, -omega };
            } else if (paramName.equals(CY)) {
                return new double[] { omega, 0 };
            }  else {
                return new double[] { cy - y[1], y[0] - cx };
            }
        }

        public double[] exactY(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            double dx0 = y0[0] - cx;
            double dy0 = y0[1] - cy;
            return new double[] {
                cx + cos * dx0 - sin * dy0,
                cy + sin * dx0 + cos * dy0
            };
        }

        public double[][] exactDyDy0(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            return new double[][] {
                { cos, -sin },
                { sin,  cos }
            };
        }

        public double[] exactDyDcx(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            return new double[] {1 - cos, -sin};
        }

        public double[] exactDyDcy(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            return new double[] {sin, 1 - cos};
        }

        public double[] exactDyDom(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            double dx0 = y0[0] - cx;
            double dy0 = y0[1] - cy;
            return new double[] { -t * (sin * dx0 + cos * dy0) , t * (cos * dx0 - sin * dy0) };
        }

    }

    /** ODE representing a point moving on a circle with provided center and angular rate. */
    private static class ParameterizedCircle extends AbstractParameterizable
        implements OrdinaryDifferentialEquation, ParametersController {

        public static final String CX = "cx";
        public static final String CY = "cy";
        public static final String OMEGA = "omega";

        private final double[] y0;
        private double cx;
        private double cy;
        private double omega;

        public ParameterizedCircle(double[] y0, double cx, double cy, double omega) {
            super(CX,CY,OMEGA);
            this.y0    = y0.clone();
            this.cx    = cx;
            this.cy    = cy;
            this.omega = omega;
        }

        @Override
        public int getDimension() {
            return 2;
        }

        @Override
        public double[] computeDerivatives(double t, double[] y) {
            return new double[] {
                omega * (cy - y[1]),
                omega * (y[0] - cx)
            };
        }

        @Override
        public double getParameter(final String name)
            throws MathIllegalArgumentException {
            if (name.equals(CX)) {
                return cx;
            } else if (name.equals(CY)) {
                    return cy;
            } else if (name.equals(OMEGA)) {
                return omega;
            } else {
                throw new MathIllegalArgumentException(LocalizedODEFormats.UNKNOWN_PARAMETER, name);
            }
        }

        @Override
        public void setParameter(final String name, final double value)
            throws MathIllegalArgumentException {
            if (name.equals(CX)) {
                cx = value;
            } else if (name.equals(CY)) {
                cy = value;
            } else if (name.equals(OMEGA)) {
                omega = value;
            } else {
                throw new MathIllegalArgumentException(LocalizedODEFormats.UNKNOWN_PARAMETER, name);
            }
        }

        public double[] exactY(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            double dx0 = y0[0] - cx;
            double dy0 = y0[1] - cy;
            return new double[] {
                cx + cos * dx0 - sin * dy0,
                cy + sin * dx0 + cos * dy0
            };
        }

        public double[][] exactDyDy0(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            return new double[][] {
                { cos, -sin },
                { sin,  cos }
            };
        }

        public double[] exactDyDcx(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            return new double[] {1 - cos, -sin};
        }

        public double[] exactDyDcy(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            return new double[] {sin, 1 - cos};
        }

        public double[] exactDyDom(double t) {
            double cos = FastMath.cos(omega * t);
            double sin = FastMath.sin(omega * t);
            double dx0 = y0[0] - cx;
            double dy0 = y0[1] - cy;
            return new double[] { -t * (sin * dx0 + cos * dy0) , t * (cos * dx0 - sin * dy0) };
        }

    }

}
