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

package org.hipparchus.ode.nonstiff;


import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Random;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.DenseOutputModel;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.ode.sampling.StepInterpolatorTestUtils;
import org.junit.Assert;
import org.junit.Test;

public class DormandPrince853StateInterpolatorTest {

    @Test
    public void derivativesConsistency()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem3 pb = new TestProblem3(0.1);
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = scalAbsoluteTolerance;
        DormandPrince853Integrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                                          scalAbsoluteTolerance,
                                                                          scalRelativeTolerance);
        StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 1.8e-12);
    }

    @Test
    public void serialization()
                    throws IOException, ClassNotFoundException,
                    MathIllegalArgumentException, MathIllegalStateException {

        TestProblem3 pb = new TestProblem3(0.9);
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = scalAbsoluteTolerance;
        DormandPrince853Integrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                                          scalAbsoluteTolerance,
                                                                          scalRelativeTolerance);
        integ.addStepHandler(new DenseOutputModel());
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        ObjectOutputStream    oos = new ObjectOutputStream(bos);
        for (ODEStepHandler handler : integ.getStepHandlers()) {
            oos.writeObject(handler);
        }

        Assert.assertTrue(bos.size () > 90000);
        Assert.assertTrue(bos.size () < 100000);

        ByteArrayInputStream  bis = new ByteArrayInputStream(bos.toByteArray());
        ObjectInputStream     ois = new ObjectInputStream(bis);
        DenseOutputModel cm  = (DenseOutputModel) ois.readObject();

        Random random = new Random(347588535632l);
        double maxError = 0.0;
        for (int i = 0; i < 1000; ++i) {
            double r = random.nextDouble();
            double time = r * pb.getInitialTime() + (1.0 - r) * pb.getFinalTime();
            double[] interpolatedY = cm.getInterpolatedState(time).getPrimaryState();
            double[] theoreticalY  = pb.computeTheoreticalState(time);
            double dx = interpolatedY[0] - theoreticalY[0];
            double dy = interpolatedY[1] - theoreticalY[1];
            double error = dx * dx + dy * dy;
            if (error > maxError) {
                maxError = error;
            }
        }

        Assert.assertTrue(maxError < 2.4e-10);

    }

}
