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
package org.hipparchus.special.elliptic.jacobi;

import org.hipparchus.complex.Complex;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.special.elliptic.jacobi.CopolarC;
import org.hipparchus.special.elliptic.jacobi.CopolarD;
import org.hipparchus.special.elliptic.jacobi.CopolarN;
import org.hipparchus.special.elliptic.jacobi.JacobiElliptic;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.Assert;
import org.junit.Test;

public class JacobiThetaTest {

    @Test
    public void testNoConvergence() {
        try {
            new JacobiTheta(Double.NaN).values(Complex.ZERO);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testRealZero() {
        final LegendreEllipticIntegral ei     = new LegendreEllipticIntegral(0.675);
        final double           q      = ei.getNome();
        final double           t3Ref  = 1 + 2 * (q + FastMath.pow(q, 4) + FastMath.pow(q, 9) + FastMath.pow(q, 16));
        final double           theta3 = new JacobiTheta(q).values(Complex.ZERO).theta3().getRealPart();
        Assert.assertEquals(t3Ref, theta3, 1.0e-12);
    }

    @Test
    public void testQuarterPeriod() {
        final LegendreEllipticIntegral ei     = new LegendreEllipticIntegral(0.675);
        final double           q      = ei.getNome();
        final double           theta3 = new JacobiTheta(q).values(Complex.ZERO).theta3().getRealPart();
        Assert.assertEquals(ei.getBigK(), MathUtils.SEMI_PI * theta3 * theta3, 1.0e-12);
    }

    @Test
    public void testEllipticFunctions() {

        final double z = 1.3;
        final LegendreEllipticIntegral ei = new LegendreEllipticIntegral(0.675);
        final double zeta = MathUtils.SEMI_PI * z / ei.getBigK();
        final JacobiTheta      jt = new JacobiTheta(ei.getNome());
        final Theta            theta0 = jt.values(Complex.ZERO);
        final Theta            thetaZ = jt.values(new Complex(zeta));

        // the theta functions are related to the elliptic functions
        // see https://dlmf.nist.gov/22.2
        final JacobiElliptic je = JacobiEllipticBuilder.build(ei.getK() * ei.getK());
        final CopolarN valuesN = je.valuesN(z);
        final CopolarD valuesD = je.valuesD(z);
        final CopolarC valuesC = je.valuesC(z);
        final double t02 = theta0.theta2().getRealPart();
        final double t03 = theta0.theta3().getRealPart();
        final double t04 = theta0.theta4().getRealPart();
        final double tz1 = thetaZ.theta1().getRealPart();
        final double tz2 = thetaZ.theta2().getRealPart();
        final double tz3 = thetaZ.theta3().getRealPart();
        final double tz4 = thetaZ.theta4().getRealPart();
        Assert.assertEquals(valuesN.sn(), t03 * tz1       / (t02 * tz4),       1.0e-15);
        Assert.assertEquals(valuesN.cn(), t04 * tz2       / (t02 * tz4),       1.0e-15);
        Assert.assertEquals(valuesN.dn(), t04 * tz3       / (t03 * tz4),       1.0e-15);
        Assert.assertEquals(valuesD.sd(), t03 * t03 * tz1 / (t02 * t04 * tz3), 1.0e-15);
        Assert.assertEquals(valuesD.cd(), t03 * tz2       / (t02 * tz3),       1.0e-15);
        Assert.assertEquals(valuesC.sc(), t03 * tz1       / (t04 * tz2),       1.0e-15);

    }

}
