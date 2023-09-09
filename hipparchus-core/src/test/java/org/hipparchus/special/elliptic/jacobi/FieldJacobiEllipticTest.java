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
package org.hipparchus.special.elliptic.jacobi;

import java.util.function.DoubleFunction;
import java.util.function.Function;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.analysis.differentiation.UnivariateDerivative1;
import org.hipparchus.dfp.Dfp;
import org.hipparchus.dfp.DfpField;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class FieldJacobiEllipticTest {

    @Test
    public void testCircular() {
        doTestCircular(Binary64Field.getInstance());
    }

    @Test
    public void testHyperbolic() {
        doTestHyperbolic(Binary64Field.getInstance());
    }

    @Test
    public void testNoConvergence() {
        doTestNoConvergence(Binary64Field.getInstance());
    }

    @Test
    public void testNegativeParameter() {
        doTestNegativeParameter(Binary64Field.getInstance());
    }

    @Test
    public void testAbramowitzStegunExample1() {
        doTestAbramowitzStegunExample1(Binary64Field.getInstance());
    }

    @Test
    public void testAbramowitzStegunExample2() {
        doTestAbramowitzStegunExample2(Binary64Field.getInstance());
    }

    @Test
    public void testAbramowitzStegunExample3() {
        doTestAbramowitzStegunExample3(Binary64Field.getInstance());
    }

    @Test
    public void testAbramowitzStegunExample4() {
        doTestAbramowitzStegunExample4(Binary64Field.getInstance());
    }

    @Test
    public void testAbramowitzStegunExample5() {
        doTestAbramowitzStegunExample5(Binary64Field.getInstance());
    }

    @Test
    public void testAbramowitzStegunExample7() {
        doTestAbramowitzStegunExample7(Binary64Field.getInstance());
    }

    @Test
    public void testAbramowitzStegunExample8() {
        doTestAbramowitzStegunExample8(Binary64Field.getInstance());
    }

    @Test
    public void testAbramowitzStegunExample9() {
        doTestAbramowitzStegunExample9(Binary64Field.getInstance());
    }

    @Test
    public void testAllFunctions() {
        doTestAllFunctions(Binary64Field.getInstance());
    }

    @Test
    public void testHighAccuracy() {
        final DfpField field = new DfpField(60);
        FieldJacobiElliptic<Dfp> je = JacobiEllipticBuilder.build(field.newDfp("0.75"));
        Dfp sn        = je.valuesN(field.newDfp("1.3")).sn();
        // this value was computed using Wolfram Alpha
        Dfp reference = field.newDfp("0.8929235150418389265984488063926925504375953835259703680383");
        Assert.assertTrue(sn.subtract(reference).abs().getReal() < 5.0e-58);
    }

    @Test
    public void testInverseCopolarN() {
        doTestInverseCopolarN(Binary64Field.getInstance());
    }

    @Test
    public void testInverseCopolarS() {
        doTestInverseCopolarS(Binary64Field.getInstance());
    }

    @Test
    public void testInverseCopolarC() {
        doTestInverseCopolarC(Binary64Field.getInstance());
    }

    @Test
    public void testInverseCopolarD() {
        doTestInverseCopolarD(Binary64Field.getInstance());
    }

    @Test
    public void testDerivatives() {

        final double m  = 0.75;
        final double m1 = 1 - m;
        final double u  = 1.3;

        JacobiElliptic jeD = JacobiEllipticBuilder.build(m);
        CopolarN valuesND = jeD.valuesN(u);
        CopolarC valuesCD = jeD.valuesC(u);
        CopolarS valuesSD = jeD.valuesS(u);
        CopolarD valuesDD = jeD.valuesD(u);

        FieldJacobiElliptic<UnivariateDerivative1> jeU = JacobiEllipticBuilder.build(new UnivariateDerivative1(m, 0.0));
        FieldCopolarN<UnivariateDerivative1> valuesNU = jeU.valuesN(new UnivariateDerivative1(u, 1.0));
        FieldCopolarC<UnivariateDerivative1> valuesCU = jeU.valuesC(new UnivariateDerivative1(u, 1.0));
        FieldCopolarS<UnivariateDerivative1> valuesSU = jeU.valuesS(new UnivariateDerivative1(u, 1.0));
        FieldCopolarD<UnivariateDerivative1> valuesDU = jeU.valuesD(new UnivariateDerivative1(u, 1.0));

        // see Abramowitz and Stegun section 16.16
        Assert.assertEquals(      valuesND.cn() * valuesND.dn(), valuesNU.sn().getFirstDerivative(), 2.0e-15);
        Assert.assertEquals(-1  * valuesND.sn() * valuesND.dn(), valuesNU.cn().getFirstDerivative(), 2.0e-15);
        Assert.assertEquals(-m  * valuesND.sn() * valuesND.cn(), valuesNU.dn().getFirstDerivative(), 2.0e-15);

        Assert.assertEquals(-m1 * valuesDD.sd() * valuesDD.nd(), valuesDU.cd().getFirstDerivative(), 2.0e-15);
        Assert.assertEquals(      valuesDD.cd() * valuesDD.nd(), valuesDU.sd().getFirstDerivative(), 2.0e-15);
        Assert.assertEquals( m  * valuesDD.sd() * valuesDD.cd(), valuesDU.nd().getFirstDerivative(), 2.0e-15);

        Assert.assertEquals( m1 * valuesCD.sc() * valuesCD.nc(), valuesCU.dc().getFirstDerivative(), 2.0e-15);
        Assert.assertEquals(      valuesCD.sc() * valuesCD.dc(), valuesCU.nc().getFirstDerivative(), 2.0e-15);
        Assert.assertEquals(      valuesCD.dc() * valuesCD.nc(), valuesCU.sc().getFirstDerivative(), 2.0e-15);

        Assert.assertEquals(-1  * valuesSD.ds() * valuesSD.cs(), valuesSU.ns().getFirstDerivative(), 2.0e-15);
        Assert.assertEquals(-1  * valuesSD.cs() * valuesSD.ns(), valuesSU.ds().getFirstDerivative(), 2.0e-15);
        Assert.assertEquals(-1  * valuesSD.ns() * valuesSD.ds(), valuesSU.cs().getFirstDerivative(), 2.0e-15);

    }

    private <T extends CalculusFieldElement<T>> void doTestCircular(final Field<T> field) {
        for (double m : new double[] { -1.0e-10, 0.0, 1.0e-10 }) {
         final double eps = 3 * FastMath.max(1.0e-14, FastMath.abs(m));
            final FieldJacobiElliptic<T> je = build(field, m);
            for (double t = -10; t < 10; t += 0.01) {
                final FieldCopolarN<T> n = je.valuesN(t);
                Assert.assertEquals(FastMath.sin(t), n.sn().getReal(), eps);
                Assert.assertEquals(FastMath.cos(t), n.cn().getReal(), eps);
                Assert.assertEquals(1.0,             n.dn().getReal(), eps);
            }
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestHyperbolic(final Field<T> field) {
        for (double m1 : new double[] { -1.0e-12, 0.0, 1.0e-12 }) {
            final double eps = 3 * FastMath.max(1.0e-14, FastMath.abs(m1));
            final FieldJacobiElliptic<T> je = build(field, 1.0 - m1);
            for (double t = -3; t < 3; t += 0.01) {
                final FieldCopolarN<T> n = je.valuesN(t);
                Assert.assertEquals(FastMath.tanh(t),       n.sn().getReal(), eps);
                Assert.assertEquals(1.0 / FastMath.cosh(t), n.cn().getReal(), eps);
                Assert.assertEquals(1.0 / FastMath.cosh(t), n.dn().getReal(), eps);
            }
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestNoConvergence(final Field<T> field) {
        Assert.assertTrue(build(field, Double.NaN).valuesS(0.0).cs().isNaN());
    }

    private <T extends CalculusFieldElement<T>> void doTestNegativeParameter(final Field<T> field) {
        Assert.assertEquals(0.49781366219021166315, build(field, -4.5).valuesN(8.3).sn().getReal(), 1.5e-10);
        Assert.assertEquals(0.86728401215332559984, build(field, -4.5).valuesN(8.3).cn().getReal(), 1.5e-10);
        Assert.assertEquals(1.45436686918553524215, build(field, -4.5).valuesN(8.3).dn().getReal(), 1.5e-10);
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample1(final Field<T> field) {
        // Abramowitz and Stegun give a result of -1667, but Wolfram Alpha gives the following value
        Assert.assertEquals(-1392.11114434139393839735, build(field, 0.64).valuesC(1.99650).nc().getReal(), 2.8e-10);
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample2(final Field<T> field) {
        Assert.assertEquals(0.996253, build(field, 0.19).valuesN(0.20).dn().getReal(), 1.0e-6);
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample3(final Field<T> field) {
        Assert.assertEquals(0.984056, build(field, 0.81).valuesN(0.20).dn().getReal(), 1.0e-6);
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample4(final Field<T> field) {
        Assert.assertEquals(0.980278, build(field, 0.81).valuesN(0.20).cn().getReal(), 1.0e-6);
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample5(final Field<T> field) {
        Assert.assertEquals(0.60952, build(field, 0.36).valuesN(0.672).sn().getReal(), 1.0e-5);
        Assert.assertEquals(1.1740, build(field, 0.36).valuesC(0.672).dc().getReal(), 1.0e-4);
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample7(final Field<T> field) {
        Assert.assertEquals(1.6918083, build(field, 0.09).valuesS(0.5360162).cs().getReal(), 1.0e-7);
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample8(final Field<T> field) {
        Assert.assertEquals(0.56458, build(field, 0.5).valuesN(0.61802).sn().getReal(), 1.0e-5);
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample9(final Field<T> field) {
        Assert.assertEquals(0.68402, build(field, 0.5).valuesC(0.61802).sc().getReal(), 1.0e-5);
    }

    private <T extends CalculusFieldElement<T>> void doTestAllFunctions(final Field<T> field) {
        // reference was computed from Wolfram Alpha, using the square relations
        // from Abramowitz and Stegun section 16.9 for the functions Wolfram Alpha
        // did not understood (i.e. for the sake of validation we did *not* use the
        // relations from section 16.3 which are used in the Hipparchus implementation)
        final double u = 1.4;
        final double m = 0.7;
        final double[] reference = {
              0.92516138673582827365, 0.37957398289798418747, 0.63312991237590996850,
              0.41027866958131945027, 0.68434537093007175683, 1.08089249544689572795,
              1.66800134071905681841, 2.63453251554589286796, 2.43736775548306830513,
              1.57945467502452678756, 1.46125047743207819361, 0.59951990180590090343
        };
        final FieldJacobiElliptic<T> je = build(field, m);
        Assert.assertEquals(reference[ 0], je.valuesN(u).sn().getReal(), 4 * FastMath.ulp(reference[ 0]));
        Assert.assertEquals(reference[ 1], je.valuesN(u).cn().getReal(), 4 * FastMath.ulp(reference[ 1]));
        Assert.assertEquals(reference[ 2], je.valuesN(u).dn().getReal(), 4 * FastMath.ulp(reference[ 2]));
        Assert.assertEquals(reference[ 3], je.valuesS(u).cs().getReal(), 4 * FastMath.ulp(reference[ 3]));
        Assert.assertEquals(reference[ 4], je.valuesS(u).ds().getReal(), 4 * FastMath.ulp(reference[ 4]));
        Assert.assertEquals(reference[ 5], je.valuesS(u).ns().getReal(), 4 * FastMath.ulp(reference[ 5]));
        Assert.assertEquals(reference[ 6], je.valuesC(u).dc().getReal(), 4 * FastMath.ulp(reference[ 6]));
        Assert.assertEquals(reference[ 7], je.valuesC(u).nc().getReal(), 4 * FastMath.ulp(reference[ 7]));
        Assert.assertEquals(reference[ 8], je.valuesC(u).sc().getReal(), 4 * FastMath.ulp(reference[ 8]));
        Assert.assertEquals(reference[ 9], je.valuesD(u).nd().getReal(), 4 * FastMath.ulp(reference[ 9]));
        Assert.assertEquals(reference[10], je.valuesD(u).sd().getReal(), 4 * FastMath.ulp(reference[10]));
        Assert.assertEquals(reference[11], je.valuesD(u).cd().getReal(), 4 * FastMath.ulp(reference[11]));
    }

    private <T extends CalculusFieldElement<T>> FieldJacobiElliptic<T> build(final Field<T> field, final double m) {
        return JacobiEllipticBuilder.build(field.getZero().newInstance(m));
    }

    private <T extends CalculusFieldElement<T>> void doTestInverseCopolarN(final Field<T> field) {
        final double m = 0.7;
        final FieldJacobiElliptic<T> je = build(field, m);
        doTestInverse(-0.80,  0.80, 100, field, u -> je.valuesN(u).sn(), x -> je.arcsn(x), x -> je.arcsn(x), 1.0e-14);
        doTestInverse(-1.00,  1.00, 100, field, u -> je.valuesN(u).cn(), x -> je.arccn(x), x -> je.arccn(x), 1.0e-14);
        doTestInverse( 0.55,  1.00, 100, field, u -> je.valuesN(u).dn(), x -> je.arcdn(x), x -> je.arcdn(x), 1.0e-14);
    }

    private <T extends CalculusFieldElement<T>> void doTestInverseCopolarS(final Field<T> field) {
        final double m = 0.7;
        final FieldJacobiElliptic<T> je = build(field, m);
        doTestInverse(-2.00,  2.00, 100, field, u -> je.valuesS(u).cs(), x -> je.arccs(x), x -> je.arccs(x), 1.0e-14);
        doTestInverse( 0.55,  2.00, 100, field, u -> je.valuesS(u).ds(), x -> je.arcds(x), x -> je.arcds(x), 1.0e-14);
        doTestInverse(-2.00, -0.55, 100, field, u -> je.valuesS(u).ds(), x -> je.arcds(x), x -> je.arcds(x), 1.0e-14);
        doTestInverse( 1.00,  2.00, 100, field, u -> je.valuesS(u).ns(), x -> je.arcns(x), x -> je.arcns(x), 1.0e-11);
        doTestInverse(-2.00, -1.00, 100, field, u -> je.valuesS(u).ns(), x -> je.arcns(x), x -> je.arcns(x), 1.0e-11);
    }

    private <T extends CalculusFieldElement<T>> void doTestInverseCopolarC(final Field<T> field) {
        final double m = 0.7;
        final FieldJacobiElliptic<T> je = build(field, m);
        doTestInverse( 1.00,  2.00, 100, field, u -> je.valuesC(u).dc(), x -> je.arcdc(x), x -> je.arcdc(x), 1.0e-14);
        doTestInverse(-2.00, -1.00, 100, field, u -> je.valuesC(u).dc(), x -> je.arcdc(x), x -> je.arcdc(x), 1.0e-14);
        doTestInverse( 1.00,  2.00, 100, field, u -> je.valuesC(u).nc(), x -> je.arcnc(x), x -> je.arcnc(x), 1.0e-14);
        doTestInverse(-2.00, -1.00, 100, field, u -> je.valuesC(u).nc(), x -> je.arcnc(x), x -> je.arcnc(x), 1.0e-14);
        doTestInverse(-2.00,  2.00, 100, field, u -> je.valuesC(u).sc(), x -> je.arcsc(x), x -> je.arcsc(x), 1.0e-14);
    }

    private <T extends CalculusFieldElement<T>> void doTestInverseCopolarD(final Field<T> field) {
        final double m = 0.7;
        final FieldJacobiElliptic<T> je = build(field, m);
        doTestInverse( 1.00,  1.80, 100, field, u -> je.valuesD(u).nd(), x -> je.arcnd(x), x -> je.arcnd(x), 1.0e-14);
        doTestInverse(-1.80,  1.80, 100, field, u -> je.valuesD(u).sd(), x -> je.arcsd(x), x -> je.arcsd(x), 1.0e-14);
        doTestInverse(-1.00,  1.00, 100, field, u -> je.valuesD(u).cd(), x -> je.arccd(x), x -> je.arccd(x), 1.0e-14);
    }

    private <T extends CalculusFieldElement<T>> void doTestInverse(final double xMin, final double xMax, final int n,
                                                                   final Field<T> field,
                                                                   final Function<T, T> direct,
                                                                   final Function<T, T> inverseField,
                                                                   final DoubleFunction<T> inverseDouble,
                                                                   final double tolerance) {
        for (int i = 0; i < n; ++i) {
            final T x             = field.getZero().newInstance(xMin + i * (xMax - xMin) / (n - 1));
            final T xFieldRebuilt = direct.apply(inverseField.apply(x));
            Assert.assertEquals(x.getReal(), xFieldRebuilt.getReal(), tolerance);
            final T xDoubleRebuilt = direct.apply(inverseDouble.apply(x.getReal()));
            Assert.assertEquals(x.getReal(), xDoubleRebuilt.getReal(), tolerance);
        }
    }

}
