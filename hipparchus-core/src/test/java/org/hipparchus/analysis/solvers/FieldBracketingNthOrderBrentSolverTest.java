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

package org.hipparchus.analysis.solvers;

import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.FieldUnivariateFunction;
import org.hipparchus.analysis.RealFieldUnivariateFunction;
import org.hipparchus.dfp.Dfp;
import org.hipparchus.dfp.DfpField;
import org.hipparchus.dfp.DfpMath;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Test case for {@link FieldBracketingNthOrderBrentSolver bracketing n<sup>th</sup> order Brent} solver.
 *
 */
public final class FieldBracketingNthOrderBrentSolverTest {

    @Test(expected=MathIllegalArgumentException.class)
    public void testInsufficientOrder3() {
        new FieldBracketingNthOrderBrentSolver<Dfp>(relativeAccuracy, absoluteAccuracy,
                                                    functionValueAccuracy, 1);
    }

    @Test
    public void testConstructorOK() {
        FieldBracketingNthOrderBrentSolver<Dfp> solver =
                new FieldBracketingNthOrderBrentSolver<Dfp>(relativeAccuracy, absoluteAccuracy,
                                                            functionValueAccuracy, 2);
        Assert.assertEquals(2, solver.getMaximalOrder());
    }

    @Test
    public void testConvergenceOnFunctionAccuracy() {
        FieldBracketingNthOrderBrentSolver<Dfp> solver =
                new FieldBracketingNthOrderBrentSolver<Dfp>(relativeAccuracy, absoluteAccuracy,
                                                            field.newDfp(1.0e-20), 20);
        FieldUnivariateFunction f = new FieldUnivariateFunction() {
            public <T extends RealFieldElement<T>> T value(T x) {
                T one     = x.getField().getOne();
                T oneHalf = one.divide(2);
                T xMo     = x.subtract(one);
                T xMh     = x.subtract(oneHalf);
                T xPh     = x.add(oneHalf);
                T xPo     = x.add(one);
                return xMo.multiply(xMh).multiply(x).multiply(xPh).multiply(xPo);
            }
        };

        Dfp result = solver.solve(20, f.toRealFieldUnivariateFunction(field), field.newDfp(0.2), field.newDfp(0.9),
                                  field.newDfp(0.4), AllowedSolution.BELOW_SIDE);
        Assert.assertTrue(f.value(result).abs().lessThan(solver.getFunctionValueAccuracy()));
        Assert.assertTrue(f.value(result).negativeOrNull());
        Assert.assertTrue(result.subtract(field.newDfp(0.5)).subtract(solver.getAbsoluteAccuracy()).positiveOrNull());
        result = solver.solve(20, f.toRealFieldUnivariateFunction(field), field.newDfp(-0.9), field.newDfp(-0.2),
                              field.newDfp(-0.4), AllowedSolution.ABOVE_SIDE);
        Assert.assertTrue(f.value(result).abs().lessThan(solver.getFunctionValueAccuracy()));
        Assert.assertTrue(f.value(result).positiveOrNull());
        Assert.assertTrue(result.add(field.newDfp(0.5)).subtract(solver.getAbsoluteAccuracy()).negativeOrNull());
    }

    @Test
    public void testToleranceLessThanUlp() {
        // function that is never zero
        Dfp zero = field.getZero();
        Dfp one = field.getOne();
        RealFieldUnivariateFunction<Dfp> f = (x) -> x.getReal() <= 2.1 ? one.negate(): one;
        // tolerance less than 1 ulp(x)
        FieldBracketingNthOrderBrentSolver<Dfp> solver =
                new FieldBracketingNthOrderBrentSolver<>(zero, field.newDfp(1e-55), zero, 5);

        // make sure it doesn't throw a maxIterations exception
        Dfp result = solver.solve(200, f, zero, zero.add(5.0), AllowedSolution.LEFT_SIDE);
        double difference = field.newDfp(2.1).subtract(result).abs().getReal();
        Assert.assertTrue("difference: " + difference, difference < FastMath.ulp(2.1));
    }

    @Test
    public void testNeta() {

        // the following test functions come from Beny Neta's paper:
        // "Several New Methods for solving Equations"
        // intern J. Computer Math Vol 23 pp 265-282
        // available here: http://www.math.nps.navy.mil/~bneta/SeveralNewMethods.PDF
        for (AllowedSolution allowed : AllowedSolution.values()) {
            check(new RealFieldUnivariateFunction<Dfp>() {
                public Dfp value(Dfp x) {
                    return DfpMath.sin(x).subtract(x.divide(2));
                }
            }, 200, -2.0, 2.0, allowed);

            check(new RealFieldUnivariateFunction<Dfp>() {
                public Dfp value(Dfp x) {
                    return DfpMath.pow(x, 5).add(x).subtract(field.newDfp(10000));
                }
            }, 200, -5.0, 10.0, allowed);

            check(new RealFieldUnivariateFunction<Dfp>() {
                public Dfp value(Dfp x) {
                    return x.sqrt().subtract(field.getOne().divide(x)).subtract(field.newDfp(3));
                }
            }, 200, 0.001, 10.0, allowed);

            check(new RealFieldUnivariateFunction<Dfp>() {
                public Dfp value(Dfp x) {
                    return DfpMath.exp(x).add(x).subtract(field.newDfp(20));
                }
            }, 200, -5.0, 5.0, allowed);

            check(new RealFieldUnivariateFunction<Dfp>() {
                public Dfp value(Dfp x) {
                    return DfpMath.log(x).add(x.sqrt()).subtract(field.newDfp(5));
                }
            }, 200, 0.001, 10.0, allowed);

            check(new RealFieldUnivariateFunction<Dfp>() {
                public Dfp value(Dfp x) {
                    return x.subtract(field.getOne()).multiply(x).multiply(x).subtract(field.getOne());
                }
            }, 200, -0.5, 1.5, allowed);
        }

    }

    private void check(RealFieldUnivariateFunction<Dfp> f, int maxEval, double min, double max,
                       AllowedSolution allowedSolution) {
        FieldBracketingNthOrderBrentSolver<Dfp> solver =
                new FieldBracketingNthOrderBrentSolver<Dfp>(relativeAccuracy, absoluteAccuracy,
                                                     functionValueAccuracy, 20);
        Dfp xResult = solver.solve(maxEval, f, field.newDfp(min), field.newDfp(max),
                                   allowedSolution);
        Dfp yResult = f.value(xResult);
        switch (allowedSolution) {
        case ANY_SIDE :
            Assert.assertTrue(yResult.abs().lessThan(functionValueAccuracy.multiply(2)));
            break;
        case LEFT_SIDE : {
            boolean increasing = f.value(xResult).add(absoluteAccuracy).greaterThan(yResult);
            Assert.assertTrue(increasing ? yResult.negativeOrNull() : yResult.positiveOrNull());
            break;
        }
        case RIGHT_SIDE : {
            boolean increasing = f.value(xResult).add(absoluteAccuracy).greaterThan(yResult);
            Assert.assertTrue(increasing ? yResult.positiveOrNull() : yResult.negativeOrNull());
            break;
        }
        case BELOW_SIDE :
            Assert.assertTrue(yResult.negativeOrNull());
            break;
        case ABOVE_SIDE :
            Assert.assertTrue(yResult.positiveOrNull());
            break;
        default :
            // this should never happen
            throw new MathRuntimeException(null);
        }
    }

    @Before
    public void setUp() {
        field                 = new DfpField(50);
        absoluteAccuracy      = field.newDfp(1.0e-45);
        relativeAccuracy      = field.newDfp(1.0e-45);
        functionValueAccuracy = field.newDfp(1.0e-45);
    }

    private DfpField field;
    private Dfp      absoluteAccuracy;
    private Dfp      relativeAccuracy;
    private Dfp      functionValueAccuracy;

}
