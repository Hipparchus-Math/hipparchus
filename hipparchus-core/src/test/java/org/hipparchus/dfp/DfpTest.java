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

package org.hipparchus.dfp;

import org.hipparchus.CalculusFieldElementAbstractTest;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class DfpTest extends CalculusFieldElementAbstractTest<Dfp> {

    @Override
    protected Dfp build(final double x) {
        return field.newDfp(x);
    }

    private DfpField field;
    private Dfp pinf;
    private Dfp ninf;
    private Dfp nan;
    private Dfp snan;
    private Dfp qnan;

    @BeforeEach
    void setUp() {
        // Some basic setup.  Define some constants and clear the status flags
        field = new DfpField(20);
        pinf = field.newDfp("1").divide(field.newDfp("0"));
        ninf = field.newDfp("-1").divide(field.newDfp("0"));
        nan = field.newDfp("0").divide(field.newDfp("0"));
        snan = field.newDfp((byte)1, Dfp.SNAN);
        qnan = field.newDfp((byte)1, Dfp.QNAN);
        ninf.getField().clearIEEEFlags();
    }

    @AfterEach
    void tearDown() {
        field = null;
        pinf    = null;
        ninf    = null;
        nan     = null;
        snan    = null;
        qnan    = null;
    }

    // Generic test function.  Takes params x and y and tests them for
    // equality.  Then checks the status flags against the flags argument.
    // If the test fail, it prints the desc string
    private void test(Dfp x, Dfp y, int flags, String desc)
    {
        boolean b = x.equals(y);

        if (!x.equals(y) && !x.unequal(y))  // NaNs involved
            b = (x.toString().equals(y.toString()));

        if (x.equals(field.newDfp("0")))  // distinguish +/- zero
            b = (b && (x.toString().equals(y.toString())));

        b = (b && x.getField().getIEEEFlags() == flags);

        if (!b)
            assertTrue(b, "assersion failed "+desc+" x = "+x.toString()+" flags = "+x.getField().getIEEEFlags());

        x.getField().clearIEEEFlags();
    }

    @Test
    void testByteConstructor() {
        assertEquals("0.", new Dfp(field, (byte) 0).toString());
        assertEquals("1.", new Dfp(field, (byte) 1).toString());
        assertEquals("-1.", new Dfp(field, (byte) -1).toString());
        assertEquals("-128.", new Dfp(field, Byte.MIN_VALUE).toString());
        assertEquals("127.", new Dfp(field, Byte.MAX_VALUE).toString());
    }

    @Test
    void testIntConstructor() {
        assertEquals("0.", new Dfp(field, 0).toString());
        assertEquals("1.", new Dfp(field, 1).toString());
        assertEquals("-1.", new Dfp(field, -1).toString());
        assertEquals("1234567890.", new Dfp(field, 1234567890).toString());
        assertEquals("-1234567890.", new Dfp(field, -1234567890).toString());
        assertEquals("-2147483648.", new Dfp(field, Integer.MIN_VALUE).toString());
        assertEquals("2147483647.", new Dfp(field, Integer.MAX_VALUE).toString());
    }

    @Test
    void testLongConstructor() {
        assertEquals("0.", new Dfp(field, 0l).toString());
        assertEquals("1.", new Dfp(field, 1l).toString());
        assertEquals("-1.", new Dfp(field, -1l).toString());
        assertEquals("1234567890.", new Dfp(field, 1234567890l).toString());
        assertEquals("-1234567890.", new Dfp(field, -1234567890l).toString());
        assertEquals("-9223372036854775808.", new Dfp(field, Long.MIN_VALUE).toString());
        assertEquals("9223372036854775807.", new Dfp(field, Long.MAX_VALUE).toString());
    }

    /*
     *  Test addition
     */
    @Test
    void testAdd()
    {
        test(field.newDfp("1").add(field.newDfp("1")),      // Basic tests   1+1 = 2
             field.newDfp("2"),
             0, "Add #1");

        test(field.newDfp("1").add(field.newDfp("-1")),     // 1 + (-1) = 0
             field.newDfp("0"),
             0, "Add #2");

        test(field.newDfp("-1").add(field.newDfp("1")),     // (-1) + 1 = 0
             field.newDfp("0"),
             0, "Add #3");

        test(field.newDfp("-1").add(field.newDfp("-1")),     // (-1) + (-1) = -2
             field.newDfp("-2"),
             0, "Add #4");

        // rounding mode is round half even

        test(field.newDfp("1").add(field.newDfp("1e-16")),     // rounding on add
             field.newDfp("1.0000000000000001"),
             0, "Add #5");

        test(field.newDfp("1").add(field.newDfp("1e-17")),     // rounding on add
             field.newDfp("1"),
             DfpField.FLAG_INEXACT, "Add #6");

        test(field.newDfp("0.90999999999999999999").add(field.newDfp("0.1")),     // rounding on add
             field.newDfp("1.01"),
             DfpField.FLAG_INEXACT, "Add #7");

        test(field.newDfp(".10000000000000005000").add(field.newDfp(".9")),     // rounding on add
             field.newDfp("1."),
             DfpField.FLAG_INEXACT, "Add #8");

        test(field.newDfp(".10000000000000015000").add(field.newDfp(".9")),     // rounding on add
             field.newDfp("1.0000000000000002"),
             DfpField.FLAG_INEXACT, "Add #9");

        test(field.newDfp(".10000000000000014999").add(field.newDfp(".9")),     // rounding on add
             field.newDfp("1.0000000000000001"),
             DfpField.FLAG_INEXACT, "Add #10");

        test(field.newDfp(".10000000000000015001").add(field.newDfp(".9")),     // rounding on add
             field.newDfp("1.0000000000000002"),
             DfpField.FLAG_INEXACT, "Add #11");

        test(field.newDfp(".11111111111111111111").add(field.newDfp("11.1111111111111111")), // rounding on add
             field.newDfp("11.22222222222222222222"),
             DfpField.FLAG_INEXACT, "Add #12");

        test(field.newDfp(".11111111111111111111").add(field.newDfp("1111111111111111.1111")), // rounding on add
             field.newDfp("1111111111111111.2222"),
             DfpField.FLAG_INEXACT, "Add #13");

        test(field.newDfp(".11111111111111111111").add(field.newDfp("11111111111111111111")), // rounding on add
             field.newDfp("11111111111111111111"),
             DfpField.FLAG_INEXACT, "Add #14");

        test(field.newDfp("9.9999999999999999999e131071").add(field.newDfp("-1e131052")), // overflow on add
             field.newDfp("9.9999999999999999998e131071"),
             0, "Add #15");

        test(field.newDfp("9.9999999999999999999e131071").add(field.newDfp("1e131052")), // overflow on add
             pinf,
             DfpField.FLAG_OVERFLOW, "Add #16");

        test(field.newDfp("-9.9999999999999999999e131071").add(field.newDfp("-1e131052")), // overflow on add
             ninf,
             DfpField.FLAG_OVERFLOW, "Add #17");

        test(field.newDfp("-9.9999999999999999999e131071").add(field.newDfp("1e131052")), // overflow on add
             field.newDfp("-9.9999999999999999998e131071"),
             0, "Add #18");

        test(field.newDfp("1e-131072").add(field.newDfp("1e-131072")), // underflow on add
             field.newDfp("2e-131072"),
             0, "Add #19");

        test(field.newDfp("1.0000000000000001e-131057").add(field.newDfp("-1e-131057")), // underflow on add
             field.newDfp("1e-131073"),
             DfpField.FLAG_UNDERFLOW, "Add #20");

        test(field.newDfp("1.1e-131072").add(field.newDfp("-1e-131072")), // underflow on add
             field.newDfp("1e-131073"),
             DfpField.FLAG_UNDERFLOW, "Add #21");

        test(field.newDfp("1.0000000000000001e-131072").add(field.newDfp("-1e-131072")), // underflow on add
             field.newDfp("1e-131088"),
             DfpField.FLAG_UNDERFLOW, "Add #22");

        test(field.newDfp("1.0000000000000001e-131078").add(field.newDfp("-1e-131078")), // underflow on add
             field.newDfp("0"),
             DfpField.FLAG_UNDERFLOW, "Add #23");

        test(field.newDfp("1.0").add(field.newDfp("-1e-20")), // loss of precision on alignment?
             field.newDfp("0.99999999999999999999"),
             0, "Add #23.1");

        test(field.newDfp("-0.99999999999999999999").add(field.newDfp("1")), // proper normalization?
             field.newDfp("0.00000000000000000001"),
             0, "Add #23.2");

        test(field.newDfp("1").add(field.newDfp("0")), // adding zeros
             field.newDfp("1"),
             0, "Add #24");

        test(field.newDfp("0").add(field.newDfp("0")), // adding zeros
             field.newDfp("0"),
             0, "Add #25");

        test(field.newDfp("-0").add(field.newDfp("0")), // adding zeros
             field.newDfp("0"),
             0, "Add #26");

        test(field.newDfp("0").add(field.newDfp("-0")), // adding zeros
             field.newDfp("0"),
             0, "Add #27");

        test(field.newDfp("-0").add(field.newDfp("-0")), // adding zeros
             field.newDfp("-0"),
             0, "Add #28");

        test(field.newDfp("1e-20").add(field.newDfp("0")), // adding zeros
             field.newDfp("1e-20"),
             0, "Add #29");

        test(field.newDfp("1e-40").add(field.newDfp("0")), // adding zeros
             field.newDfp("1e-40"),
             0, "Add #30");

        test(pinf.add(ninf), // adding infinities
             nan,
             DfpField.FLAG_INVALID, "Add #31");

        test(ninf.add(pinf), // adding infinities
             nan,
             DfpField.FLAG_INVALID, "Add #32");

        test(ninf.add(ninf), // adding infinities
             ninf,
             0, "Add #33");

        test(pinf.add(pinf), // adding infinities
             pinf,
             0, "Add #34");

        test(pinf.add(field.newDfp("0")), // adding infinities
             pinf,
             0, "Add #35");

        test(pinf.add(field.newDfp("-1e131071")), // adding infinities
             pinf,
             0, "Add #36");

        test(pinf.add(field.newDfp("1e131071")), // adding infinities
             pinf,
             0, "Add #37");

        test(field.newDfp("0").add(pinf), // adding infinities
             pinf,
             0, "Add #38");

        test(field.newDfp("-1e131071").add(pinf), // adding infinities
             pinf,
             0, "Add #39");

        test(field.newDfp("1e131071").add(pinf), // adding infinities
             pinf,
             0, "Add #40");

        test(ninf.add(field.newDfp("0")), // adding infinities
             ninf,
             0, "Add #41");

        test(ninf.add(field.newDfp("-1e131071")), // adding infinities
             ninf,
             0, "Add #42");

        test(ninf.add(field.newDfp("1e131071")), // adding infinities
             ninf,
             0, "Add #43");

        test(field.newDfp("0").add(ninf), // adding infinities
             ninf,
             0, "Add #44");

        test(field.newDfp("-1e131071").add(ninf), // adding infinities
             ninf,
             0, "Add #45");

        test(field.newDfp("1e131071").add(ninf), // adding infinities
             ninf,
             0, "Add #46");

        test(field.newDfp("9.9999999999999999999e131071").add(field.newDfp("5e131051")),  // overflow
             pinf,
             DfpField.FLAG_OVERFLOW, "Add #47");

        test(field.newDfp("9.9999999999999999999e131071").add(field.newDfp("4.9999999999999999999e131051")),  // overflow
             field.newDfp("9.9999999999999999999e131071"),
             DfpField.FLAG_INEXACT, "Add #48");

        test(nan.add(field.newDfp("1")),
             nan,
             0, "Add #49");

        test(field.newDfp("1").add(nan),
             nan,
             0, "Add #50");

        test(field.newDfp("12345678123456781234").add(field.newDfp("0.12345678123456781234")),
             field.newDfp("12345678123456781234"),
             DfpField.FLAG_INEXACT, "Add #51");

        test(field.newDfp("12345678123456781234").add(field.newDfp("123.45678123456781234")),
             field.newDfp("12345678123456781357"),
             DfpField.FLAG_INEXACT, "Add #52");

        test(field.newDfp("123.45678123456781234").add(field.newDfp("12345678123456781234")),
             field.newDfp("12345678123456781357"),
             DfpField.FLAG_INEXACT, "Add #53");

        test(field.newDfp("12345678123456781234").add(field.newDfp(".00001234567812345678")),
             field.newDfp("12345678123456781234"),
             DfpField.FLAG_INEXACT, "Add #54");

        test(field.newDfp("12345678123456781234").add(field.newDfp(".00000000123456781234")),
             field.newDfp("12345678123456781234"),
             DfpField.FLAG_INEXACT, "Add #55");

        test(field.newDfp("-0").add(field.newDfp("-0")),
             field.newDfp("-0"),
             0, "Add #56");

        test(field.newDfp("0").add(field.newDfp("-0")),
             field.newDfp("0"),
             0, "Add #57");

        test(field.newDfp("-0").add(field.newDfp("0")),
             field.newDfp("0"),
             0, "Add #58");

        test(field.newDfp("0").add(field.newDfp("0")),
             field.newDfp("0"),
             0, "Add #59");
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Test comparisons

    // utility function to help test comparisons
    private void cmptst(Dfp a, Dfp b, String op, boolean result, double num)
    {
        if (op == "equal")
            if (a.equals(b) != result)
                fail("assersion failed.  "+op+" compare #"+num);

        if (op == "unequal")
            if (a.unequal(b) != result)
                fail("assersion failed.  "+op+" compare #"+num);

        if (op == "lessThan")
            if (a.lessThan(b) != result)
                fail("assersion failed.  "+op+" compare #"+num);

        if (op == "greaterThan")
            if (a.greaterThan(b) != result)
                fail("assersion failed.  "+op+" compare #"+num);
    }

    @Test
    void testCompare()
    {
        // test equal() comparison
        // check zero vs. zero
        field.clearIEEEFlags();

        cmptst(field.newDfp("0"), field.newDfp("0"), "equal", true, 1);         // 0 == 0
        cmptst(field.newDfp("0"), field.newDfp("-0"), "equal", true, 2);        // 0 == -0
        cmptst(field.newDfp("-0"), field.newDfp("-0"), "equal", true, 3);       // -0 == -0
        cmptst(field.newDfp("-0"), field.newDfp("0"), "equal", true, 4);        // -0 == 0

        // check zero vs normal numbers

        cmptst(field.newDfp("0"), field.newDfp("1"), "equal", false, 5);         // 0 == 1
        cmptst(field.newDfp("1"), field.newDfp("0"), "equal", false, 6);         // 1 == 0
        cmptst(field.newDfp("-1"), field.newDfp("0"), "equal", false, 7);        // -1 == 0
        cmptst(field.newDfp("0"), field.newDfp("-1"), "equal", false, 8);        // 0 == -1
        cmptst(field.newDfp("0"), field.newDfp("1e-131072"), "equal", false, 9); // 0 == 1e-131072
        // check flags
        if (field.getIEEEFlags() != 0)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        cmptst(field.newDfp("0"), field.newDfp("1e-131078"), "equal", false, 10); // 0 == 1e-131078

        // check flags  -- underflow should be set
        if (field.getIEEEFlags() != DfpField.FLAG_UNDERFLOW)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        field.clearIEEEFlags();

        cmptst(field.newDfp("0"), field.newDfp("1e+131071"), "equal", false, 11); // 0 == 1e+131071

        // check zero vs infinities

        cmptst(field.newDfp("0"), pinf, "equal", false, 12);    // 0 == pinf
        cmptst(field.newDfp("0"), ninf, "equal", false, 13);    // 0 == ninf
        cmptst(field.newDfp("-0"), pinf, "equal", false, 14);   // -0 == pinf
        cmptst(field.newDfp("-0"), ninf, "equal", false, 15);   // -0 == ninf
        cmptst(pinf, field.newDfp("0"), "equal", false, 16);    // pinf == 0
        cmptst(ninf, field.newDfp("0"), "equal", false, 17);    // ninf == 0
        cmptst(pinf, field.newDfp("-0"), "equal", false, 18);   // pinf == -0
        cmptst(ninf, field.newDfp("-0"), "equal", false, 19);   // ninf == -0
        cmptst(ninf, pinf, "equal", false, 19.10);     // ninf == pinf
        cmptst(pinf, ninf, "equal", false, 19.11);     // pinf == ninf
        cmptst(pinf, pinf, "equal", true, 19.12);     // pinf == pinf
        cmptst(ninf, ninf, "equal", true, 19.13);     // ninf == ninf

        // check some normal numbers
        cmptst(field.newDfp("1"), field.newDfp("1"), "equal", true, 20);   // 1 == 1
        cmptst(field.newDfp("1"), field.newDfp("-1"), "equal", false, 21);   // 1 == -1
        cmptst(field.newDfp("-1"), field.newDfp("-1"), "equal", true, 22);   // -1 == -1
        cmptst(field.newDfp("1"), field.newDfp("1.0000000000000001"), "equal", false, 23);   // 1 == 1.0000000000000001

        // The tests below checks to ensure that comparisons don't set FLAG_INEXACT
        // 100000 == 1.0000000000000001
        cmptst(field.newDfp("1e20"), field.newDfp("1.0000000000000001"), "equal", false, 24);
        if (field.getIEEEFlags() != 0)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        cmptst(field.newDfp("0.000001"), field.newDfp("1e-6"), "equal", true, 25);

        // check some nans -- nans shouldnt equal anything

        cmptst(snan, snan, "equal", false, 27);
        cmptst(qnan, qnan, "equal", false, 28);
        cmptst(snan, qnan, "equal", false, 29);
        cmptst(qnan, snan, "equal", false, 30);
        cmptst(qnan, field.newDfp("0"), "equal", false, 31);
        cmptst(snan, field.newDfp("0"), "equal", false, 32);
        cmptst(field.newDfp("0"), snan, "equal", false, 33);
        cmptst(field.newDfp("0"), qnan, "equal", false, 34);
        cmptst(qnan, pinf, "equal", false, 35);
        cmptst(snan, pinf, "equal", false, 36);
        cmptst(pinf, snan, "equal", false, 37);
        cmptst(pinf, qnan, "equal", false, 38);
        cmptst(qnan, ninf, "equal", false, 39);
        cmptst(snan, ninf, "equal", false, 40);
        cmptst(ninf, snan, "equal", false, 41);
        cmptst(ninf, qnan, "equal", false, 42);
        cmptst(qnan, field.newDfp("-1"), "equal", false, 43);
        cmptst(snan, field.newDfp("-1"), "equal", false, 44);
        cmptst(field.newDfp("-1"), snan, "equal", false, 45);
        cmptst(field.newDfp("-1"), qnan, "equal", false, 46);
        cmptst(qnan, field.newDfp("1"), "equal", false, 47);
        cmptst(snan, field.newDfp("1"), "equal", false, 48);
        cmptst(field.newDfp("1"), snan, "equal", false, 49);
        cmptst(field.newDfp("1"), qnan, "equal", false, 50);
        cmptst(snan.negate(), snan, "equal", false, 51);
        cmptst(qnan.negate(), qnan, "equal", false, 52);

        //
        // Tests for un equal  -- do it all over again
        //

        cmptst(field.newDfp("0"), field.newDfp("0"), "unequal", false, 1);         // 0 == 0
        cmptst(field.newDfp("0"), field.newDfp("-0"), "unequal", false, 2);        // 0 == -0
        cmptst(field.newDfp("-0"), field.newDfp("-0"), "unequal", false, 3);       // -0 == -0
        cmptst(field.newDfp("-0"), field.newDfp("0"), "unequal", false, 4);        // -0 == 0

        // check zero vs normal numbers

        cmptst(field.newDfp("0"), field.newDfp("1"), "unequal", true, 5);         // 0 == 1
        cmptst(field.newDfp("1"), field.newDfp("0"), "unequal", true, 6);         // 1 == 0
        cmptst(field.newDfp("-1"), field.newDfp("0"), "unequal", true, 7);        // -1 == 0
        cmptst(field.newDfp("0"), field.newDfp("-1"), "unequal", true, 8);        // 0 == -1
        cmptst(field.newDfp("0"), field.newDfp("1e-131072"), "unequal", true, 9); // 0 == 1e-131072
        // check flags
        if (field.getIEEEFlags() != 0)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        cmptst(field.newDfp("0"), field.newDfp("1e-131078"), "unequal", true, 10); // 0 == 1e-131078

        // check flags  -- underflow should be set
        if (field.getIEEEFlags() != DfpField.FLAG_UNDERFLOW)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        field.clearIEEEFlags();

        cmptst(field.newDfp("0"), field.newDfp("1e+131071"), "unequal", true, 11); // 0 == 1e+131071

        // check zero vs infinities

        cmptst(field.newDfp("0"), pinf, "unequal", true, 12);    // 0 == pinf
        cmptst(field.newDfp("0"), ninf, "unequal", true, 13);    // 0 == ninf
        cmptst(field.newDfp("-0"), pinf, "unequal", true, 14);   // -0 == pinf
        cmptst(field.newDfp("-0"), ninf, "unequal", true, 15);   // -0 == ninf
        cmptst(pinf, field.newDfp("0"), "unequal", true, 16);    // pinf == 0
        cmptst(ninf, field.newDfp("0"), "unequal", true, 17);    // ninf == 0
        cmptst(pinf, field.newDfp("-0"), "unequal", true, 18);   // pinf == -0
        cmptst(ninf, field.newDfp("-0"), "unequal", true, 19);   // ninf == -0
        cmptst(ninf, pinf, "unequal", true, 19.10);     // ninf == pinf
        cmptst(pinf, ninf, "unequal", true, 19.11);     // pinf == ninf
        cmptst(pinf, pinf, "unequal", false, 19.12);     // pinf == pinf
        cmptst(ninf, ninf, "unequal", false, 19.13);     // ninf == ninf

        // check some normal numbers
        cmptst(field.newDfp("1"), field.newDfp("1"), "unequal", false, 20);   // 1 == 1
        cmptst(field.newDfp("1"), field.newDfp("-1"), "unequal", true, 21);   // 1 == -1
        cmptst(field.newDfp("-1"), field.newDfp("-1"), "unequal", false, 22);   // -1 == -1
        cmptst(field.newDfp("1"), field.newDfp("1.0000000000000001"), "unequal", true, 23);   // 1 == 1.0000000000000001

        // The tests below checks to ensure that comparisons don't set FLAG_INEXACT
        // 100000 == 1.0000000000000001
        cmptst(field.newDfp("1e20"), field.newDfp("1.0000000000000001"), "unequal", true, 24);
        if (field.getIEEEFlags() != 0)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        cmptst(field.newDfp("0.000001"), field.newDfp("1e-6"), "unequal", false, 25);

        // check some nans -- nans shouldnt be unequal to anything

        cmptst(snan, snan, "unequal", false, 27);
        cmptst(qnan, qnan, "unequal", false, 28);
        cmptst(snan, qnan, "unequal", false, 29);
        cmptst(qnan, snan, "unequal", false, 30);
        cmptst(qnan, field.newDfp("0"), "unequal", false, 31);
        cmptst(snan, field.newDfp("0"), "unequal", false, 32);
        cmptst(field.newDfp("0"), snan, "unequal", false, 33);
        cmptst(field.newDfp("0"), qnan, "unequal", false, 34);
        cmptst(qnan, pinf, "unequal", false, 35);
        cmptst(snan, pinf, "unequal", false, 36);
        cmptst(pinf, snan, "unequal", false, 37);
        cmptst(pinf, qnan, "unequal", false, 38);
        cmptst(qnan, ninf, "unequal", false, 39);
        cmptst(snan, ninf, "unequal", false, 40);
        cmptst(ninf, snan, "unequal", false, 41);
        cmptst(ninf, qnan, "unequal", false, 42);
        cmptst(qnan, field.newDfp("-1"), "unequal", false, 43);
        cmptst(snan, field.newDfp("-1"), "unequal", false, 44);
        cmptst(field.newDfp("-1"), snan, "unequal", false, 45);
        cmptst(field.newDfp("-1"), qnan, "unequal", false, 46);
        cmptst(qnan, field.newDfp("1"), "unequal", false, 47);
        cmptst(snan, field.newDfp("1"), "unequal", false, 48);
        cmptst(field.newDfp("1"), snan, "unequal", false, 49);
        cmptst(field.newDfp("1"), qnan, "unequal", false, 50);
        cmptst(snan.negate(), snan, "unequal", false, 51);
        cmptst(qnan.negate(), qnan, "unequal", false, 52);

        if (field.getIEEEFlags() != 0)
            fail("assersion failed.  compare unequal flags = "+field.getIEEEFlags());

        //
        // Tests for lessThan  -- do it all over again
        //

        cmptst(field.newDfp("0"), field.newDfp("0"), "lessThan", false, 1);         // 0 < 0
        cmptst(field.newDfp("0"), field.newDfp("-0"), "lessThan", false, 2);        // 0 < -0
        cmptst(field.newDfp("-0"), field.newDfp("-0"), "lessThan", false, 3);       // -0 < -0
        cmptst(field.newDfp("-0"), field.newDfp("0"), "lessThan", false, 4);        // -0 < 0

        // check zero vs normal numbers

        cmptst(field.newDfp("0"), field.newDfp("1"), "lessThan", true, 5);         // 0 < 1
        cmptst(field.newDfp("1"), field.newDfp("0"), "lessThan", false, 6);         // 1 < 0
        cmptst(field.newDfp("-1"), field.newDfp("0"), "lessThan", true, 7);        // -1 < 0
        cmptst(field.newDfp("0"), field.newDfp("-1"), "lessThan", false, 8);        // 0 < -1
        cmptst(field.newDfp("0"), field.newDfp("1e-131072"), "lessThan", true, 9); // 0 < 1e-131072
        // check flags
        if (field.getIEEEFlags() != 0)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        cmptst(field.newDfp("0"), field.newDfp("1e-131078"), "lessThan", true, 10); // 0 < 1e-131078

        // check flags  -- underflow should be set
        if (field.getIEEEFlags() != DfpField.FLAG_UNDERFLOW)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());
        field.clearIEEEFlags();

        cmptst(field.newDfp("0"), field.newDfp("1e+131071"), "lessThan", true, 11); // 0 < 1e+131071

        // check zero vs infinities

        cmptst(field.newDfp("0"), pinf, "lessThan", true, 12);    // 0 < pinf
        cmptst(field.newDfp("0"), ninf, "lessThan", false, 13);    // 0 < ninf
        cmptst(field.newDfp("-0"), pinf, "lessThan", true, 14);   // -0 < pinf
        cmptst(field.newDfp("-0"), ninf, "lessThan", false, 15);   // -0 < ninf
        cmptst(pinf, field.newDfp("0"), "lessThan", false, 16);    // pinf < 0
        cmptst(ninf, field.newDfp("0"), "lessThan", true, 17);    // ninf < 0
        cmptst(pinf, field.newDfp("-0"), "lessThan", false, 18);   // pinf < -0
        cmptst(ninf, field.newDfp("-0"), "lessThan", true, 19);   // ninf < -0
        cmptst(ninf, pinf, "lessThan", true, 19.10);     // ninf < pinf
        cmptst(pinf, ninf, "lessThan", false, 19.11);     // pinf < ninf
        cmptst(pinf, pinf, "lessThan", false, 19.12);     // pinf < pinf
        cmptst(ninf, ninf, "lessThan", false, 19.13);     // ninf < ninf

        // check some normal numbers
        cmptst(field.newDfp("1"), field.newDfp("1"), "lessThan", false, 20);   // 1 < 1
        cmptst(field.newDfp("1"), field.newDfp("-1"), "lessThan", false, 21);   // 1 < -1
        cmptst(field.newDfp("-1"), field.newDfp("-1"), "lessThan", false, 22);   // -1 < -1
        cmptst(field.newDfp("1"), field.newDfp("1.0000000000000001"), "lessThan", true, 23);   // 1 < 1.0000000000000001

        // The tests below checks to ensure that comparisons don't set FLAG_INEXACT
        // 100000 < 1.0000000000000001
        cmptst(field.newDfp("1e20"), field.newDfp("1.0000000000000001"), "lessThan", false, 24);
        if (field.getIEEEFlags() != 0)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        cmptst(field.newDfp("0.000001"), field.newDfp("1e-6"), "lessThan", false, 25);

        // check some nans -- nans shouldnt be lessThan to anything
        cmptst(snan, snan, "lessThan", false, 27);
        cmptst(qnan, qnan, "lessThan", false, 28);
        cmptst(snan, qnan, "lessThan", false, 29);
        cmptst(qnan, snan, "lessThan", false, 30);
        cmptst(qnan, field.newDfp("0"), "lessThan", false, 31);
        cmptst(snan, field.newDfp("0"), "lessThan", false, 32);
        cmptst(field.newDfp("0"), snan, "lessThan", false, 33);
        cmptst(field.newDfp("0"), qnan, "lessThan", false, 34);
        cmptst(qnan, pinf, "lessThan", false, 35);
        cmptst(snan, pinf, "lessThan", false, 36);
        cmptst(pinf, snan, "lessThan", false, 37);
        cmptst(pinf, qnan, "lessThan", false, 38);
        cmptst(qnan, ninf, "lessThan", false, 39);
        cmptst(snan, ninf, "lessThan", false, 40);
        cmptst(ninf, snan, "lessThan", false, 41);
        cmptst(ninf, qnan, "lessThan", false, 42);
        cmptst(qnan, field.newDfp("-1"), "lessThan", false, 43);
        cmptst(snan, field.newDfp("-1"), "lessThan", false, 44);
        cmptst(field.newDfp("-1"), snan, "lessThan", false, 45);
        cmptst(field.newDfp("-1"), qnan, "lessThan", false, 46);
        cmptst(qnan, field.newDfp("1"), "lessThan", false, 47);
        cmptst(snan, field.newDfp("1"), "lessThan", false, 48);
        cmptst(field.newDfp("1"), snan, "lessThan", false, 49);
        cmptst(field.newDfp("1"), qnan, "lessThan", false, 50);
        cmptst(snan.negate(), snan, "lessThan", false, 51);
        cmptst(qnan.negate(), qnan, "lessThan", false, 52);

        //lessThan compares with nans should raise FLAG_INVALID
        if (field.getIEEEFlags() != DfpField.FLAG_INVALID)
            fail("assersion failed.  compare lessThan flags = "+field.getIEEEFlags());
        field.clearIEEEFlags();

        //
        // Tests for greaterThan  -- do it all over again
        //

        cmptst(field.newDfp("0"), field.newDfp("0"), "greaterThan", false, 1);         // 0 > 0
        cmptst(field.newDfp("0"), field.newDfp("-0"), "greaterThan", false, 2);        // 0 > -0
        cmptst(field.newDfp("-0"), field.newDfp("-0"), "greaterThan", false, 3);       // -0 > -0
        cmptst(field.newDfp("-0"), field.newDfp("0"), "greaterThan", false, 4);        // -0 > 0

        // check zero vs normal numbers

        cmptst(field.newDfp("0"), field.newDfp("1"), "greaterThan", false, 5);         // 0 > 1
        cmptst(field.newDfp("1"), field.newDfp("0"), "greaterThan", true, 6);         // 1 > 0
        cmptst(field.newDfp("-1"), field.newDfp("0"), "greaterThan", false, 7);        // -1 > 0
        cmptst(field.newDfp("0"), field.newDfp("-1"), "greaterThan", true, 8);        // 0 > -1
        cmptst(field.newDfp("0"), field.newDfp("1e-131072"), "greaterThan", false, 9); // 0 > 1e-131072
        // check flags
        if (field.getIEEEFlags() != 0)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        cmptst(field.newDfp("0"), field.newDfp("1e-131078"), "greaterThan", false, 10); // 0 > 1e-131078

        // check flags  -- underflow should be set
        if (field.getIEEEFlags() != DfpField.FLAG_UNDERFLOW)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());
        field.clearIEEEFlags();

        cmptst(field.newDfp("0"), field.newDfp("1e+131071"), "greaterThan", false, 11); // 0 > 1e+131071

        // check zero vs infinities

        cmptst(field.newDfp("0"), pinf, "greaterThan", false, 12);    // 0 > pinf
        cmptst(field.newDfp("0"), ninf, "greaterThan", true, 13);    // 0 > ninf
        cmptst(field.newDfp("-0"), pinf, "greaterThan", false, 14);   // -0 > pinf
        cmptst(field.newDfp("-0"), ninf, "greaterThan", true, 15);   // -0 > ninf
        cmptst(pinf, field.newDfp("0"), "greaterThan", true, 16);    // pinf > 0
        cmptst(ninf, field.newDfp("0"), "greaterThan", false, 17);    // ninf > 0
        cmptst(pinf, field.newDfp("-0"), "greaterThan", true, 18);   // pinf > -0
        cmptst(ninf, field.newDfp("-0"), "greaterThan", false, 19);   // ninf > -0
        cmptst(ninf, pinf, "greaterThan", false, 19.10);     // ninf > pinf
        cmptst(pinf, ninf, "greaterThan", true, 19.11);     // pinf > ninf
        cmptst(pinf, pinf, "greaterThan", false, 19.12);     // pinf > pinf
        cmptst(ninf, ninf, "greaterThan", false, 19.13);     // ninf > ninf

        // check some normal numbers
        cmptst(field.newDfp("1"), field.newDfp("1"), "greaterThan", false, 20);   // 1 > 1
        cmptst(field.newDfp("1"), field.newDfp("-1"), "greaterThan", true, 21);   // 1 > -1
        cmptst(field.newDfp("-1"), field.newDfp("-1"), "greaterThan", false, 22);   // -1 > -1
        cmptst(field.newDfp("1"), field.newDfp("1.0000000000000001"), "greaterThan", false, 23);   // 1 > 1.0000000000000001

        // The tests below checks to ensure that comparisons don't set FLAG_INEXACT
        // 100000 > 1.0000000000000001
        cmptst(field.newDfp("1e20"), field.newDfp("1.0000000000000001"), "greaterThan", true, 24);
        if (field.getIEEEFlags() != 0)
            fail("assersion failed.  compare flags = "+field.getIEEEFlags());

        cmptst(field.newDfp("0.000001"), field.newDfp("1e-6"), "greaterThan", false, 25);

        // check some nans -- nans shouldnt be greaterThan to anything
        cmptst(snan, snan, "greaterThan", false, 27);
        cmptst(qnan, qnan, "greaterThan", false, 28);
        cmptst(snan, qnan, "greaterThan", false, 29);
        cmptst(qnan, snan, "greaterThan", false, 30);
        cmptst(qnan, field.newDfp("0"), "greaterThan", false, 31);
        cmptst(snan, field.newDfp("0"), "greaterThan", false, 32);
        cmptst(field.newDfp("0"), snan, "greaterThan", false, 33);
        cmptst(field.newDfp("0"), qnan, "greaterThan", false, 34);
        cmptst(qnan, pinf, "greaterThan", false, 35);
        cmptst(snan, pinf, "greaterThan", false, 36);
        cmptst(pinf, snan, "greaterThan", false, 37);
        cmptst(pinf, qnan, "greaterThan", false, 38);
        cmptst(qnan, ninf, "greaterThan", false, 39);
        cmptst(snan, ninf, "greaterThan", false, 40);
        cmptst(ninf, snan, "greaterThan", false, 41);
        cmptst(ninf, qnan, "greaterThan", false, 42);
        cmptst(qnan, field.newDfp("-1"), "greaterThan", false, 43);
        cmptst(snan, field.newDfp("-1"), "greaterThan", false, 44);
        cmptst(field.newDfp("-1"), snan, "greaterThan", false, 45);
        cmptst(field.newDfp("-1"), qnan, "greaterThan", false, 46);
        cmptst(qnan, field.newDfp("1"), "greaterThan", false, 47);
        cmptst(snan, field.newDfp("1"), "greaterThan", false, 48);
        cmptst(field.newDfp("1"), snan, "greaterThan", false, 49);
        cmptst(field.newDfp("1"), qnan, "greaterThan", false, 50);
        cmptst(snan.negate(), snan, "greaterThan", false, 51);
        cmptst(qnan.negate(), qnan, "greaterThan", false, 52);

        //greaterThan compares with nans should raise FLAG_INVALID
        if (field.getIEEEFlags() != DfpField.FLAG_INVALID)
            fail("assersion failed.  compare greaterThan flags = "+field.getIEEEFlags());
        field.clearIEEEFlags();
    }

    //
    // Test multiplication
    //
    @Test
    void testMultiply()
    {
        test(field.newDfp("1").multiply(field.newDfp("1")),      // Basic tests   1*1 = 1
             field.newDfp("1"),
             0, "Multiply #1");

        test(field.newDfp("1").multiply(1),             // Basic tests   1*1 = 1
             field.newDfp("1"),
             0, "Multiply #2");

        test(field.newDfp("-1").multiply(field.newDfp("1")),     // Basic tests   -1*1 = -1
             field.newDfp("-1"),
             0, "Multiply #3");

        test(field.newDfp("-1").multiply(1),            // Basic tests   -1*1 = -1
             field.newDfp("-1"),
             0, "Multiply #4");

        // basic tests with integers
        test(field.newDfp("2").multiply(field.newDfp("3")),
             field.newDfp("6"),
             0, "Multiply #5");

        test(field.newDfp("2").multiply(3),
             field.newDfp("6"),
             0, "Multiply #6");

        test(field.newDfp("-2").multiply(field.newDfp("3")),
             field.newDfp("-6"),
             0, "Multiply #7");

        test(field.newDfp("-2").multiply(3),
             field.newDfp("-6"),
             0, "Multiply #8");

        test(field.newDfp("2").multiply(field.newDfp("-3")),
             field.newDfp("-6"),
             0, "Multiply #9");

        test(field.newDfp("-2").multiply(field.newDfp("-3")),
             field.newDfp("6"),
             0, "Multiply #10");

        //multiply by zero

        test(field.newDfp("-2").multiply(field.newDfp("0")),
             field.newDfp("-0"),
             0, "Multiply #11");

        test(field.newDfp("-2").multiply(0),
             field.newDfp("-0"),
             0, "Multiply #12");

        test(field.newDfp("2").multiply(field.newDfp("0")),
             field.newDfp("0"),
             0, "Multiply #13");

        test(field.newDfp("2").multiply(0),
             field.newDfp("0"),
             0, "Multiply #14");

        test(field.newDfp("2").multiply(pinf),
             pinf,
             0, "Multiply #15");

        test(field.newDfp("2").multiply(ninf),
             ninf,
             0, "Multiply #16");

        test(field.newDfp("-2").multiply(pinf),
             ninf,
             0, "Multiply #17");

        test(field.newDfp("-2").multiply(ninf),
             pinf,
             0, "Multiply #18");

        test(ninf.multiply(field.newDfp("-2")),
             pinf,
             0, "Multiply #18.1");

        test(field.newDfp("5e131071").twice(),
             pinf,
             DfpField.FLAG_OVERFLOW, "Multiply #19");

        test(field.newDfp("5e131071").multiply(field.newDfp("1.999999999999999")),
             field.newDfp("9.9999999999999950000e131071"),
             0, "Multiply #20");

        test(field.newDfp("-5e131071").twice(),
             ninf,
             DfpField.FLAG_OVERFLOW, "Multiply #22");

        test(field.newDfp("-5e131071").multiply(field.newDfp("1.999999999999999")),
             field.newDfp("-9.9999999999999950000e131071"),
             0, "Multiply #23");

        test(field.newDfp("1e-65539").multiply(field.newDfp("1e-65539")),
             field.newDfp("1e-131078"),
             DfpField.FLAG_UNDERFLOW, "Multiply #24");

        test(field.newDfp("1").multiply(nan),
             nan,
             0, "Multiply #25");

        test(nan.multiply(field.newDfp("1")),
             nan,
             0, "Multiply #26");

        test(nan.multiply(pinf),
             nan,
             0, "Multiply #27");

        test(pinf.multiply(nan),
             nan,
             0, "Multiply #27");

        test(pinf.multiply(field.newDfp("0")),
             nan,
             DfpField.FLAG_INVALID, "Multiply #28");

        test(field.newDfp("0").multiply(pinf),
             nan,
             DfpField.FLAG_INVALID, "Multiply #29");

        test(pinf.multiply(pinf),
             pinf,
             0, "Multiply #30");

        test(ninf.multiply(pinf),
             ninf,
             0, "Multiply #31");

        test(pinf.multiply(ninf),
             ninf,
             0, "Multiply #32");

        test(ninf.multiply(ninf),
             pinf,
             0, "Multiply #33");

        test(pinf.multiply(1),
             pinf,
             0, "Multiply #34");

        test(pinf.multiply(0),
             nan,
             DfpField.FLAG_INVALID, "Multiply #35");

        test(nan.multiply(1),
             nan,
             0, "Multiply #36");

        test(field.newDfp("1").multiply(10000),
             field.newDfp("10000"),
             0, "Multiply #37");

        test(field.newDfp("2").multiply(1000000),
             field.newDfp("2000000"),
             0, "Multiply #38");

        test(field.newDfp("1").multiply(-1),
             field.newDfp("-1"),
             0, "Multiply #39");
    }

    @Test
    void testDivide()
    {
        test(field.newDfp("1").divide(nan),      // divide by NaN = NaN
             nan,
             0, "Divide #1");

        test(nan.divide(field.newDfp("1")),      // NaN / number = NaN
             nan,
             0, "Divide #2");

        test(pinf.divide(field.newDfp("1")),
             pinf,
             0, "Divide #3");

        test(pinf.divide(field.newDfp("-1")),
             ninf,
             0, "Divide #4");

        test(pinf.divide(pinf),
             nan,
             DfpField.FLAG_INVALID, "Divide #5");

        test(ninf.divide(pinf),
             nan,
             DfpField.FLAG_INVALID, "Divide #6");

        test(pinf.divide(ninf),
             nan,
             DfpField.FLAG_INVALID, "Divide #7");

        test(ninf.divide(ninf),
             nan,
             DfpField.FLAG_INVALID, "Divide #8");

        test(field.newDfp("0").divide(field.newDfp("0")),
             nan,
             DfpField.FLAG_DIV_ZERO, "Divide #9");

        test(field.newDfp("1").divide(field.newDfp("0")),
             pinf,
             DfpField.FLAG_DIV_ZERO, "Divide #10");

        test(field.newDfp("1").divide(field.newDfp("-0")),
             ninf,
             DfpField.FLAG_DIV_ZERO, "Divide #11");

        test(field.newDfp("-1").divide(field.newDfp("0")),
             ninf,
             DfpField.FLAG_DIV_ZERO, "Divide #12");

        test(field.newDfp("-1").divide(field.newDfp("-0")),
             pinf,
             DfpField.FLAG_DIV_ZERO, "Divide #13");

        test(field.newDfp("1").divide(field.newDfp("3")),
             field.newDfp("0.33333333333333333333"),
             DfpField.FLAG_INEXACT, "Divide #14");

        test(field.newDfp("1").divide(field.newDfp("6")),
             field.newDfp("0.16666666666666666667"),
             DfpField.FLAG_INEXACT, "Divide #15");

        test(field.newDfp("10").divide(field.newDfp("6")),
             field.newDfp("1.6666666666666667"),
             DfpField.FLAG_INEXACT, "Divide #16");

        test(field.newDfp("100").divide(field.newDfp("6")),
             field.newDfp("16.6666666666666667"),
             DfpField.FLAG_INEXACT, "Divide #17");

        test(field.newDfp("1000").divide(field.newDfp("6")),
             field.newDfp("166.6666666666666667"),
             DfpField.FLAG_INEXACT, "Divide #18");

        test(field.newDfp("10000").divide(field.newDfp("6")),
             field.newDfp("1666.6666666666666667"),
             DfpField.FLAG_INEXACT, "Divide #19");

        test(field.newDfp("1").divide(field.newDfp("1")),
             field.newDfp("1"),
             0, "Divide #20");

        test(field.newDfp("1").divide(field.newDfp("-1")),
             field.newDfp("-1"),
             0, "Divide #21");

        test(field.newDfp("-1").divide(field.newDfp("1")),
             field.newDfp("-1"),
             0, "Divide #22");

        test(field.newDfp("-1").divide(field.newDfp("-1")),
             field.newDfp("1"),
             0, "Divide #23");

        test(field.newDfp("1e-65539").divide(field.newDfp("1e65539")),
             field.newDfp("1e-131078"),
             DfpField.FLAG_UNDERFLOW, "Divide #24");

        test(field.newDfp("1e65539").divide(field.newDfp("1e-65539")),
             pinf,
             DfpField.FLAG_OVERFLOW, "Divide #24");

        test(field.newDfp("2").divide(field.newDfp("1.5")),     // test trial-divisor too high
             field.newDfp("1.3333333333333333"),
             DfpField.FLAG_INEXACT, "Divide #25");

        test(field.newDfp("2").divide(pinf),
             field.newDfp("0"),
             0, "Divide #26");

        test(field.newDfp("2").divide(ninf),
             field.newDfp("-0"),
             0, "Divide #27");

        test(field.newDfp("0").divide(field.newDfp("1")),
             field.newDfp("0"),
             0, "Divide #28");
    }

    @Test
    void testReciprocal()
    {
        test(nan.reciprocal(),
             nan,
             0, "Reciprocal #1");

        test(field.newDfp("0").reciprocal(),
             pinf,
             DfpField.FLAG_DIV_ZERO, "Reciprocal #2");

        test(field.newDfp("-0").reciprocal(),
             ninf,
             DfpField.FLAG_DIV_ZERO, "Reciprocal #3");

        test(field.newDfp("3").reciprocal(),
             field.newDfp("0.33333333333333333333"),
             DfpField.FLAG_INEXACT, "Reciprocal #4");

        test(field.newDfp("6").reciprocal(),
             field.newDfp("0.16666666666666666667"),
             DfpField.FLAG_INEXACT, "Reciprocal #5");

        test(field.newDfp("1").reciprocal(),
             field.newDfp("1"),
             0, "Reciprocal #6");

        test(field.newDfp("-1").reciprocal(),
             field.newDfp("-1"),
             0, "Reciprocal #7");

        test(pinf.reciprocal(),
             field.newDfp("0"),
             0, "Reciprocal #8");

        test(ninf.reciprocal(),
             field.newDfp("-0"),
             0, "Reciprocal #9");
    }

    @Test
    void testDivideInt()
    {
        test(nan.divide(1),      // NaN / number = NaN
             nan,
             0, "DivideInt #1");

        test(pinf.divide(1),
             pinf,
             0, "DivideInt #2");

        test(field.newDfp("0").divide(0),
             nan,
             DfpField.FLAG_DIV_ZERO, "DivideInt #3");

        test(field.newDfp("1").divide(0),
             pinf,
             DfpField.FLAG_DIV_ZERO, "DivideInt #4");

        test(field.newDfp("-1").divide(0),
             ninf,
             DfpField.FLAG_DIV_ZERO, "DivideInt #5");

        test(field.newDfp("1").divide(3),
             field.newDfp("0.33333333333333333333"),
             DfpField.FLAG_INEXACT, "DivideInt #6");

        test(field.newDfp("1").divide(6),
             field.newDfp("0.16666666666666666667"),
             DfpField.FLAG_INEXACT, "DivideInt #7");

        test(field.newDfp("10").divide(6),
             field.newDfp("1.6666666666666667"),
             DfpField.FLAG_INEXACT, "DivideInt #8");

        test(field.newDfp("100").divide(6),
             field.newDfp("16.6666666666666667"),
             DfpField.FLAG_INEXACT, "DivideInt #9");

        test(field.newDfp("1000").divide(6),
             field.newDfp("166.6666666666666667"),
             DfpField.FLAG_INEXACT, "DivideInt #10");

        test(field.newDfp("10000").divide(6),
             field.newDfp("1666.6666666666666667"),
             DfpField.FLAG_INEXACT, "DivideInt #20");

        test(field.newDfp("1").divide(1),
             field.newDfp("1"),
             0, "DivideInt #21");

        test(field.newDfp("1e-131077").divide(10),
             field.newDfp("1e-131078"),
             DfpField.FLAG_UNDERFLOW, "DivideInt #22");

        test(field.newDfp("0").divide(1),
             field.newDfp("0"),
             0, "DivideInt #23");

        test(field.newDfp("1").divide(10000),
             nan,
             DfpField.FLAG_INVALID, "DivideInt #24");

        test(field.newDfp("1").divide(-1),
             nan,
             DfpField.FLAG_INVALID, "DivideInt #25");
    }

    @Test
    void testNextAfter()
    {
        test(field.newDfp("1").nextAfter(pinf),
             field.newDfp("1.0000000000000001"),
             0, "NextAfter #1");

        test(field.newDfp("1.0000000000000001").nextAfter(ninf),
             field.newDfp("1"),
             0, "NextAfter #1.5");

        test(field.newDfp("1").nextAfter(ninf),
             field.newDfp("0.99999999999999999999"),
             0, "NextAfter #2");

        test(field.newDfp("0.99999999999999999999").nextAfter(field.newDfp("2")),
             field.newDfp("1"),
             0, "NextAfter #3");

        test(field.newDfp("-1").nextAfter(ninf),
             field.newDfp("-1.0000000000000001"),
             0, "NextAfter #4");

        test(field.newDfp("-1").nextAfter(pinf),
             field.newDfp("-0.99999999999999999999"),
             0, "NextAfter #5");

        test(field.newDfp("-0.99999999999999999999").nextAfter(field.newDfp("-2")),
             field.newDfp("-1"),
             0, "NextAfter #6");

        test(field.newDfp("2").nextAfter(field.newDfp("2")),
             field.newDfp("2"),
             0, "NextAfter #7");

        test(field.newDfp("0").nextAfter(field.newDfp("0")),
             field.newDfp("0"),
             0, "NextAfter #8");

        test(field.newDfp("-2").nextAfter(field.newDfp("-2")),
             field.newDfp("-2"),
             0, "NextAfter #9");

        test(field.newDfp("0").nextAfter(field.newDfp("1")),
             field.newDfp("1e-131092"),
             DfpField.FLAG_UNDERFLOW, "NextAfter #10");

        test(field.newDfp("0").nextAfter(field.newDfp("-1")),
             field.newDfp("-1e-131092"),
             DfpField.FLAG_UNDERFLOW, "NextAfter #11");

        test(field.newDfp("-1e-131092").nextAfter(pinf),
             field.newDfp("-0"),
             DfpField.FLAG_UNDERFLOW|DfpField.FLAG_INEXACT, "Next After #12");

        test(field.newDfp("1e-131092").nextAfter(ninf),
             field.newDfp("0"),
             DfpField.FLAG_UNDERFLOW|DfpField.FLAG_INEXACT, "Next After #13");

        test(field.newDfp("9.9999999999999999999e131078").nextAfter(pinf),
             pinf,
             DfpField.FLAG_OVERFLOW|DfpField.FLAG_INEXACT, "Next After #14");
    }

    @Test
    void testToString()
    {
        assertEquals("Infinity", pinf.toString(), "toString #1");
        assertEquals("-Infinity", ninf.toString(), "toString #2");
        assertEquals("NaN", nan.toString(), "toString #3");
        assertEquals("NaN", field.newDfp((byte) 1, Dfp.QNAN).toString(), "toString #4");
        assertEquals("NaN", field.newDfp((byte) 1, Dfp.SNAN).toString(), "toString #5");
        assertEquals("1.2300000000000000e100", field.newDfp("1.23e100").toString(), "toString #6");
        assertEquals("-1.2300000000000000e100", field.newDfp("-1.23e100").toString(), "toString #7");
        assertEquals("12345678.1234", field.newDfp("12345678.1234").toString(), "toString #8");
        assertEquals("0.00001234", field.newDfp("0.00001234").toString(), "toString #9");
    }

    @Override
    @Test
    public void testRound()
    {
        field.setRoundingMode(DfpField.RoundingMode.ROUND_DOWN);

        // Round down
        test(field.newDfp("12345678901234567890").add(field.newDfp("0.9")),
             field.newDfp("12345678901234567890"),
             DfpField.FLAG_INEXACT, "Round #1");

        test(field.newDfp("12345678901234567890").add(field.newDfp("0.99999999")),
             field.newDfp("12345678901234567890"),
             DfpField.FLAG_INEXACT, "Round #2");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.99999999")),
             field.newDfp("-12345678901234567890"),
             DfpField.FLAG_INEXACT, "Round #3");

        field.setRoundingMode(DfpField.RoundingMode.ROUND_UP);

        // Round up
        test(field.newDfp("12345678901234567890").add(field.newDfp("0.1")),
             field.newDfp("12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #4");

        test(field.newDfp("12345678901234567890").add(field.newDfp("0.0001")),
             field.newDfp("12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #5");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.1")),
             field.newDfp("-12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #6");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.0001")),
             field.newDfp("-12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #7");

        field.setRoundingMode(DfpField.RoundingMode.ROUND_HALF_UP);

        // Round half up
        test(field.newDfp("12345678901234567890").add(field.newDfp("0.4999")),
             field.newDfp("12345678901234567890"),
             DfpField.FLAG_INEXACT, "Round #8");

        test(field.newDfp("12345678901234567890").add(field.newDfp("0.5000")),
             field.newDfp("12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #9");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.4999")),
             field.newDfp("-12345678901234567890"),
             DfpField.FLAG_INEXACT, "Round #10");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.5000")),
             field.newDfp("-12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #11");

        field.setRoundingMode(DfpField.RoundingMode.ROUND_HALF_DOWN);

        // Round half down
        test(field.newDfp("12345678901234567890").add(field.newDfp("0.5001")),
             field.newDfp("12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #12");

        test(field.newDfp("12345678901234567890").add(field.newDfp("0.5000")),
             field.newDfp("12345678901234567890"),
             DfpField.FLAG_INEXACT, "Round #13");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.5001")),
             field.newDfp("-12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #14");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.5000")),
             field.newDfp("-12345678901234567890"),
             DfpField.FLAG_INEXACT, "Round #15");

        field.setRoundingMode(DfpField.RoundingMode.ROUND_HALF_ODD);

        // Round half odd
        test(field.newDfp("12345678901234567890").add(field.newDfp("0.5000")),
             field.newDfp("12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #16");

        test(field.newDfp("12345678901234567891").add(field.newDfp("0.5000")),
             field.newDfp("12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #17");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.5000")),
             field.newDfp("-12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #18");

        test(field.newDfp("-12345678901234567891").add(field.newDfp("-0.5000")),
             field.newDfp("-12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #19");

        field.setRoundingMode(DfpField.RoundingMode.ROUND_CEIL);

        // Round ceil
        test(field.newDfp("12345678901234567890").add(field.newDfp("0.0001")),
             field.newDfp("12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #20");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.9999")),
             field.newDfp("-12345678901234567890"),
             DfpField.FLAG_INEXACT, "Round #21");

        field.setRoundingMode(DfpField.RoundingMode.ROUND_FLOOR);

        // Round floor
        test(field.newDfp("12345678901234567890").add(field.newDfp("0.9999")),
             field.newDfp("12345678901234567890"),
             DfpField.FLAG_INEXACT, "Round #22");

        test(field.newDfp("-12345678901234567890").add(field.newDfp("-0.0001")),
             field.newDfp("-12345678901234567891"),
             DfpField.FLAG_INEXACT, "Round #23");

        field.setRoundingMode(DfpField.RoundingMode.ROUND_HALF_EVEN);  // reset
    }

    @Override
    @Test
    public void testCeil()
    {
        test(field.newDfp("1234.0000000000000001").ceil(),
             field.newDfp("1235"),
             DfpField.FLAG_INEXACT, "Ceil #1");
    }

    @Test
    void testCeilSmallNegative()
    {
        test(field.newDfp("-0.00009999").ceil(),
             field.newDfp("0"),
             DfpField.FLAG_INEXACT, "Ceil small positive");
    }

    @Test
    void testCeilSmallPositive()
    {
        test(field.newDfp("+0.00009999").ceil(),
             field.newDfp("+1"),
             DfpField.FLAG_INEXACT, "Ceil small positive");
    }

    @Override
    @Test
    public void testFloor()
    {
        test(field.newDfp("1234.9999999999999999").floor(),
             field.newDfp("1234"),
             DfpField.FLAG_INEXACT, "Floor #1");
    }

    @Test
    void testFloorSmallNegative()
    {
        test(field.newDfp("-0.00009999").floor(),
             field.newDfp("-1"),
             DfpField.FLAG_INEXACT, "Floor small negative");
    }

    @Test
    void testFloorSmallPositive()
    {
        test(field.newDfp("+0.00009999").floor(),
             field.newDfp("0"),
             DfpField.FLAG_INEXACT, "Floor small positive");
    }

    @Override
    @Test
    public void testRint()
    {
        test(field.newDfp("1234.50000000001").rint(),
             field.newDfp("1235"),
             DfpField.FLAG_INEXACT, "Rint #1");

        test(field.newDfp("1234.5000").rint(),
             field.newDfp("1234"),
             DfpField.FLAG_INEXACT, "Rint #2");

        test(field.newDfp("1235.5000").rint(),
             field.newDfp("1236"),
             DfpField.FLAG_INEXACT, "Rint #3");
    }

    @Test
    void testCopySign()
    {
        test(Dfp.copysign(field.newDfp("1234."), field.newDfp("-1")),
             field.newDfp("-1234"),
             0, "CopySign #1");

        test(Dfp.copysign(field.newDfp("-1234."), field.newDfp("-1")),
             field.newDfp("-1234"),
             0, "CopySign #2");

        test(Dfp.copysign(field.newDfp("-1234."), field.newDfp("1")),
             field.newDfp("1234"),
             0, "CopySign #3");

        test(Dfp.copysign(field.newDfp("1234."), field.newDfp("1")),
             field.newDfp("1234"),
             0, "CopySign #4");
    }

    @Test
    void testIntValue()
    {
        assertEquals(1234, field.newDfp("1234").intValue(), "intValue #1");
        assertEquals(-1234, field.newDfp("-1234").intValue(), "intValue #2");
        assertEquals(1234, field.newDfp("1234.5").intValue(), "intValue #3");
        assertEquals(1235, field.newDfp("1234.500001").intValue(), "intValue #4");
        assertEquals(2147483647, field.newDfp("1e1000").intValue(), "intValue #5");
        assertEquals(-2147483648, field.newDfp("-1e1000").intValue(), "intValue #6");
    }

    @Test
    void testLog10K()
    {
        assertEquals(1, field.newDfp("123456").log10K(), "log10K #1");
        assertEquals(2, field.newDfp("123456789").log10K(), "log10K #2");
        assertEquals(0, field.newDfp("2").log10K(), "log10K #3");
        assertEquals(0, field.newDfp("1").log10K(), "log10K #3");
        assertEquals(-1, field.newDfp("0.1").log10K(), "log10K #4");
    }

    @Test
    void testPower10K()
    {
        Dfp d = field.newDfp();

        test(d.power10K(0), field.newDfp("1"), 0, "Power10 #1");
        test(d.power10K(1), field.newDfp("10000"), 0, "Power10 #2");
        test(d.power10K(2), field.newDfp("100000000"), 0, "Power10 #3");

        test(d.power10K(-1), field.newDfp("0.0001"), 0, "Power10 #4");
        test(d.power10K(-2), field.newDfp("0.00000001"), 0, "Power10 #5");
        test(d.power10K(-3), field.newDfp("0.000000000001"), 0, "Power10 #6");
    }

    @Test
    public void testLog10()
    {

        assertEquals(1, field.newDfp("12").intLog10(), "log10 #1");
        assertEquals(2, field.newDfp("123").intLog10(), "log10 #2");
        assertEquals(3, field.newDfp("1234").intLog10(), "log10 #3");
        assertEquals(4, field.newDfp("12345").intLog10(), "log10 #4");
        assertEquals(5, field.newDfp("123456").intLog10(), "log10 #5");
        assertEquals(6, field.newDfp("1234567").intLog10(), "log10 #6");
        assertEquals(7, field.newDfp("12345678").intLog10(), "log10 #6");
        assertEquals(8, field.newDfp("123456789").intLog10(), "log10 #7");
        assertEquals(9, field.newDfp("1234567890").intLog10(), "log10 #8");
        assertEquals(10, field.newDfp("12345678901").intLog10(), "log10 #9");
        assertEquals(11, field.newDfp("123456789012").intLog10(), "log10 #10");
        assertEquals(12, field.newDfp("1234567890123").intLog10(), "log10 #11");

        assertEquals(0, field.newDfp("2").intLog10(), "log10 #12");
        assertEquals(0, field.newDfp("1").intLog10(), "log10 #13");
        assertEquals(-1, field.newDfp("0.12").intLog10(), "log10 #14");
        assertEquals(-2, field.newDfp("0.012").intLog10(), "log10 #15");
    }

    @Test
    void testPower10()
    {
        Dfp d = field.newDfp();

        test(d.power10(0), field.newDfp("1"), 0, "Power10 #1");
        test(d.power10(1), field.newDfp("10"), 0, "Power10 #2");
        test(d.power10(2), field.newDfp("100"), 0, "Power10 #3");
        test(d.power10(3), field.newDfp("1000"), 0, "Power10 #4");
        test(d.power10(4), field.newDfp("10000"), 0, "Power10 #5");
        test(d.power10(5), field.newDfp("100000"), 0, "Power10 #6");
        test(d.power10(6), field.newDfp("1000000"), 0, "Power10 #7");
        test(d.power10(7), field.newDfp("10000000"), 0, "Power10 #8");
        test(d.power10(8), field.newDfp("100000000"), 0, "Power10 #9");
        test(d.power10(9), field.newDfp("1000000000"), 0, "Power10 #10");

        test(d.power10(-1), field.newDfp(".1"), 0, "Power10 #11");
        test(d.power10(-2), field.newDfp(".01"), 0, "Power10 #12");
        test(d.power10(-3), field.newDfp(".001"), 0, "Power10 #13");
        test(d.power10(-4), field.newDfp(".0001"), 0, "Power10 #14");
        test(d.power10(-5), field.newDfp(".00001"), 0, "Power10 #15");
        test(d.power10(-6), field.newDfp(".000001"), 0, "Power10 #16");
        test(d.power10(-7), field.newDfp(".0000001"), 0, "Power10 #17");
        test(d.power10(-8), field.newDfp(".00000001"), 0, "Power10 #18");
        test(d.power10(-9), field.newDfp(".000000001"), 0, "Power10 #19");
        test(d.power10(-10), field.newDfp(".0000000001"), 0, "Power10 #20");
    }

    @Test
    void testRemainder()
    {
        test(field.newDfp("10").remainder(field.newDfp("3")),
             field.newDfp("1"),
             DfpField.FLAG_INEXACT, "Remainder #1");

        test(field.newDfp("9").remainder(field.newDfp("3")),
             field.newDfp("0"),
             0, "Remainder #2");

        test(field.newDfp("-9").remainder(field.newDfp("3")),
             field.newDfp("-0"),
             0, "Remainder #3");
    }

    @Override
    @Test
    public void testSqrt()
    {
        test(field.newDfp("0").sqrt(),
             field.newDfp("0"),
             0, "Sqrt #1");

        test(field.newDfp("-0").sqrt(),
             field.newDfp("-0"),
             0, "Sqrt #2");

        test(field.newDfp("1").sqrt(),
             field.newDfp("1"),
             0, "Sqrt #3");

        test(field.newDfp("2").sqrt(),
             field.newDfp("1.4142135623730950"),
             DfpField.FLAG_INEXACT, "Sqrt #4");

        test(field.newDfp("3").sqrt(),
             field.newDfp("1.7320508075688773"),
             DfpField.FLAG_INEXACT, "Sqrt #5");

        test(field.newDfp("5").sqrt(),
             field.newDfp("2.2360679774997897"),
             DfpField.FLAG_INEXACT, "Sqrt #6");

        test(field.newDfp("500").sqrt(),
             field.newDfp("22.3606797749978970"),
             DfpField.FLAG_INEXACT, "Sqrt #6.2");

        test(field.newDfp("50000").sqrt(),
             field.newDfp("223.6067977499789696"),
             DfpField.FLAG_INEXACT, "Sqrt #6.3");

        test(field.newDfp("-1").sqrt(),
             nan,
             DfpField.FLAG_INVALID, "Sqrt #7");

        test(pinf.sqrt(),
             pinf,
             0, "Sqrt #8");

        test(field.newDfp((byte) 1, Dfp.QNAN).sqrt(),
             nan,
             0, "Sqrt #9");

        test(field.newDfp((byte) 1, Dfp.SNAN).sqrt(),
             nan,
             DfpField.FLAG_INVALID, "Sqrt #9");
    }

    @Test
    void testIssue567() {
        DfpField field = new DfpField(100);
        assertEquals(0.0, field.getZero().toDouble(), Precision.SAFE_MIN);
        assertEquals(0.0, field.newDfp(0.0).toDouble(), Precision.SAFE_MIN);
        assertEquals(-1, FastMath.copySign(1, field.newDfp(-0.0).toDouble()), Precision.EPSILON);
        assertEquals(+1, FastMath.copySign(1, field.newDfp(+0.0).toDouble()), Precision.EPSILON);
    }

    @Test
    void testIsZero() {
        assertTrue(field.getZero().isZero());
        assertTrue(field.getZero().negate().isZero());
        assertTrue(field.newDfp(+0.0).isZero());
        assertTrue(field.newDfp(-0.0).isZero());
        assertFalse(field.newDfp(1.0e-90).isZero());
        assertFalse(nan.isZero());
        assertFalse(nan.negate().isZero());
        assertFalse(pinf.isZero());
        assertFalse(pinf.negate().isZero());
        assertFalse(ninf.isZero());
        assertFalse(ninf.negate().isZero());
    }

    @Test
    void testSignPredicates() {

        assertTrue(field.getZero().negativeOrNull());
        assertTrue(field.getZero().positiveOrNull());
        assertFalse(field.getZero().strictlyNegative());
        assertFalse(field.getZero().strictlyPositive());

        assertTrue(field.getZero().negate().negativeOrNull());
        assertTrue(field.getZero().negate().positiveOrNull());
        assertFalse(field.getZero().negate().strictlyNegative());
        assertFalse(field.getZero().negate().strictlyPositive());

        assertFalse(field.getOne().negativeOrNull());
        assertTrue(field.getOne().positiveOrNull());
        assertFalse(field.getOne().strictlyNegative());
        assertTrue(field.getOne().strictlyPositive());

        assertTrue(field.getOne().negate().negativeOrNull());
        assertFalse(field.getOne().negate().positiveOrNull());
        assertTrue(field.getOne().negate().strictlyNegative());
        assertFalse(field.getOne().negate().strictlyPositive());

        assertFalse(nan.negativeOrNull());
        assertFalse(nan.positiveOrNull());
        assertFalse(nan.strictlyNegative());
        assertFalse(nan.strictlyPositive());

        assertFalse(nan.negate().negativeOrNull());
        assertFalse(nan.negate().positiveOrNull());
        assertFalse(nan.negate().strictlyNegative());
        assertFalse(nan.negate().strictlyPositive());

        assertFalse(pinf.negativeOrNull());
        assertTrue(pinf.positiveOrNull());
        assertFalse(pinf.strictlyNegative());
        assertTrue(pinf.strictlyPositive());

        assertTrue(pinf.negate().negativeOrNull());
        assertFalse(pinf.negate().positiveOrNull());
        assertTrue(pinf.negate().strictlyNegative());
        assertFalse(pinf.negate().strictlyPositive());

        assertTrue(ninf.negativeOrNull());
        assertFalse(ninf.positiveOrNull());
        assertTrue(ninf.strictlyNegative());
        assertFalse(ninf.strictlyPositive());

        assertFalse(ninf.negate().negativeOrNull());
        assertTrue(ninf.negate().positiveOrNull());
        assertFalse(ninf.negate().strictlyNegative());
        assertTrue(ninf.negate().strictlyPositive());

    }

    @Test
    void testSpecialConstructors() {
        assertEquals(ninf, field.newDfp(Double.NEGATIVE_INFINITY));
        assertEquals(ninf, field.newDfp("-Infinity"));
        assertEquals(pinf, field.newDfp(Double.POSITIVE_INFINITY));
        assertEquals(pinf, field.newDfp("Infinity"));
        assertTrue(field.newDfp(Double.NaN).isNaN());
        assertTrue(field.newDfp("NaN").isNaN());
    }

    @Test
    void testHypotNoOverflow() {
        Dfp x = field.newDfp(+3);
        Dfp y = field.newDfp(-4);
        Dfp h = field.newDfp(+5);
        for (int i = 0; i < 70000; ++i) {
            x = x.multiply(10);
            y = y.multiply(10);
            h = h.multiply(10);
        }
        assertEquals(h, x.hypot(y));
    }

    @Test
    void testGetExponentVsDouble() {
        for (int i = -1000; i < 1000; ++i) {
            final double x      = FastMath.scalb(1.0, i);
            final double xMinus = 0.99 * x;
            final double xPlus  = 1.01 * x;
            final Dfp dfpMinus  = field.newDfp(xMinus);
            final Dfp dfpPlus   = field.newDfp(xPlus);
            assertEquals(FastMath.getExponent(xMinus), dfpMinus.getExponent());
            assertEquals(FastMath.getExponent(xPlus),  dfpPlus.getExponent());
        }
    }

    @Test
    void testGetExponentSpecialCases() {
        assertEquals(-435412, field.newDfp(0).getExponent());
        assertEquals(0, field.newDfp(1).getExponent());
        assertEquals(1, field.newDfp(2).getExponent());
        assertEquals(435411, field.newDfp(Double.NaN).getExponent());
        assertEquals(435411, field.newDfp(Double.POSITIVE_INFINITY).getExponent());
        assertEquals(435411, field.newDfp(Double.NEGATIVE_INFINITY).getExponent());
    }

    @Test
    void testGetExponentAutonomous() {
        for (int i = -435411; i < 435411; i += 217) {
            final Dfp x = field.newDfp(2).pow(i).multiply(1.1);
            assertEquals(i, x.getExponent());
        }
    }

    @Test
    void testEqualsHashcodeContract() {
        DfpField var1 = new DfpField(1);
        Dfp var6 = var1.newDfp(-0.0d);
        Dfp var5 = var1.newDfp(0L);

        // Checks the contract:  equals-hashcode on var5 and var6
        assertTrue(var5.equals(var6) ? var5.hashCode() == var6.hashCode() : true);
    }

    @Test
    void testZero() {
        Dfp zero = new DfpField(15).getZero();
        assertEquals(0.0, zero.toDouble(), 1.0e-15);
    }

    @Test
    void testOne() {
        Dfp one = new DfpField(15).getOne();
        assertEquals(1.0, one.toDouble(), 1.0e-15);
    }

    @Test
    void testToDegreesDefinition() {
        double epsilon = 1.0e-14;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.0; x += 0.001) {
                Dfp value = new Dfp(field, x);
                assertEquals(FastMath.toDegrees(x), value.toDegrees().getReal(), epsilon);
            }
        }
    }

    @Test
    void testToRadiansDefinition() {
        double epsilon = 1.0e-15;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.0; x += 0.001) {
                Dfp value = new Dfp(field, x);
                assertEquals(FastMath.toRadians(x), value.toRadians().getReal(), epsilon);
            }
        }
    }

    @Test
    void testDegRad() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            Dfp value = field.newDfp("x");
            Dfp rebuilt = value.toDegrees().toRadians();
            Dfp zero = rebuilt.subtract(value);
            assertEquals(0, zero.getReal(), 3.0e-16);
        }
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testMap() {
        int[] decimalDigits = new int[] { 10, 50, 100 }; 
        Map<DfpField, Integer> map = new HashMap<>();
        for (int i = 0; i < 1000; ++i) {
            // create a brand new DfpField for each derivative
            map.put(new DfpField(decimalDigits[i % decimalDigits.length]), 0);
        }

        // despite we have created numerous DfpField instances,
        // there should be only one field for each precision
        assertEquals(decimalDigits.length, map.size());
        DfpField first = map.entrySet().iterator().next().getKey();
        assertEquals(first, first);
        assertNotEquals(first, Binary64Field.getInstance());

    }

    @Test
    void testPi() {
        assertEquals("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117",
                            new DfpField(100).newDfp(1.0).getPi().toString());
    }

    @Test
    void testEquals() {
        DfpField f10A = new DfpField(10);
        DfpField f10B = new DfpField(10);
        DfpField f50  = new DfpField(50);
        assertNotEquals(f10A, f50);
        assertEquals(f10A, f10B);
        f10B.setRoundingMode(DfpField.RoundingMode.ROUND_DOWN);
        assertNotEquals(f10A, f10B);
        f10B.setRoundingMode(DfpField.RoundingMode.ROUND_HALF_EVEN);
        assertEquals(f10A, f10B);
        f10B.setIEEEFlags(DfpField.FLAG_UNDERFLOW | DfpField.FLAG_OVERFLOW);
        assertNotEquals(f10A.getIEEEFlags(), f10B.getIEEEFlags());
        assertEquals(f10A, f10B);
    }

    @Test
    void testRunTimeClass() {
        DfpField field = new DfpField(15);
        assertEquals(Dfp.class, field.getRuntimeClass());
    }

    @Override
    @Test
    public void testLinearCombinationReference() {
        final DfpField field25 = new DfpField(25);
        doTestLinearCombinationReference(field25::newDfp, 4.15e-9, 4.21e-9);
    }

    @Test
    void testConvertToSameAccuracy() {
        DfpField field13 = new DfpField(13);
        DfpField field16 = new DfpField(16); // in fact 13, 14, 15 and 16 decimal digits are all similar to 4 digits in radix 10000
        Dfp dfp = field13.newDfp(1.25);
        assertSame(dfp, dfp.newInstance(field16, DfpField.RoundingMode.ROUND_HALF_EVEN));
    }

    @Test
    void testConvertToHigherAccuracy() {

        DfpField field16 = new DfpField(16);
        DfpField field24 = new DfpField(24);

        checkConvert(field16, "1.25", field24, "1.25", null);

        assertTrue(field16.newDfp(-1).sqrt().newInstance(field24, null).isNaN());
        assertTrue(field16.newDfp().reciprocal().newInstance(field24, null).isInfinite());
    }

    @Test
    void testUlpdDfpA() {
        assertEquals(1.0e-36, new DfpField(40).getOne().ulp().getReal(), 1.0e-51);
        assertEquals(1.0e-40, new DfpField(41).getOne().ulp().getReal(), 1.0e-55);
    }

    @Test
    void testUlpdDfpB() {
        DfpField field = new DfpField(41);
        Dfp one = field.getOne();
        Dfp ulp = one.ulp();
        assertTrue(one.add(ulp).greaterThan(one));
        assertFalse(one.add(ulp.half()).greaterThan(one));
    }

    @Test
    void testConvertToLowerAccuracy() {
        DfpField field16 = new DfpField(16);
        DfpField field24 = new DfpField(24);

        checkConvert(field24, "1234.56789012345678901234", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "1234.56789012345678901234", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "1234.56789012345678901234", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "1234.56789012345678901234", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "1234.56789012345678901234", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "1234.56789012345678901234", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "1234.56789012345678901234", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "1234.56789012345678901234", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "1234.56789012345600000001", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "1234.56789012345600000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "1234.56789012345600000001", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "1234.56789012345600000001", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "1234.56789012345600000001", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "1234.56789012345600000001", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "1234.56789012345600000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "1234.56789012345600000001", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "1234.56789012345650000000", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "1234.56789012345650000000", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "1234.56789012345650000000", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "1234.56789012345650000000", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "1234.56789012345650000000", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "1234.56789012345650000000", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "1234.56789012345650000000", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "1234.56789012345650000000", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "1234.56789012345650000001", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "1234.56789012345650000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "1234.56789012345650000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "1234.56789012345650000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "1234.56789012345650000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "1234.56789012345650000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "1234.56789012345650000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "1234.56789012345650000001", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "1234.56789012345750000000", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "1234.56789012345750000000", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "1234.56789012345750000000", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "1234.56789012345750000000", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "1234.56789012345750000000", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "1234.56789012345750000000", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "1234.56789012345750000000", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "1234.56789012345750000000", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "1234.56789012345750000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "1234.56789012345750000001", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "1234.56789012345750000001", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "1234.56789012345750000001", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "1234.56789012345750000001", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "1234.56789012345750000001", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "1234.56789012345750000001", field16, "1234.567890123458", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "1234.56789012345750000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "-1234.56789012345678901234", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "-1234.56789012345678901234", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "-1234.56789012345678901234", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "-1234.56789012345678901234", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "-1234.56789012345678901234", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "-1234.56789012345678901234", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "-1234.56789012345678901234", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "-1234.56789012345678901234", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "-1234.56789012345600000001", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "-1234.56789012345600000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "-1234.56789012345600000001", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "-1234.56789012345600000001", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "-1234.56789012345600000001", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "-1234.56789012345600000001", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "-1234.56789012345600000001", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "-1234.56789012345600000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "-1234.56789012345650000000", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "-1234.56789012345650000000", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "-1234.56789012345650000000", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "-1234.56789012345650000000", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "-1234.56789012345650000000", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "-1234.56789012345650000000", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "-1234.56789012345650000000", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "-1234.56789012345650000000", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "-1234.56789012345650000001", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "-1234.56789012345650000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "-1234.56789012345650000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "-1234.56789012345650000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "-1234.56789012345650000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "-1234.56789012345650000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "-1234.56789012345650000001", field16, "-1234.567890123456", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "-1234.56789012345650000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "-1234.56789012345750000000", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "-1234.56789012345750000000", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "-1234.56789012345750000000", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "-1234.56789012345750000000", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "-1234.56789012345750000000", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "-1234.56789012345750000000", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "-1234.56789012345750000000", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "-1234.56789012345750000000", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "-1234.56789012345750000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_DOWN);
        checkConvert(field24, "-1234.56789012345750000001", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_UP);
        checkConvert(field24, "-1234.56789012345750000001", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_HALF_UP);
        checkConvert(field24, "-1234.56789012345750000001", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_HALF_DOWN);
        checkConvert(field24, "-1234.56789012345750000001", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "-1234.56789012345750000001", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "-1234.56789012345750000001", field16, "-1234.567890123457", DfpField.RoundingMode.ROUND_CEIL);
        checkConvert(field24, "-1234.56789012345750000001", field16, "-1234.567890123458", DfpField.RoundingMode.ROUND_FLOOR);

        checkConvert(field24, "1234.56789012345650000001", field16, "1234.567890123457", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "1234.56789012345600000000", field16, "1234.567890123456", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "9999.99999999999950000000", field16, "10000.", DfpField.RoundingMode.ROUND_HALF_EVEN);
        checkConvert(field24, "9999.99999999999950000000", field16, "9999.999999999999", DfpField.RoundingMode.ROUND_HALF_ODD);
        checkConvert(field24, "9999.99999999999950000001", field16, "10000.", DfpField.RoundingMode.ROUND_HALF_ODD);

        assertTrue(field24.newDfp(-1).sqrt().newInstance(field16, DfpField.RoundingMode.ROUND_HALF_EVEN).isNaN());
        assertTrue(field24.newDfp().reciprocal().newInstance(field16, DfpField.RoundingMode.ROUND_HALF_EVEN).isInfinite());

    }

    private void checkConvert(DfpField originalField, String originalValue,
                            DfpField targetField, String targetValue,
                            DfpField.RoundingMode rmode) {
        Dfp original  = originalField.newDfp(originalValue);
        Dfp converted = original.newInstance(targetField, rmode);
        assertEquals(targetValue, converted.toString());
    }

}
