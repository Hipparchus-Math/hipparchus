<!--
 Licensed to the Hipparchus project under one or more
 contributor license agreements.  See the NOTICE file distributed with
 this work for additional information regarding copyright ownership.
 The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 (the "License"); you may not use this file except in compliance with
 the License.  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->
# Fractions
## Overview
The fraction packages provides a fraction number type as well as
fraction number formatting.

## Fraction Numbers
[Fraction](../apidocs/org/hipparchus/fraction/Fraction.html)
and [BigFraction](../apidocs/org/hipparchus/fraction/BigFraction.html)
provide fraction number type that forms the basis for
the fraction functionality found in Hipparchus. The former one can be
used for fractions whose numerators and denominators are small enough
to fit in an int (taking care of intermediate values) while the second
class should be used when there is a risk the numerator and denominator
grow very large.

A fraction number, can be built from two integer arguments representing numerator
and denominator or from a double which will be approximated:

    Fraction f = new Fraction(1, 3); // 1 / 3
    Fraction g = new Fraction(0.25); // 1 / 4

Of special note with fraction construction, when a fraction is created it is always reduced to lowest terms.

The `Fraction` class provides many unary and binary
fraction operations.  These operations provide the means to add,
subtract, multiple and, divide fractions along with other functions similar to the real number functions found in
`java.math.BigDecimal`:

    Fraction lhs = new Fraction(1, 3);
    Fraction rhs = new Fraction(2, 5);
    
    Fraction answer = lhs.add(rhs);     // add two fractions
            answer = lhs.subtract(rhs); // subtract two fractions
            answer = lhs.abs();         // absolute value
            answer = lhs.reciprocal();  // reciprocal of lhs

Like fraction construction, for each of the fraction functions, the resulting fraction is reduced to lowest terms.

## Fraction Formatting and Parsing
`Fraction` instances can be converted to and from strings using the
[FractionFormat](../apidocs/org/hipparchus/fraction/FractionFormat.html) class.
`FractionFormat` is a `java.text.Format` extension and, as such,
is used like other formatting objects (e.g. `java.text.SimpleDateFormat`):

    FractionFormat format = new FractionFormat(); // default format
    Fraction f = new Fraction(2, 4);
    String s = format.format(f); // s contains "1 / 2", note the reduced fraction

To customize the formatting output, one or two
`java.text.NumberFormat` instances can be used to construct
a `FractionFormat`.  These number formats control the
formatting of the numerator and denominator of the fraction:

    NumberFormat nf = NumberFormat.getInstance(Locale.FRANCE);
    // create fraction format with custom number format
    // when one number format is used, both numerator and
    // denominator are formatted the same
    FractionFormat format = new FractionFormat(nf);
    Fraction f = new Fraction(2000, 3333);
    String s = format.format(c); // s contains "2.000 / 3.333"
    
    NumberFormat nf2 = NumberFormat.getInstance(Locale.US);
    // create fraction format with custom number formats
    format = new FractionFormat(nf, nf2);
    s = format.format(f); // s contains "2.000 / 3,333"

Formatting's inverse operation, parsing, can also be performed by
`FractionFormat`.  To parse a fraction from a string,
simply call the `parse` method:

    FractionFormat ff = new FractionFormat();
    Fraction f = ff.parse("-10 / 21");
