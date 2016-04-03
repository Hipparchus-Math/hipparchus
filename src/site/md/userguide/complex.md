# 7 Complex Numbers
## 7.1 Overview
The complex packages provides a complex number type as well as complex
versions of common transcendental functions and complex number
formatting.


## 7.2 Complex Numbers
[          Complex](../apidocs/org.hipparchus/complex/Complex.html)
provides a complex number type that forms the basis for
the complex functionality found in commons-math.

Complex functions and arithmetic operations are implemented in
commons-math by applying standard computational formulas and
following the rules for `java.lang.Double` arithmetic in
handling infinite and `NaN` values.  No attempt is made
to comply with ANSII/IEC C99x Annex G or any other standard for
Complex arithmetic.  See the class and method javadocs for the
[           Complex](../apidocs/org.hipparchus/complex/Complex.html)
and
[           ComplexUtils](../apidocs/org.hipparchus/complex/ComplexUtils.html)
classes for details on computing formulas.

To create a complex number, simply call the constructor passing in two
floating-point arguments, the first being the real part of the
complex number and the second being the imaginary part:
`Complex c = new Complex(1.0, 3.0); // 1 + 3i`

Complex numbers may also be created from polar representations
using the `polar2Complex` method in
`ComplexUtils`.

The `Complex` class provides basic unary and binary
complex number operations.  These operations provide the means to add,
subtract, multiply and divide complex numbers along with other
complex number functions similar to the real number functions found in
`java.math.BigDecimal`:

    Complex lhs = new Complex(1.0, 3.0);
    Complex rhs = new Complex(2.0, 5.0);
    
    Complex answer = lhs.add(rhs);       // add two complex numbers
            answer = lhs.subtract(rhs);  // subtract two complex numbers
            answer = lhs.abs();          // absolute value
            answer = lhs.conjugate(rhs); // complex conjugate


## 7.3 Complex Transcendental Functions
[          Complex](../apidocs/org.hipparchus/complex/Complex.html)
also provides implementations of serveral transcendental
functions involving complex number arguments.
These operations provide the means to compute the log, sine, tangent,
and other complex values :

    Complex first  = new Complex(1.0, 3.0);
    Complex second = new Complex(2.0, 5.0);
    
    Complex answer = first.log();        // natural logarithm.
            answer = first.cos();        // cosine
            answer = first.pow(second);  // first raised to the power of second


## 7.4 Complex Formatting and Parsing
`Complex` instances can be converted to and from strings
using the[          ComplexFormat](../apidocs/org.hipparchus/complex/ComplexFormat.html)
class.
`ComplexFormat` is a `java.text.Format`
extension and, as such, is used like other formatting objects (e.g.
`java.text.SimpleDateFormat`):

    ComplexFormat format = new ComplexFormat(); // default format
    Complex c = new Complex(1.1111, 2.2222);
    String s = format.format(c); // s contains "1.11 + 2.22i"

To customize the formatting output, one or two
`java.text.NumberFormat` instances can be used to construct
a `ComplexFormat`.  These number formats control the
formatting of the real and imaginary values of the complex number:

    NumberFormat nf = NumberFormat.getInstance();
    nf.setMinimumFractionDigits(3);
    nf.setMaximumFractionDigits(3);
    
    // create complex format with custom number format
    // when one number format is used, both real and
    // imaginary parts are formatted the same
    ComplexFormat cf = new ComplexFormat(nf);
    Complex c = new Complex(1.11, 2.2222);
    String s = format.format(c); // s contains "1.110 + 2.222i"
    
    NumberFormat nf2 = NumberFormat.getInstance();
    nf.setMinimumFractionDigits(1);
    nf.setMaximumFractionDigits(1);
    
    // create complex format with custom number formats
    cf = new ComplexFormat(nf, nf2);
    s = format.format(c); // s contains "1.110 + 2.2i"

Another formatting customization provided by
`ComplexFormat` is the text used for the imaginary
designation.  By default, the imaginary notation is "i" but, it can be
manipulated using the `setImaginaryCharacter` method.

Formatting inverse operation, parsing, can also be performed by
`ComplexFormat`.  Parse a complex number from a string,
simply call the `parse` method:

    ComplexFormat cf = new ComplexFormat();
    Complex c = cf.parse("1.110 + 2.222i");


