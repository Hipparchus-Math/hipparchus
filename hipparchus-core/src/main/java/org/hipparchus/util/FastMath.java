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
package org.hipparchus.util;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathRuntimeException;

/**
 * Faster, more accurate, portable alternative to {@link Math} and
 * {@link StrictMath} for large scale computation.
 * <p>
 * FastMath is a drop-in replacement for both Math and StrictMath. This
 * means that for any method in Math (say {@code Math.sin(x)} or
 * {@code Math.cbrt(y)}), user can directly change the class and use the
 * methods as is (using {@code FastMath.sin(x)} or {@code FastMath.cbrt(y)}
 * in the previous example).
 * <p>
 * FastMath's speed is achieved by relying heavily on JIT compiler optimization
 * to native code present in many JVMs today, and use of large tables.
 * The larger tables are lazily initialized on first use, so that the setup
 * time does not penalize methods that don't need them.
 *
 * TODO: where are the tables being lazily initialized?
 *
 * <p>
 * Note that FastMath is
 * extensively used inside Hipparchus, so by calling some algorithms,
 * the overhead when the the tables need to be initialized will occur
 * regardless of the end-user calling FastMath methods directly or not.
 * Performance figures for a specific JVM and hardware can be evaluated by
 * running the FastMathTestPerformance tests in the test directory of the source
 * distribution.
 * <p>
 * FastMath accuracy should be mostly independent of the JVM as it relies only
 * on IEEE-754 basic operations and on embedded tables. Almost all operations
 * are accurate to about 0.5 ulp throughout the domain range. This statement,
 * of course is only a rough global observed behavior, it is <em>not</em> a
 * guarantee for <em>every</em> double numbers input (see William Kahan's <a
 * href="http://en.wikipedia.org/wiki/Rounding#The_table-maker.27s_dilemma">Table
 * Maker's Dilemma</a>).
 * <p>
 * FastMath additionally implements the following methods not found in Math/StrictMath:
 * <ul>
 * <li>{@link #asinh(double)}</li>
 * <li>{@link #acosh(double)}</li>
 * <li>{@link #atanh(double)}</li>
 * </ul>
 * The following methods are found in Math/StrictMath since 1.6 only, they are provided
 * by FastMath even in 1.5 Java virtual machines
 * <ul>
 * <li>{@link #copySign(double, double)}</li>
 * <li>{@link #getExponent(double)}</li>
 * <li>{@link #nextAfter(double,double)}</li>
 * <li>{@link #nextUp(double)}</li>
 * <li>{@link #scalb(double, int)}</li>
 * <li>{@link #copySign(float, float)}</li>
 * <li>{@link #getExponent(float)}</li>
 * <li>{@link #nextAfter(float,double)}</li>
 * <li>{@link #nextUp(float)}</li>
 * <li>{@link #scalb(float, int)}</li>
 * </ul>
 */
public class FastMath {
    /** Archimede's constant PI, ratio of circle circumference to diameter. */
    public static final double PI = 105414357.0 / 33554432.0 + 1.984187159361080883e-9;

    /** Napier's constant e, base of the natural logarithm. */
    public static final double E = 2850325.0 / 1048576.0 + 8.254840070411028747e-8;

    /** Index of exp(0) in the array of integer exponentials. */
    static final int EXP_INT_TABLE_MAX_INDEX = 750;
    /** Length of the array of integer exponentials. */
    static final int EXP_INT_TABLE_LEN = EXP_INT_TABLE_MAX_INDEX * 2;
    /** Logarithm table length. */
    static final int LN_MANT_LEN = 1024;
    /** Exponential fractions table length. */
    static final int EXP_FRAC_TABLE_LEN = 1025; // 0, 1/1024, ... 1024/1024

    /** StrictMath.log(Double.MAX_VALUE): {@value} */
    private static final double LOG_MAX_VALUE = StrictMath.log(Double.MAX_VALUE);

    /** Indicator for tables initialization.
     * <p>
     * This compile-time constant should be set to true only if one explicitly
     * wants to compute the tables at class loading time instead of using the
     * already computed ones provided as literal arrays below.
     * </p>
     */
    private static final boolean RECOMPUTE_TABLES_AT_RUNTIME = false;

    /** Sine, Cosine, Tangent tables are for 0, 1/8, 2/8, ... 13/8 = PI/2 approx. */

    /** Tangent table, used by atan() (high bits). */
    private static final double[] TANGENT_TABLE_A =
        {
        +0.0d,
        +0.1256551444530487d,
        +0.25534194707870483d,
        +0.3936265707015991d,
        +0.5463024377822876d,
        +0.7214844226837158d,
        +0.9315965175628662d,
        +1.1974215507507324d,
        +1.5574076175689697d,
        +2.092571258544922d,
        +3.0095696449279785d,
        +5.041914939880371d,
        +14.101419448852539d,
        -18.430862426757812d,
    };

    /** Tangent table, used by atan() (low bits). */
    private static final double[] TANGENT_TABLE_B =
        {
        +0.0d,
        -7.877917738262007E-9d,
        -2.5857668567479893E-8d,
        +5.2240336371356666E-9d,
        +5.206150291559893E-8d,
        +1.8307188599677033E-8d,
        -5.7618793749770706E-8d,
        +7.848361555046424E-8d,
        +1.0708593250394448E-7d,
        +1.7827257129423813E-8d,
        +2.893485277253286E-8d,
        +3.1660099222737955E-7d,
        +4.983191803254889E-7d,
        -3.356118100840571E-7d,
    };

    /** Eighths.
     * This is used by sinQ, because its faster to do a table lookup than
     * a multiply in this time-critical routine
     */
    private static final double[]
                    EIGHTHS = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625};

    /**
     * 0x40000000 - used to split a double into two parts, both with the low order bits cleared.
     * Equivalent to 2^30.
     */
    private static final long HEX_40000000 = 0x40000000L; // 1073741824L

    /** Mask used to clear low order 30 bits */
    private static final long MASK_30BITS = -1L - (HEX_40000000 -1); // 0xFFFFFFFFC0000000L;

    /** Mask used to clear the non-sign part of an int. */
    private static final int MASK_NON_SIGN_INT = 0x7fffffff;

    /** Constant: {@value}. */
    private static final double F_1_3 = 1d / 3d;
    /** Constant: {@value}. */
    private static final double F_1_5 = 1d / 5d;
    /** Constant: {@value}. */
    private static final double F_1_7 = 1d / 7d;
    /** Constant: {@value}. */
    private static final double F_1_9 = 1d / 9d;
    /** Constant: {@value}. */
    private static final double F_1_11 = 1d / 11d;
    /** Constant: {@value}. */
    private static final double F_1_13 = 1d / 13d;
    /** Constant: {@value}. */
    private static final double F_1_15 = 1d / 15d;
    /** Constant: {@value}. */
    private static final double F_1_17 = 1d / 17d;
    /** Constant: {@value}. */
    private static final double F_3_4 = 3d / 4d;
    /** Constant: {@value}. */
    private static final double F_15_16 = 15d / 16d;
    /** Constant: {@value}. */
    private static final double F_13_14 = 13d / 14d;
    /** Constant: {@value}. */
    private static final double F_11_12 = 11d / 12d;
    /** Constant: {@value}. */
    private static final double F_9_10 = 9d / 10d;
    /** Constant: {@value}. */
    private static final double F_7_8 = 7d / 8d;
    /** Constant: {@value}. */
    private static final double F_5_6 = 5d / 6d;
    /** Constant: {@value}. */
    private static final double F_1_2 = 1d / 2d;
    /** Constant: {@value}. */
    private static final double F_1_4 = 1d / 4d;

    /**
     * Private Constructor
     */
    private FastMath() {}

    // Generic helper methods

    /**
     * Get the high order bits from the mantissa.
     * Equivalent to adding and subtracting HEX_40000 but also works for very large numbers
     *
     * @param d the value to split
     * @return the high order part of the mantissa
     */
    private static double doubleHighPart(double d) {
        if (d > -Precision.SAFE_MIN && d < Precision.SAFE_MIN){
            return d; // These are un-normalised - don't try to convert
        }
        long xl = Double.doubleToRawLongBits(d); // can take raw bits because just gonna convert it back
        xl &= MASK_30BITS; // Drop low order bits
        return Double.longBitsToDouble(xl);
    }

    /** Compute the square root of a number.
     * <p><b>Note:</b> this implementation currently delegates to {@link Math#sqrt}
     * @param a number on which evaluation is done
     * @return square root of a
     */
    public static double sqrt(final double a) {
        return Math.sqrt(a);
    }

    /** Compute the hyperbolic cosine of a number.
     * @param x number on which evaluation is done
     * @return hyperbolic cosine of x
     */
    public static double cosh(double x) {
      if (Double.isNaN(x)) {
          return x;
      }

      // cosh[z] = (exp(z) + exp(-z))/2

      // for numbers with magnitude 20 or so,
      // exp(-z) can be ignored in comparison with exp(z)

      if (x > 20) {
          if (x >= LOG_MAX_VALUE) {
              // Avoid overflow (MATH-905).
              final double t = exp(0.5 * x);
              return (0.5 * t) * t;
          } else {
              return 0.5 * exp(x);
          }
      } else if (x < -20) {
          if (x <= -LOG_MAX_VALUE) {
              // Avoid overflow (MATH-905).
              final double t = exp(-0.5 * x);
              return (0.5 * t) * t;
          } else {
              return 0.5 * exp(-x);
          }
      }

      final double[] hiPrec = new double[2];
      if (x < 0.0) {
          x = -x;
      }
      exp(x, 0.0, hiPrec);

      double ya = hiPrec[0] + hiPrec[1];
      double yb = -(ya - hiPrec[0] - hiPrec[1]);

      double temp = ya * HEX_40000000;
      double yaa = ya + temp - temp;
      double yab = ya - yaa;

      // recip = 1/y
      double recip = 1.0/ya;
      temp = recip * HEX_40000000;
      double recipa = recip + temp - temp;
      double recipb = recip - recipa;

      // Correct for rounding in division
      recipb += (1.0 - yaa*recipa - yaa*recipb - yab*recipa - yab*recipb) * recip;
      // Account for yb
      recipb += -yb * recip * recip;

      // y = y + 1/y
      temp = ya + recipa;
      yb += -(temp - ya - recipa);
      ya = temp;
      temp = ya + recipb;
      yb += -(temp - ya - recipb);
      ya = temp;

      double result = ya + yb;
      result *= 0.5;
      return result;
    }

    /** Compute the hyperbolic sine of a number.
     * @param x number on which evaluation is done
     * @return hyperbolic sine of x
     */
    public static double sinh(double x) {
      boolean negate = false;
      if (Double.isNaN(x)) {
          return x;
      }

      // sinh[z] = (exp(z) - exp(-z) / 2

      // for values of z larger than about 20,
      // exp(-z) can be ignored in comparison with exp(z)

      if (x > 20) {
          if (x >= LOG_MAX_VALUE) {
              // Avoid overflow (MATH-905).
              final double t = exp(0.5 * x);
              return (0.5 * t) * t;
          } else {
              return 0.5 * exp(x);
          }
      } else if (x < -20) {
          if (x <= -LOG_MAX_VALUE) {
              // Avoid overflow (MATH-905).
              final double t = exp(-0.5 * x);
              return (-0.5 * t) * t;
          } else {
              return -0.5 * exp(-x);
          }
      }

      if (x == 0) {
          return x;
      }

      if (x < 0.0) {
          x = -x;
          negate = true;
      }

      double[] hiPrec = new double[2];
      double result;

      if (x > 0.25) {
          exp(x, 0.0, hiPrec);

          double ya = hiPrec[0] + hiPrec[1];
          double yb = -(ya - hiPrec[0] - hiPrec[1]);

          double temp = ya * HEX_40000000;
          double yaa = ya + temp - temp;
          double yab = ya - yaa;

          // recip = 1/y
          double recip = 1.0/ya;
          temp = recip * HEX_40000000;
          double recipa = recip + temp - temp;
          double recipb = recip - recipa;

          // Correct for rounding in division
          recipb += (1.0 - yaa*recipa - yaa*recipb - yab*recipa - yab*recipb) * recip;
          // Account for yb
          recipb += -yb * recip * recip;

          recipa = -recipa;
          recipb = -recipb;

          // y = y - 1/y
          temp = ya + recipa;
          yb += -(temp - ya - recipa);
          ya = temp;
          temp = ya + recipb;
          yb += -(temp - ya - recipb);
          ya = temp;

          result = ya + yb;
          result *= 0.5;
      } else {
          expm1(x, hiPrec);

          double ya = hiPrec[0] + hiPrec[1];
          double yb = -(ya - hiPrec[0] - hiPrec[1]);

          /* Compute expm1(-x) = -expm1(x) / (expm1(x) + 1) */
          double denom = 1.0 + ya;
          double denomr = 1.0 / denom;
          double denomb = -(denom - 1.0 - ya) + yb;
          double ratio = ya * denomr;
          double temp = ratio * HEX_40000000;
          double ra = ratio + temp - temp;
          double rb = ratio - ra;

          temp = denom * HEX_40000000;
          double za = denom + temp - temp;
          double zb = denom - za;

          rb += (ya - za*ra - za*rb - zb*ra - zb*rb) * denomr;

          // Adjust for yb
          rb += yb*denomr;                        // numerator
          rb += -ya * denomb * denomr * denomr;   // denominator

          // y = y - 1/y
          temp = ya + ra;
          yb += -(temp - ya - ra);
          ya = temp;
          temp = ya + rb;
          yb += -(temp - ya - rb);
          ya = temp;

          result = ya + yb;
          result *= 0.5;
      }

      if (negate) {
          result = -result;
      }

      return result;
    }

    /**
     * Combined hyperbolic sine and hyperbolic cosine function.
     *
     * @param x Argument.
     * @return [sinh(x), cosh(x)]
     */
    public static SinhCosh sinhCosh(double x) {
      boolean negate = false;
      if (Double.isNaN(x)) {
          return new SinhCosh(x, x);
      }

      // sinh[z] = (exp(z) - exp(-z) / 2
      // cosh[z] = (exp(z) + exp(-z))/2

      // for values of z larger than about 20,
      // exp(-z) can be ignored in comparison with exp(z)

      if (x > 20) {
          final double e;
          if (x >= LOG_MAX_VALUE) {
              // Avoid overflow (MATH-905).
              final double t = exp(0.5 * x);
              e = (0.5 * t) * t;
          } else {
              e = 0.5 * exp(x);
          }
          return new SinhCosh(e, e);
      } else if (x < -20) {
          final double e;
          if (x <= -LOG_MAX_VALUE) {
              // Avoid overflow (MATH-905).
              final double t = exp(-0.5 * x);
              e = (-0.5 * t) * t;
          } else {
              e = -0.5 * exp(-x);
          }
          return new SinhCosh(e, -e);
      }

      if (x == 0) {
          return new SinhCosh(x, 1.0);
      }

      if (x < 0.0) {
          x = -x;
          negate = true;
      }

      double[] hiPrec = new double[2];
      double resultM;
      double resultP;

      if (x > 0.25) {
          exp(x, 0.0, hiPrec);

          final double ya = hiPrec[0] + hiPrec[1];
          final double yb = -(ya - hiPrec[0] - hiPrec[1]);

          double temp = ya * HEX_40000000;
          double yaa = ya + temp - temp;
          double yab = ya - yaa;

          // recip = 1/y
          double recip = 1.0/ya;
          temp = recip * HEX_40000000;
          double recipa = recip + temp - temp;
          double recipb = recip - recipa;

          // Correct for rounding in division
          recipb += (1.0 - yaa*recipa - yaa*recipb - yab*recipa - yab*recipb) * recip;
          // Account for yb
          recipb += -yb * recip * recip;

          // y = y - 1/y
          temp = ya - recipa;
          double ybM = yb - (temp - ya + recipa);
          double yaM = temp;
          temp = yaM - recipb;
          ybM += -(temp - yaM + recipb);
          yaM = temp;
          resultM = yaM + ybM;
          resultM *= 0.5;

          // y = y + 1/y
          temp = ya + recipa;
          double ybP = yb - (temp - ya - recipa);
          double yaP = temp;
          temp = yaP + recipb;
          ybP += -(temp - yaP - recipb);
          yaP = temp;
          resultP = yaP + ybP;
          resultP *= 0.5;

      } else {
          expm1(x, hiPrec);

          final double ya = hiPrec[0] + hiPrec[1];
          final double yb = -(ya - hiPrec[0] - hiPrec[1]);

          /* Compute expm1(-x) = -expm1(x) / (expm1(x) + 1) */
          double denom = 1.0 + ya;
          double denomr = 1.0 / denom;
          double denomb = -(denom - 1.0 - ya) + yb;
          double ratio = ya * denomr;
          double temp = ratio * HEX_40000000;
          double ra = ratio + temp - temp;
          double rb = ratio - ra;

          temp = denom * HEX_40000000;
          double za = denom + temp - temp;
          double zb = denom - za;

          rb += (ya - za*ra - za*rb - zb*ra - zb*rb) * denomr;

          // Adjust for yb
          rb += yb*denomr;                        // numerator
          rb += -ya * denomb * denomr * denomr;   // denominator

          // y = y - 1/y
          temp = ya + ra;
          double ybM = yb - (temp - ya - ra);
          double yaM = temp;
          temp = yaM + rb;
          ybM += -(temp - yaM - rb);
          yaM = temp;
          resultM = yaM + ybM;
          resultM *= 0.5;

          // y = y + 1/y + 2
          temp = ya - ra;
          double ybP = yb - (temp - ya + ra);
          double yaP = temp;
          temp = yaP - rb;
          ybP += -(temp - yaP + rb);
          yaP = temp;
          resultP = yaP + ybP + 2;
          resultP *= 0.5;
      }

      if (negate) {
          resultM = -resultM;
      }

      return new SinhCosh(resultM, resultP);

    }

    /**
     * Combined hyperbolic sine and hyperbolic cosine function.
     *
     * @param x Argument.
     * @param <T> the type of the field element
     * @return [sinh(x), cosh(x)]
     */
    public static <T extends CalculusFieldElement<T>> FieldSinhCosh<T> sinhCosh(T x) {
        return x.sinhCosh();
    }

    /** Compute the hyperbolic tangent of a number.
     * @param x number on which evaluation is done
     * @return hyperbolic tangent of x
     */
    public static double tanh(double x) {
      boolean negate = false;

      if (Double.isNaN(x)) {
          return x;
      }

      // tanh[z] = sinh[z] / cosh[z]
      // = (exp(z) - exp(-z)) / (exp(z) + exp(-z))
      // = (exp(2x) - 1) / (exp(2x) + 1)

      // for magnitude > 20, sinh[z] == cosh[z] in double precision

      if (x > 20.0) {
          return 1.0;
      }

      if (x < -20) {
          return -1.0;
      }

      if (x == 0) {
          return x;
      }

      if (x < 0.0) {
          x = -x;
          negate = true;
      }

      double result;
      if (x >= 0.5) {
          double[] hiPrec = new double[2];
          // tanh(x) = (exp(2x) - 1) / (exp(2x) + 1)
          exp(x*2.0, 0.0, hiPrec);

          double ya = hiPrec[0] + hiPrec[1];
          double yb = -(ya - hiPrec[0] - hiPrec[1]);

          /* Numerator */
          double na = -1.0 + ya;
          double nb = -(na + 1.0 - ya);
          double temp = na + yb;
          nb += -(temp - na - yb);
          na = temp;

          /* Denominator */
          double da = 1.0 + ya;
          double db = -(da - 1.0 - ya);
          temp = da + yb;
          db += -(temp - da - yb);
          da = temp;

          temp = da * HEX_40000000;
          double daa = da + temp - temp;
          double dab = da - daa;

          // ratio = na/da
          double ratio = na/da;
          temp = ratio * HEX_40000000;
          double ratioa = ratio + temp - temp;
          double ratiob = ratio - ratioa;

          // Correct for rounding in division
          ratiob += (na - daa*ratioa - daa*ratiob - dab*ratioa - dab*ratiob) / da;

          // Account for nb
          ratiob += nb / da;
          // Account for db
          ratiob += -db * na / da / da;

          result = ratioa + ratiob;
      }
      else {
          double[] hiPrec = new double[2];
          // tanh(x) = expm1(2x) / (expm1(2x) + 2)
          expm1(x*2.0, hiPrec);

          double ya = hiPrec[0] + hiPrec[1];
          double yb = -(ya - hiPrec[0] - hiPrec[1]);

          /* Numerator */
          double na = ya;
          double nb = yb;

          /* Denominator */
          double da = 2.0 + ya;
          double db = -(da - 2.0 - ya);
          double temp = da + yb;
          db += -(temp - da - yb);
          da = temp;

          temp = da * HEX_40000000;
          double daa = da + temp - temp;
          double dab = da - daa;

          // ratio = na/da
          double ratio = na/da;
          temp = ratio * HEX_40000000;
          double ratioa = ratio + temp - temp;
          double ratiob = ratio - ratioa;

          // Correct for rounding in division
          ratiob += (na - daa*ratioa - daa*ratiob - dab*ratioa - dab*ratiob) / da;

          // Account for nb
          ratiob += nb / da;
          // Account for db
          ratiob += -db * na / da / da;

          result = ratioa + ratiob;
      }

      if (negate) {
          result = -result;
      }

      return result;
    }

    /** Compute the inverse hyperbolic cosine of a number.
     * @param a number on which evaluation is done
     * @return inverse hyperbolic cosine of a
     */
    public static double acosh(final double a) {
        return log(a + sqrt(a * a - 1));
    }

    /** Compute the inverse hyperbolic sine of a number.
     * @param a number on which evaluation is done
     * @return inverse hyperbolic sine of a
     */
    public static double asinh(double a) {
        boolean negative = false;
        if (a < 0) {
            negative = true;
            a = -a;
        }

        double absAsinh;
        if (a > 0.167) {
            absAsinh = log(sqrt(a * a + 1) + a);
        } else {
            final double a2 = a * a;
            if (a > 0.097) {
                absAsinh = a * (1 - a2 * (F_1_3 - a2 * (F_1_5 - a2 * (F_1_7 - a2 * (F_1_9 - a2 * (F_1_11 - a2 * (F_1_13 - a2 * (F_1_15 - a2 * F_1_17 * F_15_16) * F_13_14) * F_11_12) * F_9_10) * F_7_8) * F_5_6) * F_3_4) * F_1_2);
            } else if (a > 0.036) {
                absAsinh = a * (1 - a2 * (F_1_3 - a2 * (F_1_5 - a2 * (F_1_7 - a2 * (F_1_9 - a2 * (F_1_11 - a2 * F_1_13 * F_11_12) * F_9_10) * F_7_8) * F_5_6) * F_3_4) * F_1_2);
            } else if (a > 0.0036) {
                absAsinh = a * (1 - a2 * (F_1_3 - a2 * (F_1_5 - a2 * (F_1_7 - a2 * F_1_9 * F_7_8) * F_5_6) * F_3_4) * F_1_2);
            } else {
                absAsinh = a * (1 - a2 * (F_1_3 - a2 * F_1_5 * F_3_4) * F_1_2);
            }
        }

        return negative ? -absAsinh : absAsinh;
    }

    /** Compute the inverse hyperbolic tangent of a number.
     * @param a number on which evaluation is done
     * @return inverse hyperbolic tangent of a
     */
    public static double atanh(double a) {
        boolean negative = false;
        if (a < 0) {
            negative = true;
            a = -a;
        }

        double absAtanh;
        if (a > 0.15) {
            absAtanh = 0.5 * log((1 + a) / (1 - a));
        } else {
            final double a2 = a * a;
            if (a > 0.087) {
                absAtanh = a * (1 + a2 * (F_1_3 + a2 * (F_1_5 + a2 * (F_1_7 + a2 * (F_1_9 + a2 * (F_1_11 + a2 * (F_1_13 + a2 * (F_1_15 + a2 * F_1_17))))))));
            } else if (a > 0.031) {
                absAtanh = a * (1 + a2 * (F_1_3 + a2 * (F_1_5 + a2 * (F_1_7 + a2 * (F_1_9 + a2 * (F_1_11 + a2 * F_1_13))))));
            } else if (a > 0.003) {
                absAtanh = a * (1 + a2 * (F_1_3 + a2 * (F_1_5 + a2 * (F_1_7 + a2 * F_1_9))));
            } else {
                absAtanh = a * (1 + a2 * (F_1_3 + a2 * F_1_5));
            }
        }

        return negative ? -absAtanh : absAtanh;
    }

    /** Compute the signum of a number.
     * The signum is -1 for negative numbers, +1 for positive numbers and 0 otherwise
     * @param a number on which evaluation is done
     * @return -1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a
     */
    public static double signum(final double a) {
        return (a < 0.0) ? -1.0 : ((a > 0.0) ? 1.0 : a); // return +0.0/-0.0/NaN depending on a
    }

    /** Compute the signum of a number.
     * The signum is -1 for negative numbers, +1 for positive numbers and 0 otherwise
     * @param a number on which evaluation is done
     * @return -1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a
     */
    public static float signum(final float a) {
        return (a < 0.0f) ? -1.0f : ((a > 0.0f) ? 1.0f : a); // return +0.0/-0.0/NaN depending on a
    }

    /** Compute next number towards positive infinity.
     * @param a number to which neighbor should be computed
     * @return neighbor of a towards positive infinity
     */
    public static double nextUp(final double a) {
        return nextAfter(a, Double.POSITIVE_INFINITY);
    }

    /** Compute next number towards positive infinity.
     * @param a number to which neighbor should be computed
     * @return neighbor of a towards positive infinity
     */
    public static float nextUp(final float a) {
        return nextAfter(a, Float.POSITIVE_INFINITY);
    }

    /** Compute next number towards negative infinity.
     * @param a number to which neighbor should be computed
     * @return neighbor of a towards negative infinity
     */
    public static double nextDown(final double a) {
        return nextAfter(a, Double.NEGATIVE_INFINITY);
    }

    /** Compute next number towards negative infinity.
     * @param a number to which neighbor should be computed
     * @return neighbor of a towards negative infinity
     */
    public static float nextDown(final float a) {
        return nextAfter(a, Float.NEGATIVE_INFINITY);
    }

    /** Clamp a value within an interval.
     * @param value value to clamp
     * @param inf lower bound of the clamping interval
     * @param sup upper bound of the clamping interval
     * @return value clamped within [inf; sup], or value if already within bounds.
     * @since 3.0
     */
    public static int clamp(final int value, final int inf, final int sup) {
        return max(inf, min(value, sup));
    }

    /** Clamp a value within an interval.
     * @param value value to clamp
     * @param inf lower bound of the clamping interval
     * @param sup upper bound of the clamping interval
     * @return value clamped within [inf; sup], or value if already within bounds.
     * @since 3.0
     */
    public static long clamp(final long value, final long inf, final long sup) {
        return max(inf, min(value, sup));
    }

    /** Clamp a value within an interval.
     * @param value value to clamp
     * @param inf lower bound of the clamping interval
     * @param sup upper bound of the clamping interval
     * @return value clamped within [inf; sup], or value if already within bounds.
     * @since 3.0
     */
    public static int clamp(final long value, final int inf, final int sup) {
        return (int) max(inf, min(value, sup));
    }

    /** Clamp a value within an interval.
     * <p>
     * This method assumes -0.0 is below +0.0
     * </p>
     * @param value value to clamp
     * @param inf lower bound of the clamping interval
     * @param sup upper bound of the clamping interval
     * @return value clamped within [inf; sup], or value if already within bounds.
     * @since 3.0
     */
    public static float clamp(final float value, final float inf, final float sup) {
        return max(inf, min(value, sup));
    }

    /** Clamp a value within an interval.
     * <p>
     * This method assumes -0.0 is below +0.0
     * </p>
     * @param value value to clamp
     * @param inf lower bound of the clamping interval
     * @param sup upper bound of the clamping interval
     * @return value clamped within [inf; sup], or value if already within bounds.
     * @since 3.0
     */
    public static double clamp(final double value, final double inf, final double sup) {
        return max(inf, min(value, sup));
    }

    /** Returns a pseudo-random number between 0.0 and 1.0.
     * <p><b>Note:</b> this implementation currently delegates to {@link Math#random}
     * @return a random number between 0.0 and 1.0
     */
    public static double random() {
        return Math.random();
    }

    /**
     * Exponential function.<br>
     * <br>
     *
     * Delegates to {@link Math#exp(double)}
     *
     * @param x   a double
     * @return double e<sup>x</sup>
     */
    public static double exp(double x) {
        return StrictMath.exp(x);
    }

    /**
     * Internal helper method for computing a high-precision exponential in
     *
     * <ul>
     * <li>{@link #sinh(double)},</li>
     * <li>{@link #cosh(double)},</li>
     * <li>{@link #tanh(double)}, and</li>
     * <li>{@link #sinhCosh(double)}</li>
     * </ul>
     *
     * @param x original argument of the exponential function
     * @param extra extra bits of precision on input (To Be Confirmed)
     * @param hiPrec extra bits of precision on output (To Be Confirmed)
     * @return exp(x)
     */
    private static double exp(double x, double extra, double[] hiPrec) {
        double intPartA;
        double intPartB;
        int intVal = (int) x;

        /* Lookup exp(floor(x)).
         * intPartA will have the upper 22 bits, intPartB will have the lower
         * 52 bits.
         */
        if (x < 0.0) {

            // We don't check against intVal here as conversion of large negative double values
            // may be affected by a JIT bug. Subsequent comparisons can safely use intVal
            if (x < -746d) {
                if (hiPrec != null) {
                    hiPrec[0] = 0.0;
                    hiPrec[1] = 0.0;
                }
                return 0.0;
            }

            if (intVal < -709) {
                /* This will produce a subnormal output */
                final double result = exp(x+40.19140625, extra, hiPrec) / 285040095144011776.0;
                if (hiPrec != null) {
                    hiPrec[0] /= 285040095144011776.0;
                    hiPrec[1] /= 285040095144011776.0;
                }
                return result;
            }

            if (intVal == -709) {
                /* exp(1.494140625) is nearly a machine number... */
                final double result = exp(x+1.494140625, extra, hiPrec) / 4.455505956692756620;
                if (hiPrec != null) {
                    hiPrec[0] /= 4.455505956692756620;
                    hiPrec[1] /= 4.455505956692756620;
                }
                return result;
            }

            intVal--;

        } else {
            if (intVal > 709) {
                if (hiPrec != null) {
                    hiPrec[0] = Double.POSITIVE_INFINITY;
                    hiPrec[1] = 0.0;
                }
                return Double.POSITIVE_INFINITY;
            }

        }

        intPartA = ExpIntTable.EXP_INT_TABLE_A[EXP_INT_TABLE_MAX_INDEX+intVal];
        intPartB = ExpIntTable.EXP_INT_TABLE_B[EXP_INT_TABLE_MAX_INDEX+intVal];

        /* Get the fractional part of x, find the greatest multiple of 2^-10 less than
         * x and look up the exp function of it.
         * fracPartA will have the upper 22 bits, fracPartB the lower 52 bits.
         */
        final int intFrac = (int) ((x - intVal) * 1024.0);
        final double fracPartA = ExpFracTable.EXP_FRAC_TABLE_A[intFrac];
        final double fracPartB = ExpFracTable.EXP_FRAC_TABLE_B[intFrac];

        /* epsilon is the difference in x from the nearest multiple of 2^-10.  It
         * has a value in the range 0 <= epsilon < 2^-10.
         * Do the subtraction from x as the last step to avoid possible loss of precision.
         */
        final double epsilon = x - (intVal + intFrac / 1024.0);

        /* Compute z = exp(epsilon) - 1.0 via a minimax polynomial.  z has
       full double precision (52 bits).  Since z < 2^-10, we will have
       62 bits of precision when combined with the constant 1.  This will be
       used in the last addition below to get proper rounding. */

        /* Remez generated polynomial.  Converges on the interval [0, 2^-10], error
       is less than 0.5 ULP */
        double z = 0.04168701738764507;
        z = z * epsilon + 0.1666666505023083;
        z = z * epsilon + 0.5000000000042687;
        z = z * epsilon + 1.0;
        z = z * epsilon + -3.940510424527919E-20;

        /* Compute (intPartA+intPartB) * (fracPartA+fracPartB) by binomial
       expansion.
       tempA is exact since intPartA and intPartB only have 22 bits each.
       tempB will have 52 bits of precision.
         */
        double tempA = intPartA * fracPartA;
        double tempB = intPartA * fracPartB + intPartB * fracPartA + intPartB * fracPartB;

        /* Compute the result.  (1+z)(tempA+tempB).  Order of operations is
       important.  For accuracy add by increasing size.  tempA is exact and
       much larger than the others.  If there are extra bits specified from the
       pow() function, use them. */
        final double tempC = tempB + tempA;

        // If tempC is positive infinite, the evaluation below could result in NaN,
        // because z could be negative at the same time.
        if (tempC == Double.POSITIVE_INFINITY) {
            if (hiPrec != null) {
                hiPrec[0] = Double.POSITIVE_INFINITY;
                hiPrec[1] = 0.0;
            }
            return Double.POSITIVE_INFINITY;
        }

        final double result;
        if (extra != 0.0) {
            result = tempC*extra*z + tempC*extra + tempC*z + tempB + tempA;
        } else {
            result = tempC*z + tempB + tempA;
        }

        if (hiPrec != null) {
            // If requesting high precision
            hiPrec[0] = tempA;
            hiPrec[1] = tempC*extra*z + tempC*extra + tempC*z + tempB;
        }

        return result;
    }

    /** Compute exp(x) - 1
     * @param x number to compute shifted exponential
     * @return exp(x) - 1
     */
    public static double expm1(double x) {
      return expm1(x, null);
    }

    /** Internal helper method for expm1
     * @param x number to compute shifted exponential
     * @param hiPrecOut receive high precision result for -1.0 < x < 1.0
     * @return exp(x) - 1
     */
    private static double expm1(double x, double[] hiPrecOut) {
        if (Double.isNaN(x) || x == 0.0) { // NaN or zero
            return x;
        }

        if (x <= -1.0 || x >= 1.0) {
            // If not between +/- 1.0
            //return exp(x) - 1.0;
            double[] hiPrec = new double[2];
            exp(x, 0.0, hiPrec);
            if (x > 0.0) {
                return -1.0 + hiPrec[0] + hiPrec[1];
            } else {
                final double ra = -1.0 + hiPrec[0];
                double rb = -(ra + 1.0 - hiPrec[0]);
                rb += hiPrec[1];
                return ra + rb;
            }
        }

        double baseA;
        double baseB;
        double epsilon;
        boolean negative = false;

        if (x < 0.0) {
            x = -x;
            negative = true;
        }

        {
            int intFrac = (int) (x * 1024.0);
            double tempA = ExpFracTable.EXP_FRAC_TABLE_A[intFrac] - 1.0;
            double tempB = ExpFracTable.EXP_FRAC_TABLE_B[intFrac];

            double temp = tempA + tempB;
            tempB = -(temp - tempA - tempB);
            tempA = temp;

            temp = tempA * HEX_40000000;
            baseA = tempA + temp - temp;
            baseB = tempB + (tempA - baseA);

            epsilon = x - intFrac/1024.0;
        }


        /* Compute expm1(epsilon) */
        double zb = 0.008336750013465571;
        zb = zb * epsilon + 0.041666663879186654;
        zb = zb * epsilon + 0.16666666666745392;
        zb = zb * epsilon + 0.49999999999999994;
        zb *= epsilon;
        zb *= epsilon;

        double za = epsilon;
        double temp = za + zb;
        zb = -(temp - za - zb);
        za = temp;

        temp = za * HEX_40000000;
        temp = za + temp - temp;
        zb += za - temp;
        za = temp;

        /* Combine the parts.   expm1(a+b) = expm1(a) + expm1(b) + expm1(a)*expm1(b) */
        double ya = za * baseA;
        //double yb = za*baseB + zb*baseA + zb*baseB;
        temp = ya + za * baseB;
        double yb = -(temp - ya - za * baseB);
        ya = temp;

        temp = ya + zb * baseA;
        yb += -(temp - ya - zb * baseA);
        ya = temp;

        temp = ya + zb * baseB;
        yb += -(temp - ya - zb*baseB);
        ya = temp;

        //ya = ya + za + baseA;
        //yb = yb + zb + baseB;
        temp = ya + baseA;
        yb += -(temp - baseA - ya);
        ya = temp;

        temp = ya + za;
        //yb += (ya > za) ? -(temp - ya - za) : -(temp - za - ya);
        yb += -(temp - ya - za);
        ya = temp;

        temp = ya + baseB;
        //yb += (ya > baseB) ? -(temp - ya - baseB) : -(temp - baseB - ya);
        yb += -(temp - ya - baseB);
        ya = temp;

        temp = ya + zb;
        //yb += (ya > zb) ? -(temp - ya - zb) : -(temp - zb - ya);
        yb += -(temp - ya - zb);
        ya = temp;

        if (negative) {
            /* Compute expm1(-x) = -expm1(x) / (expm1(x) + 1) */
            double denom = 1.0 + ya;
            double denomr = 1.0 / denom;
            double denomb = -(denom - 1.0 - ya) + yb;
            double ratio = ya * denomr;
            temp = ratio * HEX_40000000;
            final double ra = ratio + temp - temp;
            double rb = ratio - ra;

            temp = denom * HEX_40000000;
            za = denom + temp - temp;
            zb = denom - za;

            rb += (ya - za * ra - za * rb - zb * ra - zb * rb) * denomr;

            // f(x) = x/1+x
            // Compute f'(x)
            // Product rule:  d(uv) = du*v + u*dv
            // Chain rule:  d(f(g(x)) = f'(g(x))*f(g'(x))
            // d(1/x) = -1/(x*x)
            // d(1/1+x) = -1/( (1+x)^2) *  1 =  -1/((1+x)*(1+x))
            // d(x/1+x) = -x/((1+x)(1+x)) + 1/1+x = 1 / ((1+x)(1+x))

            // Adjust for yb
            rb += yb * denomr;                      // numerator
            rb += -ya * denomb * denomr * denomr;   // denominator

            // negate
            ya = -ra;
            yb = -rb;
        }

        if (hiPrecOut != null) {
            hiPrecOut[0] = ya;
            hiPrecOut[1] = yb;
        }

        return ya + yb;
    }

    /**
     * Natural logarithm.<br>
     * <br>
     *
     * Delegates to {@link Math#log(double)}
     *
     * @param x   a double
     * @return log(x)
     */
    public static double log(final double x) {
        return StrictMath.log(x);
    }

    /**
     * Computes log(1 + x).<br>
     * <br>
     *
     * Delegates to {@link Math#log1p(double)}
     *
     * @param x Number.
     * @return {@code log(1 + x)}.
     */
    public static double log1p(final double x) {
        return StrictMath.log1p(x);
    }

    /** Compute the base 10 logarithm.<br>
     * <br>
     *
     * Delegates to {@link Math#log10(double)}
     *
     * @param x a number
     * @return log10(x)
     */
    public static double log10(final double x) {
        return StrictMath.log10(x);
    }

    /**
     * Computes the <a href="http://mathworld.wolfram.com/Logarithm.html">
     * logarithm</a> in a given base.
     *
     * Returns {@code NaN} if either argument is negative.
     * If {@code base} is 0 and {@code x} is positive, 0 is returned.
     * If {@code base} is positive and {@code x} is 0,
     * {@code Double.NEGATIVE_INFINITY} is returned.
     * If both arguments are 0, the result is {@code NaN}.
     *
     * @param base Base of the logarithm, must be greater than 0.
     * @param x Argument, must be greater than 0.
     * @return the value of the logarithm, i.e. the number {@code y} such that
     * <code>base<sup>y</sup> = x</code>.
     */
    public static double log(double base, double x) {
        return log(x) / log(base);
    }

    /**
     * Power function.  Compute x^y.<br>
     * <br>
     *
     * Delegates to {@link Math#pow(double, double)}
     *
     * @param x   a double
     * @param y   a double
     * @return double
     */
    public static double pow(final double x, final double y) {
        return StrictMath.pow(x, y);
    }

    /**
     * Raise a double to an int power.
     *
     * @param d Number to raise.
     * @param e Exponent.
     * @return d<sup>e</sup>
     */
    public static double pow(double d, int e) {
        return pow(d, (long) e);
    }

    /**
     * Raise a double to a long power.
     *
     * @param d Number to raise.
     * @param e Exponent.
     * @return d<sup>e</sup>
     */
    public static double pow(double d, long e) {
        if (e == 0) {
            return 1.0;
        } else if (e > 0) {
            return new Split(d).pow(e).full;
        } else {
            return new Split(d).reciprocal().pow(-e).full;
        }
    }

    /** Class operator on double numbers split into one 26 bits number and one 27 bits number. */
    private static class Split {

        /** Split version of NaN. */
        public static final Split NAN = new Split(Double.NaN, 0);

        /** Split version of positive infinity. */
        public static final Split POSITIVE_INFINITY = new Split(Double.POSITIVE_INFINITY, 0);

        /** Split version of negative infinity. */
        public static final Split NEGATIVE_INFINITY = new Split(Double.NEGATIVE_INFINITY, 0);

        /** Full number. */
        private final double full;

        /** High order bits. */
        private final double high;

        /** Low order bits. */
        private final double low;

        /** Simple constructor.
         * @param x number to split
         */
        Split(final double x) {
            full = x;
            high = Double.longBitsToDouble(Double.doubleToRawLongBits(x) & ((-1L) << 27));
            low  = x - high;
        }

        /** Simple constructor.
         * @param high high order bits
         * @param low low order bits
         */
        Split(final double high, final double low) {
            this(high == 0.0 ? (low == 0.0 && Double.doubleToRawLongBits(high) == Long.MIN_VALUE /* negative zero */ ? -0.0 : low) : high + low, high, low);
        }

        /** Simple constructor.
         * @param full full number
         * @param high high order bits
         * @param low low order bits
         */
        Split(final double full, final double high, final double low) {
            this.full = full;
            this.high = high;
            this.low  = low;
        }

        /** Multiply the instance by another one.
         * @param b other instance to multiply by
         * @return product
         */
        public Split multiply(final Split b) {
            // beware the following expressions must NOT be simplified, they rely on floating point arithmetic properties
            final Split  mulBasic  = new Split(full * b.full);
            final double mulError  = low * b.low - (((mulBasic.full - high * b.high) - low * b.high) - high * b.low);
            return new Split(mulBasic.high, mulBasic.low + mulError);
        }

        /** Compute the reciprocal of the instance.
         * @return reciprocal of the instance
         */
        public Split reciprocal() {

            final double approximateInv = 1.0 / full;
            final Split  splitInv       = new Split(approximateInv);

            // if 1.0/d were computed perfectly, remultiplying it by d should give 1.0
            // we want to estimate the error so we can fix the low order bits of approximateInvLow
            // beware the following expressions must NOT be simplified, they rely on floating point arithmetic properties
            final Split product = multiply(splitInv);
            final double error  = (product.high - 1) + product.low;

            // better accuracy estimate of reciprocal
            return Double.isNaN(error) ? splitInv : new Split(splitInv.high, splitInv.low - error / full);

        }

        /** Computes this^e.
         * @param e exponent (beware, here it MUST be > 0; the only exclusion is Long.MIN_VALUE)
         * @return d^e, split in high and low bits
         */
        private Split pow(final long e) {

            // prepare result
            Split result = new Split(1);

            // d^(2p)
            Split d2p = new Split(full, high, low);

            for (long p = e; p != 0; p >>>= 1) {

                if ((p & 0x1) != 0) {
                    // accurate multiplication result = result * d^(2p) using Veltkamp TwoProduct algorithm
                    result = result.multiply(d2p);
                }

                // accurate squaring d^(2(p+1)) = d^(2p) * d^(2p) using Veltkamp TwoProduct algorithm
                d2p = d2p.multiply(d2p);

            }

            if (Double.isNaN(result.full)) {
                if (Double.isNaN(full)) {
                    return NAN;
                } else {
                    // some intermediate numbers exceeded capacity,
                    // and the low order bits became NaN (because infinity - infinity = NaN)
                    if (abs(full) < 1) {
                        return new Split(copySign(0.0, full), 0.0);
                    } else if (full < 0 && (e & 0x1) == 1) {
                        return NEGATIVE_INFINITY;
                    } else {
                        return POSITIVE_INFINITY;
                    }
                }
            } else {
                return result;
            }

        }

    }

    /**
     * Sine function.<br>
     * <br>
     *
     * Delegates to {@link Math#sin(double)}
     *
     * @param x Argument.
     * @return sin(x)
     */
    public static double sin(double x) {
        return StrictMath.sin(x);
    }

    /**
     * Cosine function.<br>
     * <br>
     *
     * Delegates to {@link Math#cos(double)}
     *
     * @param x Argument.
     * @return cos(x)
     */
    public static double cos(double x) {
        return StrictMath.cos(x);
    }

    /**
     * Combined Sine and Cosine function.<br>
     * <br>
     *
     * Delegates to {@link Math#sin(double)} and {@link Math#cos(double)}
     *
     * @param x Argument.
     * @return [sin(x), cos(x)]
     */
    public static SinCos sinCos(double x) {
        return new SinCos(sin(x), cos(x));
    }

    /**
     * Combined Sine and Cosine function.
     *
     * @param x Argument.
     * @param <T> the type of the field element
     * @return [sin(x), cos(x)]
     * @since 1.4
     */
    public static <T extends CalculusFieldElement<T>> FieldSinCos<T> sinCos(T x) {
        return x.sinCos();
    }

    /**
     * Tangent function.
     *
     * @param x Argument.
     * @return tan(x)
     */
    public static double tan(double x) {
        return StrictMath.tan(x);
    }

    /**
     * Arctangent function
     *  @param x a number
     *  @return atan(x)
     */
    public static double atan(double x) {
        return atan(x, 0.0, false);
    }

    /** Internal helper function to compute arctangent.
     * @param xa number from which arctangent is requested
     * @param xb extra bits for x (may be 0.0)
     * @param leftPlane if true, result angle must be put in the left half plane
     * @return atan(xa + xb) (or angle shifted by {@code PI} if leftPlane is true)
     */
    private static double atan(double xa, double xb, boolean leftPlane) {
        if (xa == 0.0) { // Matches +/- 0.0; return correct sign
            return leftPlane ? copySign(Math.PI, xa) : xa;
        }

        final boolean negate;
        if (xa < 0) {
            // negative
            xa = -xa;
            xb = -xb;
            negate = true;
        } else {
            negate = false;
        }

        if (xa > 1.633123935319537E16) { // Very large input
            return (negate ^ leftPlane) ? (-Math.PI * F_1_2) : (Math.PI * F_1_2);
        }

        /* Estimate the closest tabulated arctan value, compute eps = xa-tangentTable */
        final int idx;
        if (xa < 1) {
            idx = (int) (((-1.7168146928204136 * xa * xa + 8.0) * xa) + 0.5);
        } else {
            final double oneOverXa = 1 / xa;
            idx = (int) (-((-1.7168146928204136 * oneOverXa * oneOverXa + 8.0) * oneOverXa) + 13.07);
        }

        final double ttA = TANGENT_TABLE_A[idx];
        final double ttB = TANGENT_TABLE_B[idx];

        double epsA = xa - ttA;
        double epsB = -(epsA - xa + ttA);
        epsB += xb - ttB;

        double temp = epsA + epsB;
        epsB = -(temp - epsA - epsB);
        epsA = temp;

        /* Compute eps = eps / (1.0 + xa*tangent) */
        temp = xa * HEX_40000000;
        double ya = xa + temp - temp;
        double yb = xb + xa - ya;
        xa = ya;
        xb += yb;

        //if (idx > 8 || idx == 0)
        if (idx == 0) {
            /* If the slope of the arctan is gentle enough (< 0.45), this approximation will suffice */
            //double denom = 1.0 / (1.0 + xa*tangentTableA[idx] + xb*tangentTableA[idx] + xa*tangentTableB[idx] + xb*tangentTableB[idx]);
            final double denom = 1d / (1d + (xa + xb) * (ttA + ttB));
            //double denom = 1.0 / (1.0 + xa*tangentTableA[idx]);
            ya = epsA * denom;
            yb = epsB * denom;
        } else {
            double temp2 = xa * ttA;
            double za = 1d + temp2;
            double zb = -(za - 1d - temp2);
            temp2 = xb * ttA + xa * ttB;
            temp = za + temp2;
            zb += -(temp - za - temp2);
            za = temp;

            zb += xb * ttB;
            ya = epsA / za;

            temp = ya * HEX_40000000;
            final double yaa = (ya + temp) - temp;
            final double yab = ya - yaa;

            temp = za * HEX_40000000;
            final double zaa = (za + temp) - temp;
            final double zab = za - zaa;

            /* Correct for rounding in division */
            yb = (epsA - yaa * zaa - yaa * zab - yab * zaa - yab * zab) / za;

            yb += -epsA * zb / za / za;
            yb += epsB / za;
        }


        epsA = ya;
        epsB = yb;

        /* Evaluate polynomial */
        final double epsA2 = epsA * epsA;

        /*
    yb = -0.09001346640161823;
    yb = yb * epsA2 + 0.11110718400605211;
    yb = yb * epsA2 + -0.1428571349122913;
    yb = yb * epsA2 + 0.19999999999273194;
    yb = yb * epsA2 + -0.33333333333333093;
    yb = yb * epsA2 * epsA;
         */

        yb = 0.07490822288864472;
        yb = yb * epsA2 - 0.09088450866185192;
        yb = yb * epsA2 + 0.11111095942313305;
        yb = yb * epsA2 - 0.1428571423679182;
        yb = yb * epsA2 + 0.19999999999923582;
        yb = yb * epsA2 - 0.33333333333333287;
        yb = yb * epsA2 * epsA;


        ya = epsA;

        temp = ya + yb;
        yb = -(temp - ya - yb);
        ya = temp;

        /* Add in effect of epsB.   atan'(x) = 1/(1+x^2) */
        yb += epsB / (1d + epsA * epsA);

        final double eighths = EIGHTHS[idx];

        //result = yb + eighths[idx] + ya;
        double za = eighths + ya;
        double zb = -(za - eighths - ya);
        temp = za + yb;
        zb += -(temp - za - yb);
        za = temp;

        double result = za + zb;

        if (leftPlane) {
            // Result is in the left plane
            final double resultb = -(result - za - zb);
            final double pia = 1.5707963267948966 * 2;
            final double pib = 6.123233995736766E-17 * 2;

            za = pia - result;
            zb = -(za - pia + result);
            zb += pib - resultb;

            result = za + zb;
        }


        if (negate ^ leftPlane) {
            result = -result;
        }

        return result;
    }

    /**
     * Two arguments arctangent function
     * @param y ordinate
     * @param x abscissa
     * @return phase angle of point (x,y) between {@code -PI} and {@code PI}
     */
    public static double atan2(double y, double x) {
        if (Double.isNaN(x) || Double.isNaN(y)) {
            return Double.NaN;
        }

        if (y == 0) {
            final double result = x * y;
            final double invx = 1d / x;

            if (invx == 0) { // X is infinite
                if (x > 0) {
                    return y; // return +/- 0.0
                } else {
                    return copySign(Math.PI, y);
                }
            }

            if (x < 0 || invx < 0) {
                return copySign(Math.PI, y);
            } else {
                return result;
            }
        }

        // y cannot now be zero

        if (y == Double.POSITIVE_INFINITY) {
            if (x == Double.POSITIVE_INFINITY) {
                return Math.PI * F_1_4;
            }

            if (x == Double.NEGATIVE_INFINITY) {
                return Math.PI * F_3_4;
            }

            return Math.PI * F_1_2;
        }

        if (y == Double.NEGATIVE_INFINITY) {
            if (x == Double.POSITIVE_INFINITY) {
                return -Math.PI * F_1_4;
            }

            if (x == Double.NEGATIVE_INFINITY) {
                return -Math.PI * F_3_4;
            }

            return -Math.PI * F_1_2;
        }

        if (x == Double.POSITIVE_INFINITY) {
            return copySign(0d, y);
        }

        if (x == Double.NEGATIVE_INFINITY)
        {
            return copySign(Math.PI, y);
        }

        // Neither y nor x can be infinite or NAN here

        if (x == 0) {
            return copySign(Math.PI * F_1_2, y);
        }

        // Compute ratio r = y/x
        final double r = y / x;
        if (Double.isInfinite(r)) { // bypass calculations that can create NaN
            return atan(r, 0, x < 0);
        }

        double ra = doubleHighPart(r);
        double rb = r - ra;

        // Split x
        final double xa = doubleHighPart(x);
        final double xb = x - xa;

        rb += (y - ra * xa - ra * xb - rb * xa - rb * xb) / x;

        final double temp = ra + rb;
        rb = -(temp - ra - rb);
        ra = temp;

        if (ra == 0) { // Fix up the sign so atan works correctly
            ra = copySign(0d, y);
        }

        // Call atan
        return atan(ra, rb, x < 0);

    }

    /** Compute the arc sine of a number.
     * @param x number on which evaluation is done
     * @return arc sine of x
     */
    public static double asin(double x) {
      if (Double.isNaN(x)) {
          return Double.NaN;
      }

      if (x > 1.0 || x < -1.0) {
          return Double.NaN;
      }

      if (x == 1.0) {
          return Math.PI/2.0;
      }

      if (x == -1.0) {
          return -Math.PI/2.0;
      }

      if (x == 0.0) { // Matches +/- 0.0; return correct sign
          return x;
      }

      /* Compute asin(x) = atan(x/sqrt(1-x*x)) */

      /* Split x */
      double temp = x * HEX_40000000;
      final double xa = x + temp - temp;
      final double xb = x - xa;

      /* Square it */
      double ya = xa*xa;
      double yb = xa*xb*2.0 + xb*xb;

      /* Subtract from 1 */
      ya = -ya;
      yb = -yb;

      double za = 1.0 + ya;
      double zb = -(za - 1.0 - ya);

      temp = za + yb;
      zb += -(temp - za - yb);
      za = temp;

      /* Square root */
      double y;
      y = sqrt(za);
      temp = y * HEX_40000000;
      ya = y + temp - temp;
      yb = y - ya;

      /* Extend precision of sqrt */
      yb += (za - ya*ya - 2*ya*yb - yb*yb) / (2.0*y);

      /* Contribution of zb to sqrt */
      double dx = zb / (2.0*y);

      // Compute ratio r = x/y
      double r = x/y;
      temp = r * HEX_40000000;
      double ra = r + temp - temp;
      double rb = r - ra;

      rb += (x - ra*ya - ra*yb - rb*ya - rb*yb) / y;  // Correct for rounding in division
      rb += -x * dx / y / y;  // Add in effect additional bits of sqrt.

      temp = ra + rb;
      rb = -(temp - ra - rb);
      ra = temp;

      return atan(ra, rb, false);
    }

    /** Compute the arc cosine of a number.
     * @param x number on which evaluation is done
     * @return arc cosine of x
     */
    public static double acos(double x) {
      if (Double.isNaN(x)) {
          return Double.NaN;
      }

      if (x > 1.0 || x < -1.0) {
          return Double.NaN;
      }

      if (x == -1.0) {
          return Math.PI;
      }

      if (x == 1.0) {
          return 0.0;
      }

      if (x == 0) {
          return Math.PI/2.0;
      }

      /* Compute acos(x) = atan(sqrt(1-x*x)/x) */

      /* Split x */
      double temp = x * HEX_40000000;
      final double xa = x + temp - temp;
      final double xb = x - xa;

      /* Square it */
      double ya = xa*xa;
      double yb = xa*xb*2.0 + xb*xb;

      /* Subtract from 1 */
      ya = -ya;
      yb = -yb;

      double za = 1.0 + ya;
      double zb = -(za - 1.0 - ya);

      temp = za + yb;
      zb += -(temp - za - yb);
      za = temp;

      /* Square root */
      double y = sqrt(za);
      temp = y * HEX_40000000;
      ya = y + temp - temp;
      yb = y - ya;

      /* Extend precision of sqrt */
      yb += (za - ya*ya - 2*ya*yb - yb*yb) / (2.0*y);

      /* Contribution of zb to sqrt */
      yb += zb / (2.0*y);
      y = ya+yb;
      yb = -(y - ya - yb);

      // Compute ratio r = y/x
      double r = y/x;

      // Did r overflow?
      if (Double.isInfinite(r)) { // x is effectively zero
          return Math.PI/2; // so return the appropriate value
      }

      double ra = doubleHighPart(r);
      double rb = r - ra;

      rb += (y - ra*xa - ra*xb - rb*xa - rb*xb) / x;  // Correct for rounding in division
      rb += yb / x;  // Add in effect additional bits of sqrt.

      temp = ra + rb;
      rb = -(temp - ra - rb);
      ra = temp;

      return atan(ra, rb, x<0);
    }

    /** Compute the cubic root of a number.<br>
     * <br>
     *
     * Delegates to {@link Math#cbrt(double)}
     *
     * @param x number on which evaluation is done
     * @return cubic root of x
     */
    public static double cbrt(double x) {
        return StrictMath.cbrt(x);
    }

    /**
     *  Convert degrees to radians, with error of less than 0.5 ULP
     *  @param x angle in degrees
     *  @return x converted into radians
     */
    public static double toRadians(double x)
    {
        if (Double.isInfinite(x) || x == 0.0) { // Matches +/- 0.0; return correct sign
            return x;
        }

        // These are PI/180 split into high and low order bits
        final double facta = 0.01745329052209854;
        final double factb = 1.997844754509471E-9;

        double xa = doubleHighPart(x);
        double xb = x - xa;

        double result = xb * factb + xb * facta + xa * factb + xa * facta;
        if (result == 0) {
            result *= x; // ensure correct sign if calculation underflows
        }
        return result;
    }

    /**
     *  Convert radians to degrees, with error of less than 0.5 ULP
     *  @param x angle in radians
     *  @return x converted into degrees
     */
    public static double toDegrees(double x)
    {
        if (Double.isInfinite(x) || x == 0.0) { // Matches +/- 0.0; return correct sign
            return x;
        }

        // These are 180/PI split into high and low order bits
        final double facta = 57.2957763671875;
        final double factb = 3.145894820876798E-6;

        double xa = doubleHighPart(x);
        double xb = x - xa;

        return xb * factb + xb * facta + xa * factb + xa * facta;
    }

    /**
     * Absolute value.
     * @param x number from which absolute value is requested
     * @return abs(x)
     */
    public static int abs(final int x) {
        final int i = x >>> 31;
        return (x ^ (~i + 1)) + i;
    }

    /**
     * Absolute value.
     * @param x number from which absolute value is requested
     * @return abs(x)
     */
    public static long abs(final long x) {
        final long l = x >>> 63;
        // l is one if x negative zero else
        // ~l+1 is zero if x is positive, -1 if x is negative
        // x^(~l+1) is x is x is positive, ~x if x is negative
        // add around
        return (x ^ (~l + 1)) + l;
    }

    /**
     * Absolute value.
     * @param x number from which absolute value is requested
     * @return abs(x), or throws an exception for {@code Integer.MIN_VALUE}
     * @since 2.0
     */
    public static int absExact(final int x) {
        if (x == Integer.MIN_VALUE) {
            throw new ArithmeticException();
        }
        return abs(x);
    }

    /**
     * Absolute value.
     * @param x number from which absolute value is requested
     * @return abs(x), or throws an exception for {@code Long.MIN_VALUE}
     * @since 2.0
     */
    public static long absExact(final long x) {
        if (x == Long.MIN_VALUE) {
            throw new ArithmeticException();
        }
        return abs(x);
    }

    /**
     * Absolute value.
     * @param x number from which absolute value is requested
     * @return abs(x)
     * @since 2.0
     */
    public static float abs(final float x) {
        return Float.intBitsToFloat(MASK_NON_SIGN_INT & Float.floatToRawIntBits(x));
    }

    /**
     * Absolute value.<br>
     * <br>
     *
     * Delegates to {@link Math#abs(double)}
     *
     * @param x number from which absolute value is requested
     * @return abs(x)
     */
    public static double abs(double x) {
        return StrictMath.abs(x);
    }

    /**
     * Negates the argument.
     * @param x number from which opposite value is requested
     * @return -x, or throws an exception for {@code Integer.MIN_VALUE}
     * @since 2.0
     */
    public static int negateExact(final int x) {
        if (x == Integer.MIN_VALUE) {
            throw new ArithmeticException();
        }
        return -x;
    }

    /**
     * Negates the argument.
     * @param x number from which opposite value is requested
     * @return -x, or throws an exception for {@code Long.MIN_VALUE}
     * @since 2.0
     */
    public static long negateExact(final long x) {
        if (x == Long.MIN_VALUE) {
            throw new ArithmeticException();
        }
        return -x;
    }

    /**
     * Compute least significant bit (Unit in Last Position) for a number.
     * @param x number from which ulp is requested
     * @return ulp(x)
     */
    public static double ulp(double x) {
        if (Double.isInfinite(x)) {
            return Double.POSITIVE_INFINITY;
        }
        return abs(x - Double.longBitsToDouble(Double.doubleToRawLongBits(x) ^ 1));
    }

    /**
     * Compute least significant bit (Unit in Last Position) for a number.
     * @param x number from which ulp is requested
     * @return ulp(x)
     */
    public static float ulp(float x) {
        if (Float.isInfinite(x)) {
            return Float.POSITIVE_INFINITY;
        }
        return abs(x - Float.intBitsToFloat(Float.floatToRawIntBits(x) ^ 1));
    }

    /**
     * Multiply a double number by a power of 2.
     * @param d number to multiply
     * @param n power of 2
     * @return d &times; 2<sup>n</sup>
     */
    public static double scalb(final double d, final int n) {

        // first simple and fast handling when 2^n can be represented using normal numbers
        if ((n > -1023) && (n < 1024)) {
            return d * Double.longBitsToDouble(((long) (n + 1023)) << 52);
        }

        // handle special cases
        if (Double.isNaN(d) || Double.isInfinite(d) || (d == 0)) {
            return d;
        }
        if (n < -2098) {
            return (d > 0) ? 0.0 : -0.0;
        }
        if (n > 2097) {
            return (d > 0) ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
        }

        // decompose d
        final long bits = Double.doubleToRawLongBits(d);
        final long sign = bits & 0x8000000000000000L;
        int  exponent   = ((int) (bits >>> 52)) & 0x7ff;
        long mantissa   = bits & 0x000fffffffffffffL;

        // compute scaled exponent
        int scaledExponent = exponent + n;

        if (n < 0) {
            // we are really in the case n <= -1023
            if (scaledExponent > 0) {
                // both the input and the result are normal numbers, we only adjust the exponent
                return Double.longBitsToDouble(sign | (((long) scaledExponent) << 52) | mantissa);
            } else if (scaledExponent > -53) {
                // the input is a normal number and the result is a subnormal number

                // recover the hidden mantissa bit
                mantissa |= 1L << 52;

                // scales down complete mantissa, hence losing least significant bits
                final long mostSignificantLostBit = mantissa & (1L << (-scaledExponent));
                mantissa >>>= 1 - scaledExponent;
                if (mostSignificantLostBit != 0) {
                    // we need to add 1 bit to round up the result
                    mantissa++;
                }
                return Double.longBitsToDouble(sign | mantissa);

            } else {
                // no need to compute the mantissa, the number scales down to 0
                return (sign == 0L) ? 0.0 : -0.0;
            }
        } else {
            // we are really in the case n >= 1024
            if (exponent == 0) {

                // the input number is subnormal, normalize it
                while ((mantissa >>> 52) != 1) {
                    mantissa <<= 1;
                    --scaledExponent;
                }
                ++scaledExponent;
                mantissa &= 0x000fffffffffffffL;

                if (scaledExponent < 2047) {
                    return Double.longBitsToDouble(sign | (((long) scaledExponent) << 52) | mantissa);
                } else {
                    return (sign == 0L) ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
                }

            } else if (scaledExponent < 2047) {
                return Double.longBitsToDouble(sign | (((long) scaledExponent) << 52) | mantissa);
            } else {
                return (sign == 0L) ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
            }
        }

    }

    /**
     * Multiply a float number by a power of 2.
     * @param f number to multiply
     * @param n power of 2
     * @return f &times; 2<sup>n</sup>
     */
    public static float scalb(final float f, final int n) {

        // first simple and fast handling when 2^n can be represented using normal numbers
        if ((n > -127) && (n < 128)) {
            return f * Float.intBitsToFloat((n + 127) << 23);
        }

        // handle special cases
        if (Float.isNaN(f) || Float.isInfinite(f) || (f == 0f)) {
            return f;
        }
        if (n < -277) {
            return (f > 0) ? 0.0f : -0.0f;
        }
        if (n > 276) {
            return (f > 0) ? Float.POSITIVE_INFINITY : Float.NEGATIVE_INFINITY;
        }

        // decompose f
        final int bits = Float.floatToIntBits(f);
        final int sign = bits & 0x80000000;
        int  exponent  = (bits >>> 23) & 0xff;
        int mantissa   = bits & 0x007fffff;

        // compute scaled exponent
        int scaledExponent = exponent + n;

        if (n < 0) {
            // we are really in the case n <= -127
            if (scaledExponent > 0) {
                // both the input and the result are normal numbers, we only adjust the exponent
                return Float.intBitsToFloat(sign | (scaledExponent << 23) | mantissa);
            } else if (scaledExponent > -24) {
                // the input is a normal number and the result is a subnormal number

                // recover the hidden mantissa bit
                mantissa |= 1 << 23;

                // scales down complete mantissa, hence losing least significant bits
                final int mostSignificantLostBit = mantissa & (1 << (-scaledExponent));
                mantissa >>>= 1 - scaledExponent;
                if (mostSignificantLostBit != 0) {
                    // we need to add 1 bit to round up the result
                    mantissa++;
                }
                return Float.intBitsToFloat(sign | mantissa);

            } else {
                // no need to compute the mantissa, the number scales down to 0
                return (sign == 0) ? 0.0f : -0.0f;
            }
        } else {
            // we are really in the case n >= 128
            if (exponent == 0) {

                // the input number is subnormal, normalize it
                while ((mantissa >>> 23) != 1) {
                    mantissa <<= 1;
                    --scaledExponent;
                }
                ++scaledExponent;
                mantissa &= 0x007fffff;

                if (scaledExponent < 255) {
                    return Float.intBitsToFloat(sign | (scaledExponent << 23) | mantissa);
                } else {
                    return (sign == 0) ? Float.POSITIVE_INFINITY : Float.NEGATIVE_INFINITY;
                }

            } else if (scaledExponent < 255) {
                return Float.intBitsToFloat(sign | (scaledExponent << 23) | mantissa);
            } else {
                return (sign == 0) ? Float.POSITIVE_INFINITY : Float.NEGATIVE_INFINITY;
            }
        }

    }

    /**
     * Get the next machine representable number after a number, moving
     * in the direction of another number.
     * <p>
     * The ordering is as follows (increasing):
     * </p>
     * <ul>
     * <li>-INFINITY</li>
     * <li>-MAX_VALUE</li>
     * <li>-MIN_VALUE</li>
     * <li>-0.0</li>
     * <li>+0.0</li>
     * <li>+MIN_VALUE</li>
     * <li>+MAX_VALUE</li>
     * <li>+INFINITY</li>
     * </ul>
     * <p>
     * If arguments compare equal, then the second argument is returned.
     * </p>
     * <p>
     * If {@code direction} is greater than {@code d},
     * the smallest machine representable number strictly greater than
     * {@code d} is returned; if less, then the largest representable number
     * strictly less than {@code d} is returned.
     * </p>
     * <p>
     * If {@code d} is infinite and direction does not
     * bring it back to finite numbers, it is returned unchanged.
     * </p>
     *
     * @param d base number
     * @param direction (the only important thing is whether
     * {@code direction} is greater or smaller than {@code d})
     * @return the next machine representable number in the specified direction
     */
    public static double nextAfter(double d, double direction) {

        // handling of some important special cases
        if (Double.isNaN(d) || Double.isNaN(direction)) {
            return Double.NaN;
        } else if (d == direction) {
            return direction;
        } else if (Double.isInfinite(d)) {
            return (d < 0) ? -Double.MAX_VALUE : Double.MAX_VALUE;
        } else if (d == 0) {
            return (direction < 0) ? -Double.MIN_VALUE : Double.MIN_VALUE;
        }
        // special cases MAX_VALUE to infinity and  MIN_VALUE to 0
        // are handled just as normal numbers
        // can use raw bits since already dealt with infinity and NaN
        final long bits = Double.doubleToRawLongBits(d);
        final long sign = bits & 0x8000000000000000L;
        if ((direction < d) ^ (sign == 0L)) {
            return Double.longBitsToDouble(sign | ((bits & 0x7fffffffffffffffL) + 1));
        } else {
            return Double.longBitsToDouble(sign | ((bits & 0x7fffffffffffffffL) - 1));
        }

    }

    /**
     * Get the next machine representable number after a number, moving
     * in the direction of another number.
     * <p>* The ordering is as follows (increasing):</p>
     * <ul>
     * <li>-INFINITY</li>
     * <li>-MAX_VALUE</li>
     * <li>-MIN_VALUE</li>
     * <li>-0.0</li>
     * <li>+0.0</li>
     * <li>+MIN_VALUE</li>
     * <li>+MAX_VALUE</li>
     * <li>+INFINITY</li>
     * </ul>
     * <p>
     * If arguments compare equal, then the second argument is returned.
     * </p>
     * <p>
     * If {@code direction} is greater than {@code f},
     * the smallest machine representable number strictly greater than
     * {@code f} is returned; if less, then the largest representable number
     * strictly less than {@code f} is returned.
     * </p>
     * <p>
     * If {@code f} is infinite and direction does not
     * bring it back to finite numbers, it is returned unchanged.
     * </p>
     *
     * @param f base number
     * @param direction (the only important thing is whether
     * {@code direction} is greater or smaller than {@code f})
     * @return the next machine representable number in the specified direction
     */
    public static float nextAfter(final float f, final double direction) {

        // handling of some important special cases
        if (Double.isNaN(f) || Double.isNaN(direction)) {
            return Float.NaN;
        } else if (f == direction) {
            return (float) direction;
        } else if (Float.isInfinite(f)) {
            return (f < 0f) ? -Float.MAX_VALUE : Float.MAX_VALUE;
        } else if (f == 0f) {
            return (direction < 0) ? -Float.MIN_VALUE : Float.MIN_VALUE;
        }
        // special cases MAX_VALUE to infinity and  MIN_VALUE to 0
        // are handled just as normal numbers

        final int bits = Float.floatToIntBits(f);
        final int sign = bits & 0x80000000;
        if ((direction < f) ^ (sign == 0)) {
            return Float.intBitsToFloat(sign | ((bits & 0x7fffffff) + 1));
        } else {
            return Float.intBitsToFloat(sign | ((bits & 0x7fffffff) - 1));
        }

    }

    /** Get the largest whole number smaller than x.<br>
     * <br>
     *
     * Delegates to {@link Math#floor(double)}
     *
     * @param x number from which floor is requested
     * @return a double number f such that f is an integer f &lt;= x &lt; f + 1.0
     */
    public static double floor(double x) {
        return StrictMath.floor(x);
    }

    /** Get the smallest whole number larger than x.<br>
     * <br>
     *
     * Delegates to {@link Math#ceil(double)}
     *
     * @param x number from which ceil is requested
     * @return a double number c such that c is an integer c - 1.0 &lt; x &lt;= c
     */
    public static double ceil(double x) {
        return StrictMath.ceil(x);
    }

    /** Get the whole number that is the nearest to x, or the even one if x is exactly half way between two integers.
     * @param x number from which nearest whole number is requested
     * @return a double number r such that r is an integer r - 0.5 &lt;= x &lt;= r + 0.5
     */
    public static double rint(double x) {
        double y = floor(x);
        double d = x - y;

        if (d > 0.5) {
            if (y == -1.0) {
                return -0.0; // Preserve sign of operand
            }
            return y+1.0;
        }
        if (d < 0.5) {
            return y;
        }

        /* half way, round to even */
        long z = (long) y;
        return (z & 1) == 0 ? y : y + 1.0;
    }

    /** Get the closest long to x.
     * @param x number from which closest long is requested
     * @return closest long to x
     */
    public static long round(double x) {
        final long bits = Double.doubleToRawLongBits(x);
        final int biasedExp = ((int)(bits>>52)) & 0x7ff;
        // Shift to get rid of bits past comma except first one: will need to
        // 1-shift to the right to end up with correct magnitude.
        final int shift = (52 - 1 + Double.MAX_EXPONENT) - biasedExp;
        if ((shift & -64) == 0) {
            // shift in [0,63], so unbiased exp in [-12,51].
            long extendedMantissa = 0x0010000000000000L | (bits & 0x000fffffffffffffL);
            if (bits < 0) {
                extendedMantissa = -extendedMantissa;
            }
            // If value is positive and first bit past comma is 0, rounding
            // to lower integer, else to upper one, which is what "+1" and
            // then ">>1" do.
            return ((extendedMantissa >> shift) + 1L) >> 1;
        } else {
            // +-Infinity, NaN, or a mathematical integer.
            return (long) x;
        }
    }

    /** Get the closest int to x.
     * @param x number from which closest int is requested
     * @return closest int to x
     */
    public static int round(final float x) {
        final int bits = Float.floatToRawIntBits(x);
        final int biasedExp = (bits>>23) & 0xff;
        // Shift to get rid of bits past comma except first one: will need to
        // 1-shift to the right to end up with correct magnitude.
        final int shift = (23 - 1 + Float.MAX_EXPONENT) - biasedExp;
        if ((shift & -32) == 0) {
            // shift in [0,31], so unbiased exp in [-9,22].
            int extendedMantissa = 0x00800000 | (bits & 0x007fffff);
            if (bits < 0) {
                extendedMantissa = -extendedMantissa;
            }
            // If value is positive and first bit past comma is 0, rounding
            // to lower integer, else to upper one, which is what "+1" and
            // then ">>1" do.
            return ((extendedMantissa >> shift) + 1) >> 1;
        } else {
            // +-Infinity, NaN, or a mathematical integer.
            return (int) x;
        }
    }

    /** Compute the minimum of two values
     * @param a first value
     * @param b second value
     * @return a if a is lesser or equal to b, b otherwise
     */
    public static int min(final int a, final int b) {
        return (a <= b) ? a : b;
    }

    /** Compute the minimum of two values
     * @param a first value
     * @param b second value
     * @return a if a is lesser or equal to b, b otherwise
     */
    public static long min(final long a, final long b) {
        return (a <= b) ? a : b;
    }

    /** Compute the minimum of two values
     * @param a first value
     * @param b second value
     * @return a if a is lesser or equal to b, b otherwise
     */
    public static float min(final float a, final float b) {
        if (a > b) {
            return b;
        }
        if (a < b) {
            return a;
        }
        /* if either arg is NaN, return NaN */
        if (a != b) {
            return Float.NaN;
        }
        /* min(+0.0,-0.0) == -0.0 */
        /* 0x80000000 == Float.floatToRawIntBits(-0.0d) */
        int bits = Float.floatToRawIntBits(a);
        if (bits == 0x80000000) {
            return a;
        }
        return b;
    }

    /** Compute the minimum of two values
     * @param a first value
     * @param b second value
     * @return a if a is lesser or equal to b, b otherwise
     */
    public static double min(final double a, final double b) {
        if (a > b) {
            return b;
        }
        if (a < b) {
            return a;
        }
        /* if either arg is NaN, return NaN */
        if (a != b) {
            return Double.NaN;
        }
        /* min(+0.0,-0.0) == -0.0 */
        /* 0x8000000000000000L == Double.doubleToRawLongBits(-0.0d) */
        long bits = Double.doubleToRawLongBits(a);
        if (bits == 0x8000000000000000L) {
            return a;
        }
        return b;
    }

    /** Compute the maximum of two values
     * @param a first value
     * @param b second value
     * @return b if a is lesser or equal to b, a otherwise
     */
    public static int max(final int a, final int b) {
        return (a <= b) ? b : a;
    }

    /** Compute the maximum of two values
     * @param a first value
     * @param b second value
     * @return b if a is lesser or equal to b, a otherwise
     */
    public static long max(final long a, final long b) {
        return (a <= b) ? b : a;
    }

    /** Compute the maximum of two values
     * @param a first value
     * @param b second value
     * @return b if a is lesser or equal to b, a otherwise
     */
    public static float max(final float a, final float b) {
        if (a > b) {
            return a;
        }
        if (a < b) {
            return b;
        }
        /* if either arg is NaN, return NaN */
        if (a != b) {
            return Float.NaN;
        }
        /* min(+0.0,-0.0) == -0.0 */
        /* 0x80000000 == Float.floatToRawIntBits(-0.0d) */
        int bits = Float.floatToRawIntBits(a);
        if (bits == 0x80000000) {
            return b;
        }
        return a;
    }

    /** Compute the maximum of two values
     * @param a first value
     * @param b second value
     * @return b if a is lesser or equal to b, a otherwise
     */
    public static double max(final double a, final double b) {
        if (a > b) {
            return a;
        }
        if (a < b) {
            return b;
        }
        /* if either arg is NaN, return NaN */
        if (a != b) {
            return Double.NaN;
        }
        /* min(+0.0,-0.0) == -0.0 */
        /* 0x8000000000000000L == Double.doubleToRawLongBits(-0.0d) */
        long bits = Double.doubleToRawLongBits(a);
        if (bits == 0x8000000000000000L) {
            return b;
        }
        return a;
    }

    /**
     * Returns the hypotenuse of a triangle with sides {@code x} and {@code y}
     * - sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)<br>
     * avoiding intermediate overflow or underflow.
     *
     * <ul>
     * <li> If either argument is infinite, then the result is positive infinity.</li>
     * <li> else, if either argument is NaN then the result is NaN.</li>
     * </ul>
     *
     * Delegates to {@link Math#hypot(double, double)}
     *
     * @param x a value
     * @param y a value
     * @return sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     */
    public static double hypot(final double x, final double y) {
        return StrictMath.hypot(x, y);
    }

    /**
     * Computes the remainder as prescribed by the IEEE 754 standard.
     * <p>
     * The remainder value is mathematically equal to {@code x - y*n}
     * where {@code n} is the mathematical integer closest to the exact mathematical value
     * of the quotient {@code x/y}.
     * If two mathematical integers are equally close to {@code x/y} then
     * {@code n} is the integer that is even.
     * </p>
     * <ul>
     * <li>If either operand is NaN, the result is NaN.</li>
     * <li>If the result is not NaN, the sign of the result equals the sign of the dividend.</li>
     * <li>If the dividend is an infinity, or the divisor is a zero, or both, the result is NaN.</li>
     * <li>If the dividend is finite and the divisor is an infinity, the result equals the dividend.</li>
     * <li>If the dividend is a zero and the divisor is finite, the result equals the dividend.</li>
     * </ul>
     * @param dividend the number to be divided
     * @param divisor the number by which to divide
     * @return the remainder, rounded
     */
    public static double IEEEremainder(final double dividend, final double divisor) {
        if (getExponent(dividend) == 1024 || getExponent(divisor) == 1024 || divisor == 0.0) {
            // we are in one of the special cases
            if (Double.isInfinite(divisor) && !Double.isInfinite(dividend)) {
                return dividend;
            } else {
                return Double.NaN;
            }
        } else {
            // we are in the general case
            final double n         = rint(dividend / divisor);
            final double remainder = Double.isInfinite(n) ? 0.0 : dividend - divisor * n;
            return (remainder == 0) ? copySign(remainder, dividend) : remainder;
        }
    }

    /** Convert a long to interger, detecting overflows
     * @param n number to convert to int
     * @return integer with same valie as n if no overflows occur
     * @exception MathRuntimeException if n cannot fit into an int
     */
    public static int toIntExact(final long n) throws MathRuntimeException {
        if (n < Integer.MIN_VALUE || n > Integer.MAX_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW);
        }
        return (int) n;
    }

    /** Increment a number, detecting overflows.
     * @param n number to increment
     * @return n+1 if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static int incrementExact(final int n) throws MathRuntimeException {

        if (n == Integer.MAX_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_ADDITION, n, 1);
        }

        return n + 1;

    }

    /** Increment a number, detecting overflows.
     * @param n number to increment
     * @return n+1 if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static long incrementExact(final long n) throws MathRuntimeException {

        if (n == Long.MAX_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_ADDITION, n, 1);
        }

        return n + 1;

    }

    /** Decrement a number, detecting overflows.
     * @param n number to decrement
     * @return n-1 if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static int decrementExact(final int n) throws MathRuntimeException {

        if (n == Integer.MIN_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_SUBTRACTION, n, 1);
        }

        return n - 1;

    }

    /** Decrement a number, detecting overflows.
     * @param n number to decrement
     * @return n-1 if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static long decrementExact(final long n) throws MathRuntimeException {

        if (n == Long.MIN_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_SUBTRACTION, n, 1);
        }

        return n - 1;

    }

    /** Add two numbers, detecting overflows.
     * @param a first number to add
     * @param b second number to add
     * @return a+b if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static int addExact(final int a, final int b) throws MathRuntimeException {

        // compute sum
        final int sum = a + b;

        // check for overflow
        if ((a ^ b) >= 0 && (sum ^ b) < 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_ADDITION, a, b);
        }

        return sum;

    }

    /** Add two numbers, detecting overflows.
     * @param a first number to add
     * @param b second number to add
     * @return a+b if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static long addExact(final long a, final long b) throws MathRuntimeException {

        // compute sum
        final long sum = a + b;

        // check for overflow
        if ((a ^ b) >= 0 && (sum ^ b) < 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_ADDITION, a, b);
        }

        return sum;

    }

    /** Subtract two numbers, detecting overflows.
     * @param a first number
     * @param b second number to subtract from a
     * @return a-b if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static int subtractExact(final int a, final int b) {

        // compute subtraction
        final int sub = a - b;

        // check for overflow
        if ((a ^ b) < 0 && (sub ^ b) >= 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_SUBTRACTION, a, b);
        }

        return sub;

    }

    /** Subtract two numbers, detecting overflows.
     * @param a first number
     * @param b second number to subtract from a
     * @return a-b if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static long subtractExact(final long a, final long b) {

        // compute subtraction
        final long sub = a - b;

        // check for overflow
        if ((a ^ b) < 0 && (sub ^ b) >= 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_SUBTRACTION, a, b);
        }

        return sub;

    }

    /** Multiply two numbers, detecting overflows.
     * @param a first number to multiply
     * @param b second number to multiply
     * @return a*b if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static int multiplyExact(final int a, final int b) {
        if (((b  >  0)  && (a > Integer.MAX_VALUE / b || a < Integer.MIN_VALUE / b)) ||
            ((b  < -1)  && (a > Integer.MIN_VALUE / b || a < Integer.MAX_VALUE / b)) ||
            ((b == -1)  && (a == Integer.MIN_VALUE))) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_MULTIPLICATION, a, b);
        }
        return a * b;
    }

    /** Multiply two numbers, detecting overflows.
     * @param a first number to multiply
     * @param b second number to multiply
     * @return a*b if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     * @since 1.3
     */
    public static long multiplyExact(final long a, final int b) {
        return multiplyExact(a, (long) b);
    }

    /** Multiply two numbers, detecting overflows.
     * @param a first number to multiply
     * @param b second number to multiply
     * @return a*b if no overflows occur
     * @exception MathRuntimeException if an overflow occurs
     */
    public static long multiplyExact(final long a, final long b) {
        if (((b  >  0l)  && (a > Long.MAX_VALUE / b || a < Long.MIN_VALUE / b)) ||
            ((b  < -1l)  && (a > Long.MIN_VALUE / b || a < Long.MAX_VALUE / b)) ||
            ((b == -1l)  && (a == Long.MIN_VALUE))) {
                throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_MULTIPLICATION, a, b);
            }
            return a * b;
    }

    /** Multiply two integers and give an exact result without overflow.
     * @param a first factor
     * @param b second factor
     * @return a * b exactly
     * @since 1.3
     */
    public static long multiplyFull(final int a, final int b) {
        return ((long) a) * ((long) b); // NOPMD - casts are intentional here
    }

    /** Multiply two long integers and give the 64 most significant bits of the result.
     * <p>
     * Beware that as Java primitive long are always considered to be signed, there are some
     * intermediate values {@code a} and {@code b} for which {@code a * b} exceeds {@code Long.MAX_VALUE}
     * but this method will still return 0l. This happens for example for {@code a = 2³¹} and
     * {@code b = 2³²} as {@code a * b = 2⁶³ = Long.MAX_VALUE + 1}, so it exceeds the max value
     * for a long, but still fits in 64 bits, so this method correctly returns 0l in this case,
     * but multiplication result would be considered negative (and in fact equal to {@code Long.MIN_VALUE}
     * </p>
     * @param a first factor
     * @param b second factor
     * @return a * b / 2<sup>64</sup>
     * @since 1.3
     */
    public static long multiplyHigh(final long a, final long b) {

        // all computations below are performed on unsigned numbers because we start
        // by using logical shifts (and not arithmetic shifts). We will therefore
        // need to take care of sign before returning
        // a negative long n between -2⁶³ and -1, interpreted as an unsigned long
        // corresponds to 2⁶⁴ + n (which is between 2⁶³ and 2⁶⁴-1)
        // so if this number is multiplied by p, what we really compute
        // is (2⁶⁴ + n) * p = 2⁶⁴ * p + n * p, therefore the part above 64 bits
        // will have an extra term p that we will need to remove
        final long tobeRemoved = ((a < 0) ? b : 0) + ((b < 0) ? a : 0);

        return unsignedMultiplyHigh(a, b) - tobeRemoved;

    }

    /** Multiply two long unsigned integers and give the 64 most significant bits of the unsigned result.
     * <p>
     * Beware that as Java primitive long are always considered to be signed, there are some
     * intermediate values {@code a} and {@code b} for which {@code a * b} exceeds {@code Long.MAX_VALUE}
     * but this method will still return 0l. This happens for example for {@code a = 2³¹} and
     * {@code b = 2³²} as {@code a * b = 2⁶³ = Long.MAX_VALUE + 1}, so it exceeds the max value
     * for a long, but still fits in 64 bits, so this method correctly returns 0l in this case,
     * but multiplication result would be considered negative (and in fact equal to {@code Long.MIN_VALUE}
     * </p>
     * @param a first factor
     * @param b second factor
     * @return a * b / 2<sup>64</sup>
     * @since 3.0
     */
    public static long unsignedMultiplyHigh(final long a, final long b) {

        // split numbers in two 32 bits parts
        final long aHigh  = a >>> 32;
        final long aLow   = a & 0xFFFFFFFFl;
        final long bHigh  = b >>> 32;
        final long bLow   = b & 0xFFFFFFFFl;

        // ab = aHigh * bHigh * 2⁶⁴ + (aHigh * bLow + aLow * bHigh) * 2³² + aLow * bLow
        final long hh     = aHigh * bHigh;
        final long hl1    = aHigh * bLow;
        final long hl2    = aLow  * bHigh;
        final long ll     = aLow  * bLow;

        // adds up everything in the above 64 bit part, taking care to avoid overflow
        final long hlHigh = (hl1 >>> 32) + (hl2 >>> 32);
        final long hlLow  = (hl1 & 0xFFFFFFFFl) + (hl2 & 0xFFFFFFFFl);
        final long carry  = (hlLow + (ll >>> 32)) >>> 32;

        return hh + hlHigh + carry;

    }

    /** Divide two integers, checking for overflow.
     * @param x dividend
     * @param y divisor
     * @return x / y
     * @exception MathRuntimeException if an overflow occurs
     * @since 3.0
     */
    public static int divideExact(final int x, final int y) {
        if (y == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }
        if (y == -1 && x == Integer.MIN_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_FRACTION, x, y);
        }
        return x / y;
    }

    /** Divide two long integers, checking for overflow.
     * @param x dividend
     * @param y divisor
     * @return x / y
     * @exception MathRuntimeException if an overflow occurs
     * @since 3.0
     */
    public static long divideExact(final long x, final long y) {
        if (y == 0l) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }
        if (y == -1l && x == Long.MIN_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_FRACTION, x, y);
        }
        return x / y;
    }

    /** Finds q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are opposite signs, but returns a different value when
     * they are same (i.e. q is positive).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0
     * @since 3.0
     */
    public static int ceilDiv(final int a, final int b) throws MathRuntimeException {

        if (b == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }

        final int m = a % b;
        if ((a ^ b) < 0 || m == 0) {
            // a and b have opposite signs, or division is exact
            return a / b;
        } else {
            // a and b have same signs and division is not exact
            return (a / b) + 1;
        }

    }

    /** Finds q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are opposite signs, but returns a different value when
     * they are same (i.e. q is positive).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0 or if a == {@code Integer.MIN_VALUE} and b = -1
     * @since 3.0
     */
    public static int ceilDivExact(final int a, final int b) throws MathRuntimeException {

        if (a == Integer.MIN_VALUE && b == -1) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_FRACTION, a, b);
        }

        return ceilDiv(a, b);

    }

    /** Finds q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are opposite signs, but returns a different value when
     * they are same (i.e. q is positive).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0
     * @since 3.0
     */
    public static long ceilDiv(final long a, final long b) throws MathRuntimeException {

        if (b == 0l) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }

        final long m = a % b;
        if ((a ^ b) < 0 || m == 0l) {
            // a and b have opposite signs, or division is exact
            return a / b;
        } else {
            // a and b have same signs and division is not exact
            return (a / b) + 1l;
        }

    }

    /** Finds q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are opposite signs, but returns a different value when
     * they are same (i.e. q is positive).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0 or if a == {@code Long.MIN_VALUE} and b = -1
     * @since 3.0
     */
    public static long ceilDivExact(final long a, final long b) throws MathRuntimeException {

        if (a == Long.MIN_VALUE && b == -1l) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_FRACTION, a, b);
        }

        return ceilDiv(a, b);

    }

    /** Finds q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are opposite signs, but returns a different value when
     * they are same (i.e. q is positive).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0
     * @since 3.0
     */
    public static long ceilDiv(final long a, final int b) throws MathRuntimeException {
        return ceilDiv(a, (long) b);
    }

    /** Finds r such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer modulo when
     * a and b are opposite signs, but returns a different value when
     * they are same (i.e. q is positive).
     *
     * @param a dividend
     * @param b divisor
     * @return r such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0
     * @since 3.0
     */
    public static int ceilMod(final int a, final int b) throws MathRuntimeException {

        if (b == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }

        final int m = a % b;
        if ((a ^ b) < 0 || m == 0) {
            // a and b have opposite signs, or division is exact
            return m;
        } else {
            // a and b have same signs and division is not exact
            return m - b;
        }

    }

    /** Finds r such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer modulo when
     * a and b are opposite signs, but returns a different value when
     * they are same (i.e. q is positive).
     *
     * @param a dividend
     * @param b divisor
     * @return r such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0
     * @since 3.0
     */
    public static int ceilMod(final long a, final int b) throws MathRuntimeException {
        return (int) ceilMod(a, (long) b);
    }

    /** Finds r such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer modulo when
     * a and b are opposite signs, but returns a different value when
     * they are same (i.e. q is positive).
     *
     * @param a dividend
     * @param b divisor
     * @return r such that {@code a = q b + r} with {@code b < r <= 0} if {@code b > 0} and {@code 0 <= r < b} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0
     * @since 3.0
     */
    public static long ceilMod(final long a, final long b) throws MathRuntimeException {

        if (b == 0l) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }

        final long m = a % b;
        if ((a ^ b) < 0l || m == 0l) {
            // a and b have opposite signs, or division is exact
            return m;
        } else {
            // a and b have same signs and division is not exact
            return m - b;
        }

    }

    /** Finds q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are same signs, but returns a different value when
     * they are opposite (i.e. q is negative).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0
     * @see #floorMod(int, int)
     */
    public static int floorDiv(final int a, final int b) throws MathRuntimeException {

        if (b == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }

        final int m = a % b;
        if ((a ^ b) >= 0 || m == 0) {
            // a and b have same sign, or division is exact
            return a / b;
        } else {
            // a and b have opposite signs and division is not exact
            return (a / b) - 1;
        }

    }

    /** Finds q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code  b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are same signs, but returns a different value when
     * they are opposite (i.e. q is negative).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code  b < 0}
     * @exception MathRuntimeException if b == 0 or if a == {@code Integer.MIN_VALUE} and b = -1
     * @see #floorMod(int, int)
     * @since 3.0
     */
    public static int floorDivExact(final int a, final int b) throws MathRuntimeException {

        if (a == Integer.MIN_VALUE && b == -1) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_FRACTION, a, b);
        }

        return floorDiv(a, b);

    }

    /** Finds q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are same signs, but returns a different value when
     * they are opposite (i.e. q is negative).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
     * @exception MathRuntimeException if b == 0
     * @see #floorMod(long, int)
     * @since 1.3
     */
    public static long floorDiv(final long a, final int b) throws MathRuntimeException {
        return floorDiv(a, (long) b);
    }

    /** Finds q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are same signs, but returns a different value when
     * they are opposite (i.e. q is negative).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
     * @exception MathRuntimeException if b == 0
     * @see #floorMod(long, long)
     */
    public static long floorDiv(final long a, final long b) throws MathRuntimeException {

        if (b == 0l) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }

        final long m = a % b;
        if ((a ^ b) >= 0l || m == 0l) {
            // a and b have same sign, or division is exact
            return a / b;
        } else {
            // a and b have opposite signs and division is not exact
            return (a / b) - 1l;
        }

    }

    /** Finds q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
     * <p>
     * This methods returns the same value as integer division when
     * a and b are same signs, but returns a different value when
     * they are opposite (i.e. q is negative).
     *
     * @param a dividend
     * @param b divisor
     * @return q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
     * @exception MathRuntimeException if b == 0 or if a == {@code Long.MIN_VALUE} and b = -1
     * @see #floorMod(long, long)
     * @since 3.0
     */
    public static long floorDivExact(final long a, final long b) throws MathRuntimeException {

        if (a == Long.MIN_VALUE && b == -1l) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_FRACTION, a, b);
        }

        return floorDiv(a, b);

    }

    /** Finds r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
     * <p>
     * This methods returns the same value as integer modulo when
     * a and b are same signs, but returns a different value when
     * they are opposite (i.e. q is negative).
     * </p>
     * @param a dividend
     * @param b divisor
     * @return r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
     * @exception MathRuntimeException if b == 0
     * @see #floorDiv(int, int)
     */
    public static int floorMod(final int a, final int b) throws MathRuntimeException {

        if (b == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }

        final int m = a % b;
        if ((a ^ b) >= 0 || m == 0) {
            // a and b have same sign, or division is exact
            return m;
        } else {
            // a and b have opposite signs and division is not exact
            return b + m;
        }

    }

    /** Finds r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
     * <p>
     * This methods returns the same value as integer modulo when
     * a and b are same signs, but returns a different value when
     * they are opposite (i.e. q is negative).
     * </p>
     * @param a dividend
     * @param b divisor
     * @return r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
     * @exception MathRuntimeException if b == 0
     * @see #floorDiv(long, int)
     * @since 1.3
     */
    public static int floorMod(final long a, final int b) {
        return (int) floorMod(a, (long) b);
    }

    /** Finds r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
     * <p>
     * This methods returns the same value as integer modulo when
     * a and b are same signs, but returns a different value when
     * they are opposite (i.e. q is negative).
     * </p>
     * @param a dividend
     * @param b divisor
     * @return r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
     * @exception MathRuntimeException if b == 0
     * @see #floorDiv(long, long)
     */
    public static long floorMod(final long a, final long b) {

        if (b == 0l) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }

        final long m = a % b;
        if ((a ^ b) >= 0l || m == 0l) {
            // a and b have same sign, or division is exact
            return m;
        } else {
            // a and b have opposite signs and division is not exact
            return b + m;
        }

    }

    /**
     * Returns the first argument with the sign of the second argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param magnitude the value to return
     * @param sign the sign for the returned value
     * @return the magnitude with the same sign as the {@code sign} argument
     */
    public static double copySign(double magnitude, double sign){
        // The highest order bit is going to be zero if the
        // highest order bit of m and s is the same and one otherwise.
        // So (m^s) will be positive if both m and s have the same sign
        // and negative otherwise.
        final long m = Double.doubleToRawLongBits(magnitude); // don't care about NaN
        final long s = Double.doubleToRawLongBits(sign);
        if ((m ^ s) >= 0) {
            return magnitude;
        }
        return -magnitude; // flip sign
    }

    /**
     * Returns the first argument with the sign of the second argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param magnitude the value to return
     * @param sign the sign for the returned value
     * @return the magnitude with the same sign as the {@code sign} argument
     */
    public static float copySign(float magnitude, float sign){
        // The highest order bit is going to be zero if the
        // highest order bit of m and s is the same and one otherwise.
        // So (m^s) will be positive if both m and s have the same sign
        // and negative otherwise.
        final int m = Float.floatToRawIntBits(magnitude);
        final int s = Float.floatToRawIntBits(sign);
        if ((m ^ s) >= 0) {
            return magnitude;
        }
        return -magnitude; // flip sign
    }

    /**
     * Return the exponent of a double number, removing the bias.
     * <p>
     * For double numbers of the form 2<sup>x</sup>, the unbiased
     * exponent is exactly x.
     * </p>
     * @param d number from which exponent is requested
     * @return exponent for d in IEEE754 representation, without bias
     */
    public static int getExponent(final double d) {
        // NaN and Infinite will return 1024 anyhow so can use raw bits
        return (int) ((Double.doubleToRawLongBits(d) >>> 52) & 0x7ff) - 1023;
    }

    /**
     * Return the exponent of a float number, removing the bias.
     * <p>
     * For float numbers of the form 2<sup>x</sup>, the unbiased
     * exponent is exactly x.
     * </p>
     * @param f number from which exponent is requested
     * @return exponent for d in IEEE754 representation, without bias
     */
    public static int getExponent(final float f) {
        // NaN and Infinite will return the same exponent anyhow so can use raw bits
        return ((Float.floatToRawIntBits(f) >>> 23) & 0xff) - 127;
    }

    /** Compute Fused-multiply-add operation a * b + c.
     * <p>
     * This method was introduced in the regular {@code Math} and {@code StrictMath}
     * methods with Java 9, and then added to Hipparchus for consistency. However,
     * a more general method was available in Hipparchus that also allow to repeat
     * this computation across several terms: {@link MathArrays#linearCombination(double[], double[])}.
     * The linear combination method should probably be preferred in most cases.
     * </p>
     * @param a first factor
     * @param b second factor
     * @param c additive term
     * @return a * b + c, using extended precision in the multiplication
     * @see MathArrays#linearCombination(double[], double[])
     * @see MathArrays#linearCombination(double, double, double, double)
     * @see MathArrays#linearCombination(double, double, double, double, double, double)
     * @see MathArrays#linearCombination(double, double, double, double, double, double, double, double)
     * @since 1.3
     */
    public static double fma(final double a, final double b, final double c) {
        return MathArrays.linearCombination(a, b, 1.0, c);
    }

    /** Compute Fused-multiply-add operation a * b + c.
     * <p>
     * This method was introduced in the regular {@code Math} and {@code StrictMath}
     * methods with Java 9, and then added to Hipparchus for consistency. However,
     * a more general method was available in Hipparchus that also allow to repeat
     * this computation across several terms: {@link MathArrays#linearCombination(double[], double[])}.
     * The linear combination method should probably be preferred in most cases.
     * </p>
     * @param a first factor
     * @param b second factor
     * @param c additive term
     * @return a * b + c, using extended precision in the multiplication
     * @see MathArrays#linearCombination(double[], double[])
     * @see MathArrays#linearCombination(double, double, double, double)
     * @see MathArrays#linearCombination(double, double, double, double, double, double)
     * @see MathArrays#linearCombination(double, double, double, double, double, double, double, double)
     */
    public static float fma(final float a, final float b, final float c) {
        return (float) MathArrays.linearCombination(a, b, 1.0, c);
    }

    /** Compute the square root of a number.
     * @param a number on which evaluation is done
     * @param <T> the type of the field element
     * @return square root of a
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T sqrt(final T a) {
        return a.sqrt();
    }

    /** Compute the hyperbolic cosine of a number.
     * @param x number on which evaluation is done
     * @param <T> the type of the field element
     * @return hyperbolic cosine of x
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T cosh(final T x) {
        return x.cosh();
    }

    /** Compute the hyperbolic sine of a number.
     * @param x number on which evaluation is done
     * @param <T> the type of the field element
     * @return hyperbolic sine of x
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T sinh(final T x) {
        return x.sinh();
    }

    /** Compute the hyperbolic tangent of a number.
     * @param x number on which evaluation is done
     * @param <T> the type of the field element
     * @return hyperbolic tangent of x
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T tanh(final T x) {
        return x.tanh();
    }

    /** Compute the inverse hyperbolic cosine of a number.
     * @param a number on which evaluation is done
     * @param <T> the type of the field element
     * @return inverse hyperbolic cosine of a
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T acosh(final T a) {
        return a.acosh();
    }

    /** Compute the inverse hyperbolic sine of a number.
     * @param a number on which evaluation is done
     * @param <T> the type of the field element
     * @return inverse hyperbolic sine of a
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T asinh(final T a) {
        return a.asinh();
    }

    /** Compute the inverse hyperbolic tangent of a number.
     * @param a number on which evaluation is done
     * @param <T> the type of the field element
     * @return inverse hyperbolic tangent of a
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T atanh(final T a) {
        return a.atanh();
    }

    /** Compute the sign of a number.
     * The sign is -1 for negative numbers, +1 for positive numbers and 0 otherwise,
     * for Complex number, it is extended on the unit circle (equivalent to z/|z|,
     * with special handling for 0 and NaN)
     * @param a number on which evaluation is done
     * @param <T> the type of the field element
     * @return -1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a
     * @since 2.0
     */
    public static <T extends CalculusFieldElement<T>> T sign(final T a) {
        return a.sign();
    }

    /**
     * Exponential function.
     *
     * Computes exp(x), function result is nearly rounded.   It will be correctly
     * rounded to the theoretical value for 99.9% of input values, otherwise it will
     * have a 1 ULP error.
     *
     * Method:
     *    Lookup intVal = exp(int(x))
     *    Lookup fracVal = exp(int(x-int(x) / 1024.0) * 1024.0 );
     *    Compute z as the exponential of the remaining bits by a polynomial minus one
     *    exp(x) = intVal * fracVal * (1 + z)
     *
     * Accuracy:
     *    Calculation is done with 63 bits of precision, so result should be correctly
     *    rounded for 99.9% of input values, with less than 1 ULP error otherwise.
     *
     * @param x   a double
     * @param <T> the type of the field element
     * @return double e<sup>x</sup>
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T exp(final T x) {
        return x.exp();
    }

    /** Compute exp(x) - 1
     * @param x number to compute shifted exponential
     * @param <T> the type of the field element
     * @return exp(x) - 1
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T expm1(final T x) {
        return x.expm1();
    }

    /**
     * Natural logarithm.
     *
     * @param x   a double
     * @param <T> the type of the field element
     * @return log(x)
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T log(final T x) {
        return x.log();
    }

    /**
     * Computes log(1 + x).
     *
     * @param x Number.
     * @param <T> the type of the field element
     * @return {@code log(1 + x)}.
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T log1p(final T x) {
        return x.log1p();
    }

    /** Compute the base 10 logarithm.
     * @param x a number
     * @param <T> the type of the field element
     * @return log10(x)
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T log10(final T x) {
        return x.log10();
    }

    /**
     * Power function.  Compute x<sup>y</sup>.
     *
     * @param x   a double
     * @param y   a double
     * @param <T> the type of the field element
     * @return x<sup>y</sup>
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T pow(final T x, final T y) {
        return x.pow(y);
    }

    /**
     * Power function.  Compute x<sup>y</sup>.
     *
     * @param x   a double
     * @param y   a double
     * @param <T> the type of the field element
     * @return x<sup>y</sup>
     * @since 1.7
     */
    public static <T extends CalculusFieldElement<T>> T pow(final T x, final double y) {
        return x.pow(y);
    }

    /**
     * Raise a double to an int power.
     *
     * @param d Number to raise.
     * @param e Exponent.
     * @param <T> the type of the field element
     * @return d<sup>e</sup>
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T pow(T d, int e) {
        return d.pow(e);
    }

    /**
     * Sine function.
     *
     * @param x Argument.
     * @param <T> the type of the field element
     * @return sin(x)
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T sin(final T x) {
        return x.sin();
    }

    /**
     * Cosine function.
     *
     * @param x Argument.
     * @param <T> the type of the field element
     * @return cos(x)
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T cos(final T x) {
        return x.cos();
    }

    /**
     * Tangent function.
     *
     * @param x Argument.
     * @param <T> the type of the field element
     * @return tan(x)
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T tan(final T x) {
        return x.tan();
    }

    /**
     * Arctangent function
     *  @param x a number
     * @param <T> the type of the field element
     *  @return atan(x)
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T atan(final T x) {
        return x.atan();
    }

    /**
     * Two arguments arctangent function
     * @param y ordinate
     * @param x abscissa
     * @param <T> the type of the field element
     * @return phase angle of point (x,y) between {@code -PI} and {@code PI}
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T atan2(final T y, final T x) {
        return y.atan2(x);
    }

    /** Compute the arc sine of a number.
     * @param x number on which evaluation is done
     * @param <T> the type of the field element
     * @return arc sine of x
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T asin(final T x) {
        return x.asin();
    }

    /** Compute the arc cosine of a number.
     * @param x number on which evaluation is done
     * @param <T> the type of the field element
     * @return arc cosine of x
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T acos(final T x) {
        return x.acos();
    }

    /** Compute the cubic root of a number.
     * @param x number on which evaluation is done
     * @param <T> the type of the field element
     * @return cubic root of x
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T cbrt(final T x) {
        return x.cbrt();
    }

    /**
     * Norm.
     * @param x number from which norm is requested
     * @param <T> the type of the field element
     * @return norm(x)
     * @since 2.0
     */
    public static <T extends CalculusFieldElement<T>> double norm(final T x) {
        return x.norm();
    }

    /**
     * Absolute value.
     * @param x number from which absolute value is requested
     * @param <T> the type of the field element
     * @return abs(x)
     * @since 2.0
     */
    public static <T extends CalculusFieldElement<T>> T abs(final T x) {
        return x.abs();
    }

    /**
     *  Convert degrees to radians, with error of less than 0.5 ULP
     *  @param x angle in degrees
     *  @param <T> the type of the field element
     *  @return x converted into radians
     */
    public static <T extends CalculusFieldElement<T>> T toRadians(T x) {
        return x.toRadians();
    }

    /**
     *  Convert radians to degrees, with error of less than 0.5 ULP
     *  @param x angle in radians
     *  @param <T> the type of the field element
     *  @return x converted into degrees
     */
    public static <T extends CalculusFieldElement<T>> T toDegrees(T x) {
        return x.toDegrees();
    }

    /**
     * Multiply a double number by a power of 2.
     * @param d number to multiply
     * @param n power of 2
     * @param <T> the type of the field element
     * @return d &times; 2<sup>n</sup>
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T scalb(final T d, final int n) {
        return d.scalb(n);
    }

    /**
     * Compute least significant bit (Unit in Last Position) for a number.
     * @param x number from which ulp is requested
     * @param <T> the type of the field element
     * @return ulp(x)
     * @since 2.0
     */
    public static <T extends CalculusFieldElement<T>> T ulp(final T x) {
        if (Double.isInfinite(x.getReal())) {
            return x.newInstance(Double.POSITIVE_INFINITY);
        }
        return x.ulp();
    }

    /** Get the largest whole number smaller than x.
     * @param x number from which floor is requested
     * @param <T> the type of the field element
     * @return a double number f such that f is an integer f &lt;= x &lt; f + 1.0
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T floor(final T x) {
        return x.floor();
    }

    /** Get the smallest whole number larger than x.
     * @param x number from which ceil is requested
     * @param <T> the type of the field element
     * @return a double number c such that c is an integer c - 1.0 &lt; x &lt;= c
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T ceil(final T x) {
        return x.ceil();
    }

    /** Get the whole number that is the nearest to x, or the even one if x is exactly half way between two integers.
     * @param x number from which nearest whole number is requested
     * @param <T> the type of the field element
     * @return a double number r such that r is an integer r - 0.5 &lt;= x &lt;= r + 0.5
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T rint(final T x) {
        return x.rint();
    }

    /** Get the closest long to x.
     * @param x number from which closest long is requested
     * @param <T> the type of the field element
     * @return closest long to x
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> long round(final T x) {
        return x.round();
    }

    /** Compute the minimum of two values
     * @param a first value
     * @param b second value
     * @param <T> the type of the field element
     * @return a if a is lesser or equal to b, b otherwise
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T min(final T a, final T b) {
        final double aR = a.getReal();
        final double bR = b.getReal();
        if (aR < bR) {
            return a;
        } else if (bR < aR) {
            return b;
        } else {
            // either the numbers are equal, or one of them is a NaN
            return Double.isNaN(aR) ? a : b;
        }
    }

    /** Compute the minimum of two values
     * @param a first value
     * @param b second value
     * @param <T> the type of the field element
     * @return a if a is lesser or equal to b, b otherwise
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T min(final T a, final double b) {
        final double aR = a.getReal();
        if (aR < b) {
            return a;
        } else if (b < aR) {
            return a.getField().getZero().add(b);
        } else {
            // either the numbers are equal, or one of them is a NaN
            return Double.isNaN(aR) ? a : a.getField().getZero().add(b);
        }
    }

    /** Compute the maximum of two values
     * @param a first value
     * @param b second value
     * @param <T> the type of the field element
     * @return b if a is lesser or equal to b, a otherwise
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T max(final T a, final T b) {
        final double aR = a.getReal();
        final double bR = b.getReal();
        if (aR < bR) {
            return b;
        } else if (bR < aR) {
            return a;
        } else {
            // either the numbers are equal, or one of them is a NaN
            return Double.isNaN(aR) ? a : b;
        }
    }

    /** Compute the maximum of two values
     * @param a first value
     * @param b second value
     * @param <T> the type of the field element
     * @return b if a is lesser or equal to b, a otherwise
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T max(final T a, final double b) {
        final double aR = a.getReal();
        if (aR < b) {
            return a.getField().getZero().add(b);
        } else if (b < aR) {
            return a;
        } else {
            // either the numbers are equal, or one of them is a NaN
            return Double.isNaN(aR) ? a : a.getField().getZero().add(b);
        }
    }

    /**
     * Returns the hypotenuse of a triangle with sides {@code x} and {@code y}
     * - sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)<br>
     * avoiding intermediate overflow or underflow.
     *
     * <ul>
     * <li> If either argument is infinite, then the result is positive infinity.</li>
     * <li> else, if either argument is NaN then the result is NaN.</li>
     * </ul>
     *
     * @param x a value
     * @param y a value
     * @param <T> the type of the field element
     * @return sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T hypot(final T x, final T y) {
        return x.hypot(y);
    }

    /**
     * Computes the remainder as prescribed by the IEEE 754 standard.
     * <p>
     * The remainder value is mathematically equal to {@code x - y*n}
     * where {@code n} is the mathematical integer closest to the exact mathematical value
     * of the quotient {@code x/y}.
     * If two mathematical integers are equally close to {@code x/y} then
     * {@code n} is the integer that is even.
     * </p>
     * <ul>
     * <li>If either operand is NaN, the result is NaN.</li>
     * <li>If the result is not NaN, the sign of the result equals the sign of the dividend.</li>
     * <li>If the dividend is an infinity, or the divisor is a zero, or both, the result is NaN.</li>
     * <li>If the dividend is finite and the divisor is an infinity, the result equals the dividend.</li>
     * <li>If the dividend is a zero and the divisor is finite, the result equals the dividend.</li>
     * </ul>
     * @param dividend the number to be divided
     * @param divisor the number by which to divide
     * @param <T> the type of the field element
     * @return the remainder, rounded
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T IEEEremainder(final T dividend, final double divisor) {
        return dividend.remainder(divisor);
    }

    /**
     * Computes the remainder as prescribed by the IEEE 754 standard.
     * <p>
     * The remainder value is mathematically equal to {@code x - y*n}
     * where {@code n} is the mathematical integer closest to the exact mathematical value
     * of the quotient {@code x/y}.
     * If two mathematical integers are equally close to {@code x/y} then
     * {@code n} is the integer that is even.
     * </p>
     * <ul>
     * <li>If either operand is NaN, the result is NaN.</li>
     * <li>If the result is not NaN, the sign of the result equals the sign of the dividend.</li>
     * <li>If the dividend is an infinity, or the divisor is a zero, or both, the result is NaN.</li>
     * <li>If the dividend is finite and the divisor is an infinity, the result equals the dividend.</li>
     * <li>If the dividend is a zero and the divisor is finite, the result equals the dividend.</li>
     * </ul>
     * @param dividend the number to be divided
     * @param divisor the number by which to divide
     * @param <T> the type of the field element
     * @return the remainder, rounded
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T IEEEremainder(final T dividend, final T divisor) {
        return dividend.remainder(divisor);
    }

    /**
     * Returns the first argument with the sign of the second argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param magnitude the value to return
     * @param sign the sign for the returned value
     * @param <T> the type of the field element
     * @return the magnitude with the same sign as the {@code sign} argument
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T copySign(T magnitude, T sign) {
        return magnitude.copySign(sign);
    }

    /**
     * Returns the first argument with the sign of the second argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param magnitude the value to return
     * @param sign the sign for the returned value
     * @param <T> the type of the field element
     * @return the magnitude with the same sign as the {@code sign} argument
     * @since 1.3
     */
    public static <T extends CalculusFieldElement<T>> T copySign(T magnitude, double sign) {
        return magnitude.copySign(sign);
    }

//    /**
//     * Print out contents of arrays, and check the length.
//     * <p>used to generate the preset arrays originally.</p>
//     * @param a unused
//     */
//    public static void main(String[] a) {
//        FastMathCalc.printarray(System.out, "EXP_INT_TABLE_A", EXP_INT_TABLE_LEN, ExpIntTable.EXP_INT_TABLE_A);
//        FastMathCalc.printarray(System.out, "EXP_INT_TABLE_B", EXP_INT_TABLE_LEN, ExpIntTable.EXP_INT_TABLE_B);
//        FastMathCalc.printarray(System.out, "EXP_FRAC_TABLE_A", EXP_FRAC_TABLE_LEN, ExpFracTable.EXP_FRAC_TABLE_A);
//        FastMathCalc.printarray(System.out, "EXP_FRAC_TABLE_B", EXP_FRAC_TABLE_LEN, ExpFracTable.EXP_FRAC_TABLE_B);
//        FastMathCalc.printarray(System.out, "LN_MANT",LN_MANT_LEN, lnMant.LN_MANT);
//        FastMathCalc.printarray(System.out, "SINE_TABLE_A", SINE_TABLE_LEN, SINE_TABLE_A);
//        FastMathCalc.printarray(System.out, "SINE_TABLE_B", SINE_TABLE_LEN, SINE_TABLE_B);
//        FastMathCalc.printarray(System.out, "COSINE_TABLE_A", SINE_TABLE_LEN, COSINE_TABLE_A);
//        FastMathCalc.printarray(System.out, "COSINE_TABLE_B", SINE_TABLE_LEN, COSINE_TABLE_B);
//        FastMathCalc.printarray(System.out, "TANGENT_TABLE_A", SINE_TABLE_LEN, TANGENT_TABLE_A);
//        FastMathCalc.printarray(System.out, "TANGENT_TABLE_B", SINE_TABLE_LEN, TANGENT_TABLE_B);
//    }

    /** Enclose large data table in nested static class so it's only loaded on first access. */
    private static class ExpIntTable {
        /** Exponential evaluated at integer values,
         * exp(x) =  expIntTableA[x + EXP_INT_TABLE_MAX_INDEX] + expIntTableB[x+EXP_INT_TABLE_MAX_INDEX].
         */
        private static final double[] EXP_INT_TABLE_A;
        /** Exponential evaluated at integer values,
         * exp(x) =  expIntTableA[x + EXP_INT_TABLE_MAX_INDEX] + expIntTableB[x+EXP_INT_TABLE_MAX_INDEX]
         */
        private static final double[] EXP_INT_TABLE_B;

        static {
            if (RECOMPUTE_TABLES_AT_RUNTIME) {
                EXP_INT_TABLE_A = new double[EXP_INT_TABLE_LEN];
                EXP_INT_TABLE_B = new double[EXP_INT_TABLE_LEN];

                final double[] tmp = new double[2];
                final double[] recip = new double[2];

                // Populate expIntTable
                for (int i = 0; i < EXP_INT_TABLE_MAX_INDEX; i++) {
                    FastMathCalc.expint(i, tmp);
                    EXP_INT_TABLE_A[i + EXP_INT_TABLE_MAX_INDEX] = tmp[0];
                    EXP_INT_TABLE_B[i + EXP_INT_TABLE_MAX_INDEX] = tmp[1];

                    if (i != 0) {
                        // Negative integer powers
                        FastMathCalc.splitReciprocal(tmp, recip);
                        EXP_INT_TABLE_A[EXP_INT_TABLE_MAX_INDEX - i] = recip[0];
                        EXP_INT_TABLE_B[EXP_INT_TABLE_MAX_INDEX - i] = recip[1];
                    }
                }
            } else {
                EXP_INT_TABLE_A = FastMathLiteralArrays.loadExpIntA();
                EXP_INT_TABLE_B = FastMathLiteralArrays.loadExpIntB();
            }
        }
    }

    /** Enclose large data table in nested static class so it's only loaded on first access. */
    private static class ExpFracTable {
        /** Exponential over the range of 0 - 1 in increments of 2^-10
         * exp(x/1024) =  expFracTableA[x] + expFracTableB[x].
         * 1024 = 2^10
         */
        private static final double[] EXP_FRAC_TABLE_A;
        /** Exponential over the range of 0 - 1 in increments of 2^-10
         * exp(x/1024) =  expFracTableA[x] + expFracTableB[x].
         */
        private static final double[] EXP_FRAC_TABLE_B;

        static {
            if (RECOMPUTE_TABLES_AT_RUNTIME) {
                EXP_FRAC_TABLE_A = new double[EXP_FRAC_TABLE_LEN];
                EXP_FRAC_TABLE_B = new double[EXP_FRAC_TABLE_LEN];

                final double[] tmp = new double[2];

                // Populate expFracTable
                final double factor = 1d / (EXP_FRAC_TABLE_LEN - 1);
                for (int i = 0; i < EXP_FRAC_TABLE_A.length; i++) {
                    FastMathCalc.slowexp(i * factor, tmp);
                    EXP_FRAC_TABLE_A[i] = tmp[0];
                    EXP_FRAC_TABLE_B[i] = tmp[1];
                }
            } else {
                EXP_FRAC_TABLE_A = FastMathLiteralArrays.loadExpFracA();
                EXP_FRAC_TABLE_B = FastMathLiteralArrays.loadExpFracB();
            }
        }
    }
}
