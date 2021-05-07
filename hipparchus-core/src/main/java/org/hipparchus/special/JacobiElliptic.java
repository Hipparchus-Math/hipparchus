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
package org.hipparchus.special;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.SinCos;

/** Computation of Jacobi elliptic functions.
 * The Jacobi elliptic functions are related to elliptic integrals.
 * We use hare the notations from <a
 * href="https://en.wikipedia.org/wiki/Abramowitz_and_Stegun">Abramowitz and
 * Stegun</a> (Ch. 16) with a parameter {@code m}. The notations from
 * <a href="https://dlmf.nist.gov/22">Digital Library of Mathematical Functions (Ch. 22)</a>
 * are different as they use modulus {@code k} instead of parameter {@code m},
 * with {@code k² = m}.
 * @since 2.0
 */
public class JacobiElliptic {

    /** Threshold near 0 for using specialized algorithm. */
    private static final double NEAR_ZERO = 1.0e-9;

    /** Threshold near 1 for using specialized algorithm. */
    private static final double NEAR_ONE = 1.0 - NEAR_ZERO;

    /** Algorithm to use for evaluating the functions. */
    private final Algorithm algorithm;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function
     */
    public JacobiElliptic(final double m) {
        algorithm = selectAlgorithm(m);
    }

    /** Evaluate the three principal Jacobi elliptic functions with pole at point n in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three principal Jacobi
     * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
     */
    public CopolarN valuesN(final double u) {
        return algorithm.values(u);
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point s in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
     */
    public CopolarS valuesS(final double u) {
        return new CopolarS(algorithm.values(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point c in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
     */
    public CopolarC valuesC(final double u) {
        return new CopolarC(algorithm.values(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point d in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
     */
    public CopolarD valuesD(final double u) {
        return new CopolarD(algorithm.values(u));
    }

    /** Select an algorithm for computing Jacobi elliptic functions.
     * @param m parameter of the Jacobi elliptic function
     * @return selected algorithm
     */
    private static Algorithm selectAlgorithm(final double m) {
        if (m < 0) {
            return new Negative(m);
        } else if (m > 1) {
            return new Big(m);
        } else if (m < NEAR_ZERO) {
            return new NearZero(m);
        } else if (m > NEAR_ONE) {
            return new NearOne(m);
        } else {
            return new Bounded(m);
        }
    }

    /** Interface for computing the principal Jacobi functions. */
    private interface Algorithm {
        /** Evaluate the three principal Jacobi elliptic functions.
         * @param u argument of the functions
         * @return copolar trio containing the three principal Jacobi
         * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
         */
        CopolarN values(double u);
    }

    /** Algorithm for computing the principal Jacobi functions for parameter m in [0; 1]. */
    private static class Bounded implements Algorithm {

        /** Max number of iterations of the AGM scale. */
        private static final int N_MAX = 16;

        /** Parameter of the Jacobi elliptic function. */
        private final double m;

        /** Initial value for arithmetic-geometric mean. */
        private final double b0;

        /** Initial value for arithmetic-geometric mean. */
        private final double c0;

        /** Simple constructor.
         * @param m parameter of the Jacobi elliptic function
         */
        public Bounded(final double m) {
            this.m  = m;
            this.b0 = FastMath.sqrt(1.0 - m);
            this.c0 = FastMath.sqrt(m);
        }

        /** {@inheritDoc}
         * <p>
         * The algorithm for evaluating the functions is based on arithmetic-geometric
         * mean. It is given in Abramowitz and Stegun, sections 16.4 and 17.6.
         * </p>
         * @param u argument of the functions, must be close far from boundaries here
         * @return copolar trio containing the three principal Jacobi
         * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
         */
        public CopolarN values(double u) {

            // initialize scale
            final double[] a = new double[N_MAX];
            final double[] c = new double[N_MAX];
            a[0]      = 1.0;
            double bi = b0;
            c[0]      = c0;

            // iterate down
            double phi = u;
            for (int i = 1; i < N_MAX; ++i) {

                // 2ⁿ u
                phi += phi;
                c[i] = 0.5 * (a[i - 1] - bi);
                
                // arithmetic mean
                a[i] = 0.5 * (a[i - 1] + bi);

                // geometric mean
                bi = FastMath.sqrt(a[i - 1] * bi);

                // convergence (by the inequality of arithmetic and geometric means, this is non-negative)
               if (c[i] <= FastMath.ulp(a[i])) {
                    // convergence has been reached

                    // iterate up
                    phi *= a[i];
                    for (int j = i; j > 0; --j) {
                        // equation 16.4.3 in Abramowitz and Stegun
                        phi = 0.5 * (phi + FastMath.asin(c[j] * FastMath.sin(phi) / a[j]));
                    }
                    // using 16.1.5 rather than 16.4.4 to avoid computing another cosine
                    final SinCos scPhi0 = FastMath.sinCos(phi);
                    return new CopolarN(scPhi0.sin(), scPhi0.cos(),
                                        FastMath.sqrt(1 - m * scPhi0.sin() * scPhi0.sin()));

                }

            }

            // we were not able to compute the value
            throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

        }

    }
    
    /** Algorithm for computing the principal Jacobi functions for parameters slightly above zero.
     * <p>
     * The algorithm for evaluating the functions is based on approximation
     * in terms of circular functions. It is given in Abramowitz and Stegun,
     * sections 16.13.
     * </p>
     */
    private static class NearZero implements Algorithm {

        /** Parameter of the Jacobi elliptic function. */
        private final double m;

        /** Simple constructor.
         * @param m parameter of the Jacobi elliptic function (must be negative here)
         */
        public NearZero(final double m) {
            this.m = m;
        }

        /** {@inheritDoc} */
        @Override
        public CopolarN values(final double u) {
            final SinCos sc     = FastMath.sinCos(u);
            final double factor = 0.25 * m * (u - sc.sin() * sc.cos());
            return new CopolarN(sc.sin() - factor * sc.cos(),       // 16.13.1
                                sc.cos() + factor * sc.sin(),       // 16.13.2
                                1 - 0.5 * m * sc.sin() * sc.sin()); // 16.13.3
        }

    }
    
    /** Algorithm for computing the principal Jacobi functions for parameters slightly below one.
     * <p>
     * The algorithm for evaluating the functions is based on approximation
     * in terms of hyperbolic functions. It is given in Abramowitz and Stegun,
     * sections 16.15.
     * </p>
     */
    private static class NearOne implements Algorithm {

        /** Complementary parameter of the Jacobi elliptic function. */
        private final double m1;

        /** Simple constructor.
         * @param m parameter of the Jacobi elliptic function (must be negative here)
         */
        public NearOne(final double m) {
            this.m1 = 1.0 - m;
        }

        /** {@inheritDoc} */
        @Override
        public CopolarN values(final double u) {
            final double s      = FastMath.sinh(u);
            final double c      = FastMath.cosh(u);
            final double sech   =  1.0 / c;
            final double t      = s * sech;
            final double factor = 0.25 * m1 * (s * c  - u) * sech;
            return new CopolarN(t + factor * sech,  // 16.15.1
                                sech - factor * t,  // 16.15.2
                                sech + factor * t); // 16.15.3
        }

    }
    
    /** Algorithm for computing the principal Jacobi functions for negative parameter m.
     * <p>
     * The rules for negative parameter change are given in Abramowitz and Stegun, section 16.10.
     * </p>
     */
    private static class Negative implements Algorithm {

        /** Algorithm to use for the positive parameter. */
        private final Algorithm algorithm;

        /** Input scaling factor. */
        private final double inputScale;

        /** output scaling factor. */
        private final double outputScale;

        /** Simple constructor.
         * @param m parameter of the Jacobi elliptic function (must be negative here)
         */
        public Negative(final double m) {
            final double omM = 1.0 - m;
            algorithm        = selectAlgorithm(-m / omM);
            inputScale       = FastMath.sqrt(omM);
            outputScale      = 1.0 / inputScale;
        }

        /** {@inheritDoc} */
        @Override
        public CopolarN values(final double u) {
            final CopolarD trioD = new CopolarD(algorithm.values(u * inputScale));
            return new CopolarN(outputScale * trioD.sd(), trioD.cd(), trioD.nd());
        }

    }
    
    /** Algorithm for computing the principal Jacobi functions for parameter m greater than 1.
     * <p>
     * The rules for reciprocal parameter change are given in Abramowitz and Stegun, section 16.11.
     * </p>
     */
    private static class Big implements Algorithm {

        /** Algorithm to use for the positive parameter. */
        private final Algorithm algorithm;

        /** Input scaling factor. */
        private final double inputScale;

        /** output scaling factor. */
        private final double outputScale;

        /** Simple constructor.
         * @param m parameter of the Jacobi elliptic function (must be greater than 1 here)
         */
        public Big(final double m) {
            algorithm   = selectAlgorithm(1.0 / m);
            inputScale  = FastMath.sqrt(m);
            outputScale = 1.0 / inputScale;
        }

        /** {@inheritDoc} */
        @Override
        public CopolarN values(final double u) {
            final CopolarN trioN = algorithm.values(u * inputScale);
            return new CopolarN(outputScale * trioN.sn(), trioN.dn(), trioN.cn());
        }

    }
    
    /** Copolar trio with pole at point n in Glaisher’s Notation.
     * <p>
     * This is a container for the three principal Jacobi elliptic functions
     * {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
     * </p>
     */
    public static class CopolarN {

        /** Value of the sn function. */
        private final double sn;

        /** Value of the cn function. */
        private final double cn;

        /** Value of the dn function. */
        private final double dn;

        /** Simple constructor.
         * @param sn value of the sn function
         * @param cn value of the cn function
         * @param dn value of the dn function
         */
        private CopolarN(final double sn, final double cn, final double dn) {
            this.sn = sn;
            this.cn = cn;
            this.dn = dn;
        }

        /** Get the value of the sn function.
         * @return sn(u|m)
         */
        public double sn() {
            return sn;
        }

        /** Get the value of the cn function.
         * @return cn(u|m)
         */
        public double cn() {
            return cn;
        }

        /** Get the value of the dn function.
         * @return dn(u|m)
         */
        public double dn() {
            return dn;
        }

    }

    /** Copolar trio with pole at point s in Glaisher’s Notation.
     * <p>
     * This is a container for the three subsidiary Jacobi elliptic functions
     * {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
     * </p>
     */
    public static class CopolarS {

        /** Value of the cs function. */
        private final double cs;

        /** Value of the dn function. */
        private final double ds;

        /** Value of the ns function. */
        private final double ns;

        /** Simple constructor.
         * @param trioN copolar trio with pole at point n in Glaisher’s Notation
         */
        private CopolarS(final CopolarN trioN) {
            this.ns = 1.0 / trioN.sn();
            this.cs = ns  * trioN.cn();
            this.ds = ns  * trioN.dn();
        }

        /** Get the value of the cs function.
         * @return cs(u|m)
         */
        public double cs() {
            return cs;
        }

        /** Get the value of the ds function.
         * @return ds(u|m)
         */
        public double ds() {
            return ds;
        }

        /** Get the value of the ns function.
         * @return ns(u|m)
         */
        public double ns() {
            return ns;
        }

    }

    /** Copolar trio with pole at point c in Glaisher’s Notation.
     * <p>
     * This is a container for the three subsidiary Jacobi elliptic functions
     * {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
     * </p>
     */
    public static class CopolarC {

        /** Value of the dc function. */
        private final double dc;

        /** Value of the nc function. */
        private final double nc;

        /** Value of the sc function. */
        private final double sc;

        /** Simple constructor.
         * @param trioN copolar trio with pole at point n in Glaisher’s Notation
         */
        private CopolarC(final CopolarN trioN) {
            this.nc = 1.0 / trioN.cn();
            this.sc = nc  * trioN.sn();
            this.dc = nc  * trioN.dn();
        }

        /** Get the value of the dc function.
         * @return dc(u|m)
         */
        public double dc() {
            return dc;
        }

        /** Get the value of the nc function.
         * @return nc(u|m)
         */
        public double nc() {
            return nc;
        }

        /** Get the value of the sc function.
         * @return sc(u|m)
         */
        public double sc() {
            return sc;
        }

    }

    /** Copolar trio with pole at point d in Glaisher’s Notation.
     * <p>
     * This is a container for the three subsidiary Jacobi elliptic functions
     * {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
     * </p>
     */
    public static class CopolarD {

        /** Value of the nd function. */
        private final double nd;

        /** Value of the sd function. */
        private final double sd;

        /** Value of the cd function. */
        private final double cd;

        /** Simple constructor.
         * @param trioN copolar trio with pole at point n in Glaisher’s Notation
         */
        private CopolarD(final CopolarN trioN) {
            this.nd = 1.0 / trioN.dn();
            this.sd = nd  * trioN.sn();
            this.cd = nd  * trioN.cn();
        }

        /** Get the value of the nd function.
         * @return nd(u|m)
         */
        public double nd() {
            return nd;
        }

        /** Get the value of the sd function.
         * @return sd(u|m)
         */
        public double sd() {
            return sd;
        }

        /** Get the value of the cd function.
         * @return cd(u|m)
         */
        public double cd() {
            return cd;
        }

    }

}
