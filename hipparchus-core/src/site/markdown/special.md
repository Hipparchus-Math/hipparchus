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
# Special Functions

## Overview

The `special` package of Hipparchus gathers several useful special
functions not provided by `java.lang.Math`.

## Erf functions

[Erf](../apidocs/org.hipparchus/special/Erf.html) contains several useful functions
involving the Error Function, Erf.

|  Function      |  Method |  Reference                                             |
|----------------|---------|--------------------------------------------------------|
| Error Function |  erf    | See [MathWorld](http://mathworld.wolfram.com/Erf.html) |

## Gamma functions

Class [`Gamma`](../apidocs/org.hipparchus/special/Gamma.html)
contains several useful functions involving the Gamma Function.

### Gamma

`Gamma.gamma(x)` computes the Gamma function, `\(\Gamma(x)\)`
(see [MathWorld](http://mathworld.wolfram.com/GammaFunction.html),
[DLMF](http://dlmf.nist.gov/5)\). The accuracy of the Hipparchus
implementation is assessed by comparison with high precision values computed
with the [Maxima](http://maxima.sourceforge.net/) Computer Algebra System.

|  Interval           |  Values tested                                |  Average error  |  Standard deviation  |  Maximum error |
|---------------------|-----------------------------------------------|-----------------|----------------------|----------------|
| `\(-5  < x  < -4\)` | `x[i] = i / 1024, i = -5119, ..., -4097`      | 0.49 ulps       | 0.57 ulps            | 3.0 ulps       |
| `\(-4  < x  < -3\)` | `x[i] = i / 1024, i = -4095, ..., -3073`      | 0.36 ulps       | 0.51 ulps            | 2.0 ulps       |
| `\(-3  < x  < -2\)` | `x[i] = i / 1024, i = -3071, ..., -2049`      | 0.41 ulps       | 0.53 ulps            | 2.0 ulps       |
| `\(-2  < x  < -1\)` | `x[i] = i / 1024, i = -2047, ..., -1025`      | 0.37 ulps       | 0.50 ulps            | 2.0 ulps       |
| `\(-1  < x  < 0\)`  | `x[i] = i / 1024, i = -1023, ..., -1`         | 0.46 ulps       | 0.54 ulps            | 2.0 ulps       |
| `\(0  < x ≤ 8\)`    | `x[i] = i / 1024, i = 1, ..., 8192`           | 0.33 ulps       | 0.48 ulps            | 2.0 ulps       |
| `\(8  < x ≤ 141\)`  | `x[i] = i / 64, i = 513, ..., 9024`           | 1.32 ulps       | 1.19 ulps            | 7.0 ulps       |

### Log Gamma

`Gamma.logGamma(x)` computes the natural logarithm of the Gamma function,
`\(\log \Gamma(x)\)`, for `\(x > 0\)`
(see [MathWorld](http://mathworld.wolfram.com/LogGammaFunction.html),
[DLMF](http://dlmf.nist.gov/5)\). The accuracy of the Hipparchus
implementation is assessed by comparison with high precision values computed
with the [Maxima](http://maxima.sourceforge.net/) Computer Algebra System.

|  Interval                                               |  Values tested                                |  Average error  |  Standard deviation  |  Maximum error |
|---------------------------------------------------------|-----------------------------------------------|-----------------|----------------------|----------------|
| `\(0  < x \le 8\)`                                      | `x[i] = i / 1024, i = 1, ..., 8192`           | 0.32 ulps       | 0.50 ulps            | 4.0 ulps       |
| `\(8  < x \le 1024\)`                                   | `x[i] = i / 8, i = 65, ..., 8192`             | 0.43 ulps       | 0.53 ulps            | 3.0 ulps       |
| `\(1024  < x \le 8192\)`                                | `x[i], i = 1025, ..., 8192`                   | 0.53 ulps       | 0.56 ulps            | 3.0 ulps       |
| `\(8933.439345993791 \le x \le 1.75555970201398e+305\)` | `x[i] = 2**(i / 8), i = 105, ..., 8112`       | 0.35 ulps       | 0.49 ulps            | 2.0 ulps       |
                                                                                                      
### Regularized Gamma

`Gamma.regularizedGammaP(a, x)` computes the value of the regularized
Gamma function, P(a, x)
(see [MathWorld](http://mathworld.wolfram.com/RegularizedGammaFunction.html)\).

## Beta functions

[Beta](../apidocs/org.hipparchus/special/Beta.html) contains
several useful functions involving the Beta Function.

### Log Beta

`Beta.logBeta(a, b)` computes the value of the natural logarithm of the
Beta function, log B(a, b).
(see [MathWorld](http://mathworld.wolfram.com/BetaFunction.html),
[DLMF](http://dlmf.nist.gov/5.12)\). The accuracy of the Hipparchus
implementation is assessed by comparison with high precision values computed
with the [Maxima](http://maxima.sourceforge.net/) Computer Algebra System.

|  Interval             |  Values tested                          |  Average error  |  Standard deviation  |  Maximum error  |
|-----------------------|-----------------------------------------|-----------------|----------------------|-----------------|
| `\(0  < x \le 8\)`<br/>`\(0  < y \le 8\)`  | `x[i] = i / 32, i = 1, ..., 256`<br/>`y[j] = j / 32, j = 1, ..., 256`   | 1.80 ulps       | 81.08 ulps           | 14031.0 ulps    |
| `\(0  < x \le 8\)`<br/>`\(8  < y \le 16\)` | `x[i] = i / 32, i = 1, ..., 256`<br/>`y[j] = j / 32, j = 257, ..., 512` | 0.50 ulps       | 3.64 ulps            | 694.0 ulps      |
| `\(0  < x \le 8\)`<br/>`\(16  < y \le 256\)` | `x[i] = i / 32, i = 1, ..., 256`<br/>`y[j] = j, j = 17, ..., 256` | 1.04 ulps       | 139.32 ulps          | 34509.0 ulps    |
| `\(8  < x \le 16\)`<br/>`\(8  < y \le 16\)`  | `x[i] = i / 32, i = 257, ..., 512`<br/>`y[j] = j / 32, j = 257, ..., 512` | 0.35 ulps       | 0.48 ulps            | 2.0 ulps        |
| `\(8  < x \le 16\)`<br/>`\(16  < y \le 256\)`   | `x[i] = i / 32, i = 257, ..., 512`<br/>`y[j] = j, j = 17, ..., 256` | 0.31 ulps       | 0.47 ulps            | 2.0 ulps        |
| `\(16  < x \le 256\)`<br/>`\(16  < y \le 256\)` | `x[i] = i, i = 17, ..., 256`<br/>`y[j] = j, j = 17, ..., 256`       | 0.35 ulps       | 0.49 ulps            | 2.0 ulps        |

### Regularized Beta

(see [MathWorld](http://mathworld.wolfram.com/RegularizedBetaFunction.html))

### Elliptic functions and integrals

Notations in the domain of elliptic functions and integrals is often confusing and inconsistent
across text books. Hipparchus implementation uses the `parameter m` to define both Jacobi elliptic
functions and Legendre elliptic integrals and uses the `nome q` to define Jacobi theta functions.
The elliptic modulus `k` (which is the square of parameter m) is not used at all in Hipparchus. All these
parameters are linked together.

See in MathWorld [elliptic integrals of the first kind](https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html),
[elliptic integrals of the second kind](https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html),
[elliptic integrals of the third kind](https://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html),
[Jacobi elliptic functions](https://mathworld.wolfram.com/JacobiEllipticFunctions.html) and
[Jacobi theta functions](https://mathworld.wolfram.com/JacobiThetaFunctions.html).

`JacobiEllipticBuilder.build(m)` builds a `JacobiElliptic` (or `FieldJacobiElliptic`) implementation for the parameter
`m` that computes the twelve elliptic functions `\(sn(u|m)\)`, `\(cn(u|m)\)`, `\(dn(u|m)\)`,
`\(cs(u|m)\)`, `\(ds(u|m)\)`, `\(ns(u|m)\)`, `\(dc(u|m)\)`, `\(nc(u|m)\)`, `\(sc(u|m)\)`, `\(nd(u|m)\)`,
`\(sd(u|m)\)`, and `\(cd(u|m)\)`. The functions are computed as copolar triplets as when one function is needed
in an expression, the two other are often also needed. The inverse functions `\(arcsn(u|m)\)`, `\(arccn(u|m)\)`, `\(arcdn(u|m)\)`,
`\(arccs(u|m)\)`, `\(arcds(u|m)\)`, `\(arcns(u|m)\)`, `\(arcdc(u|m)\)`, `\(arcnc(u|m)\)`, `\(arcsc(u|m)\)`, `\(arcnd(u|m)\)`,
`\(arcsd(u|m)\)`, and `\(arccd(u|m)\)` are also available.

`JacobiTheta` (and `FieldJacobiTheta`) computes the four Jacobi theta functions `\(\theta_1(z|\tau)\)`, `\(\theta_2(z|\tau)\)`,
`\(\theta_3(z|\tau)\)`, and `\(\theta_4(z|\tau)\)`. The half-period ratio `\(\tau\)` is linked to the nome `q`:
`\(q = e^{i\pi\tau}\)`. Here again, the four functions are computed at once and a quadruplet is returned.

`CarlsonEllipticIntegrals` is a utility class that computes the following integrals in Carlson symmetric form,
for both primitive double, `CalculusFieldElement`, `Complex` and `FieldComplex`:

| Name               |  Definition                                                                                                         |
|--------------------|---------------------------------------------------------------------------------------------------------------------|
| `\(R_F(x,y,z)\)`   | `\(\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)}\)`                                                          |
| `\(R_J(x,y,z,p)\)` | `\(\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)(t+p)}\)`                                                     |
| `\(R_G(x,y,z)\)`   | `\(\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}\left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t \)` |
| `\(R_D(x,y,z)\)`   | `\(R_J(x,y,z,z)\)`                                                                                                  |
| `\(R_C(x,y)\)`     | `\(R_F(x,y,y)\)`                                                                                                    |

where `\(s(t) = \sqrt{t+x}\sqrt{t+y}\sqrt{t+z}\)`.

`LegendreEllipticIntegrals` is a utility class that computes the following integrals,
for both primitive double, `CalculusFieldElement`, `Complex` and `FieldComplex`.
(the implementation uses `CarlsonEllipticIntegrals` internally):

| Name                  | Type       |  Definition                                                                            |
|-----------------------|------------|----------------------------------------------------------------------------------------|
| `\(K(m)\)`            |  complete  | `\(\int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}\)`                   |
| `\(K'(m)\)`           |  complete  | `\(\int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-(1-m) \sin^2\theta}}\)`               |
| `\(E(m) \)`           |  complete  | `\(\int_0^{\frac{\pi}{2}} \sqrt{1-m \sin^2\theta} d\theta\)`                           |
| `\(D(m) \)`           |  complete  | `\(\frac{K(m) - E(m)}{m}\)`                                                            |
| `\(\Pi(n, m)\)`       |  complete  | `\(\int_0^{\frac{\pi}{2}} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}\)` |
| `\(F(\phi, m)\)`      | incomplete | `\(\int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}}\)`                            |
| `\(E(\phi, m)\)`      | incomplete | `\(\int_0^{\phi} \sqrt{1-m \sin^2\theta} d\theta\)`                                    |
| `\(D(\phi, m)\)`      | incomplete | `\(\frac{K(\phi, m) - E(\phi, m)}{m}\)`                                                |
| `\(\Pi(n, \phi, m)\)` | incomplete | `\(\int_0^{\phi} \frac{d\theta}{\sqrt{1-m \sin^2\theta}(1-n \sin^2\theta)}\)`          |


