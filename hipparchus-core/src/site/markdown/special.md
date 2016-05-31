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

(see [MathWorld](http://mathworld.wolfram.com/RegularizedBetaFunction.html)\)
