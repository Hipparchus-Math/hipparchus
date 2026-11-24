<!--
 Licensed to the Hipparchus project under one or more
 contributor license agreements.  See the NOTICE file distributed with
 this work for additional information regarding copyright ownership.
 The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 (the "License"); you may not use this file except in compliance with
 the License.  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->
## Hipparchus: a mathematics Library

Hipparchus is a library of lightweight, self-contained
 mathematics and statistics components addressing the most common
 problems not available in the Java programming language.

## Modularized

Hipparchus is a large library, so it has been split into several
modules. Users can select only the modules they need
and developers can update modules more easily and more
often than would be possible in a large monolithic library.

---

## Guiding principles

  * Real-world application use cases determine development priority.
  * This package emphasizes small, easily integrated components
     rather than large libraries with complex dependencies and
     configurations.
  * All algorithms are fully documented and follow generally
     accepted best practices.
  * In situations where multiple standard algorithms exist, a
     Strategy pattern is used to support multiple implementations.
  * Limited dependencies. No external dependencies beyond the
     core Java platform (at least Java 21 starting with version
     5.0 of the library).

## Community

There are two channels for community discussions. The first one is a
forum, there are three categories [Hipparchus
annoucements](https://forum.orekit.org/c/hipparchus-announcements),
[Hipparchus usage](https://forum.orekit.org/c/hipparchus-usage) and
[Hipparchus
development](https://forum.orekit.org/c/hipparchus-development) in the
forum. Regular users should subscribe to the forum and configure their
preferences to add these categories to the "Watched" list so they are
notified when a new message appears there. The second one is a
[GitHub discussions](https://github.com/Hipparchus-Math/hipparchus/discussions)
page. Users are free to select the channel they are more comfortable
with.

## Fork
 
Hipparchus started as a fork of [Apache Commons Math](https://commons.apache.org/math/).
The fork was initiated by most of the main developers and a few contributors of
Apache Commons Math.

Version 1.0 of Hipparchus is therefore very similar to
Apache Commons Math 3.6.1 with some elements of the development version
that would ultimately lead to Apache Commons Math 4.0 which was not released
at fork time.

 ---

## 5.0 Release is out!

Hipparchus 5.0 is now available for download from the
[Hipparchus download page](downloads.html) or on Maven central and its
mirrors under the groupId **org.hipparchus**.

Hipparchus version 5.0. The main changes are:

   * Fixed log normal distribution inverse CDF.
   * Removed deprecated AbstractUnivariateDifferentiableSolver.computeObjectiveValueAndDerivative(double), replaced by computeObjectiveValueAndDerivative1(double).
   * Removed deprecated SQPOption.getConvCriteria, setConvCriteria and DEFAULT_CONV_CRITERIA.
   * Removed SQPOption qpMaxLoop and useFunHessian settings (getters, setters and default constants), which were used only by the removed SQPOptimizerS and SQPOptimizerGM.
   * Changed ADMMQPOption.DEFAULT_SCALING from true to false.
   * Removed deprecated classes and methods in geometry and optimizers.
   * Fixed Complex isFinite method.
   * Fixed wrong advertised order of Adams-Moulton integrator (double and Field versions).
   * Added Software Build Of Materials generation.
   * Fixed Wilcoxon signed rank test whenever the positive and negative rank sums are equal.
   * Moved DistanceMeasure from clustering module to core module.
   * Added DataPoint interface and replaced Clusterable and Fittable interfaces.
   * Added compute(DataPoint, DataPoint) to DistanceMeasure.
   * Fixed documentation.
