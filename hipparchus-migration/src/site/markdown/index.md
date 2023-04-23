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

## Hipparchus Migration

This module is devoted to migrating from [Apache Commons Math](https://commons.apache.org/math/)
to [Hipparchus](../index.html).

The module provides a transition script and a temporary migration jar. The migration
process is composed of several phases:

  * make sure you don't use any class that has been already deprecated in Apache Commons Math
  * run the [migration script](script.html) on the application source
  * modify the build configuration to use Hipparchus modules instead of Apache Commons Math jar
  * fix the (hopefully few) compilation errors that were not handled by the script
  * check everything runs, ignoring the warnings about deprecated interfaces that are
    provided by the intermediate hipparchus-migration jar
  * progressively get rid of the deprecated interfaces, following the deprecation statements hints
  * remove the dependency to the intermediate hipparchus-migration jar

### Modifying the build configuration

As your application depended on Apache Commons Math which as only one jar, your build
system knows about this jar. It has to be replaced by the various jars provided by
Hipparchus. All uses will at least need the hipparchus-core jar. Complex uses may
need all the jars.

The simplest way to identify which jars you need is to add them one by one, looking
at the compilation errors. This is easy and fast as there are only a dozen jars
and you will certainly not depend on some of them (hipparchus-perf and hipparchus-samples
are typically not used). The first jar that will be mandatory is hipparchus-core. If this
one is not sufficient, the second one to add is hipparchus-migration because it adds some
glue code for migration that is often needed. Then add other jars depending on the compilation
errors you get, looking only at the package names of the missing classes.

### Fixing compilation errors

There may be some problems not handled by the script. The known problems
are the following ones:

  * specific exceptions not thrown anymore. The huge number of specialized exception
     has been drastically reduced with only `MathIllegalArgumentException`,
     `MathIllegalStateException`, both of which extending `MathRuntimeException`. The
     specialize `NullArgumentException` which extends the Java language standard
     `NullPointerException` and should almost never been caught. The recommended
     way to catch Hipparchus exceptions is therefore to catch _only_ `MathRuntimeException`
     and to call the new `getSpecifier` method if you are really looking for a specific
     exception, which is a very rare need.
  * `getContext()` method, `ExceptionContextProvider` interface and `ExceptionContext`
     class not found. The exception context feature that was present in Apache Commons
     Math has been merged into the `MathRuntimeException` base class. Its `getValue(key)`
     method is essentially replaced by `getParts[index]`.
  * classes from packages `org.apache.commons.math[34].geometry.partitioning.utilities`
     or `org.apache.commons.math[34].optimization` not found
     almost all the classes that were already deprecated in Apache Commons Math have been
     removed from Hipparchus, users must  replace them before performing the migration.
  * method `derivative()` in class `PolynomialFunction` not found.
     Despite this method was not deprecated it was removed as superseded by the
     already existing `polynomialDerivative` method which is fully compatible with it
     as it returns a more specialized type
  * `getDimension()` and `getFirstIndex() not found in `EquationsMapper` and
     `getSecondaryMappers()` not present in `ExpandableODE`.
     There is now only one equations mapper for both the primary equation and
     the various secondary equations in an ODE, so the way to map the composite
     state into primary and secondary states is different, and in fact much simpler.
     Users must retrieve the single mapper by calling `getMapper()` and then
     get directly the primary and secondary states from this single mapper.
   * `Action` enumerate not found. The `Action` enumerate used as a return
     value in events handlers `eventOccurred` method was an internal enumerate
     defined in the `EventHandler` interface, so no dedicated `import` statement
     was needed. It is now a separate enumerate, shared by both the double-based
     and field-based events handlers, and needs to be imported specifically.
   * conflicts between application class names and Apache Commons Math
     class names. As the script is based on text substitution without knowledge
     about the Java language syntax, names from the application may be confused
     by the script with similar names in the Apache Commons Math library even if
     they are in different packages. Such names may be converted despite they should
     not be touched, so the change must be reverted manually for these classes.

### Checking everything runs

Once everything compiles, you should run your application tests to see if everything
runs. If not, ask to the [developers mailing lists](../mail-lists.html).

Deprecation warnings at this step are normal. They correspond to APIs that were present
in Apache Commons Math and have been modified in Hipparchus. In order to help users
switch, we have moved these APIs in the hipparchus-migration jar and they act as
facades between user code and the new Hipparchus APIs. They are marked as deprecated
so users can easily spot them and change to the newer APIs progressively while still
having a running application.

### Getting rid of the deprecated interfaces

When the application runs in a full Hipparchus environment, you can start updating
your code to use the new APIs provided by Hipparchus. This is just a traditional
maintenance task. The deprecation hints in the hipparchus-migration sources will
help you find the new interfaces. They are often very close to the previous
ones and modification is expected to be straightforward.

If you want to get rid of all problems at once, you can also remove the migration
jar from your build configuration, so instead of deprecation warnings, you will
get compilation errors which are easier to spot in integrated development
environments.

### Removing the dependency to the intermediate hipparchus-migration jar

When all deprecation warnings have been fixed, you can safely remove the
hipparchus-migration jar from your build system, everything should work.
