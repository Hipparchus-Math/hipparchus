<!---
 Licensed to the Hipparchus project under one or more
 contributor license agreements.  See the NOTICE file distributed with
 this work for additional information regarding copyright ownership.
 The Hipparchus project this file to You under the Apache License, Version 2.0
 (the "License"); you may not use this file except in compliance with
 the License.  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->
Hipparchus
===================

[![License](http://img.shields.io/:license-apache-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0.html)

The Hipparchus project is a library of lightweight, self-contained
mathematics and statistics components addressing the most common
problems not available in the Java programming language.

Documentation
-------------

More information can be found on the [homepage](https://hipparchus.org/).
The [JavaDoc](https://hipparchus.org/apidocs) can be browsed.
Questions related to the usage of Hipparchus should be posted to the [users mailing list](mailto:users@hipparchus.org).

Where can I get the latest release?
-----------------------------------
You can download source and binaries from our [download page](https://hipparchus.org/downloads.html).

Alternatively you can pull it from the central Maven repositories:

```xml
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-core</artifactId>
  <version>1.1</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-clustering</artifactId>
  <version>1.1</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-fft</artifactId>
  <version>1.1</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-fitting</artifactId>
  <version>1.1</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-geometry</artifactId>
  <version>1.1</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-ode</artifactId>
  <version>1.1</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-optim</artifactId>
  <version>1.1</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-stat</artifactId>
  <version>1.1</version>
</dependency>

```

If your project previously depended on [Apache Commons Math](http://commons.apache.org/commons-math/)
and you want to switch to Hipparchus, you can also add the temporary migration jar

```xml
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-migration</artifactId>
  <version>1.1</version>
</dependency>
```

Contributing
------------

We accept PRs via github. The [developers mailing list](mailto:developers@hipparchus.org) is the main channel of communication for contributors.
There are some guidelines which will make applying PRs easier for us:

+ No tabs! Please use spaces for indentation.
+ Respect the code style.
+ Create minimal diffs - disable on save actions like reformat source code or organize imports. If you feel the source code should be reformatted create a separate PR for this change.
+ Provide JUnit tests for your changes and make sure your changes don't break any existing tests by running ```mvn clean test```.

If you plan to contribute on a regular basis, please consider filing a [Contributor License Agreement](https://hipparchus.org/clas/ICLA.pdf).
You can learn more about contributing to Hipparchus via our [developers guide](https://www.hipparchus.org/developers.html).
License
-------
Code is under the [Apache Licence v2](https://www.apache.org/licenses/LICENSE-2.0.txt).

Additional Resources
--------------------

+ [Project home page](https://hipparchus.org/)
+ [Contributor License Agreement](https://hipparchus.org/clas/ICLA.pdf)
+ [Mailing lists](https://hipparchus.org/mail-lists.html)
