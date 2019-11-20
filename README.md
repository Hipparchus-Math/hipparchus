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
==========

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

Alternatively you can pull it from the central Maven repositories using `pom.xml` settings:

```xml
<properties>
  <!-- change the Hipparchus version number to the one suiting your needs -->
  <myprojectname.hipparchus.version>1.6</myprojectname.hipparchus.version>
</properties>

<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-core</artifactId>
  <version>${myprojectname.hipparchus.version}</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-clustering</artifactId>
  <version>${myprojectname.hipparchus.version}</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-fft</artifactId>
  <version>${myprojectname.hipparchus.version}</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-fitting</artifactId>
  <version>${myprojectname.hipparchus.version}</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-geometry</artifactId>
  <version>${myprojectname.hipparchus.version}</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-ode</artifactId>
  <version>${myprojectname.hipparchus.version}</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-optim</artifactId>
  <version>${myprojectname.hipparchus.version}</version>
</dependency>
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-stat</artifactId>
  <version>${myprojectname.hipparchus.version}</version>
</dependency>

```

If your project previously depended on [Apache Commons Math](https://commons.apache.org/proper/commons-math/)
and you want to switch to Hipparchus, you can also add the temporary migration jar

```xml
<dependency>
  <groupId>org.hipparchus</groupId>
  <artifactId>hipparchus-migration</artifactId>
  <version>${myprojectname.hipparchus.version}</version>
</dependency>
```

Contributing
------------

There are some guidelines which
will make applying contributions easier for us. Please read through our
[contributing guidelines](https://github.com/Hipparchus-Math/hipparchus/blob/master/CONTRIBUTING.md).

To contact us, there are (temporarily) two communication channels. The recommended
one is a shared [forum](https://forum.orekit.org/categories) where several categories
are dedicated to Hipparchus. The oldedr one is a [developers mailing list](mailto:developers@hipparchus.org).
Beware that the mail-based channel is buggy (read loses messages) and unmaintained, so the forum
is probably a safer channel as of release 1.6.

License
-------

Code is under the [Apache Licence v2](https://www.apache.org/licenses/LICENSE-2.0.txt).
