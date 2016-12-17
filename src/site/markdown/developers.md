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
# Developers Guide

## Aims

Creating and maintaining a mathematical and statistical library that is
accurate requires a greater degree of communication than might be the
case for other components. It is important that developers follow
guidelines laid down by the community to ensure that the code they create
can be successfully maintained by others.

## Guidelines

Developers are asked to comply with the following development guidelines.
Code that does not comply with the guidelines including the word _must_
will not be committed.  Our aim will be to fix all of the exceptions to the
_should_ guidelines prior to a release.

## Contributing

### Getting Started

  * Start by reviewing the overall objectives stated in the
    [main page](index.html) upon which the project is founded.
  * Download the Hipparchus source code. The git url for the
    current development sources of Hipparchus
    [https://github.com/Hipparchus-Math/hipparchus](https://github.com/Hipparchus-Math/hipparchus)
  * Hipparchus uses Apache Maven as our build tool. The sources can also be built
    using any IDE. To build Hipparchus, you can follow the instructions
    from the [building](building.html) page
  * Be sure to join the users and developers [email lists](mail-lists.html)
    and use them appropriately. Make any proposals here where the group can
    comment on them.
  * Setup an account on GitHub and use it to submit pull requests and
    identify bugs. Read the [GitHub help](https://help.github.com/).
  * Generating patches: The requested format for generating patches is
    either GitHub pull requests or patches using the Unified Diff format,
    which can be easily generated using the git client or various IDEs.

### Contributing ideas and code

Follow the steps below when making suggestions for additions or
enhancements to Hipparchus. This will make it easier for the community
to comment on your ideas and for the committers to keep track of them. 
Thanks in advance!

  * Start with a post to the developers mailing list, with a good, short title
    describing the new feature or enhancement.  For example,
    "Principal Components Analysis." The body of the post should include each
    of the following items (but be _as brief as possible_):

      * A concise description of the new feature / enhancement

      * References to definitions and algorithms. Using standard
        definitions and algorithms makes communication much easier and will
        greatly increase the chances that we will accept the code / idea

      * Some indication of why the addition / enhancement is practically useful

  * Assuming a generally favorable response to the idea on developers list,
    the next step is to create a ticket for a new issue in the
    Hipparchus [issues tracker](https://github.com/Hipparchus-Math/hipparchus/issues)
  * Submit code as [pull requests](https://github.com/Hipparchus-Math/hipparchus/pulls)
    or as patchs. Please use one ticket for each feature, adding multiple pull requests
    to the ticket as necessary.  Use the git diff command to generate your patches as
    diffs. Please do not submit modified copies of existing java files. Be
    patient (but not _too_ patient) with  committers reviewing patches. Post a
    _nudge_ message to developers list with a reference to the
    issue if a patch goes more than a few days with no comment or commit.

## Coding Style

Hipparchus follows a specific coding style, similar to the ones used
in Apache Commons projects or the Orekit project. The reference
settings correspond to the `src/conf/checkstyle.xml` file found in the
Hipparchus source tree. There is also a `src/conf/hipparchus-eclipse.xml`
file that configures the eclipse formater in a mostly compatible way
(there are some parts that cannot be configured automatically in
Eclipse to match checkstyle settings, mainly in the indentation parts).
One thing that Checkstyle will complain about is tabs included in the source code.
Please make sure to set your IDE or editor to use spaces instead of tabs.

Committers should configure the `user.name`, `user.email` and `core.autocrlf`
git repository or global settings with `git config`.
The first two settings define the identity and mail of the committer.
The third setting deals with line endings to achieve consistency
in line endings. Windows users should configure this setting to
`true` (thus forcing git to convert CR/LF line endings
in the workspace while maintaining LF only line endings in the repository)
while OS X and Linux users should configure it to `input`
(thus forcing git to only strip accidental CR/LF when committing into
the repository, but never when cheking out files from the repository).
See [Customizing Git - Git Configuration](http://www.git-scm.com/book/en/Customizing-Git-Git-Configuration)
in the git book for explanation about how to configure these settings and more.

## Documentation

  * Committed code _must_ include full javadoc.
  * All component contracts _must_ be fully specified in the javadoc class,
    interface or method comments, including specification of acceptable ranges
    of values, exceptions or special return values.
  * External references or full statements of definitions for all mathematical
    terms used in component documentation _must_ be provided.
  * Hipparchus javadoc generation supports embedded LaTeX formulas via the
    [MathJax](http://www.mathjax.org) javascript display engine. To
    embed mathematical expressions formatted in LaTeX in javadoc, simply surround
    the expression to be formatted with either `\\(` and  `\\)` for inline
    formulas, or `\\[` and `\\]` to have the formula appear on a separate line.
    For example, `\\(a^2 + b^2 = c^2\\)` will render an in-line formula
    saying that (a, b, c) is Pythagorean triplet.  Using `\\[` and `\\]` on
    the ends will render the same formula on a separate line.  See the MathJax
    and LaTex documentation for details on how to represent formulas and
    escape special characters.
  * Hipparchus modules documentation for the web site is generated using
    markdown syntax. It also supports embedded LaTeX formulas via the
    [MathJax](http://www.mathjax.org) javascript display engine. To
    embed mathematical expressions formatted in LaTeX in markdown, simply surround
    the expression to be formatted with either `\`\\(` and  `\\)\`` for inline
    formulas, or `\`\\[` and `\\]\`` to have the formula appear on a separate line.
    For example, `\`\\(a^2 + b^2 = c^2\\)\`` will render an in-line formula
    saying that (a, b, c) is Pythagorean triplet.  Using `\`\\[` and `\\]\`` on
    the ends will render the same formula on a separate line.  See the MathJax,
    markdown and LaTex documentation for details on how to represent formulas and
    escape special characters.
  * Implementations _should_ use standard algorithms and references or full
    descriptions of all algorithms _should_ be provided.
  * Additions and enhancements _should_ include updates to the User Guide.

## Exceptions

  * Exceptions generated by Hipparchus are all unchecked.
  * All public methods advertise all exceptions that they can generate.
    Exceptions _must_ be documented in both javadoc and method signatures
    and the documentation in the javadoc _must_ include full description
    of the conditions under which exceptions are thrown.
  * Methods _should_ fully specify parameter preconditions required for
    successful activation.  When preconditions are violated, a
    MathIllegalArgumentException should be thrown. Exception
    messages _must_ contain sufficient information on parameter values to
    determine the exact precondition failure.

## Unit Tests

  * Committed code _must_ include unit tests.
  * Unit tests _should_ provide full path coverage.
  * Unit tests _should_ verify all boundary conditions specified in
    interface contracts, including verification that exceptions are thrown or
    special values (e.g. `Double.NaN`, `Double.Infinity`) are returned as
    expected.

## Licensing and copyright

  * All new source file submissions _must_ include the Apache Software
    License in a comment that begins the file
  * All contributions must comply with the terms of the Hipparchus
    Individual Contributor License Agreement (ICLA)
    ([Open Document Format](https://hipparchus.org/clas/ICLA.odt),
    [Portable Document Format](https://hipparchus.org/clas/ICLA.pdf))
    or Corporate Contributor License Agreement (CCLA)
    ([Open Document Format](https://hipparchus.org/clas/CCLA.odt),
    [Portable Document Format](https://hipparchus.org/clas/CCLA.pdf))
  * Patches _must_ be accompanied by a clear reference to a _source_, - if code has been
    _ported_ from another language, clearly state the source of the original implementation.
    If the _expression_ of a given algorithm is derivative, please note the original source
    (textbook, paper, etc.).
  * References to source materials covered by restrictive proprietary
    licenses should be avoided.  In particular, contributions should not
    implement or include references to algorithms in
    [Numerical Recipes (NR)](http://www.nr.com/). Any questions about copyright or
    patent issues should be raised on the developers mailing list before contributing or
    committing code.

## Recommended Readings

Here is a list of relevant materials.  Much of the discussion surrounding
the development of this component will refer to the various sources
listed below, and frequently the Javadoc for a particular class or
interface will link to a definition contained in these documents.

### Concerning floating point arithmetic.

  * [http://www.validlab.com/goldberg/paper.pdf](http://www.validlab.com/goldberg/paper.pdf)
  * [http://www.cs.berkeley.edu/~wkahan/ieee754status/IEEE754.PDF](http://www.cs.berkeley.edu/~wkahan/ieee754status/IEEE754.PDF)
  * [http://www.cs.berkeley.edu/~wkahan/JAVAhurt.pdf](http://www.cs.berkeley.edu/~wkahan/JAVAhurt.pdf)

### Numerical analysis

  * [Scientific Computing FAQ @ Mathcom](http://www.mathcom.com/corpdir/techinfo.mdir/scifaq/index.html)
  * [Bibliography of accuracy and stability of numerical algorithms](http://www.ma.man.ac.uk/~higham/asna/asna2.pdf)
  * [SUNY Stony Brook numerical methods page](http://tonic.physics.sunysb.edu/docs/num_meth.html)
  * [SIAM Journal of Numerical Analysis Online](http://epubs.siam.org/sam-bin/dbq/toclist/SINUM)

### Probability and statistics

  * [Statlib at CMU](http://lib.stat.cmu.edu/)
  * [NIST Engineering Statistics Handbook](http://www.itl.nist.gov/div898/handbook/)
  * [Online Introductory Statistics (David W. Stockburger)](http://www.psychstat.smsu.edu/sbk00.htm)
  * [Online Journal of Statistical Software](http://www.jstatsoft.org/)

### References for mathematical definitions

  * [http://rd11.web.cern.ch/RD11/rkb/titleA.html](http://rd11.web.cern.ch/RD11/rkb/titleA.html)
  * [http://mathworld.wolfram.com/](http://mathworld.wolfram.com/)
  * [http://www.itl.nist.gov/div898/handbook](http://www.itl.nist.gov/div898/handbook)
  * [Chan, T. F. and J. G. Lewis 1979, _Communications of the ACM_, vol. 22 no. 9, pp. 526-531](http://doi.acm.org/10.1145/359146.359152)
  * [http://www.openmath.org](http://www.openmath.org)
