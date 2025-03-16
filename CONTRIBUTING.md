<!---
 Licensed to the Hipparchus project under one or more
 contributor license agreements.  See the NOTICE file distributed with
 this work for additional information regarding copyright ownership.
 The Hipparchus project this file to You under the Apache License, Version 2.0
 (the "License"); you may not use this file except in compliance with
 the License.  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->
Contributing to Hipparchus
==========================

You have found a bug or you have an idea for a cool new feature? Contributing code is a great way to give something back to
the open source community. Before you dig right into the code there are a few guidelines that we need contributors to
follow so that we can have a chance of keeping on top of things.

Getting Started
---------------

+ Make sure you have a [GitHub account](https://github.com/signup/free).
+ If you're planning to implement a new feature it makes sense to discuss your changes on the
  [forum](https://forum.orekit.org/c/hipparchus-development) first. This way you can make sure
  you're not wasting your time on something that isn't considered to be in Hipparchus's scope.
+ Submit a ticket for your issue, assuming one does not already exist.
  + Clearly describe the issue including steps to reproduce when it is a bug.
  + Make sure you fill in the earliest version that you know has the issue.
+ Clone the repository on GitHub.

Making Changes
--------------

+ Create a topic branch from where you want to base your work (this is usually the main branch).
+ Make commits of logical units.
+ Respect the original code style:
  + Only use spaces for indentation.
  + Create minimal diffs - disable on save actions like reformat source code or organize imports. If you feel the source code should be reformatted create a separate PR for this change.
  + Check for unnecessary whitespace with git diff --check before committing.
+ Make sure you have added the necessary tests for your changes.
  + Run all the tests with `mvn clean verify` to assure nothing else was accidentally broken.
+ Edit the relevant `changes.xml` file(s).
+ Make sure your commit messages are in the proper format.

Making Trivial Changes
----------------------

For changes of a trivial nature to comments and documentation, it is not always necessary to create a new ticket in GitHub.
In this case, it is appropriate to start the first line of a commit with '(doc)' instead of a ticket number.

Submitting Changes
------------------

+ Sign the [Contributor License Agreement](https://hipparchus.org/licenses/#clas) if you haven't already.
+ Push your changes to a topic branch in your fork of the repository.
+ Submit a pull request to the [repository](https://github.com/Hipparchus-Math/hipparchus).
+ Update your GitHub ticket and include a link to the pull request in the ticket.

Additional Resources
--------------------

+ [Project home page](https://hipparchus.org/)
+ [General GitHub documentation](https://help.github.com/)
+ [GitHub pull request documentation](https://help.github.com/send-pull-requests/)
+ [Contributor License Agreement ](https://hipparchus.org/licenses/#clas)
+ [Developers forum](https://forum.orekit.org/c/hipparchus-development)
+ [Users forum](https://forum.orekit.org/c/hipparchus-usage)
+ [Announcements forum](https://forum.orekit.org/c/hipparchus-announcements)
