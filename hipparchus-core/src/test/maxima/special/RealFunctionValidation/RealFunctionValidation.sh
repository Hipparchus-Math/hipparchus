#!/bin/bash
#
# Licensed to the Apache Software Foundation (ASF) under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Location of the Hipparchus jar file
HIPPARCHUS_JAR=$HOME/.m2/repository/org/hipparchus/hipparchus-core/1.0-SNAPSHOT/hipparchus-core-1.0-SNAPSHOT.jar

# Location of file RealFunctionValidation.jar
APP_JAR=$HOME/Documents/workspace/hipparchus/hipparchus-core/src/test/maxima/special/RealFunctionValidation/RealFunctionValidation.jar

java -cp $HIPPARCHUS_JAR:$APP_JAR RealFunctionValidation $1
