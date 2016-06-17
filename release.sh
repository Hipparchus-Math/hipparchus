#!/bin/sh
#
# Licensed to the Hipparchus under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# Hipparchus licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------
#
# Script to create a full release candidate (all modules).
#
# Run from top-level directory in full hipparchus checkout.
# Signed and hashed .bz2 files are created in /target.
#

# Config
keyID=0xCC6CB50E6A59DD12
version=1.0
notes=false

# Temporarily drop -SNAPSHOT from version
# find . -name pom.xml -exec sed -i '' 's/\-SNAPSHOT//g' {} \;

# Release notes - only top-level, refer to the online changelogs for modules
# Set notes=true to execute this, typically just for the first RC.
# Remember to check the edited notes in after executing this.
if $notes
then
	mkdir -p hipparchus-parent/src/changes
	cp src/changes/* hipparchus-parent/src/changes/
	mvn -Dversion=$version changes:announcement-generate -pl hipparchus-parent
	cp hipparchus-parent/target/announcement/release-notes.vm RELEASE-NOTES.txt
	rm -rf hipparchus-parent/src
    # Maybe hack the release notes a little ....
    read -p "Press [Enter] key start building tarballs"
fi

# Stage the release artifacts
mvn clean
mvn deploy -Dgpg.keyname=$keyID -Dwagon.skip=true -DskipStagingRepositoryClose=true -Prelease

# Cleanup target
# Add hashes and remove pom files 
cd target
for f in *.zip
do
	md5 $f > $f.md5
done
for f in *.bz2
do
	md5 $f > $f.md5
done
rm -rf archive-tmp*
rm *.pom
rm *.pom.asc
cd ..

# Undo pom changes
#find . -name pom.xml -exec git checkout {} \;


