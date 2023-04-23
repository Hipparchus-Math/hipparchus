#!/usr/bin/env python

#   Licensed to the Hipparchus project under one or more
#   contributor license agreements.  See the NOTICE file distributed with
#   this work for additional information regarding copyright ownership.
#   The Hipparchus project licenses this file to You under the Apache License, Version 2.0
#   (the "License"); you may not use this file except in compliance with
#   the License.  You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import argparse
import atexit
import filecmp
import shutil
import re
import tempfile
import os

def is_valid_directory(parser, arg):
    if not os.path.isdir(arg):
        parser.error('directory {} does not exist!'.format(arg))
    else:
        # File exists so return the directory
        return arg

def fnmatch(fname, extensions):
	for ext in extensions:
		if fname.endswith(ext):
			return True
	return False

def log(msg, args, newline = True):
	if args.verbose == True:
		if newline == True:
			print msg
		else:
			print msg,

def processDir(rootDir, args):
	for dirName, subdirList, fileList in os.walk(rootDir):
		for fname in fileList:
			if fnmatch(fname, args.ext):
				fullName = os.path.join(dirName, fname)
				processFile(fullName, args)
		
		# do not traverse into subdirs that match an ignore
		for ignore in args.ignore:
			if ignore in subdirList:
				index = subdirList.index(ignore)
				if index >= 0:
					del subdirList[index]


def processFile(file, args):
	log('processing file ' + file, args, False)

	g = os.path.basename(file)
	tmpfile = os.path.join(args.tempdir, g)
	with open(file, 'r') as sources:
		lines = sources.readlines()
		
	modified = processLines(lines, args)

	with open(tmpfile, 'w') as sources:
		for line in modified:
			sources.write(line)

	if filecmp.cmp(file, tmpfile) == False:
		log('-> changed', args)

		if args.dry_run == True:
			if args.verbose == False:
				print file
			os.remove(tmpfile)
		else:
			if args.save == True:
				shutil.move(file, file + '.orig')
			shutil.move(tmpfile, file)
	else:
		log('-> unchanged', args)

def processLines(lines, args):
	imported_subpackages = {}
	found_classes = {}

	import_pattern = re.compile(r'^\s*import\s+([a-zA-Z0-9\.\*]+);')
	static_import_pattern = re.compile(r'^\s*import\s+static\s+([a-zA-Z0-9\.\*]+);')

	# collect all classnames that appear in the
	# source file, except if they are in import
	# statements.
	for line in lines:
		m = re.match(import_pattern, line)
		if m == None:
			if re.match(static_import_pattern, line) != None:
				continue

			for classname, subpackage in args.classnames.iteritems():
				if subpackage in imported_subpackages and line.find(classname) >= 0:
					# search rule for class names, taking care of *not* finding
					# names where the pattern is only a substring of the name in
					# the file so if we want to substitute OriginalClass with
					# ReplacementClass, we do not want to substitute MyOwnOriginalClass
					# or OriginalClassExtended.
					classname_pattern = '(^|[^A-Za-z0-9_])' + classname + '([^A-Za-z0-9_]|$)'
					if re.search(classname_pattern, line):
						subst = args.class_substitutions[classname]
						new_classname = subst[subst.rfind('.')+1:]
						found_classes[classname] = new_classname
		else:
			import_string = m.group(1)
			prefix_len = prefixLength(import_string, args)
			if prefix_len > 0:
				subpackage = import_string[prefix_len+1:import_string.rfind('.')]
				imported_subpackages[subpackage] = True

	if len(found_classes) == 0:
		return lines

	# now re-process all source lines again:
	#  * add import statements for all used classes
	#  * remove other import statements that contain
	#    either a from or to prefix
	#  * other lines are just appended
	modified = []
	first_import = True
	for line in lines:
		if re.match(import_pattern, line) != None:
			if first_import:
				modified.extend(addImports(found_classes, args))
				first_import = False
			
			if not containsAnyFromPrefix(line, args) and not args.to_prefix in line:
				modified.append(line)
		elif re.match(static_import_pattern, line) != None:
			if containsAnyFromPrefix(line, args) and not args.to_prefix in line:
				for prefix in args.from_prefix:
					if prefix in line:
						line = re.sub(prefix, args.to_prefix, line)
				modified.append(line)
			else:
				modified.append(line)
		else:
			for classname in found_classes:
				new_classname = found_classes[classname]
				if classname in line:
					classname_pattern = r'(^|[^A-Za-z0-9_])' + classname + r'([^A-Za-z0-9_]|$)'
					line = re.sub(classname_pattern, r'\1' + new_classname + r'\2', line)
			modified.append(line)
	
	return modified

def addImports(found_classes, args):
	imports = []
	for classname in found_classes:
		imports.append('import ' + args.class_substitutions[classname] + ';\n')
	
	imports.sort()
	return imports

def parseClassSubstitutionRules(args):
	with open(args.classes_subst, 'r') as sources:
		lines = sources.readlines()

	args.classnames = {}
	args.class_substitutions = {}
	for line in lines:
		tokens = line.split()
		oldname = tokens[0]
		newname = tokens[1]
		
		index = oldname.rfind('.')
		classname = oldname[index+1:]
		subpackage = oldname[len('${fromprefix}')+1:oldname.rfind('.')]
		args.classnames[classname] = subpackage

		newname = newname.replace('${toprefix}', args.to_prefix)
		args.class_substitutions[classname] = newname

def containsAnyFromPrefix(line, args):
	for prefix in args.from_prefix:
		if prefix in line:
			return True
	
	return False

def prefixLength(str, args):
	for prefix in args.from_prefix:
		if prefix in str:
			return len(prefix)

	return -1

# start main

parser = argparse.ArgumentParser(prog='migrate')
parser.add_argument('--dir', action='append', required=True, type=lambda x: is_valid_directory(parser, x))
parser.add_argument('--ignore', action='append', default=[])
parser.add_argument('--ext', action='append', required=True)
parser.add_argument('--nosave', action='store_false', dest='save', default=True)
parser.add_argument('--from-prefix', action='append', dest='from_prefix', default=['org.apache.commons.math3', 'org.apache.commons.math4'])
parser.add_argument('--to-prefix', dest='to_prefix', default='org.hipparchus')
parser.add_argument('--classes-subst', default=None)
parser.add_argument('--dry-run', action='store_true', default=False)
parser.add_argument('--verbose', '-v', action='store_true', dest='verbose', default=False)

args = parser.parse_args()

tempdir = tempfile.mkdtemp()
atexit.register(lambda dir=tempdir: shutil.rmtree(dir))
args.tempdir = tempdir

if args.classes_subst == None:
	args.classes_subst = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'classes-update-deprecated-exceptions.subst')

parseClassSubstitutionRules(args)

# process the files
for dir in args.dir:
	log('Processing files in dir ' + dir, args)
	processDir(dir, args)

