#!/usr/bin/env python

import argparse
import os
import tempfile
import filecmp
import shutil
import re
import atexit

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
			index = subdirList.index(ignore)
			if index >= 0:
				del subdirList[index]


def processFile(file, args):
	log("processing file " + file, args, False)

	g = os.path.basename(file)
	tmpfile = os.path.join(args.tempdir, g)
	with open(file, "r") as sources:
		lines = sources.readlines()
	with open(tmpfile, "w") as sources:
		for line in lines:
			line = packageReplace(line, args)
			line = classesReplace(line, args)
			sources.write(line)

	if filecmp.cmp(file, tmpfile) == False:
		log("-> changed", args)

		if args.dry_run == True:
			if args.verbose == False:
				print file
			os.remove(tmpfile)
		else:
			if args.save == True:
				shutil.move(file, file + ".orig")
			shutil.move(tmpfile, file)
	else:
		log("-> unchanged", args)

def packageReplace(line, args):
	for p in args.packages_subst_pattern:
		line = re.sub(p[1], p[2], line)
	return line

def parsePackageSubstitutionRules(args):
	with open(args.packages_subst, "r") as sources:
		lines = sources.readlines()

	for line in lines:
		pattern = re.sub(r'^\s*"([^"]*)"\s*"([^"]*)"\s*', r's:\1:\2', line)
		if pattern.startswith('s:'):
			args.packages_subst_pattern.append(pattern.split(':'))

def classesReplace(line, args):
	for p in args.classes_subst_pattern:
		line = re.sub(p[1], p[2], line)
	return line

def parseClassSubstitutionRules(args):
	with open(args.classes_subst, "r") as sources:
		lines = sources.readlines()

	# substitution rules for class names, taking care of *not* substituting
	# names where the pattern is only a substring of the name in the file
	# so if we want to substitute OriginalClass with ReplacementClass, we
	# do not want to substitute MyOwnOriginalClass or OriginalClassExtended
	for line in lines:
		pattern = re.sub(r'^\s*([A-Za-z_][A-Za-z0-9_]*)\s*([A-Za-z_][A-Za-z0-9_]*)\s*$', r's:(^|[^A-Za-z0-9_])\1([^A-Za-z0-9_]|$):\\1\2\\2', line)
		if pattern.startswith('s:'):
			args.classes_subst_pattern.append(pattern.split(':'))


# start main

parser = argparse.ArgumentParser(prog='migrate')
parser.add_argument('--dir', action='append', required=True, type=lambda x: is_valid_directory(parser, x))
parser.add_argument('--ignore', action='append', default=[])
parser.add_argument('--ext', action='append', required=True)
parser.add_argument('--nosave', action='store_false', dest='save', default=True)
parser.add_argument('--packages-subst', default=None)
parser.add_argument('--classes-subst', default=None, required=True)
parser.add_argument('--dry-run', action='store_true', default=False)
parser.add_argument('--verbose', '-v', action='store_true', dest='verbose', default=False)

args = parser.parse_args()

tempdir = tempfile.mkdtemp()
atexit.register(lambda dir=tempdir: shutil.rmtree(dir))
args.tempdir = tempdir

args.packages_subst_pattern = []
if args.packages_subst == None:
	print 'WARNING: package substitutions disabled at user request'
else:
	parsePackageSubstitutionRules(args)

args.classes_subst_pattern = []
parseClassSubstitutionRules(args)

#print args

# process the files
for dir in args.dir:
	log("Processing files in dir " + dir, args)
	processDir(dir, args)

