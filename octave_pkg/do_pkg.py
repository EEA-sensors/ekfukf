'''
 Copyright (C) 2017 - Juan Pablo Carbajal

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
'''

# Author: Juan Pablo Carbajal <ajuanpi+dev@gmail.com>

import os, sys, getopt
import shutil
import fnmatch, re

WD = os.getcwd()

PKG_STRUCT = {
              'folders':
                        ['inst',
                        'src'],
              'files':
                        ['COPYING',
                        'DESCRIPTION',
                        'INDEX',
                        'Makefile',
                        'NEWS',
                        'CITATION']
              }


def get_files_byext (files, ext):
  ptr = '.*\.{}$'.format(ext)
  return [f for f in files if re.match(ptr, f)]

def analyse_source (source):
  content = os.walk (source)
  mfiles = dict()
  for path, dirn, files in content:
      if '.' not in path:

          folder = '.'
          if path != source:
            folder = ''
            while path != source:
              path, tmp = os.path.split(path)
              folder = os.path.join (tmp,folder)

          tmp = get_files_byext (files, 'm')
          if tmp:
            mfiles[folder] = tmp

  return mfiles

def do_empty_pkg (dest):
  os.mkdir (dest)
  for f in PKG_STRUCT['folders']:
    p = os.path.join (dest, f)


def put_mfiles (mfiles, source, dest):
  print ('Copying m-files to {}'.format(dest))

  for dirn,files in mfiles.items():
    if dirn == '.':
      dirn = ''

    out_dir  = os.path.join (dest, 'inst', dirn)
    in_files = [os.path.join (source, dirn, f) for f in files]

    if not os.path.exists (out_dir):
       os.makedirs (out_dir)

    [shutil.copy(f, out_dir) for f in in_files]
  return

def copy_done_file (in_file, dest):
  if os.path.isfile (in_file):
    shutil.copy(in_file, dest)
  else:
    print ('No {} to copy over'.format(in_file))
  return

def do_pkg (argv):
  help_str = 'do_pkg.py -i <input_folder> -o <output_folder>'

  # Parse input arguments
  SOURCE = WD[:]          # Source folder
  DEST   = WD[:] + '_of'  # Destination folder
  try:
    opts, args = getopt.getopt(argv,"hi:o:",["idir=","odir="])
  except getopt.GetoptError:
    print (help_str)
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
       print (help_str)
       sys.exit()
    elif opt in ("-i", "--idir"):
       SOURCE = arg
    elif opt in ("-o", "--odir"):
       DEST = arg
  print ('Input folder is ', SOURCE)
  print ('Output folder is ', DEST)
  ###### end parse arguments

  do_empty_pkg (DEST)

  mfiles = analyse_source (SOURCE)
  put_mfiles (mfiles, SOURCE, DEST)

  # Pre-generated files
  for f in PKG_STRUCT['files']:
    copy_done_file (f, DEST)

if __name__ == '__main__':
  do_pkg(sys.argv[1:]);
