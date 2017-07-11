#!/bin/bash

# COPYING
cp ../LICENSE.txt COPYING
# INDEX
python3 do_index.py
# PACKAGE
python3 do_pkg.py -i .. -o ekfukf_of
# remove created files
rm COPYING INDEX
