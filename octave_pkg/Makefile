## Copyright 2015-2016 Mike Miller
## Copyright 2015-2016 Carnë Draug
## Copyright 2015-2016 Oliver Heimlich
## Copyright 2017 Juan Pablo Carbajal
##
## Copying and distribution of this file, with or without modification,
## are permitted in any medium without royalty provided the copyright
## notice and this notice are preserved.  This file is offered as-is,
## without any warranty.

## Makefile to simplify Octave Forge package maintenance tasks

## Some shell programs
MD5SUM    ?= md5sum
SED       ?= sed
GREP      ?= grep
TAR       ?= tar

## Helper function
TOLOWER   := $(SED) -e 'y/ABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyz/'

### Note the use of ':=' (immediate set) and not just '=' (lazy set).
### http://stackoverflow.com/a/448939/1609556
PACKAGE := $(shell $(SED) -n -e 's/^Name: *\(\w\+\)/\1/p' DESCRIPTION | $(TOLOWER))
VERSION := $(shell $(SED) -n -e 's/^Version: *\(\w\+\)/\1/p' DESCRIPTION | $(TOLOWER))
DEPENDS := $(shell $(SED) -n -e 's/^Depends[^,]*, \(.*\)/\1/p' DESCRIPTION | $(SED) 's/ *([^()]*),*/ /g')

## This are the files that will be created for the releases.
TARGET_DIR      := OF
RELEASE_DIR     := $(TARGET_DIR)/$(PACKAGE)-$(VERSION)
RELEASE_TARBALL := $(TARGET_DIR)/$(PACKAGE)-$(VERSION).tar.gz
HTML_DIR        := $(TARGET_DIR)/$(PACKAGE)-html
HTML_TARBALL    := $(TARGET_DIR)/$(PACKAGE)-html.tar.gz

## Octave binaries
MKOCTFILE ?= mkoctfile
OCTAVE    ?= octave --no-gui

## Remove if not needed, most packages do not have PKG_ADD directives.
M_SOURCES   := $(wildcard inst/*.m)

## Targets that are not filenames.
## https://www.gnu.org/software/make/manual/html_node/Phony-Targets.html
.PHONY: help dist html release install all check run clean

## make will display the command before runnning them.  Use @command
## to not display it (makes specially sense for echo).
help:
	@echo "Targets:"
	@echo "   dist             - Create $(RELEASE_TARBALL) for release"
	@echo "   html             - Create $(HTML_TARBALL) for release"
	@echo "   release          - Create both of the above and show md5sums"
	@echo
	@echo "   install          - Install the package in GNU Octave"
	@echo "   all              - Build all oct files"
	@echo "   run              - Run Octave with development in PATH (no install)"
	@echo "   check            - Execute package tests (w/o install)"
	@echo "   doctest          - Tests only the help text via the doctest package"
	@echo
	@echo "   clean            - Remove releases, html documentation, and oct files"

# dist and html targets are only PHONY/alias targets to the release
# and html tarballs.
dist: $(RELEASE_TARBALL)
html: $(HTML_TARBALL)

# An implicit rule with a recipe to build the tarballs correctly.
%.tar.gz: %
	tar -c -f - --posix -C "$(TARGET_DIR)/" "$(notdir $<)" | gzip -9n > "$@"

# Some packages are distributed outside Octave Forge in non tar.gz format.
%.zip: %
	cd "$(TARGET_DIR)" ; zip -9qr - "$(notdir $<)" > "$(notdir $@)"

# Create the unpacked package.
#
# Notes:
#    * having this recipe separate from the one that makes the tarball
#      makes it easy to have packages in alternative formats (such as zip)
#    * note that if a commands needs to be ran in a specific directory,
#      the command to "cd" needs to be on the same line.  Each line restores
#      the original working directory.
$(RELEASE_DIR):
	@echo "Creating package version $(VERSION) release ..."
	mkdir -p $@
	rsync -a -q --exclude='Makefile' * $@

# install is a prerequesite to the html directory (note that the html
# tarball will use the implicit rule for ".tar.gz" files).
$(HTML_DIR): install
	@echo "Generating HTML documentation. This may take a while ..."
	$(RM) -r "$@"
	$(OCTAVE) --no-window-system --silent \
	  --eval "pkg load generate_html; " \
	  --eval "pkg load $(PACKAGE);" \
	  --eval 'generate_package_html ("${PACKAGE}", "$@", "octave-forge");'
	chmod -R a+rX,u+w,go-w $@

# To make a release, build the distribution and html tarballs.
release: dist html
	@$(MD5SUM) $(RELEASE_TARBALL) $(HTML_TARBALL)
	@echo "Upload @ https://sourceforge.net/p/octave/package-releases/new/"

install: $(RELEASE_TARBALL)
	@echo "Installing package locally ..."
	$(OCTAVE) --silent --eval 'pkg ("install", "-verbose", "$(RELEASE_TARBALL)")'

clean:
	$(RM) -r $(RELEASE_DIR) $(RELEASE_TARBALL) $(HTML_TARBALL) $(HTML_DIR)

#
# Recipes for testing purposes
#

# Build any requires oct files.  Some packages may not need this at all.
# Other packages may require a configure file to be created and run first.
all:
	@echo "all: Nothing to be done"

# Start an Octave session with the package directories on the path for
# interactice test of development sources.
run: all
	$(OCTAVE) --silent --persist --path "inst/" --path "src/" \
	  --eval 'if(!isempty("$(DEPENDS)")); pkg load $(DEPENDS); endif;' \
	  --eval '$(PKG_ADD)'

# Test example blocks in the documentation.  Needs doctest package
#  http://octave.sourceforge.net/doctest/index.html
doctest: all
	$(OCTAVE) --path "inst/" --path "src/" \
	  --eval '${PKG_ADD}' \
	  --eval 'pkg load doctest;' \
	  --eval "targets = '$(shell (ls inst; ls src | grep .oct) | cut -f2 -d@ | cut -f1 -d.)';" \
	  --eval "targets = strsplit (targets, ' ');" \
	  --eval "doctest (targets);"

# Note "doctest" as prerequesite.  When testing the package, also check
# the documentation.
check: all doctest
	@echo "No recipies to run tests"
