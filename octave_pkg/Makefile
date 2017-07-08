## Makefile to simplify Octave Forge package maintenance tasks
OCTPKG_DIR = octave_pkg

PACKAGE = $(shell $(SED) -n -e 's/^Name: *\(\w\+\)/\1/p' $(OCTPKG_DIR)/DESCRIPTION | $(TOLOWER))
VERSION = $(shell $(SED) -n -e 's/^Version: *\(\w\+\)/\1/p' $(OCTPKG_DIR)/DESCRIPTION | $(TOLOWER))
DEPENDS = $(shell $(SED) -n -e 's/^Depends[^,]*, \(.*\)/\1/p' $(OCTPKG_DIR)/DESCRIPTION | $(SED) 's/ *([^()]*),*/ /g')

RELEASE_DIR     = $(PACKAGE)-$(VERSION)
RELEASE_TARBALL = $(PACKAGE)-$(VERSION).tar.gz
HTML_DIR        = $(PACKAGE)-html
HTML_TARBALL    = $(PACKAGE)-html.tar.gz

MD5SUM    ?= md5sum
MKOCTFILE ?= mkoctfile
#OCTAVE    ?= octave --no-window-system
OCTAVE    ?= octave --no-gui
SED       ?= sed
TAR       ?= tar

TOLOWER = $(SED) -e 'y/ABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyz/'

M_SOURCES := $(wildcard octave_pkg/*.m)
PKG_ADD   := $(shell grep -Pho '(?<=\#\# PKG_ADD: ).*' $(M_SOURCES))

SUBFOLDERS  = cov mean prior inf lik util doc
SOURCES     = util/solve_chol.c util/make.m util/lbfgsb/*
HG_EXCLUDED = --exclude $(OCTPKG_DIR) --exclude startup.m \
              --exclude test_sparse.m --exclude "util/minimize_*.m" \
              --exclude ".octaverc" --exclude "Copyright"

.PHONY: help dist html release install all check run doc clean maintainer-clean

help:
	@echo "Targets:"
	@echo "   dist             - Create $(RELEASE_TARBALL) for release"
	@echo "   html             - Create $(HTML_TARBALL) for release"
	@echo "   release          - Create both of the above and show md5sums"
	@echo
	@echo "   install          - Install the package in GNU Octave"
	@echo "   all              - Build all oct files"
	@echo "   check            - Execute package tests (w/o install)"
	@echo "   run              - Run Octave with development in PATH (no install)"
	@echo "   doc              - Build Texinfo package manual"
	@echo
	@echo "   clean            - Remove releases, html documentation, and oct files"
	@echo "   maintainer-clean - Additionally remove all generated files"

$(RELEASE_DIR): .hg/dirstate $(OCTPKG_DIR)/index.sh $(OCTPKG_DIR)/news.sh
	@echo "Creating package version $(VERSION) release ..."
	-rm -rf $@
	hg archive --exclude ".hg*" --exclude Makefile $(HG_EXCLUDED) --type files $@
	cp $(OCTPKG_DIR)/DESCRIPTION "$@"
	cp $(OCTPKG_DIR)/CITATION "$@"
	cp $(OCTPKG_DIR)/__gpml_package_register__.m "$@"
	@echo "Creating INDEX file ..."
	$(OCTPKG_DIR)/index.sh "$@" > "$@/INDEX"
	@echo "Creating NEWS file ..."
	$(OCTPKG_DIR)/news.sh > "$@/NEWS"
	@echo "Creating COPYING file ..."
	cp "Copyright" "$@/COPYING"
	@echo "Package structure ..."
	mkdir "$@/inst" && mkdir "$@/src"
	cd "$@" && mv $(SOURCES) "src/" && rm -rf "util/lbfgsb" && \
	                                    mv $(SUBFOLDERS) "inst/" && mv *.m "inst/"
	cp $(OCTPKG_DIR)/Makefile "$@/src"
	chmod -R a+rX,u+w,go-w $@

$(RELEASE_TARBALL): $(RELEASE_DIR)
	$(TAR) cf - --posix $< | gzip -9n > $@
	-rm -rf $<

$(HTML_DIR): install
	@echo "Generating HTML documentation. This may take a while ..."
	-rm -rf $@
	$(OCTAVE) --silent \
	  --eval 'graphics_toolkit ("gnuplot");' \
	  --eval 'pkg load generate_html $(PACKAGE);' \
	  --eval 'generate_package_html ("$(PACKAGE)", "$@", "octave-forge");'
	chmod -R a+rX,u+w,go-w $@

$(HTML_TARBALL): $(HTML_DIR)
	$(TAR) cf - --posix $< | gzip -9n > $@
	-rm -rf $<

dist: $(RELEASE_TARBALL)

html: $(HTML_TARBALL)

release: dist html
	@$(MD5SUM) $(RELEASE_TARBALL) $(HTML_TARBALL)
	@echo "Upload @ https://sourceforge.net/p/octave/package-releases/new/"
	@echo "Execute: hg tag \"$(VERSION)\""

install: $(RELEASE_TARBALL)
	@echo "Installing package locally ..."
	$(OCTAVE) --silent --eval 'pkg install $(RELEASE_TARBALL);'

all: $(RELEASE_DIR)
	cd "$(RELEASE_DIR)/src" && $(MAKE) $@

check: all $(RELEASE_DIR)
	cd $(RELEASE_DIR) && \
	$(OCTAVE) --silent \
	  --eval 'if(!isempty("$(DEPENDS)")); pkg load $(DEPENDS); endif;' \
	  --eval 'addpath (fullfile ([pwd filesep "inst"]));' \
	  --eval '$(PKG_ADD)' \
	  --eval 'cellfun(@(x)runtests (x), __gpml_package_register__);'

run: all $(RELEASE_DIR)
	cd $(RELEASE_DIR) && \
	$(OCTAVE) --silent --persist \
	  --eval 'if(!isempty("$(DEPENDS)")); pkg load $(DEPENDS); endif;' \
	  --eval 'addpath (fullfile ([pwd filesep "inst"]));' \
	  --eval '$(PKG_ADD)'

doc:

clean:
	-rm -rf $(RELEASE_DIR) $(RELEASE_TARBALL) $(HTML_TARBALL) $(HTML_DIR)

maintainer-clean: clean
