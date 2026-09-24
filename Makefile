# kycg — functional analysis of DNA methylation at CpG resolution.
#
# Links against YAME's libyame.a (git submodule under external/YAME, with
# htslib vendored inside it rather than as a submodule of its own), which
# carries the bit-packed CpG formats and the bit-parallel `summary` counting.
# kycg adds statistical inference on top; YAME stays a data engine.
#
# The submodule must be initialized before building:
#   git submodule update --init --recursive
# `make` drives YAME's own `make lib` target, so no manual step is needed.

# ?= does not work here: make predefines CC, so its origin is `default` rather
# than `undefined` and ?= never fires. Override only that case, leaving an
# explicit CC from the environment or the command line (which is how conda
# builds set it) untouched.
ifeq ($(origin CC),default)
  CC = cc
endif
# _GNU_SOURCE makes glibc declare strcasestr, which sits behind __USE_GNU in
# <string.h> under every glibc. Without it the call is an implicit declaration
# -- a warning on gcc <= 13 but a hard error from gcc 14 on. kycg's own last
# caller left with src/ui.c; the flag stays because YAME's headers are compiled
# into these objects and the pin below is the reason it is safe at all.
# Safe here only because conda-recipe/conda_build_config.yaml pins the build to
# the glibc 2.17 sysroot: from glibc 2.38 on, _GNU_SOURCE also switches on the
# C23 strtol redirects (strtol -> __isoc23_strtol@GLIBC_2.38), which would make
# the binary unrunnable on RHEL 9 / Ubuntu 22.04. Do not drop that pin.
CFLAGS = -W -Wall -finline-functions -std=gnu99 -Wno-unused-result -O3 -D_GNU_SOURCE
CLIB = -lpthread -lz -lm

PREFIX ?= /usr/local

OS := $(shell uname)
ifeq ($(OS), Darwin)
	CFLAGS += -Wno-unused-function
else
	CLIB += -lrt
endif

PROG = kycg

SRC_DIR = src
TEST_DIR = tests

# YAME submodule: we consume its static library and its headers directly.
# There is no installed-header split upstream, so we point at src/.
YAME_DIR = external/YAME
YAME_LIB = $(YAME_DIR)/libyame.a
YAME_HTSLIB = $(YAME_DIR)/htslib/libhts.a

CFLAGS += -I$(SRC_DIR) -I$(YAME_DIR)/src -I$(YAME_DIR)/htslib

# libcurl is no longer kycg's concern: the fetch/verify engine moved into
# libyame (src/assets.c), and YAME owns the curl link. Everything kycg needs at
# link time -- libyame, htslib, zlib/pthread/curl -- comes from `yame-config
# --libs`, which is a build product of `make -C external/YAME lib`. It is
# expanded with $$(...) in the link recipes so the shell runs it after
# $(YAME_LIB) has been built. Whether `kycg fetch` can reach the network is
# therefore a runtime question (yame_assets_have_curl()), not a compile flag.

SOURCES := $(wildcard $(SRC_DIR)/*.c)
OBJECTS := $(SOURCES:$(SRC_DIR)/%.c=$(SRC_DIR)/%.o)
# Everything except the CLI dispatcher, so tests can link the library half.
LIBOBJECTS = $(filter-out $(SRC_DIR)/main.o, $(OBJECTS))

.PHONY: all build debug clean distclean test test-bins test-docs coverage registry yame yame-bin install

all: build

build: $(PROG)

# `clean` first, deliberately: target-specific variables do not participate in
# up-to-date checks, so `make && make debug` would otherwise recompile nothing
# and leave the optimized, unsanitized binary in place -- reporting success
# while giving you no sanitizer at all.
debug: CFLAGS += -g -fsanitize=address,undefined
debug: CFLAGS := $(filter-out -O3,$(CFLAGS))
debug: clean build

#####################
##### libraries #####
#####################

# Delegate to YAME's own `lib` target; it builds htslib as a prerequisite.
yame:
	$(MAKE) -C $(YAME_DIR) lib

# The yame COMMAND, which `lib` does not build. kycg has no `info` subcommand:
# checking that a query and a knowledgebase index the same row space is
# `yame info`, and the docs tell users to run it, so a source tree has to be
# able to produce it. Packaged users get it from the separate `yame` conda
# package; this target is the source-build equivalent. Not a prerequisite of
# `build` -- kycg links the library, not the binary.
yame-bin:
	$(MAKE) -C $(YAME_DIR) build
	@echo "built $(YAME_DIR)/yame -- add it to PATH, or install the yame conda package"

$(YAME_LIB) $(YAME_HTSLIB): yame

###################
### compilation ###
###################

# Objects depend on every kycg header, not just their own .c. Without this a
# regenerated registry.h or kbinfo.h leaves a stale binary that reports the
# previous contents -- which looks like a data bug and is not one.
# YAME's headers count too: src/test.c and src/info.c include cdata.h, cfile.h,
# summary.h and wzmisc.h from the submodule. Without them, bumping the
# submodule to a YAME whose cdata_t layout changed relinks against the new
# libyame.a while keeping stale objects compiled for the old struct -- the link
# succeeds and the offsets disagree at run time.
KYCG_HEADERS := $(wildcard $(SRC_DIR)/*.h) $(wildcard $(YAME_DIR)/src/*.h)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(KYCG_HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@

###################
###   linking   ###
###################

# $(LDFLAGS) is expanded explicitly: this is an explicit recipe, so no implicit
# rule supplies it. conda-forge puts -L$$PREFIX/lib ONLY in LDFLAGS, and without
# it a linux-64 conda build fails to find -lz -- libcurl resolves by luck,
# because curl-config emits its own -L.
$(PROG): $(YAME_LIB) $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $$($(YAME_DIR)/yame-config --libs)

###################
###    tests    ###
###################

# The statistics are pure functions over counts, and argv permutation is a
# pure function over argv. store.o, however, is now a thin shim over libyame's
# yame_assets_* primitives, so anything linking it must also pull in libyame
# (and the htslib/zlib/curl behind it) -- hence the $(YAME_LIB) prereq and the
# `yame-config --libs` on the link line.
TEST_SRC := $(wildcard $(TEST_DIR)/*.c)
TEST_BIN := $(TEST_SRC:$(TEST_DIR)/%.c=$(TEST_DIR)/%)

TEST_OBJ = $(SRC_DIR)/hypergeo.o $(SRC_DIR)/enrich.o $(SRC_DIR)/args.o \
           $(SRC_DIR)/store.o $(SRC_DIR)/digest.o

$(TEST_DIR)/%: $(TEST_DIR)/%.c $(TEST_OBJ) $(YAME_LIB)
	$(CC) $(CFLAGS) -o $@ $< $(TEST_OBJ) $$($(YAME_DIR)/yame-config --libs)

# The whole suite: the C unit binaries above, then the command-level tests in
# tests/t_*.sh. `yame-bin` is a prerequisite because those tests pack their
# own fixtures with `yame pack` -- kycg reads what yame writes, so a committed
# fixture would quietly become a format the current YAME no longer produces.
test: test-bins yame-bin
	bash $(TEST_DIR)/run.sh

# The documented-workflow gate: every runnable block on the page, as a reader
# would run it, against THIS checkout's binaries. Deliberately NOT part of
# `make test` and never run in CI or a conda build -- it needs the network, a
# git clone and a few minutes. The sandbox persists between runs, so a second
# run downloads nothing. On the HPC run it under sbatch (release SOP step 4b).
#
#   make test-docs
#   KYCG_DOCS_SANDBOX=~/tmp/kycg/docs_gate make test-docs
#   KYCG_DOCS_BIN=/path/to/env/bin make test-docs     # test an installed copy
test-docs: build yame-bin
	python3 $(TEST_DIR)/docs_gate.py

# Just the unit binaries. `make` builds the CLI only, so a coverage run that
# wants the unit tests too has to ask for them by name -- without this they
# are simply absent from the scratch tree and run.sh skips them, which reads
# as a suite that passes while measuring nothing of hypergeo.c or digest.c.
test-bins: $(TEST_BIN)

###################
###  generated  ###
###################

# The two generated headers, projected from YAME's shared file table in the
# pinned submodule: src/registry.h (per-file urls, digests, sizes and row
# counts) and src/kbinfo.h (what each set is, for the browser's info pane).
#
# Regenerate with every submodule bump -- that bump IS a catalog bump, and a
# header left behind describes the previous one. Both are committed, and
# `make test` runs each generator with --check, so a forgotten regeneration
# fails the suite rather than shipping quietly. Emission is offline and
# locale-pinned, so it reproduces byte for byte wherever it runs.
registry:
	./tools/make_registry.sh -o src/registry.h
	./tools/make_kbinfo.sh -o src/kbinfo.h

# Line coverage of the suite over src/, in a scratch copy: the instrumented
# binary is built -O0 and drops .gcda beside itself, so it must never become
# the one in the repo. Rewrites docs/coverage.json, which CI checks.
coverage: build yame-bin
	scripts/coverage.sh

###################
###   install   ###
###################

# Install the binary into $(PREFIX)/bin, honouring $(DESTDIR). Used by the
# conda recipe; `make install PREFIX=/usr/local` for a manual install.
#
# Only the binary is installed. kycg carries no runtime data files -- the
# registry is compiled in and knowledgebases are fetched into the user's own
# store -- so there is nothing else to place.
install: $(PROG)
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp -f $(PROG) $(DESTDIR)$(PREFIX)/bin/$(PROG)
	chmod 755 $(DESTDIR)$(PREFIX)/bin/$(PROG)

###################
###    clean    ###
###################

clean:
	rm -f $(SRC_DIR)/*.o $(PROG) $(TEST_BIN)

distclean: clean
	$(MAKE) -C $(YAME_DIR) clean
