# kycg — functional analysis of DNA methylation at CpG resolution.
#
# Links against YAME's libyame.a (git submodule under external/YAME), which
# carries the bit-packed CpG formats and the bit-parallel `summary` counting.
# kycg adds statistical inference on top; see DESIGN.md for the division.
#
# The submodule must be initialized before building:
#   git submodule update --init --recursive
# `make` drives YAME's own `make lib` target, so no manual step is needed.

CC ?= gcc
CFLAGS = -W -Wall -finline-functions -std=gnu99 -Wno-unused-result -O3
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

# libcurl is OPTIONAL and used by `kycg fetch` alone. Everything that analyzes
# data -- test, info, list -- works without it, so a build box with no curl
# headers still produces a fully functional analysis tool; fetch then reports
# that the build has no network support instead of silently missing. Set
# CURL=0 to force it off even where curl-config exists.
CURL ?= auto
ifeq ($(CURL),auto)
  CURL_CFLAGS := $(shell curl-config --cflags 2>/dev/null)
  CURL_LIBS   := $(shell curl-config --libs 2>/dev/null)
else ifeq ($(CURL),0)
  CURL_CFLAGS :=
  CURL_LIBS   :=
else
  CURL_CFLAGS := $(shell curl-config --cflags 2>/dev/null)
  CURL_LIBS   := $(shell curl-config --libs 2>/dev/null)
endif

ifneq ($(CURL_LIBS),)
  CFLAGS += -DKYCG_HAVE_CURL $(CURL_CFLAGS)
  CLIB   += $(CURL_LIBS)
endif

SOURCES := $(wildcard $(SRC_DIR)/*.c)
OBJECTS := $(SOURCES:$(SRC_DIR)/%.c=$(SRC_DIR)/%.o)
# Everything except the CLI dispatcher, so tests can link the library half.
LIBOBJECTS = $(filter-out $(SRC_DIR)/main.o, $(OBJECTS))

.PHONY: all build debug clean test yame install

all: build

build: $(PROG)

debug: CFLAGS += -g -fsanitize=address,undefined
debug: CFLAGS := $(filter-out -O3,$(CFLAGS))
debug: build

#####################
##### libraries #####
#####################

# Delegate to YAME's own `lib` target; it builds htslib as a prerequisite.
yame:
	$(MAKE) -C $(YAME_DIR) lib

$(YAME_LIB) $(YAME_HTSLIB): yame

###################
### compilation ###
###################

# Objects depend on every kycg header, not just their own .c. Without this a
# regenerated registry.h or kbinfo.h leaves a stale binary that reports the
# previous contents -- which looks like a data bug and is not one.
KYCG_HEADERS := $(wildcard $(SRC_DIR)/*.h)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(KYCG_HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@

###################
###   linking   ###
###################

$(PROG): $(YAME_LIB) $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(YAME_LIB) $(YAME_HTSLIB) $(CLIB)

###################
###    tests    ###
###################

# The statistics are pure functions over counts, and argv permutation is a
# pure function over argv, so the tested modules all link without YAME.
TEST_SRC := $(wildcard $(TEST_DIR)/*.c)
TEST_BIN := $(TEST_SRC:$(TEST_DIR)/%.c=$(TEST_DIR)/%)

TEST_OBJ = $(SRC_DIR)/hypergeo.o $(SRC_DIR)/enrich.o $(SRC_DIR)/args.o

$(TEST_DIR)/%: $(TEST_DIR)/%.c $(TEST_OBJ)
	$(CC) $(CFLAGS) -o $@ $< $(TEST_OBJ) -lm

test: $(TEST_BIN)
	@for t in $(TEST_BIN); do echo "== $$t"; ./$$t || exit 1; done

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
