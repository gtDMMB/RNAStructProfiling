# Copyright (c) 2011 the authors listed at the following URL, and/or
# the authors of referenced articles or incorporated external code:
# http://en.literateprograms.org/Hash_table_(C)?action=history&offset=20100620072342
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 

## If the user supplies an external libgtfold.a path, use it here:
LIBGTFOLD=libgtfold.a
EXT_LIBGTFOLD=$(shell env | grep LIBGTFOLD | sed -e 's/^.*=//')
ifneq "$(EXT_LIBGTFOLD)" ""
	LIBGTFOLD=$(EXT_LIBGTFOLD)
endif


CC=g++
CFLAGS=-Wall -Wextra -g -I./include -std=c++11
STATIC_CFLAGS=-Wall -Wextra -g -I./include -std=c++11
LD=gcc
LDFLAGS=-Wall -Wextra -g -lm -lgomp -lgmp
STATIC_LDFLAGS=-Wall -Wextra -g\
	-Wl,-Bdynamic -lm -lgomp -lgmp\
	-Wl,--verbose -Wl,--as-needed

# Define special CFLAGS and LDFLAGS if we intend on having a 
# statically linked binary at the end of the make process:
BUILD_STATIC_LOCAL=$(shell env | grep STATIC_LINKAGE_ONLY | sed -e 's/^.*=//')
ifneq "$(BUILD_STATIC_LOCAL)" ""
	CFLAGS=$(STATIC_CFLAGS)
	LDFLAGS=$(STATIC_LDFLAGS)
endif

EXE=RNAprofile
ALLOBJS=hashtbl.o Options.o Profile.o helix_class.o\
	graph.o memoryDFS.o Profnode.o Set.o main.o

all: RNAprofile

RNAprofile: $(ALLOBJS)
	$(LD) -o $(EXE) $(LDFLAGS) $(ALLOBJS) $(LIBGTFOLD)

main.o: main.c Set.h hashtbl.h Options.h graph.h memoryDFS.h\
	./include/boltzmann_main.h
	$(CC) $(CFLAGS) -c main.c

hashtbl.o: hashtbl.c hashtbl.h
	$(CC) $(CFLAGS) -c hashtbl.c

Set.o: Set.c Set.h hashtbl.h helix_class.h Profile.h Options.h\
	graph.h Profnode.h
	$(CC) $(CFLAGS) -c Set.c

Profile.o: Profile.c Profile.h hashtbl.h
	$(CC) $(CFLAGS) -c Profile.c

helix_class.o: helix_class.c helix_class.h 
	$(CC) $(CFLAGS) -c helix_class.c

Options.o: Options.c Options.h
	$(CC) $(CFLAGS) -c Options.c

graph.o: graph.c graph.h Set.h
	$(CC) $(CFLAGS) -c graph.c

memoryDFS.o: memoryDFS.c memoryDFS.h graph.h hashtbl.h
	$(CC) $(CFLAGS) -c memoryDFS.c

Profnode.o: Profnode.c Profnode.h
	$(CC) $(CFLAGS) -c Profnode.c

clean:
	rm -f *.o *.code $(EXE)
