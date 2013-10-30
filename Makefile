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

INCLUDE = -I./include
CC=g++

all: RNAprofile

RNAprofile: main.o hashtbl.o Set.o profile.o helix_class.o options.o graph.o memoryDFS.o Profnode.o
	$(CC) -o RNAprofile -Wall -g $(INCLUDE) hashtbl.o main.o Set.o profile.o helix_class.o options.o graph.o memoryDFS.o profnode.o -lm  -lgomp libgtfold.a -lgmp

main.o: main.c Set.h hashtbl.h Options.h graph.h memoryDFS.h ./include/boltzmann_main.h
	$(CC) -o main.o -Wall -g $(INCLUDE) -c main.c

hashtbl.o: hashtbl.c hashtbl.h
	$(CC) -o hashtbl.o -Wall -g -c hashtbl.c

Set.o: Set.c Set.h hashtbl.h helix_class.h Profile.h Options.h graph.h Profnode.h
	$(CC) -o Set.o -Wall -g -c -lm Set.c

profile.o: Profile.c Profile.h hashtbl.h
	$(CC) -o profile.o -Wall -g -lm -c Profile.c

helix_class.o: helix_class.c helix_class.h 
	$(CC) -o helix_class.o -Wall -g -c helix_class.c

options.o: Options.c Options.h
	$(CC) -o options.o -Wall -g -c Options.c

graph.o: graph.c graph.h Set.h
	$(CC) -o graph.o -Wall -g -c graph.c

memoryDFS.o: memoryDFS.c memoryDFS.h graph.h hashtbl.h
	$(CC) -o memoryDFS.o -Wall -g -c memoryDFS.c

Profnode.o: Profnode.c Profnode.h
	$(CC) -o profnode.o -Wall -g -c Profnode.c

clean:
	rm -f *.o main RNAprofile
