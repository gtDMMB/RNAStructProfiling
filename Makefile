# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 

all: main

main: hashtbl.o main.o profile.o summaryGraph.o graph.o memoryDFS.o
	cc -o main -Wall -pedantic -g -lm hashtbl.o main.o profile.o summaryGraph.o graph.o memoryDFS.o

hashtbl.o: hashtbl.c hashtbl.h
	cc -o hashtbl.o -Wall -pedantic -g -c hashtbl.c

main.o: main.c hashtbl.h profile.h 
	cc -o main.o -Wall -pedantic -g -c main.c

profile.o: profile.c hashtbl.h profile.h
	cc -o profile.o -Wall -pedantic -g -c profile.c

summaryGraph.o: summaryGraph.c hashtbl.h summaryGraph.h memoryDFS.h memoryDFS.c
	cc -o summaryGraph.o -Wall -pedantic -g -c summaryGraph.c

graph.o : graph.c graph.h
	cc -o graph.o -Wall -pedantic -g -c graph.c

memoryDFS.o : memoryDFS.c memoryDFS.h graph.h
	cc -o memoryDFS.o -Wall -pedantic -g -c memoryDFS.c

clean:
	rm -f *.o main
