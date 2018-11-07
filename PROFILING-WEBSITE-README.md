# Steps to re-creating a sane executable for the gtDMMB profiling website

## Introduction

Note that we have built the new files on a *CentOS 7.5 / x86_64* box. We have deduced by 
upgrade timelines and when the functionality on the website began to disappear 
that this is the best possible platform ([see here]()): 
```
$ echo $(cat /etc/centos-release) "[$(uname -m)]"
CentOS Linux release 7.5.1804 (Core) [x86_64]
```
If the architecture of the local machine changes and this process needs to be repeated, it 
may be necessary to *cross compile* to obtain an executable which is consistent with the 
target (GA Tech webhosting server) architecture. 
However, since the resulting replacement libraries and executables are designed to be 
*statically linked*, i.e., the resulting execution of these binary files does not depend on 
linker references to *dynamic* libraries (which, as we have seen, can change over time), 
it **shouldn't** be necessary to repeat the steps in this README. It is nonetheless being 
created as user documentation to accompany the project. 

We also need to install the following CentOS packages with ``yum``:
```
$ sudo yum install glibc-static libstdc++-static
$ sudo yum install gmp-devel
$ sudo yum install glibc-devel.i686 libgcc.i686 libstdc++-devel.i686 ncurses-devel.i686
```
This will allow us to statically link standard C++ library functions, including the math
library with ``-lm``, into our profiling code binary. Again, the purpose of this is that if 
the dynamic libraries change due to a re-install of the server platform, as long as the 
architecture (i.e., machine type of ``uname -m``) stay the same our code should still work 
after the re-install is complete -- even, say, 5-6 years down the line. 

### Other possibilities to look into

Another solid option to look into if this process is necessary in the future is the 
[ELF statifier](http://statifier.sourceforge.net/) project. This utility takes an executable 
with possibly dynamically linked libraries and produces one *large(r)* executable which has 
no dynamic library dependencies. There are other good utilities like this, e.g., 
[Ermine](http://www.magicermine.com/), but to my knowledge these are not freeware applications. 

## Compiling GTFold as a statically linked library (*.a file)

Run the following commands:
```
$ git clone https://github.gatech.edu/akirkpatrick3/gtfold-internal.git
$ cd gtfold-internal/gtfold-mfe
$ ./configure --enable-64bit --enable-debug CFLAGS="-fPIC" CXXFLAGS="-fPIC" && make
$ cd src
$ export OBJFILES=$(ls *.o | tr "\n" " " | sed -e 's/ main.o//g')
$ g++ -o libgtfold.a -shared -Wl,-Bdynamic -lm -lgomp -lgmp -Wl,--verbose $OBJFILES -lc
$ export LIBGTFOLD=$(readlink -f ./libgtfold.a)
$ cd ~
```

## Compiling the <RNA Profiling Software> as a standalone (statically-linked) executable

First, setup the ``Makefile`` check to see if we need to compile with static linkage of 
libraries only (*we DO want this here*):
``
$ export STATIC_LINKAGE_ONLY=yes
$ cp $LIBGTFOLD ./
$ unset LIBGTFOLD
$ git clone https://github.com/gtDMMB/RNAStructProfiling.git
$ cd RNAStructProfiling
$ make
$ ldd RNAprofile
$ cp RNAprofile libgtfold.a /lib64/libm.so* /lib64/libgomp.so.1 /lib64/libgmp.so.10 /lib64/libstdc++.so.6 /lib64/libgcc_s.so.1 ProfilingWebsiteBinary/
```

## Recompiling the dot utility (from graphviz)

```
$ wget https://graphviz.gitlab.io/pub/graphviz/stable/SOURCES/graphviz.tar.gz
$ tar xvzf graphviz.tar.gz
$ cd graphviz*
$ mkdir -p /private
$ ./configure --enable-shared --with-rsvg=yes --prefix=/private CFLAGS="-fPIC" CXXFLAGS="-fPIC"
$ make && sudo make install
#$ sudo yum install graphviz.x86_64
$ cp /private/bin/dot ~/RNAStructProfiling/ProfilingWebsiteBinary/
$ ldd /private/bin/dot
$ cp /private/lib/libgvc.so.6 /lib64/libdl.so.2 /lib64/libltdl.so.7 /private/lib/libxdot.so.4 /private/lib/libcgraph.so.6 /private/lib/libpathplan.so.4 /lib64/libexpat.so.1 /lib64/libz.so.1 /lib64/libm.so.6 /private/lib/libcdt.so.5 /lib64/libc.so.6 /lib64/ld-linux-x86-64.so.2 ~/RNAStructProfiling/ProfilingWebsiteBinary/
$ mkdir -p ~/RNAStructProfiling/ProfilingWebsiteBinary/lib/graphviz/
$ cp -r /private/lib/graphviz/config6 ~/RNAStructProfiling/ProfilingWebsiteBinary/lib/graphviz/
```

## Installing the newly built utilities

### Backups (if necessary):

First, ``ssh`` into the rnaprofiling GA Tech ip address. Once you have a shell, we need to backup 
the existing compiled utilities for posterity's sake before proceeding to replace them:
```
$ ssh ourusername@OurServersIP
$ cp -rpa private private-dist-$(date +"%F-%H%M%S")
$ exit
```

### Copy over local files

```
$ scp ProfilingWebsiteBinary/* ourusername@OurServersIP:private/
$ ssh ourusername@OurServersIP
$ export LD_LIBRARY_PATH=/private
$ cd private
$ mv ../httpdocs/dot ../private-*/
$ ./RNAprofile
```

## Demonstration of point

Now head over to the [profiling website](http://rnaprofiling.gatech.edu/), 
enter some ascii string filled with "a/t/c/u/g"'s and observe that the previously 
existing errors have disappeared. :smile: :check_mark: :exclamation_point:
