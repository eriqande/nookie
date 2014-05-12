# nookie

A simple C program to simulate genotypes at unlinked
markers from specified matings.



## Terms 

As a work of the United States Government, the source code in nookie.c is in
public domain within the United States. Additionally, we waive
copyright and related rights in the work worldwide through the CC0 1.0
Universal public domain dedication.

the eca-shared libraries were written by Eric Anderson before joining the 
federal workforce, and I think they are under the MIT license.

the ranlib code is under a different license, but I can't find which.

See TERMS.md for more information.

## Disclaimer

I mostly want to say that I wrote this a long time ago in C.  It works, and it reasonable
fast, but it is breathtaking to see how much work it was to write applications in
C back in the day.  All that memory management!  I'll bet most of the functionality
of this program could be written in only about 5% of the code in R.  Oh well!

## Compiling and Using

Here is the general sequence assuming that you have a reasonable C build
system and a good command line:

```
# get it from git
git clone https://github.com/eriqande/nookie.git

# change into the nookie directory
cd nookie

# get the submodules you need
git submodule init
git submodule update

# compile the executable
./CompileNookie.sh

# get the long version of the help for the program:
./nookie --help-full

# or, if you want that formatted more nicely, you can do:
./nookie --help-nroff | nroff -man | less
```

