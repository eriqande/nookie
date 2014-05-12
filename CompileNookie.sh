# simple script to compile nookie
gcc -o nookie \
	src/nookie.c \
	eca-shared/ecalibs/ECA_Opt3.c \
	eca-shared/ecalibs/ECA_MemAlloc.c \
	eca-shared/ecalibs/MathStatRand.c \
	eca-shared/ranlib/src/ranlib.c \
	eca-shared/ranlib/src/com.c \
	eca-shared/ranlib/linpack/linpack.c \
	-Ieca-shared/ecalibs -Ieca-shared/ranlib/src
