# Makefile for compiling Frag1D
export BUILDSEQMAP := ${PWD}/..
export LD_LIBRARY_PATH := ../lib
CP = /bin/cp -f 

all:
	make -f myfunc.makefile
	make install -f myfunc.makefile
	make -f createseqmap.makefile
	make -f getseqresseq.makefile
	make -f shapestring.makefile
	make -f splitfasta.makefile
	make -f getpdbfilepath.makefile
	make -f getpdbaafilepath.makefile
	make -f getdsspfilepath.makefile
install:
	make install  -f createseqmap.makefile
	make install  -f getseqresseq.makefile
	make install  -f shapestring.makefile
	make install  -f splitfasta.makefile
	make install -f getpdbfilepath.makefile
	make install -f getpdbaafilepath.makefile
	make install -f getdsspfilepath.makefile
	cp -uf classify-seqid.awk  ../bin/
	cp -uf id2pdbid.awk        ../bin/
	cp -uf scopid2stdid.awk    ../bin/
	cp -uf createdssp.sh       ../bin/
	cp -uf createpdbaa.sh      ../bin/
	cp -uf cShapeStringList.sh ../bin/
clean:
	make clean -f myfunc.makefile
	make clean -f createseqmap.makefile
	make clean -f getseqresseq.makefile
	make clean -f shapestring.makefile
	make clean -f splitfasta.makefile
	make clean -f getpdbfilepath.makefile
	make clean -f getpdbaafilepath.makefile
	make clean -f getdsspfilepath.makefile
