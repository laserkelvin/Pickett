CC=gcc
CFLAGS=-g -O3 -Wall
#CFLAGS=-O2 -Wall
EXEQ=spfit spcat calmrg dpfit dpcat calbak
EXEA=${EXEQ} iamcalc moiam stark termval sortn
#next line for atlas blas
#BLASLIB=-lcblas -latlas
#next line for fortran blas and cblas wrappers
#BLASLIB=-lcblas -lblas
#next line for supplied routines
BLASLIB=
ifndef ($(BLASLIB))
	LBLAS=dblas.o
endif
default: ${EXEQ} 
all: ${EXEA}
install:  
	-mv ${EXEQ} /usr/local/bin 
dpfit: calfit.o subfit.o dpi.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
dpcat: calcat.o sortsub.o dpi.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
spfit: calfit.o subfit.o spinv.o spinit.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
spcat: calcat.o sortsub.o spinv.o spinit.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
calmrg: calmrg.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
calbak: calbak.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
termval: termval.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
stark: stark.o splib.a ; gcc -o $@ $^ $(BLASLIB) -lm
moiam: moiam.o ftran.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
iamcalc: iamcalc.o ftran.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
calbak: calbak.o splib.a ; gcc -o $@ $^ $(BLASLIB) -lm
cnvwn: cnvwn.o splib.a ; gcc -o $@ $^ $(BLASLIB) -lm

splib.a: ulib.o cnjj.o slibgcc.o catutil.o lsqfit.o $(LBLAS)
	ar r splib.a $^
	ranlib splib.a

calfit.o:calfit.c calpgm.h
subfit.o:subfit.c calpgm.h
lsqfit.o:lsqfit.c lsqfit.h
calcat.o:calcat.c calpgm.h
sortsub.o: sortsub.c calpgm.h
calmrg.o:calmrg.c calpgm.h
termval.o:termval.c calpgm.h
stark.o:stark.c calpgm.h
iamcalc.o:iamcalc.c calpgm.h
moiam.o:moiam.c calpgm.h
ulib.o:ulib.c calpgm.h
cnjj.o:cnjj.c cnjj.h
slibgcc.o:slibgcc.c calpgm.h
spinv.o:spinv.c calpgm.h spinit.h
spinit.o:spinit.c calpgm.h spinit.h
dpi.o:dpi.c calpgm.h
ftran.o:ftran.c
sortn.o: sortn.c
sortsub.o: sortsub.c
dblas.o: dblas.c
util.o:util.c
calbak.o:calbak.c
cnvwn.o:cnvwn.c
sortn: sortn.o sortsub.o; gcc -o $@ $^
