CC = c++
CFLAGS=-O2 -std=c++14 -g

LOCALPATH=/home/nirmalya/local/

INC = -I${LOCALPATH}/include

PROG_OPT_LIB=${LOCALPATH}/lib/libboost_program_options.a
LIBDIR=${LOCALPATH}/lib/
LIBS=${LIBDIR}/libhts.so $(PROG_OPT_LIB)

all: clean tools
	
tools:
	$(CC) $(INC) $(STXXLINC) $(CFLAGS)  UMINorm.cpp -o umi_norm $(LIBS) -lstdc++fs
	
clean:
	rm -f umi_norm 

