VLIB= -g -O0

MY_CXX_FLAGS= -std=c++11 #-Wall -Wextra -DNDEBUG
MY_CXX_OPT_FLAGS= -O3 -m64 
MY_CXX=g++

M=0
MODE = -DM=$(M) 
B=0
BIN = -DB=$(B)
LFLAGS = -lm -ldl

CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(LFLAGS) $(MODE) $(BIN)

####

LIBOBJ = \
	../external/rankbv/rankbv.o\
	../external/malloc_count/malloc_count.o

####

INPUT=../../dataset/reads.100.fastq
OUTPUT=../../output.fq

####


all: main 

%.o:%.cpp
	$(MY_CXX) -c -o $@ $<  $(CXX_FLAGS) 

main: fq_compression.cpp ${LIBOBJ} 
	$(MY_CXX) fq_compression.cpp -o fq_compression ${LIBOBJ} $(CXX_FLAGS) 

clean:
	rm -f fq_compression 

run:
	./fq_compression -e $(INPUT).bwt -q $(INPUT).bwt.qs -f $(INPUT).bwt -o $(OUTPUT) 
	
valgrind: 
	$(MY_CXX) fq_compression.cpp -o fq_compression $(VLIB)
#	valgrind --tool=memcheck --leak-check=full --track-origins=yes ./fq_compression -e $(INPUT).bwt -q $(INPUT).bwt.qs -f $(INPUT).bwt -o $(OUTPUT) 

