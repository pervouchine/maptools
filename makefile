SRC=src/
BIN=bin/

PROGS=${BIN}map_single ${BIN}map_agnostic ${BIN}best_match ${BIN}net_filter ${BIN}syntenic_filter ${BIN}getwind ${BIN}getsegm ${BIN}transf ${BIN}progressbar.o ${BIN}genutils.o ${BIN}getsegm2

all :: ${PROGS}

${BIN}progressbar.o:  ${SRC}progressbar.c ${SRC}progressbar.h
	g++ -c ${SRC}progressbar.c -o ${BIN}progressbar.o

${BIN}genutils.o:     ${SRC}genutils.c ${SRC}genutils.h
	g++ -c ${SRC}genutils.c -o ${BIN}genutils.o

${BIN}map_single : ${SRC}map_single.c ${BIN}genutils.o ${BIN}progressbar.o
	g++ -o ${BIN}map_single ${SRC}map_single.c ${BIN}genutils.o ${BIN}progressbar.o

${BIN}map_agnostic : ${SRC}map_agnostic.c ${BIN}genutils.o ${BIN}progressbar.o
	g++ -o ${BIN}map_agnostic ${SRC}map_agnostic.c ${BIN}genutils.o ${BIN}progressbar.o

${BIN}best_match : ${SRC}best_match.c ${BIN}genutils.o ${BIN}progressbar.o
	g++ -o ${BIN}best_match ${SRC}best_match.c ${BIN}genutils.o ${BIN}progressbar.o

${BIN}net_filter : ${SRC}net_filter.c ${BIN}genutils.o ${BIN}progressbar.o
	g++ -o ${BIN}net_filter ${SRC}net_filter.c ${BIN}genutils.o ${BIN}progressbar.o

${BIN}syntenic_filter : ${SRC}syntenic_filter.c ${BIN}genutils.o ${BIN}progressbar.o
	g++ -o ${BIN}syntenic_filter ${SRC}syntenic_filter.c ${BIN}genutils.o ${BIN}progressbar.o 

${BIN}getwind : ${SRC}getwind.c ${BIN}genutils.o ${BIN}progressbar.o
	g++ -o ${BIN}getwind ${SRC}getwind.c ${BIN}genutils.o ${BIN}progressbar.o

${BIN}getsegm : ${SRC}getsegm.c ${BIN}genutils.o ${BIN}progressbar.o
	g++ -o ${BIN}getsegm ${SRC}getsegm.c ${BIN}genutils.o ${BIN}progressbar.o

${BIN}getsegm2 : ${SRC}getsegm2.c ${BIN}genutils.o ${BIN}progressbar.o
	g++ -o ${BIN}getsegm2 ${SRC}getsegm2.c ${BIN}genutils.o ${BIN}progressbar.o

${BIN}transf : ${SRC}transf.c ${BIN}genutils.o ${BIN}progressbar.o
	g++ -o ${BIN}transf ${SRC}transf.c ${BIN}genutils.o ${BIN}progressbar.o

clean :: 
	rm -f ${PROGS}

