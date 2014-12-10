ERRORFILE=err.log
BINS=trans_brute
LIBS=digraph.o editing.o
RESULTS=$(ls trans_out*)
SSH_SRV_ADR=igel@mpc632.mata.uni-jena.de
#SSH_SRV_ADR=mweller@ipc710.inf-i1.uni-jena.de
#SSH_SRV_ADR=igel@login.minet.uni-jena.de
CFLAGS=-O2 -march=native
GPP=g++ ${CFLAGS}

all: trans_brute

########### library ##################

digraph.o: digraph.cpp digraph.h
	${GPP} -Wall -g -o digraph.o -c digraph.cpp 2>&1 | tee -a ${ERRORFILE}

editing.o: editing.cpp editing.h
	${GPP} -Wall -g -o editing.o -c editing.cpp 2>&1 | tee -a ${ERRORFILE}


############ solution ###################

trans_brute: trans_brute.cpp ${LIBS}
	${GPP} -Wall -g -o trans_brute trans_brute.cpp ${LIBS} 2>&1 | tee -a ${ERRORFILE}

search_dag: search_dag.cpp ${LIBS}
	${GPP} -Wall -g -o search_dag search_dag.cpp ${LIBS} 2>&1 | tee -a ${ERRORFILE}

brute_graph: brute_graph.cpp ${LIBS}
	${GPP} -Wall -g -o brute_graph brute_graph.cpp ${LIBS} 2>&1 | tee -a ${ERRORFILE}

########## general stuff #############

install: all
	cp trans_brute /usr/local/bin

clean:
	rm -f ${BINS} ${LIBS} err.log

package:
	tar -cjf trans.tar.bz2 *.h *.cpp

copy:
	scp *.h *.c* *.txt ${SSH_SRV_ADR}:~/diplomarbeit

test_results:
	scp -r ${SSH_SRV_ADR}:~/diplomarbeit/[0-9]*x* ../data/

view_test_results:
	for file in ${RESULTS}; do \
		cat ${file} | grep edits | sed -e 's/: .*(time/ time/g' | tr -d '\n';\
		echo;\
	done
check_test_results:
	cat ${RESULTS} | grep edits | sed -e 's/edits.*//g' | tr -d '\n' | sed -e 's/ /-/' | sed -e 's/$/\n/' | bc -l ; echo
