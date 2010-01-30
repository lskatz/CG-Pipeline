MANIFEST = *.pl *.sh lib Makefile mats logs utils README AUTHORS INSTALL FILES
#test-data
C_MANIFEST = mcmc-core.c mcmc-core.h Makefile
DIST_MANIFEST = ${MANIFEST} c
PROJECT = binner
HOST = topaz
AUTHOR = kislyuk

TEST_OPTS = -keep
#TEST_OPTS := ${TEST_OPTS} -plotfile=/storage2/kislyuk/tmp/plots.m

TEST_SETDIR = test-data
#TEST_SETDIR = test-data/trouble
#TEST_SETDIR = test-data/compostbin-sets
#TEST_SETDIR = test-data/srijak-tests
#TEST_SETDIR = .
#TEST_SETFILE = set1.fasta
#TEST_SETFILE = set2.fasta
#TEST_SETFILE = set3.fasta
#TEST_SETFILE = set4.fasta
#TEST_SETFILE = set5.fasta
#TEST_SETFILE = set6.fasta
#TEST_SETFILE = set7.fasta
##TEST_SETFILE = set8.fasta
#TEST_SETFILE = set9.fasta
#TEST_SETFILE = set10.fasta
#TEST_SETFILE = set11.fasta
#TEST_SETFILE = set12.fasta
#TEST_SETFILE = test.fasta
#TEST_SETFILE = amd.fasta
#TEST_SETFILE = amd-4way.fasta
TEST_SETFILE = amdx.fasta

all:
	make

up:
	scp -r ${MANIFEST} ${AUTHOR}@${HOST}:~/${PROJECT}/ ; \
	(cd c; scp -r ${C_MANIFEST} ${AUTHOR}@${HOST}:~/${PROJECT}/c)

roll-dist:
	make clean
	mkdir -p roll-dist-tmp/${PROJECT}
	-cp -R ${DIST_MANIFEST} roll-dist-tmp/${PROJECT}
	(cd roll-dist-tmp && tar -h -c --gzip ${PROJECT} > ../`date +%Y%m%d`-${PROJECT}.tar.gz)
	rm -rf roll-dist-tmp

clean:
	-rm -f bin/*
	(cd c; make clean)

test:
	(cd c; make clean; make)
	time ./test_binner.pl ${TEST_OPTS} ${TEST_SETDIR}/${TEST_SETFILE}

test-compostbin:
	(cd c; make clean; make)
	for chain_order in 3 4; do \
		for CB_SETFILE in $$(cd test-data/compostbin-sets; ls set*.fasta); do \
			time ./test_binner.pl ${TEST_OPTS} test-data/compostbin-sets/$${CB_SETFILE} \
				-chain_order=$$chain_order -plotfile=logs/$${CB_SETFILE}.$${chain_order}o.plot.m \
				> logs/$${CB_SETFILE}.$${chain_order}o.log ; \
		done; \
	done

test-all:
	(cd c; make clean; make)
	for chain_order in 2 3; do \
		for frag_num in 250 2500; do \
			for model_var in 1e-7 5e-7 1e-6 5e-6 1e-5 5e-5 1e-4 5e-4 1e-3; do \
				time ./test_binner.pl ${TEST_OPTS} ${TEST_SETDIR}/${TEST_SETFILE} \
					-frag_num=$$frag_num -model_var=$$model_var -chain_order=$$chain_order -plotfile=logs/${TEST_SETFILE}.$${frag_num}f.$${model_var}v.$${chain_order}o.plot.m \
					> logs/${TEST_SETFILE}.$${frag_num}f.$${model_var}v.$${chain_order}o.log 2>&1; \
			done; \
		done; \
	done

test-leaks:
	valgrind --leak-check=yes c/mcmc-core -o 3 -s 5000 -n 4 -m 5e-05 -f 5e-05 -h /tmp/mcmc.pl.19713.P3hT5/core_hint_model -i 0 < /tmp/mcmc.pl.19713.P3hT5/core_in > valgrind2.out 2>valgrind2.err
