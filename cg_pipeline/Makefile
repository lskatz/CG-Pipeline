MANIFEST = Makefile  README  conf  doc  lib  scripts data
DIST_MANIFEST = ${MANIFEST}
PROJECT = cg_pipeline
HOST = topaz
AUTHOR = kislyuk

TEST_OPTS = -keep

TEST_SETDIR = test-data
TEST_SETFILE = amdx.fasta

all: init_databases

init_databases:
	(cd data; ../scripts/update_uniprot)
	(cd data; ../scripts/update_vfdb)

install:
	echo "todo: make install"

up:
	scp -r ${MANIFEST} ${AUTHOR}@${HOST}:~/${PROJECT}/

roll-dist:
	make clean
	mkdir -p roll-dist-tmp/${PROJECT}
	-cp -R ${DIST_MANIFEST} roll-dist-tmp/${PROJECT}
	(cd roll-dist-tmp && tar -h -c --gzip ${PROJECT} > ../`date +%Y%m%d`-${PROJECT}.tar.gz)
	rm -rf roll-dist-tmp

clean:
	-rm -f bin/*

test:
	(cd c; make clean; make)
	time ./test_binner.pl ${TEST_OPTS} ${TEST_SETDIR}/${TEST_SETFILE}
