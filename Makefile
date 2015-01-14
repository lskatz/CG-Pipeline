MANIFEST = Makefile README conf doc lib scripts
DIST_MANIFEST = ${MANIFEST}
PROJECT = cg_pipeline
HOST = topaz
AUTHOR = kislyuk
DESTDIR = "/opt/${PROJECT}"

TEST_SETDIR = test-data
TEST_SETFILE = amdx.fasta

install:
	-echo "Installing ${PROJECT} to ${DESTDIR}..."
	make install-app
	-echo "Initializing databases in ${DESTDIR}/data..."
	make init_databases

install-app:
	mkdir -p ${DESTDIR}
	cp -R ${MANIFEST} "${DESTDIR}"
	#make revise

init_databases:
	(mkdir -p "${DESTDIR}/data"; cd "${DESTDIR}/data"; PATH="${DESTDIR}/scripts:$${PATH}"; update_cg_databases)

#up:
#	scp -r ${MANIFEST} ${AUTHOR}@${HOST}:~/${PROJECT}/

#roll-dist:
#	make clean
#	mkdir -p roll-dist-tmp/${PROJECT}
#	-cp -R ${DIST_MANIFEST} roll-dist-tmp/${PROJECT}
#	(cd roll-dist-tmp && tar -h -c --gzip ${PROJECT} > ../`date +%Y%m%d`-${PROJECT}.tar.gz)
#	rm -rf roll-dist-tmp

test:
	echo "Running test suite..."
	./scripts/run_pipeline_checkPrereqs.pl
	#./scripts/test_cg_deps


config:
	echo "Running configuration..."
	./scripts/reconfigure.pl

#revise:
#	echo "revising Makefile..."
#	(./scripts/revise_make.pl ${DESTDIR}/Makefile ${DESTDIR})
