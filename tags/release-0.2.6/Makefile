MANIFEST = Makefile README conf doc lib scripts data
DIST_MANIFEST = ${MANIFEST}
PROJECT = cg_pipeline
HOST = topaz
AUTHOR = kislyuk
DESTDIR = "/opt/${PROJECT}"

TEST_SETDIR = test-data
TEST_SETFILE = amdx.fasta

init_databases:
	(cd "${DESTDIR}/data"; PATH="${DESTDIR}/scripts:$${PATH}" update_cg_databases)

install:
	-echo "Installing ${PROJECT} to ${DESTDIR}..."
	mkdir -p ${DESTDIR}
	cp -R ${MANIFEST} "${DESTDIR}"
	-echo "Initializing databases in ${DESTDIR}/data..."
	make init_databases

up:
	scp -r ${MANIFEST} ${AUTHOR}@${HOST}:~/${PROJECT}/

roll-dist:
	make clean
	mkdir -p roll-dist-tmp/${PROJECT}
	-cp -R ${DIST_MANIFEST} roll-dist-tmp/${PROJECT}
	(cd roll-dist-tmp && tar -h -c --gzip ${PROJECT} > ../`date +%Y%m%d`-${PROJECT}.tar.gz)
	rm -rf roll-dist-tmp

test:
	echo "Running test suite..."
	./scripts/test_cg_deps
