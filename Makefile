MANIFEST = Makefile README conf doc lib scripts
DIST_MANIFEST = ${MANIFEST}
PROJECT = cg_pipeline
DESTDIR := "/opt/${PROJECT}"

TEST_SETDIR = test-data
TEST_SETFILE = amdx.fasta


install: config
	-echo "Run 'make test' to check your prerequesites. Run 'make data' to install all databases (takes a long time!)"

data:
	(mkdir -p "${DESTDIR}/data"; cd "${DESTDIR}/data"; PATH="${DESTDIR}/scripts:$${PATH}"; update_cg_databases)

test:
	./scripts/run_pipeline_checkPrereqs.pl

config:
	./scripts/reconfigure.pl

