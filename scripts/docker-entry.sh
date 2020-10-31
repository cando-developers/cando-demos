#! /bin/bash
source /opt/clasp/bin/setenv-clasp
export CLASP_FEATURES=jit-log-symbols
/opt/clasp/bin/cando -f no-auto-lparallel -f cando-jupyter -e "(ql:quickload :cando-jupyter)" -e '(cando-user:jupyterlab-fork-server "/tmp/clasp-fork-server/")' &
sleep 30
exec "$@"
