#!/bin/csh

set pdbab=$1
pisa test -analyse ${pdbab}
pisa test -detail interfaces 1

