#!/bin/bash

export PATH=$PATH:/usr_storage/software/cactus/bin
export PYTHONPATH=/usr_storage/software/cactus/submodules/:$PYTHONPATH
source /usr_storage/software/cactus/cactus_env/bin/activate
export PATH="$PATH:/usr_storage/software/bedtools2/bin/"

nohup hal2maf cactus.hal Pkor_as_ref.maf --refGenome Pkor --noAncestors --onlyOrthologs --noDupes 






