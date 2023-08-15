#!/bin/bash
##########################################################################
#BSUB	-J	np-stretching

#BSUB	-q	gpu
#BSUB	-o	stdout.%J.dat
#BSUB	-e	stderr.%J.dat

##########################################################################
nGPU=1
#BSUB	-n	1
nvidia-smi
source queryGPU ${nGPU}
export MV2_ENABLE_AFFINITY=0
##########################################################################

python2 np-stretchingHB.gala --gpu=${LSB_GPU}>a.log

##########################################################################

