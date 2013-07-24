#!/bin/bash

for i in {0..15}
do
    ssh compute-0-$i "echo -e '\nexport LD_LIBRARY_PATH=/share/apps/openmpi_intel_20130712/lib:$LD_LIBRARY_PATH\nexport LD_LIBRARY_PATH=/share/apps/intel_20130618/composer_xe_2013.3.163/compiler/lib/intel64:$LD_LIBRARY_PATH\nexport PATH=/share/apps/intel_20130618/composer_xe_2013.3.163/bin/intel64:$PATH\nexport PATH=/share/apps/openmpi_intel_20130712/bin:$PATH' >> /etc/bashrc"
done
