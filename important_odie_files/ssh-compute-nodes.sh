#!/bin/bash

export PATH=/share/apps/openmpi_intel_20130712/bin:/share/apps/intel_20130618/composer_xe_2013.3.163/bin/intel64:/share/apps/openmpi_intel_20130712/bin:/share/apps/intel_20130618/composer_xe_2013.3.163/bin/intel64:/share/apps/openmpi_intel_20130712/bin:/usr/lib64/qt-3.3/bin:/state/partition1/apps/intel_20130618/composer_xe_2013.3.163/bin/intel64:/state/partition1/apps/intel_20130618/composer_xe_2013.3.163/mpirt/bin/intel64:/state/partition1/apps/intel_20130618/composer_xe_2013.3.163/bin/intel64:/state/partition1/apps/intel_20130618/composer_xe_2013.3.163/bin/intel64_mic:/state/partition1/apps/intel_20130618/composer_xe_2013.3.163/debugger/gui/intel64:/share/apps/intel_20130618/advisor_xe_2013/bin64:/share/apps/intel_20130618/inspector_xe_2013/bin64:/share/apps/intel_20130618/vtune_amplifier_xe_2013/bin64:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/eclipse:/opt/ganglia/bin:/opt/ganglia/sbin:/usr/java/latest/bin:/opt/maven/bin:/opt/pdsh/bin:/opt/rocks/bin:/opt/rocks/sbin:/opt/condor/bin:/opt/condor/sbin:/usr/local/samba/sbin:/opt/gridengine/bin/linux-x64
export SGE_CELL=default
export SGE_ARCH=linux-x64
export SGE_EXECD_PORT=537
export SGE_QMASTER_PORT=536
export SGE_ROOT=/opt/gridengine



# Delete old files.
rm /home/jjmaldonis/OdieCode/important_odie_files/ssh-compute-nodes.log
rm /home/jjmaldonis/OdieCode/important_odie_files/qstat.log
rm /home/jjmaldonis/OdieCode/important_odie_files/email.tmp

cd /home/jjmaldonis/OdieCode/important_odie_files

# Run qstat.
/opt/gridengine/bin/linux-x64/qstat -u "*" 2>&1 | tee /home/jjmaldonis/OdieCode/important_odie_files/qstat.log

# Test compute nodes.
date >> /home/jjmaldonis/OdieCode/important_odie_files/ssh-compute-nodes.log
echo "" >> /home/jjmaldonis/OdieCode/important_odie_files/ssh-compute-nodes.log
echo "Checking compute node status..."
for i in {0..15}
do
	if ! /usr/bin/ssh -q compute-0-$i -C "exit"
		then echo "compute-0-$i is down!" >> /home/jjmaldonis/OdieCode/important_odie_files/ssh-compute-nodes.log && SHOULD_EMAIL=1
		else echo "compute-0-$i is up." >> /home/jjmaldonis/OdieCode/important_odie_files/ssh-compute-nodes.log
	fi
done
echo "Finished checking compute node status."

# Email
SHOULD_EMAIL=1
cat /home/jjmaldonis/OdieCode/important_odie_files/ssh-compute-nodes.log > /home/jjmaldonis/OdieCode/important_odie_files/email.tmp && echo "" >> /home/jjmaldonis/OdieCode/important_odie_files/email.tmp && echo "Current jobs running on Odie:" >> /home/jjmaldonis/OdieCode/important_odie_files/email.tmp && cat /home/jjmaldonis/OdieCode/important_odie_files/qstat.log >> /home/jjmaldonis/OdieCode/important_odie_files/email.tmp
if [ $SHOULD_EMAIL -eq 1 ]
	then cat /home/jjmaldonis/OdieCode/important_odie_files/email.tmp | mail -s "ssh-compute-nodes.log" -a "/home/jjmaldonis/OdieCode/important_odie_files/email.tmp" jjmaldonis@gmail.com
fi


# I linked this file to /etc/cron.daily so that it runs every day. I cant figure out how I added it for a user like I did last time.
