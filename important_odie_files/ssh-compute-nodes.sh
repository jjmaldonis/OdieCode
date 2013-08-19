#!/bin/bash

# Delete old files.
rm ssh-compute-nodes.log
rm qstat.log
rm email.tmp

# Run qstat.
qstat > qstat.log

# Test compute nodes.
date >> ssh-compute-nodes.log
echo "" >> ssh-compute-nodes.log
echo "Checking compute node status..."
for i in {0..15}
do
	if ! ssh -q compute-0-$i -C "exit"
		then echo "compute-0-$i is down!" >> ssh-compute-nodes.log && SHOULD_EMAIL=1
		else echo "compute-0-$i is up." >> ssh-compute-nodes.log
	fi
done
echo "Finished checking compute node status."

# Email
SHOULD_EMAIL=1
cat /home/jjmaldonis/OdieCode/important_odie_files/ssh-compute-nodes.log > email.tmp && echo "" >> email.tmp && echo "Current jobs running on Odie:" >> email.tmp && cat /home/jjmaldonis/OdieCode/important_odie_files/qstat.log >> email.tmp
if [ $SHOULD_EMAIL -eq 1 ]
	then cat /home/jjmaldonis/OdieCode/important_odie_files/email.tmp | mail -s "ssh-compute-nodes.log" -a "/home/jjmaldonis/OdieCode/important_odie_files/email.tmp" jjmaldonis@gmail.com
fi


# I linked this file to /etc/cron.daily so that it runs every day. I cant figure out how I added it for a user like I did last time.
