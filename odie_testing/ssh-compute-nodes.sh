#!/bin/bash
rm ssh-compute-nodes.log
date >> ssh-compute-nodes.log
SHOULD_EMAIL=0
echo "Checking compute node status..."
for i in {0..15}
do
	if ! ssh -q compute-0-$i -C "exit"
		then echo "compute-0-$i is down!" >> ssh-compute-nodes.log && SHOULD_EMAIL=1
		else echo "compute-0-$i is up." >> ssh-compute-nodes.log
	fi
done
echo "Finished checking compute node status."
if [ $SHOULD_EMAIL -eq 1 ]
	then cat /home/jjmaldonis/OdieCode/odie_testing/ssh-compute-nodes.log | mail -s "ssh-compute-nodes.log" -a "/home/jjmaldonis/OdieCode/odie_testing/ssh-compute-nodes.log" jjmaldonis@gmail.com
fi
