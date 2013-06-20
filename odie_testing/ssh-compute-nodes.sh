#!/bin/bash
echo "Checking compute node status..."
for i in {0..15}
do
	if ! ssh -q compute-0-$i -C "exit"; then echo "compute-0-$i is down!"; else echo "compute-0-$i is up."; fi
done
echo "Finished checking compute node status."
