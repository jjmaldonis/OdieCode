OdieCode
========

Below I will detail how to install ROCKS 6.1.

Base installation of ROCKS 6.1

Currently our RAID card can only handle four drives - they are connected to the card. We can buy a new RAID card to support more drives - talk to Dan at ACE (daniel@acecomputers.com).

Partitioning:
We have four 2TB hard drives that we are making into RAID 5.

Using the hardware RAID controller, configure the RAID (Alt+3 during boot) to RAID 5, set the stripe to 256kb, and make the boot disk 250 (GB). The rest will be made into another partition (5.7TB). Note: 256kb striping is the maximum. We can also do 64kb striping which may be better if we have smaller file sizes. If our files are usually over 256kb then the 256kb is what we want. However, it may only be a minor speed improvement, I don't know.

Since we still have the 2 other hard drives in the other slots which are connected to the motherboard, there will be four drives seen in 'parted' now. During the previous parted setup, sda was the RAID 5 250GB disk, sdb was the remainder of the RAID 5 config (5.7TB), and sdc and sdd were the other two hard drives. However, when I booted using ROCKS 6.1, sda and sdb were the two 2TB hard drives connected to the motherboard and the others were the RAID config with the boot as sdc and the rest as sdd. Just be careful and know which one you are messing with.

Boot using a rescue cd with parted on it.
For /dev/sdb:
parted /dev/sdb
(parted) mklabel gpt
(parted) mkpart primary 0 -0 ****** The 0 -0 arguments mean start at sector 0 and go until sector -0, i.e. the end. ******
quit
mkfs.xfs /dev/sdb1

For /dev/sda:
parted /dev/sda
mklabel msdos
quit

reboot

Install ROCKS 6.1
Put in the ROCKS 6.1 install DVD will all the rolls on it. Enter the following information when prompted:
Select the following rolls:
base, hpc, kernel, os, sge, ganglia, web-server, area51
Fully-Qualified Host Name: odie.msae.wisc.edu
Cluster Name: Odie
Leave the rest of the fields as is.
eth1 IP should be 128.104.188.162 with netmask 255.255.252.0
eth0 IP is 10.1.1.1 with netmask 255.255.0.0
Gateway 128.104.188.1
DNS Servers: 128.104.201.24,144.92.12.24

Select manual partitioning and create the following:
On the 250GB disk ONLY (uncheck the other boxs for the following 4 partitions) create:
/boot   with 256MB (ext4)
swap   with 32768MB (swap)
/var     with 32768MB (ext4)
/          with the remaining space on the 250GB disk (ext4)
On the other RAID disk DO NOT change the type, simple write in the name/label: /state/partition1
Click Next and ROCKS should download and install the OS. It takes a bit, maybe 20-30 min. It will then reboot automatically and, having ejected the disk, will boot into ROCKS.

You now how a functional CentOS 6.3 installation with ROCKS 6.1!
Login using 'root' as the username and the password you specified during installation.




Post installation (need to install in /export/apps/<name> for all the nodes to access - on the nodes, this will be /share/apps/<name>).

To install yum (only need to install on the head node): 
wget yum-3.2.29-40.el6.centos.noarch.rpm    (link taken from http://mirror.stanford.edu/yum/pub/centos/6/os/x86_64/Packages/ -- Ctrl+F "yum" to find a similar one)
rpm -Uvh yum-3.2.29-40.el6.centos.noarch.rpm

To install tw_cli (only need to install on the head node):
wget http://dl.atrpms.net/el5-x86_64/atrpms/stable/tw_cli-2.00.03.018-7.x86_64.rpm
rpm -Uvh tw_cli-2.00.03.018-7.x86_64.rpm

To install the Intel compilers, winSCP the .tgz file over and follow this guide:
http://software.intel.com/sites/default/files/article/251099/release-notes-studio-xe-2013-l.pdf
Below is a general outline of what I did.
tar -xvfz intel...tgz
cd into directory that was created
./install.sh 2>&1 | tee intel_install_<date>.out
Follow the prompts. Set the install directory to: /export/apps/intel_<date>
I got the error:
"32-bit libraries not found on this system.
This product release requires the presence of 32-bit compatibility libraries
when running on Intel(R) 64 architecture systems. One or more of these libraries
could not be found:
    libstdc++ (including libstdc++6)
    glibc
    libgcc
Without these libraries, the compiler will not function properly.  Please refer 
to Release Notes for more information."
at one point so I installed these libraries.
I searched for them using yum search <name> and then installed all the *.i686 ones with yum install <full_name> as taken from the yum search.
I then rechecked the installation dependencies and that problem was solved.
However, I still have the problem:
"The system does not use an Intel Architecture processor. The drivers for
Hardware Event-based Sampling (EBS) data collection and Power analysis require a
genuine Intel processor and will be disabled. To enable this functionality,
install the product on a system with an Intel Architecture processor."
Considering it wants us to put in a new processor, I just ignored this error and continued with the installation.
It finished successfully. I then executed the "source" commands that the output provides. Namely,
source /share/apps/intel/intel-original-2013-x84-mv0.0/vtune_amplifier_xe_2013/amplxe-vars.sh
source /share/apps/intel/intel-original-2013-x84-mv0.0/inspector_xe_2013/inspxe-vars.sh
source /share/apps/intel/intel-original-2013-x84-mv0.0/advisor_xe_2013/advixe-vars.sh
source /share/apps/intel/intel-original-2013-x84-mv0.0/bin/compilervars.sh intel64
Your compilers are now:
    For C++: icpc
    For C: icc
    For Fortran: ifort
 Now you need to make sure to install all other software that will run on the nodes with the Intel compilers (i.e. OpenMPI).

To install OpenMPI:
Download from http://www.open-mpi.org/software/ompi/v1.6/ (the .gz one)
tar xvfz openmpi-1.6.4.tar.gz
cd openmpi-1.6.4
./configure --prefix=/export/apps/openmpi_intel_20130618/ CC=icc CXX=icpc F77=ifort FC=ifort --with-sge 2>&1 | tee -a openmpi-intel_configure_20130618.log
make all install 2>&1 | tee -a openmpi-intel_make_20130618.log
 
To add users (besides for root) do:
adduser <user_name>
passwd <user_name> (then enter your password)
rocks sync users    (this takes a while)
   *** IMPORTANT NOTE: "rocks sync users" seems to push the user data onto the compute nodes (once they login and only if they ssh into the compute node I think) which makes sense, but it also seems to push the /share/apps as well. This is important and not exactly intuitive. I don't know then if each user needs to ssh into each compute node???? That doesn't make sense...***

Set up your environment variables so that everything runs.

Install the compute nodes with "insert-ethers" following the ROCKS installation guide. Our ethernet switch didn't connect last time.
To reboot the compute nodes if you have installed them already:
*** Read this first. ***
rocks run host compute


***Note*** that /share seems empty on the compute nodes but apps/ it is actually in there. You can check by cd /share/apps followed by ls to see the contents of /share/apps/.
