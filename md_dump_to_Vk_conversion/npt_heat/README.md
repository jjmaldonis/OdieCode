To submit all jobs run:
python xyz_to_femsim.py

However, some setup is required.
Every directory in this dir MUST be a directory that you want to run femsim in.
For example, I have 5 directories. Their names do not matter.
However, each must have certain files and directories in them:
vk/ directory for vk output
xyz_files/ directory containing all the xyz models you want run through femsim
param_file.in containing the correct parameters and the fem input filename
A fem input file with the correct relative path (i.e. 'dir_name/fem.txt'). Be careful with this one, the filename must be short. If you get output errors about "Attempt to use pointer K when it is not associated wtih a target" then your 'dir_name/fem.txt' is too long. Shorten it and rerun.
