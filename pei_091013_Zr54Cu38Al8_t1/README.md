The fem experimental data is for as cast Zr54Cu38Al8 with oxidation layer subtraction. 

Date: 9/9/13. I am a couple days into making this stuff now, and it wont be ready for another couple, but this is the approx date of when we did this.


Description of the models:
Zr54Cu38Al8_1000atoms.xyz is the initial model that I generated with generate_model.py (cutoff distance of 1.76, model size is in the file).
Zr54Cu38Al8_sm_tot.dump is the .dump file after the MD simulation that optimized the box size.
56500000.xyz was taken from the last timestep in Zr54Cu38Al8_sm_tot.dump. The RDF was run on this model using Jinwoos code, which turns out not to be correct. Paul says it looks like the correction for the near density has been done twice.
Paul found that with 1000 atoms, the box side length is 24.581 Angstroms (at 300K). This gives us an atom density of 1000/(24.581)^3 = 0.06732887526; and therefore for a box side length of 20*sqrt(20) Angstroms we need 1000/(24.581)^3*(20*sqrt(2))^3 = 1523.4785 ~= 1523 atoms. This is the number of atoms we need for a box side length of 20*sqrt(2)A = 28.28427A.
Zr54Cu38Al8_1523atoms_post_md.xyz is the model created by generate_model_pei.py (cutoff of 2.0 (still too low), smaller box size (see file)).

The above model (Zr54Cu38Al8_1523atoms_post_md.dat) (converted from Zr54Cu38Al8_1000atoms_post_md.xyz using 'xyz_dat_conversion Zr54Cu38Al8_1523atoms_post_md.dat') was submitted to MD using lmp_emin_exp.in in order to get it to its lowest energy state within its current local minimum. The resulting output is lmp_emin_exp.restart which was converted to a dat model named Zr54Cu38Al8_1523atoms_start.dat.
Zr54Cu38Al8_1523atoms_start.dat was converted to xyz format using 'xyz_dat_conversion Zr54Cu38Al8_1523atoms_post_md.xyz'. This will be the one that is submitted to HRMC.


Here is the general outline of approach:
    1. Make a DRP model using generate_model.py with the correct composition and 1000 atoms. Size of the box isnt that big of a deal, just make it not unusual. 
    2. Run an MD sim on that model (after converting to .dat) that will give us a .dump file and an atom density of the model.
    3. Extract the last timestep in the dump file (see dump_to_xyz.py) and get the cutoff distance, which is found by calculating a summed RDF (see Pauls Igor and note at the end). Use the beginning of the first peak as the cutoff distance in the next step.
    4. For the HRMC sim, we want a model with side lengths of 20*sqrt(2)A = 28.28427A. So modify generate_model.py with these side lengths, calculate the number of atoms in the model using this and the atom density found in step 2, and run generate_model.py. Chances are the code will run indefinitely, at which point you have to reduce the cutoff (until it works) and then do step 5.
    5. Take the generated model from step 4, convert it to .dat, and run an MD sim (its fast) to move the atoms around so that the energy of the model is at its local minimum (the energy will always go down in the MD sim). When you run this (with 'lmp_odie <lmp_emin_exp.in' after modifying lmp_emin_exp.in with the correct input filename) you have to make sure the difference in the last two energy steps is not too big. Also, make sure the energy/atom is on the order of a few eV (mine was ~ -3.8eV). This MD sim will result in a .restart file which you must convert to .dat with restart2data. Then convert the .dat to .xyz using xyz_dat_conversion. Feel free to run an RDF on this to make sure it looks good.
    6. Finally you have the model that you can submit to HRMC! This model now has the correct atomic density and is "random" as much as is possible while still keeping the atoms far enough apart to not spike the potential energy of the structure.


    Note on calculating the RDF: The input xyz file must be tab deliminated. Use Excel to load and re-save if necessary (you may also need to remove comments after lines; the first line should be a comment). Load the Igor functions from: https://github.com/paul-voyles/Igor/blob/master/STEM%20simulations/XYZ.ipf. Load the .xyz file into Igor with LoadXYZ() and select the .xyz file. Then RDF(xyz, 200) and PartialRDF(xyz,200) create the RDFs which you can display.



Peis simulation parameter settings:
Cutoffs in param_file.in = 2.272A taken from the RDF described above in step 3.
Starting temperature = 100,000 and decreasing it to temp*sqrt(0.7) every 50k moves until temp=30. (This is what Jinwoo used).
Max_move distance = 1.5A. (Same as Jinwoos).
Autoslice is not used.
Since we are using using the fem data in the HRMC, we set its weighting factor based on the the below calculation. Beta is the weighting factor.
Peis glass thickness is, on average, 27nm and the thickness gradient is 35-50nm. This means beta = 1/3*ts/te = 1/3*20*sqrt(2)/270 = 0.03492. Then we plug 1/0.03492 = 28.6378 into param_file.in.
There is another weighting factor - that between chi2 and the energy. We call this one alpha. See pg 25 in my notebook. We just set this basically by making a best guess. I ran the simulation with weighting factors of 1 for alpha and of 28.6378 for beta, and the looked at the output chi2 over energy in the output chi_squared file. I looked at the chi2/energy column, and it instabily flattens out around 1000 accpeted moves (I only ran to a bit less than 2000 accepted moves though!). So I took the "flat" part and averged the chi2/energy values. The approximate reciprocal of this number is the number I then used for alpha (390.0). However, when I ran the simulation with alpha=390.0, chi2/energy was still noticably less than 1.0, so I recalculated alpha with this data and got alpha=550.0. This is currenty what is running and I will continue to keep an eye on things.

