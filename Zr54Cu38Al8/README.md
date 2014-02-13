The fem experimental data is for Zr54Cu38Al8 with oxidation layer subtraction. 
t1 is as cast.
t2 is annealed for 10 minutes at 300C.
t3 is annealed for 60 minutes at 300C.

Date: 9/9/13. I am a couple days into making this stuff now, and it wont be ready for another couple, but this is the approx date of when we started this. 2/11/14 is about when we finished the simulations. Still need to analyze.


##Description of the models:
Zr54Cu38Al8_1000atoms.xyz is the initial model that I generated with generate_model.py (cutoff distance of 1.76, model size is in the file).
I ran two MD simulations, one with 1000 atoms and one with 10000 atoms. Both were used to determine the atom density (#atoms/volume) by letting the volume equilibrate. The 1k atom simulation is fast, and gave a final volume (see out file) that resulted in an atom density of 0.05455541 atoms/A^3. The 10k atom simulation had an atom density of 0.0546284232. The slight difference in atom density turns out to be irrelavent, because the energy/atom of the final models are both within kT of each other at room temp. Also, the pressures both fluctuate around 0 which is good. See the Zr54_sm_tot.in files for the MD specifics.
(Note: Masses should not be defined in .dat file, they are defined in the potential. This keeps things confusing. Make sure the order of the elements on the line in the lammps input file corresponds with the atomtype numbers in the .dat file. Our convention is that both should be in order of decreasing atomic number.)
Checking the RDF/PDF of the final model (converted from the restart file to dat then to xyz for RINGS) gives us the cutoff we need.
We want a box size of 20*sqrt(2), which gives us 1234 atoms in the cubic world (box).
Zr54Cu38Al8_1234atoms_post_md.xyz is the model created by generate_model_pei.py (cutoff of 2.25 was used to generate the model (still too low)).
The above model (Zr54Cu38Al8_1234atoms_post_md.dat) (converted to xyz format) was submitted to MD using lmp_emin_exp.in in order to get it to its lowest energy state within its current local minimum (conjugate gradient method). The resulting output is lmp_emin_exp.restart which was converted to a dat model named Zr54Cu38Al8_1234atoms_start.dat, then to xyz format with the same name (different extension).
Make sure the composition didnt get screwed up in the conversion. Run check_model.py on the xyz file (Zr54Cu38Al8_1234atoms_start.xyz)!
This Zr54Cu38Al8_1234atoms_start.xyz was the one that was submitted to HRMC for all three calculations.


##Here is the general outline of approach:
    1. Make a DRP model using generate_model.py with the correct composition and 1000 atoms. Size of the box isnt that big of a deal, just make it not unusual. 
    2. Run an MD sim on that model (after converting to .dat) that will equilibrate the model volumetrically so we can extract the atom density (see Zr54_sm_tot.in).
    3. Extract the final model and get the cutoff distance, which is found by calculating a summed RDF (use Rings, or see Pauls Igor and note at the end). Use the beginning of the first peak as the cutoff distance in the next step.
    4. For the HRMC sim, we want a model with side lengths of 20*sqrt(2)A = 28.28427A. So modify generate_model.py with these side lengths, calculate the number of atoms in the model using this and the atom density found in step 2, and run generate_model.py. Chances are the code will run indefinitely, at which point you have to reduce the cutoff (until it works). Then you can continue to step 5.
    5. Take the generated model from step 4, convert it to .dat, and run a CGM sim (its super fast) to move the atoms around so that the energy of the model is at its local minimum (the energy will always go down here). When you run this (with 'lmp_odie < lmp_emin_exp.in' after modifying lmp_emin_exp.in with the correct input filename) you have to make sure the difference in the last two energy steps is not too big. Also, make sure the energy/atom is on the order of a few eV (mine was ~ -3.8eV). This lammps sim will result in a .restart file which you must convert to .dat with restart2data. Then convert the .dat to .xyz using model.py. Run an RDF on this to make sure it looks good.
    5a. Run check_model.py on this xyz file and make sure everything looks good - composition especially because it could have gotten screwed up pretty easily. Check everything closely, dont overlook things.
    6. Finally you have the model that you can submit to HRMC! This model now has the correct atomic density and is "random" as much as is possible while still keeping the atoms far enough apart to not spike the potential energy of the structure.

    Note on calculating the RDF in Igor: The input xyz file must be tab deliminated. Use Excel to load and re-save if necessary (you may also need to remove comments after lines; the first line should be a comment). Load the Igor functions from: https://github.com/paul-voyles/Igor/blob/master/STEM%20simulations/XYZ.ipf. Load the .xyz file into Igor with LoadXYZ() and select the .xyz file. Then RDF(xyz, 200) and PartialRDF(xyz,200) create the RDFs which you can display.



Peis simulation parameter settings:
Cutoffs in param_file.in = 2.272A taken from the RDF described above in step 3.
Starting temperature = 200,000 and decreasing it to temp*sqrt(0.7) every 50k moves until you decide to quit.
Max_move distance = 1.5A. (Same as Jinwoos).
Autoslice is not used.
Since we are using using the fem data in the HRMC, we set its weighting factor based on the the below calculation. Beta is the weighting factor.
Peis glass thickness is, on average, 27nm and the thickness gradient is 35-50nm. This means beta = 1/3*ts/te = 1/3*20*sqrt(2)/270 = 0.03492. Then we plug 1/0.03492 = 28.6378 into param_file.in.
The other weighting factor, that between chi2 and the energy, was set to 40. This is an approximate number and based on my calculations in first group meeting research presentation (10-25-13). See that for details. Note I forgot to change alpha based on the number of atoms, but its not that big of a deal. Precision is not especially important with this factor (correct alpha should have been 31.6).


--------------------------------------------------------

###This is info on the actual jobs.
Everything is in t1, t2, and t3 directories. Data is good.
