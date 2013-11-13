This readme was created by Jason to try to give some description of the code for the reader, and for himself!

However, you should be careful. I may forget to update this if I make a change!

##Pixels
Ideally, we want every atom in the model to contribute once and only once to the intensity. For RMC this is *essential*. With round pixels we will either have overlapping atoms (square inscribed in circle idea) or we will leave out atoms (cirle inscribed in square). By the way, the model is a box (it cant be a sphere due to boundary effects) and thats where the squares come from.

So for RMC we use square pixels to remove this effect, even though it does not exactly replicate the experiment. Removing this effect for RMC is essential because if an atom moves into a space that is not in the intensity calculation, or is duplicated in the intensity calculation, the intensity is per atom can be seen to be different. This is not good. Note the best explanation...

For Femsim, we do want to exactly replicate the experiment and we are simply calculating the variance of the model. Therefore, we are not moving atoms (and have no input data by the way) and so the RMC probelms are irrelevant here. We will use overlapping pixels instead of leaving atoms out here.

In the code we want the pixel placement to be as obvious for the user as possible. Page 11 in my notebook shows pictures of the pixel setup for a 3x3 case.
