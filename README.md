# Atom-Molecules
CSE232 project6

Here is the project's description.

The distant galaxy MilkyWayx has a peculiar chemistry that is similar to ours, but is much simpler. For example, molecules in the MilkyWayx galaxy are only made up of 3 types of atoms; Hydrogenyx (Hh), Carbonyx (Cc), and Sulphuryx (Ss).
The properties of atoms and molecules are:
Every atom has a certain number of bonds that it can form with other atoms, and it can form at most one bond with any other single atom. When 2 or more atoms interact, they form a molecule. Of course, a molecule is also formed when an atom and a molecule, or 2 molecules interact with each other. The number of bonds each of the atoms can form is as follows: N(Hh)=1, N(Cc)=4, and N(Ss)=7 .
Every atom also has a particular atomic weight associated with it. When a molecule is formed, its molecular weight is simply the combined atomic weights of all the atoms in that molecule. The weights of each atom is as follows: W(Hh)=2.4, W(Cc)=5.6, and W(Ss)=10.8 .
In the MilyWayx galaxy, 2 molecules are equal, if and only if, both molecules have exactly the same number of each type of atom.
In the MilyWayx galaxy, a molecule is stable, if and only if, every atom in the molecule can find a unique atom to bond with every one of its available bonds. Fortunately, Prof. Hakimyx has discovered how molecules become stable, and this gives us an easy way to test if a particular molecule is stable. What you need to do is lay out all the atoms in the molecule in decreasing order of their available bonds. Then simply repeat the following three steps:
Take the first atom (this will have the maximum number of available bonds), and connect each of its available bonds to a unique atom (in decreasing order of available bonds).
Remove the first atom from the list.
For each of the other atoms that was just connected, remove those bonds.
At any point, if any of the above steps can't be done, then the molecule is not stable. If the above procedure stops, or there are no remaining atoms with any available bonds, then the molecule is stable.
Here's a 6 step illustration of the above process which shows the formation of a stable molecule made up of 1 Ss atom, 3 Cc atoms, and 7 Hh atoms.
//s1.jpg
//s2.jpg
//s3.jpg
//s4.jpg
//s5.jpg
//s6.jpg
Here is an illustration of how un unstable molecule made up of 4 Cc atoms and 2 Hh atoms tries and fails to become stable.
//r1.jpg  
//r2.jpg
//r3.jpg
//r4.jpg
//r5.jpg
//r6.jpg
 
Now that you have learnt how the chemistry of the MilkyWay Galaxy works, you should be able to analyse chemical equations that have been broadcast to us by the residents of that galaxy. These equations come in three different forms:
 Simple Molecules: These are simple strings in the form <atom_symbol1><number1><atom_symbol2><number2> ... <atom_symbolN><numberN>, where <number> is assumed to be 1 if it's missing, and the resulting value is just the molecule made up of the atoms in their specified quantities. e.g. Cc3Hh is a molecule made up of 3 Cc atoms and 1 Hh atom. HhSs2 is a molecule made up of 1 Hh atom and 2 Ss atoms.
Reactions: These are strings of the form <simple_molecule1> + <simple_molecule2> + ... + <simple_moleculeN>, where the resulting value is simply the molecule formed by combining N Simple Molecules.
Transformations: These are strings of the form <ReactionLHS> -> <ReactionRHS>. A transformation is valid if the molecule formed by ReactionLHS is equal to the molecule formed by ReactionRHS.
For this project, you are provided 3 files, atoms_and_molecules.h, atoms_and_molecules.cpp (which are empty since you have to write them), and test_atoms_and_molecules.cpp, for testing purposes. Use this test file to guide you through the development of your project.
Note: Implementing this project without templates will be very painful.
