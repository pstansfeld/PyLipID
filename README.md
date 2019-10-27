# PyLipID

pylipid.py: is a toolkit to calculate lipid interactions with membrane proteins. 
It calculates: 
1/ lipid interactions with the proteins in terms of their duration, occupancy, num. of lipids surrounding given residues and koff;
2/ lipid binding sites via interaction networks. 

It plots:
1/ lipid interaction with the protein as a function protein residues. 
2/ the calculated lipid koff to each protein residue. 
3/ binding site interaction network. 


Usage:

-f: Trajectories to check. Can be a series of trajectories with similar system settings. Read in by mdtraj.load().
-c: structural information of the trajectories given to -f. Read in by mdtraj.load(). Supported format include gro, pdb xyz, etc. 
-tu: time unit of all the calculations. Available options include ns and us. 
-save_dir: directory where all the results will be located. Will use current working directory if not specified. 
-cutoffs: the double cutoffs used to define lipid interactions. A continuous lipid contact with a given residue starts when the lipid
gets closer to the given residue than the smaller cutoff and ends when the lipid gets farther than the larger cutoff. 
-lipids: specify the lipid residue name 
-lipid_atoms: specify the atoms to check
-nprot: num. of proteins in the system
-resi_offset: Shift the residue index of the protein. Can be useful when a protein with missing residues at its N-terminus was martinized 
to Martini force field, as martinize.py shift the residue index of the first residue to 1 regardless of its original index. 
-plot_koff: plot koff to each protein residues and save all the fiures in a directory koff_{lipid} for each lipid species.
-plot_duration: Plot the averaged interaction duration as a funtion of residue ID for each lipid species.
-plot_occupancy: Plot the average occupancy as a function of reisude ID fror each lipid species.
-plot_lipidcount: Plot the average num. of surrounding lipid as a function of residue ID for each lipid species.


Usage example: 
python pylipid.py -f ./run_1/md.xtc ./run_2/md.xtc -c ./run_1/protein_lipids.gro ./run_2/protein_lipids.gro 
-cutoffs 0.55 1.4 -lipids POPC CHOL POP2 -nprot 1 -resi_offset 5 -plot_koff -plot_duration -plot_lipidcount -plot_occupancy


