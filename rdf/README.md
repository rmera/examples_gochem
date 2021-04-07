RDF:  A simple tool for the solute-solvent interaction analysis using the goChem library.

USE:

./rdf  [flags] pdbname xtcname solventname

Where "solventname" is the residue name for the solvent molecules in the PDB file. Will obtain the RDF or MDDF for the given solute and solvent in a Gromacs XTC trajectory.
Requires the open source xdrfile library.


Flags:

*  **-refatom int** Index of the atom to be used as reference, counting from 0. (optional, a full residue ID can be given instead, see option -refid)
*  **-refid int** Index of the residue to be used as reference, as it appears in the PDB reference
*  **-com** Use the centre of mass (if available, centroid otherwise) for obtaining the distance to a solvent molecule, instead of the closest  atom
*  **-end** float The largest distance considered, in A (default 10)
*  **-getsystem**  No MDDF/RDF is calculated. Instead, a subsystem with the solute and solvent up to _end_ A is produced for each processed frame of the traj
*  **-skip** int Every how many frames should a frame be read. (default 1)
*  **-step float** The interval for the RDF/MDDF calculation, in A (default 0.1)

One of the two first flags must be given, while the others are optional. 

This program is primarily meant as an example of goChem capabilities. Thus, the options are somewhat limited.

