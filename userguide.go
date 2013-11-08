// +build one
package main

import (
	"fmt"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/xtc"
)

func main() {
	One()
	Two()
	Three()

}

//One opens the sample.xyz file in the test directory, and pull a number of hardcoded atoms
//In the direction of a hardcoded vectos. It builds 12 files with the pulled atoms  displaced by
//different ammounts along the pulling vector
func One() {
	pulled_atoms := []int{43, 41, 42, 40, 85, 86, 87} //indexes
	pulling_vector := []int{40, 88}                   //The direction in which we pull is given by atoms 40 and 88 counting from 0
	pull_ammount := 4.0                               //how much do we want to pull, in Angstroms
	mol, err := chem.XYZRead("sample.xyz")
	if err != nil {
		panic(err.Error())
	}
	pulled := chem.ZeroVecs(7)
	pulled.SomeVecs(mol.Coords[0], pulled_atoms) //We use the pulled_atoms set of indexes to get the atoms we will be pulling
	at1 := mol.Coord(pulling_vector[0], 0)
	vector := chem.ZeroVecs(1)
	vector.Copy(mol.Coord(pulling_vector[1], 0)) //We copy to avoid altering the atom 88 coordinates
	vector.Sub(vector, at1)
	vector.Unit(vector)
	vector.Scale(pull_ammount, vector) //we started with an unit lenght bector and now multiply by the desired pull ammount
	pulled.AddRow(pulled, vector)
	if err != nil {
		panic(err.Error())
	}
	mol.Coords[0].SetVecs(pulled, pulled_atoms) //Now we put the pulled coordinates into the original molecule
	chem.XYZWrite("sample_pulled.xyz", mol.Coords[0], mol)
}

func Two() {
	mol, err := chem.PDBRead("2c9v.pdb", true)
	if err != nil {
		panic(err.Error())
	}
	//Use the reference to get the residue indexes for all histidine residues.
	//This is the function written below.
	ResSelected := SelectResidue(mol, "HIS")
	fmt.Println("The histidines in this proteins are residues number:", ResSelected)
	allowed_chains := []string{"A", "B"}
	//With the follwing we obtains all atoms that belong to the desired
	//molecules, in the allowed chains. In this case we allow both chains
	//(A and B) in the PDB file.
	HisAtoms := chem.Molecules2Atoms(mol, ResSelected, allowed_chains)
	//now you can use HisAtoms as indexes for chem.SomeCoords() and follow the first workflow.
	//I'm too lazy to do it here ;-)
	fmt.Println("Atom indexes for all histidine atoms:", HisAtoms)

}

//Uses a reference (mol) and a residue name (3-letter code) and returns a list with the ID of those residues in the molecule.
func SelectResidue(mol chem.Ref, residue string) []int {
	selected_residues := []int{} //we still dont select anything, so empty slice
	prevmol := -1                //This is the index of the molecular ID of the last atom read (a meaningless negative number initially)
	//it helps to make sure we process each molecule only once.
	for j := 0; j < mol.Len(); j++ {
		if mol.Atom(j).Molid != prevmol && mol.Atom(j).Molname == residue {
			selected_residues = append(selected_residues, mol.Atom(j).Molid) //if the residue match and we have not processed
			//this residue before, we add it to the list.
		}
		prevmol = mol.Atom(j).Molid
	}
	return selected_residues
}

func Three() {
	traj, err := xtc.New("test.xtc")
	if err != nil {
		panic(err.Error())
	}
	ProcessTraj(traj)
}

//Obtains and prints the distance between atoms 2 and 10 (counting from zero) for each frame of trajectory
//traj.
func ProcessTraj(traj chem.Traj) {
	coords := chem.ZeroVecs(traj.Len())
	for i := 0; ; i++ { //infinite loop, we only break out of it by using "break"
		err := traj.Next(coords) //Obtain the next frame of the trajectory.
		if err != nil && err.Error() != "No more frames" {
			panic(err.Error())
			break
		} else if err == nil {
			atom10 := coords.VecView(10)
			atom2 := coords.VecView(2)
			fmt.Println("Distance between the third and tenth atoms in the ", i+1, " frame: ", Distance(atom10, atom2), "A")
		} else {
			break //leave the loop if the trajectory is over (i.e. if we get a "No more frames" error when trying to read a new fram)
		}
	}
}

//Calculates and returns the distance between two atoms.
func Distance(atom1, atom2 *chem.VecMatrix) float64 {
	res := chem.ZeroVecs(1)
	res.Sub(atom1, atom2)
	return res.Norm(0)
}
