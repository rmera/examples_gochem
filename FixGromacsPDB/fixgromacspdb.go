package main

import (
	"os"
	"strings"

	chem "github.com/rmera/gochem"
)

func main() {
	mol, err := chem.PDBFileRead(os.Args[1], true)
	if err != nil {
		panic(err.Error())
	}
	coord := mol.Coords[0]
	chem.FixGromacsPDB(mol)
	newname := strings.Replace(os.Args[1], ".pdb", "_fixed.pdb", 1)
	chem.PDBFileWrite(newname, coord, mol, nil)

}
