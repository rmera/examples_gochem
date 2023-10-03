package main

import (
	"os"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
)

func main() {
	xyz, err := chem.XYZFileRead(os.Args[1])
	if err != nil {
		panic(err.Error())
	}

	for i := 0; i < xyz.Len(); i++ {
		a := xyz.Atom(i)
		if len(os.Args) > 2 && os.Args[2] == "-res" {
			a.MolID = i
		} else {
			a.MolID = 1
		}
		a.ID = i
		a.Chain = "A"
		a.MolName = "UNK"
		a.Name = strings.ToUpper(a.Symbol) + strconv.Itoa(i)
	}
	name := strings.Replace(strings.ToLower(os.Args[1]), ".xyz", ".pdb", 1)
	chem.PDBFileWrite(name, xyz.Coords[0], xyz, nil)
}
