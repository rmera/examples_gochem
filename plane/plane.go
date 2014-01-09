package main

import (
	"fmt"
	"github.com/rmera/gochem"
	"github.com/rmera/scu"
	"os"
)

//This program will align the best plane passing through a set of atoms in a molecule with the XY-plane.
//Usage:
func main() {
	if len(os.Args) < 2 {
		fmt.Printf("Usage:\n%s file.xyz [indexes.dat]\nindexes.dat is a file containing one single line, with all the atoms defining the plane separated by spaces. If it is not given, all the atoms of the molecule will be taken to define the plane.\n", os.Args[0])
		os.Exit(1)
	}
	mol, err := chem.XYZRead(os.Args[1])
	if err != nil {
		panic(err.Error())
	}
	var indexes []int
	//if no file with indexes given, will just use all the atoms.
	if len(os.Args) < 3 {
		indexes = make([]int, mol.Len())
		for k, v := range indexes {
			indexes[k] = v
		}
	} else {
		indexes, err = scu.IndexFileParse(os.Args[2])
		if err != nil {
			panic(err.Error())
		}
	}
	some := chem.ZeroVecs(len(indexes)) //will contain the atoms selected to define the plane.
	some.SomeVecs(mol.Coords[0], indexes)
	//for most rotation things it is good to have the molecule centered on its mean.
	mol.Coords[0], _, _ = chem.MassCentrate(mol.Coords[0], some, nil)
	//As we changed the atomic positions, must extract the plane-defining atoms again.
	some.SomeVecs(mol.Coords[0], indexes)
	//The strategy is: Take the normal to the plane of the molecule (os molecular subset), and rotate it until it matches the Z-axis
	//This will mean that the plane of the molecule will now match the XY-plane.
	best, err := chem.BestPlane(some, nil)
	if err != nil {
		panic(err.Error())
	}
	z, _ := chem.NewVecs([]float64{0, 0, 1})
	zero, _ := chem.NewVecs([]float64{0, 0, 0})
	fmt.Fprintln(os.Stderr, "Best  Plane", best, z, indexes)
	axis := chem.ZeroVecs(1)
	axis.Cross(best, z)
	fmt.Fprintln(os.Stderr, "axis", axis)
	//The main part of the program, where the rotation actually happens. Note that we rotate the whole
	//molecule, not just the planar subset, this is only used to calculate the rotation angle.
	mol.Coords[0], err = chem.RotateAbout(mol.Coords[0], zero, axis, chem.Angle(best, z))
	if err != nil {
		panic(err.Error())
	}
	//Now we write the rotated result.
	final,err:=chem.XYZStringWrite(mol.Coords[0], mol)
	fmt.Print(final)
	fmt.Fprintln(os.Stderr,err)
}
