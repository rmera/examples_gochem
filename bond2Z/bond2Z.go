package main

import (
	"os"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)


//cross Takes 2 3-len column or row vectors and returns a column or a row
//vector, respectively, with the Cross product of them.
//should panic
func cross(a, b *v3.Matrix) *v3.Matrix {
	c := v3.Zeros(1)
	c.Cross(a, b)
	return c
}

func main () {
	//runtime.GOMAXPROCS(2) ///////////////////////////
	mol, err := chem.XYZFileRead(os.Args[1])
	if err != nil {
		panic(err.Error())
	}
	//The selection thing
	orient_atoms,err := scu.IndexStringParse(os.Args[2])
	if err!= nil{
		panic(err.Error())
	}
	//Get the axis of rotation
	//ov1:=mol.Coord(orient_atoms[0], 0)
	ov2 := mol.Coord(orient_atoms[1], 0)
	//now we center the thing in the beta carbon of D124
	mol.Coords[0].SubVec(mol.Coords[0], ov2)
	//Now the rotation
	ov1 := mol.Coord(orient_atoms[0], 0) //make sure we have the correct versions
	ov2 = mol.Coord(orient_atoms[1], 0)  //same
	orient := v3.Zeros(ov2.NVecs())
	orient.Sub(ov2, ov1)
	//	PDBWrite(mol,"test/2c9v-124centered.pdb")
	Z, _ := v3.NewMatrix([]float64{0, 0, 1})
	axis := cross(orient, Z)
	angle := chem.Angle(orient, Z)
	oldcoords := v3.Zeros(mol.Coords[0].NVecs())
	oldcoords.Copy(mol.Coords[0])
	mol.Coords[0] = chem.Rotate(oldcoords, mol.Coords[0], axis, angle)
	if err != nil {
		panic(err.Error())
	}
	chem.XYZFileWrite("aligned.xyz", mol.Coords[0], mol)
}
