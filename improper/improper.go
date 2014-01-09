package main

import (
	"fmt"
	"github.com/rmera/gochem"
	"github.com/rmera/scu"
	"os"
)


//improper.go Takes an XYZ file and the zero-based indexes of 4 atoms (a,b,c and d) It calculates the improper dihedral between 
//a plane defined by the vectors ab and ac and the vector ad. The atoms b and c must be given in an order such that the 
//cross produc ab x ac points in the same general direction as the vector ad.

func main() {
	mol, err := chem.XYZRead(os.Args[1])
	if err != nil {
		panic(err.Error())
	}
	ndx,_:=scu.IndexStringParse(os.Args[2])
	fmt.Println(ndx)
	coord := mol.Coords[0]
	piv:=coord.VecView(ndx[0])
	v1:=coord.VecView(ndx[1])
	v2:=coord.VecView(ndx[2])
	v1.Sub(v1,piv)
	v2.Sub(v2,piv)
	plane:=chem.ZeroVecs(1)
	plane.Cross(v1,v2)
	v3:=coord.VecView(ndx[3])
	v3.Sub(v3,piv)
	angle:=chem.Angle(v3,plane)
	fmt.Println(90-angle*chem.Rad2Deg, 90+angle*chem.Rad2Deg)
}
