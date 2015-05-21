package main

import (
	"fmt"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
	"os"
)

//improper.go Takes an XYZ file and the zero-based indexes of 4 atoms (a,b,c and d) It calculates the improper
//dihedral between a plane defined by the vectors ab and ac and the vector ad. The atoms b and c must be given
//in an order such that the cross produc ab x ac points in the same general direction as the vector ad. The
//indexes a,b,c and d have to be given after the XYZ file name as a quoted string separated by spaces.

//example:

// ./improper file.xyz "0 1 2 3"

func main() {
	mol, err := chem.XYZFileRead(os.Args[1])
	if err != nil {
		panic(err.Error())
	}
	ndx, _ := scu.IndexStringParse(os.Args[2])
	fmt.Println(ndx)
	coord := mol.Coords[0]
	piv := coord.VecView(ndx[0])
	v1 := coord.VecView(ndx[1])
	v2 := coord.VecView(ndx[2])
	v1.Sub(v1, piv)
	v2.Sub(v2, piv)
	plane := v3.Zeros(1)
	plane.Cross(v1, v2)
	vfin := coord.VecView(ndx[3])
	vfin.Sub(vfin, piv)
	angle := chem.Angle(vfin, plane)
	fmt.Println(90-angle*chem.Rad2Deg, 90+angle*chem.Rad2Deg)
}
