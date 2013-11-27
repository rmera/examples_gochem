package main

import (
	"fmt"
	"github.com/rmera/gochem"
)


//A very basic analysis of the moment of inertia tensor.
func main() {
	mol, err := chem.XYZRead("../sample.xyz")
	if err != nil {
		panic(err.Error())
	}
	mass, err := mol.Masses()
	if err != nil {
		mass = nil //MomentTensor will simply assign 1 to all masses
	}
	coords := mol.Coords[0]
	moment, err := chem.MomentTensor(coords, mass)
	if err != nil {
		panic(err.Error())
	}
	eigvectors, eigvalues, err := chem.EigenWrap(moment,-1)
	if err != nil {
		panic(err.Error())
	}
	main := eigvectors.VecView(2) //The last is the widest.
	fmt.Printf("Widest eigenvector: %v  Corresponding eigenvalue: %4.1f\n", main, eigvalues[2])
}
