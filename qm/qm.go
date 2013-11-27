package main

import (
	"fmt"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/qm"
	"math"
)


//This is a throw-away mini program, so the data is hardcoded.
func main() {
	mol, err := chem.XYZRead("../sample.xyz")
	if err != nil {
		panic(err.Error())
	}
	mol.SetCharge(1)
	mol.SetMulti(1)
	//Setting up the calculations
	//We will mosstly go with the defaults.
	calc := new(qm.Calc)
	calc.SCFTightness = 1 //rather demanding.
	calc.Method = "TPSS"
	calc.Dielectric = 80  //COSMO with epsilon=80. Just delete this line to avoid COSMO usage.
	calc.Basis = "def2-TZVP"
	calc.RI = true //RI approximation (also called charge-density fitting, in some programs).
	calc.Disperssion = "D3"
	orca := qm.NewOrcaHandle() //Change this line to use MOPAC2012 or Turbomole
	//Now we play with a bond and make orca inputs
	// to calculate a SP energy for each geometry.
	//*******************************************
	//I just hard-coded these ones, they make sense
	//for my test system but you will have to change them for yours.
	axis1 := mol.Coords[0].VecView(6) //the 2 atoms defining the rotation axis
	axis2 := mol.Coords[0].VecView(7)
	torotate_indexes := []int{8, 9, 10, 11, 70, 83, 69}
	torotate:=chem.ZeroVecs(len(torotate_indexes))
	torotate.SomeVecs(mol.Coords[0],torotate_indexes)                     //The atoms that will rotate
	angles := []float64{0, 0.2 * math.Pi, 0.4 * math.Pi, 0.5 * math.Pi} //The rotation angles in radians.
	for _, angle := range angles {
		rotated, err := chem.RotateAbout(torotate, axis1, axis2, angle)
		if err != nil {
			panic(err.Error())
		}
		//now we put the rotated coords in the molecule
		mol.Coords[0].SetVecs(rotated, torotate_indexes)
		orca.SetName(fmt.Sprintf("angle%1.1f", angle))
		//We first write the QM input and then an XYZ witht he non optimized geometry.
		if err := orca.BuildInput(mol.Coords[0],mol, calc); err != nil {
			panic(err.Error())
		}
		chem.XYZWrite(fmt.Sprintf("angle%1.1f.xyz", angle),mol.Coords[0],mol)
	}
}

