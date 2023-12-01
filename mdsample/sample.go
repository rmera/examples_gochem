package main

import (
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/qm"
	"github.com/rmera/gochem/traj/dcd"
	v3 "github.com/rmera/gochem/v3"
)

func Err(err error) {
	if err != nil {
		log.Fatal(err.Error())
	}
}

func MDErr(err error) bool {
	if _, ok := err.(chem.LastFrameError); ok {
		return true
	} else {
		log.Fatal(err.Error())
	}
	return true
}

func main() {
	name := os.Args[1]
	namext := strings.Split(name, ".")
	var mol *chem.Molecule
	var err error
	if strings.ToLower(namext[len(namext)-1]) == "pdb" {
		mol, err = chem.PDBFileRead(name)
	} else {
		mol, err = chem.XYZFileRead(name)
	}
	Err(err)
	trjname := os.Args[2]
	trjext := strings.Split(trjname, ".")
	var traj chem.Traj

	switch trjext[len(trjext)-1] { //we will be supporting these 3 formats, althoug we could add pDynamo's crd and Gromacs' xtc
	case "xyz":
		_, traj, err = chem.XYZFileAsTraj(trjname)
	case "pdb":
		traj, err = chem.PDBFileRead(trjname, false)
	case "dcd":
		traj, err = dcd.New(trjname)
	}

	start, err := strconv.Atoi(os.Args[3])
	Err(err)
	skip, err := strconv.Atoi(os.Args[4])
	Err(err)
	XTB := qm.NewXTBHandle()
	basename := strings.Replace(os.Args[1], ".xyz", "_hess", 1)
	Calc := new(qm.Calc)
	Calc.Method = "gfn2"
	Calc.Job = &qm.Job{Forces: true}
	epsilon, err := strconv.Atoi(os.Args[5])
	Calc.Dielectric = float64(epsilon)
	Err(err)
	N := 0
	totalG := 0.0
	totalGsq := 0.0
	coords := v3.Zeros(traj.Len())
	for i := 0; ; i++ {
		if i%skip != 0 || i < start {
			err := traj.Next(nil)
			if err != nil {
				MDErr(err)
				break
			}
			continue
		}
		err := traj.Next(coords)
		if err != nil {
			MDErr(err)
			break //MDErr will kill the program if the error is not EOF, so we can confidently just "break" here
		}
		XTB.SetName(fmt.Sprintf("%s_G_%d", basename, i))
		XTB.BuildInput(coords, mol, Calc)
		XTB.Run(true)
		_, err = XTB.OptimizedGeometry(mol)
		Err(err)
		img, err := XTB.LargestImaginary()
		if err != nil || img > 0.1 { //I'm leaving something for float/rounding issues. Small imaginary modes are probably OK, so you can raise this value if you want.
			log.Printf("Skipping frame %d because of negative frequency. Wavenumber: %5.3f", i, img*-1)
			os.Rename("vibspectrum", fmt.Sprintf("vibspectrum_%d", i))
			continue //we simply skip non-minima.
		}
		G, err := XTB.FreeEnergy()
		Err(err)
		fmt.Printf("G: %5.3f Frame: %d\n", G, i)
		N++
		totalG += G
		totalGsq += (G * G)
	}
	n := float64(N)
	variance := totalGsq/n - math.Pow(totalG/(n), 2)
	stdev := math.Sqrt(variance)
	fmt.Printf("Average G: %5.3f , stdev: %5.3f Frames: %d\n", totalG/float64(N), stdev, N)
}
