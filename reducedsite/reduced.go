package main

import (
	"fmt"
	"log"
	"os"
	"slices"
	"strings"

	chem "github.com/rmera/gochem"

	"github.com/rmera/gochem/qm"
	"github.com/rmera/gochem/traj/xtc"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

func LastFrame(errt error) error {
	if errt != nil {
		if _, ok := errt.(chem.LastFrameError); ok {
			return errt
		} else {
			log.Fatal("Error reading trajectory: ", errt.Error())
		}
	}
	return nil
}

func ParseIndexesFile(fname string) [][]int {
	ret := make([][]int, 0, 3)
	f, err := scu.NewMustReadFile(fname)
	scu.QErr(err)
	defer f.Close()
	for l := f.Next(); l != "EOF"; l = f.Next() {
		s, err := scu.IndexStringParse(strings.Replace(l, "\n", "", -1))
		scu.QErr(err)
		ret = append(ret, s)
	}
	return ret
}

func main() {
	mol, err := chem.PDBFileRead(os.Args[1], false)
	scu.QErr(err)
	subseqs := ParseIndexesFile(os.Args[3])
	ch := os.Args[4]
	chains := make([]string, 0, len(subseqs))
	for range subseqs {
		chains = append(chains, ch)
	}
	indexes, err := chem.CutBackRef(mol, chains, subseqs)
	scu.QErr(err)
	qmc := v3.Zeros(len(indexes))
	qmm := chem.NewTopology(1, 1)
	qmm.SomeAtoms(mol, indexes)
	qmc.SomeVecs(mol.Coords[0], indexes)
	basename := strings.Replace(os.Args[1], ".xyz", "", -1)
	chem.XYZFileWrite(basename+"first.xyz", qmc, qmm)
	fixed := make([]int, 0, 80)
	fixednames := []string{"CTZ", "NTZ", "C", "CA", "N"}
	//we will optimize only the Hs
	for i := 0; i < len(qmm.Atoms); i++ {
		at := qmm.Atom(i)
		if slices.Contains(fixednames, at.Name) {
			fixed = append(fixed, i)
		}
	}
	calc := &qm.Calc{}
	calc.Method = "gfn2"
	calc.CConstraints = fixed
	calc.Job = &qm.Job{Opti: true}
	xtb := qm.NewXTBHandle()
	coord := mol.Coords[0]
	//Show time!
	traj, err := xtc.New(os.Args[2])
	scu.QErr(err)
	skip := 10
	var errt error
	for i := 0; ; i++ {
		if i%skip == 0 {
			errt = traj.Next(coord)
			if LastFrame(errt) != nil {
				break
			}
		} else {
			errt = traj.Next(nil)
			if LastFrame(errt) != nil {
				break
			}
			continue
		}
		qmc.SomeVecs(coord, indexes)
		xtb.SetName(fmt.Sprintf("%s%d", basename, i))
		err := xtb.BuildInput(qmc, qmm, calc)
		scu.QErr(err)
		err = xtb.Run(true)
		scu.QErr(err)
		energy, err := xtb.Energy()
		scu.QErr(err)
		geo, err := xtb.OptimizedGeometry(qmm)
		scu.QErr(err)
		chem.XYZFileWrite(fmt.Sprintf("%s%d.xyz", basename, i), geo, qmm)
		fmt.Printf("Frame %d, Energy: %5.3f\n", i, energy)
	}
}
