package main

import (
	"flag"
	//	"fmt"
	"os"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/amberold"
	"github.com/rmera/gochem/dcd"
	v3 "github.com/rmera/gochem/v3"

	//	"github.com/rmera/gochem/xtc"
	//	"github.com/rmera/scu"
	//	"gonum.org/v1/gonum/mat"
	//	"math"
	//	"sort"
	//	"strconv"
	"strings"
)

////use:  program [-skip=number -begin=number2] pdbfile trajname
func main() {
	//The skip options
	skip := flag.Int("skip", 0, "How many frames to skip between reads.")
	begin := flag.Int("begin", 1, "The frame from where to start reading.")
	format := flag.Int("format", 0, "0 for OldAmber (crd, default), 2 for dcd (NAMD)")
	end := flag.Int("end", 100000, "The last frame")
	flag.Parse()
	args := flag.Args()
	//	println("SKIP", *skip, *begin, args) ///////////////////////////
	mol, err := chem.PDBFileRead(args[0], false)
	if err != nil {
		panic(err.Error())
	}
	var traj chem.Traj
	switch *format {
	//	case 0:
	//		traj, err = xtc.New(args[2])
	//		if err != nil {
	//			panic(err.Error())
	//		}
	case 0:
		traj, err = amberold.New(args[1], mol.Len(), false)
		if err != nil {
			panic(err.Error())
		}
	case 2:
		traj, err = dcd.New(args[1])
		if err != nil {
			panic(err.Error())
		}
	case 3:
		traj, err = chem.PDBFileRead(args[1], false)
		if err != nil {
			panic(err.Error())
		}
	}

	Coords := make([]*v3.Matrix, 0, 0)
	var coords *v3.Matrix
	lastread := -1
	for i := 0; i < *end; i++ { //infinite loop, we only break out of it by using "break"  //modified for profiling
		if lastread < 0 || (i >= lastread+(*skip) && i >= (*begin)-1) {
			coords = v3.Zeros(traj.Len())
		}
		err := traj.Next(coords) //Obtain the next frame of the trajectory.
		if err != nil {
			_, ok := err.(chem.LastFrameError)
			if ok {
				break //We processed all frames and are ready, not a real error.

			} else {
				panic(err.Error())
			}
		}
		if (lastread >= 0 && i < lastread+(*skip)) || i < (*begin)-1 { //not so nice check for this twice
			continue
		}
		lastread = i
		Coords = append(Coords, coords)
		coords = nil // Not sure this works
	}
	pdbname := strings.Replace(args[1], ".crd", ".pdb", 1)
	fout, err := os.Create(pdbname)
	if err != nil {
		panic(err.Error())
	}
	defer fout.Close()
	err = chem.MultiPDBWrite(fout, Coords, mol, nil)
	if err != nil {
		panic("Couldn't write multipdb: " + err.Error())
	}
}
