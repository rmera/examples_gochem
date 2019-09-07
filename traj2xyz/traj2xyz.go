
package main

import (
	"flag"
//	"fmt"
	"os"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/amberold"
	"github.com/rmera/gochem/dcd"
	"github.com/rmera/gochem/v3"
//	"github.com/rmera/gochem/xtc"
//	"github.com/rmera/scu"
//	"gonum.org/v1/gonum/mat"
//	"math"
	//	"sort"
//	"strconv"
//	"strings"
)

////use:  program [-skip=number -begin=number2] xyzname trajname multixyzname
func main() {
	//The skip options
	skip := flag.Int("skip", 0, "How many frames to skip between reads.")
	begin := flag.Int("begin", 1, "The frame from where to start reading.")
	format := flag.Int("format", 0, "0 for OldAmber (crd, default), 2 for dcd (NAMD)")
	flag.Parse()
	args := flag.Args()
	//	println("SKIP", *skip, *begin, args) ///////////////////////////
	mol, err := chem.XYZFileRead(args[0])
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
		traj, err = amberold.New(args[1],mol.Len(),false)
		if err != nil {
			panic(err.Error())
		}
	case 2:
		traj, err = dcd.New(args[1])
		if err != nil {
			panic(err.Error())
		}
	 case 3:
		 traj, err =  chem.PDBFileRead(args[1], false)
		if err != nil {
			panic(err.Error())
		}
	}
	out, err := os.Create(args[2])
	if err != nil {
		panic(err.Error())
	}
	defer out.Close()
	var coords *v3.Matrix
	lastread := -1
	for i := 0; ; i++ { //infinite loop, we only break out of it by using "break"  //modified for profiling
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
		err = chem.XYZWrite(out, coords, mol)
		if err!=nil{
			panic(err.Error())
		}
		coords = nil // Not sure this works
	}
}

