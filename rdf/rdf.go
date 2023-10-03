/*
 * rdf.go

 * Copyright 2020 Raul Mera A. (raulpuntomeraatusachpuntocl)
 *
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 2.1 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *
*/
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package main

import (
	"flag"
	"fmt"
	"log"
	"time"

	chem "github.com/rmera/gochem"

	"github.com/rmera/gochem/solv"
	"github.com/rmera/gochem/traj/xtc"
	v3 "github.com/rmera/gochem/v3"
)

// lazy wrapping function for critical errors.
func Error(err error) {
	if err != nil {
		panic(err.Error())
	}
}

// Justa silly little error wrapping for non critical errors
func Warning(err error) {
	if err != nil {
		fmt.Println("Warning!,", err.Error())
	}
}

func main() {
	getsystem := flag.Bool("getsystem", false, "No MDDF/RDF is calculated. Instead, a subsystem with the solute and solvent up to _end_ A is produced for each processed frame of the traj")
	com := flag.Bool("com", false, " Use the centre of mass (if available, centroid otherwise) for obtaining the distance to a solvent molecule, instead of the closest  atom")
	cpus := flag.Int("cpus", -1, "How many gorutines to use. If <0, as many as the logical CPUs will be used.")
	skip := flag.Int("skip", 1, "Every how many frames should a frame be read.")
	step := flag.Float64("step", 0.1, "The interval for the MDDF calculation, in A")
	end := flag.Float64("end", 10, "The larger distance considered")
	refatom := flag.Int("refatom", -1, "Index of the atom to be used as reference, counting from 0. (optional, a full residue ID can be given instead, see option -refid)")
	refID := flag.Int("refid", -1, "Index of the residue to be used as reference, as it appears in the PDB reference")
	flag.Parse()
	args := flag.Args()
	var start, fin time.Time
	var elapsed time.Duration
	trajname := args[1]
	pdbname := args[0]
	solname := []string{args[2]}
	var solute []int
	//	Error(err)
	mol, err := chem.PDBFileRead(pdbname, false)
	if *refID >= 0 {
		solute = chem.Molecules2Atoms(mol, []int{*refID}, []string{})
	} else if *refatom >= 0 {
		solute = []int{*refatom}
	} else {
		log.Fatal("Either -refid or -refatom needs to be specified")
	}
	if *skip < 1 {
		*skip = 1
	}
	Error(err)
	fmt.Println("skip", *skip, "step", *step, "end", *end, "refatom", *refatom, "refid", *refID, "trajname", trajname, "pdbname", pdbname, "solname", solname, "solute", solute, args) ////////////////////////////////

	options := solv.DefaultOptions()
	options.Step(*step)
	options.End(*end)
	options.COM(*com)
	options.Cpus(*cpus)

	traj1, err := xtc.New(trajname)
	if *getsystem {
		GetSystem(traj1, mol, solute, solname, options)
		return

	}
	Error(err)
	start = time.Now()
	rdf1, solv1, err := solv.ConcMolRDF(traj1, mol, solute, solname, options)
	Error(err)
	fin = time.Now()
	elapsed = fin.Sub(start)
	for i, v := range rdf1 {
		fmt.Printf("%7.4f\t%7.4f\n", v, solv1[i])
	}
	fmt.Println(elapsed)

}

// GetSystem, For each read frame, creates XYZ file with the solute and the solvent molecues(molecules with residue names contained in residues) within end A.
// Either any atom from a solvent molecule, or the COM is used to define whether a molecule iswhithin the range, or not.
func GetSystem(traj chem.Traj, mol chem.Atomer, refindexes []int, residues []string, options *solv.Options) {
	topol := chem.NewTopology(0, 1)
	coords := v3.Zeros(mol.Len())
	var err error
reading:

	for i := 0; ; i++ {
		//	fmt.Println(refindexes)
		if i > 0 && i%options.Skip() != 0 && err == nil {
			err = traj.Next(nil) //if this err is not nil, the next traj.Next() will not be excecuted, whether it's a skip or a read. Instead, we'll go directly to error processing.
			continue
		} else if err == nil { //in case the frame we skipped before gave an error
			err = traj.Next(coords)
		}
		if err != nil {
			switch err := err.(type) {
			default:
				panic(err.Error()) ///what can you do.
			case chem.LastFrameError:
				break reading
			}
		}
		whithin := solv.DistRank(coords, mol, refindexes, residues, options)
		indexes := make([]int, 0, len(refindexes)+whithin.Len())
		indexes = append(indexes, refindexes...)
		fmt.Println(indexes)
		//	indexes = append(indexes, whithin.AtomIDs(mol)...)
		basename := fmt.Sprintf("System_%dA_Frame%d_Natoms%d", int(options.End()), i, len(indexes))
		chem.PDBFileWrite("Full"+basename+".pdb", coords, mol, nil)
		topol.SomeAtoms(mol, indexes)
		c := v3.Zeros(len(indexes))
		c.SomeVecs(coords, indexes)
		fmt.Printf("%d\t%d\n", i, len(indexes))
		chem.XYZFileWrite(basename+".xyz", c, topol) //If there was an error we just moveto the next one.

	}
	return

}
