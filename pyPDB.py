#!/usr/bin/python

import os
import sys

class Atom(object):

    """Atom Class"""

    def __init__(self, rectype="ATOM", id=-1, name="", altLoc=" ", resName=" ",
                 chainID=-1, resSeq=-1, iCode=" ",
                 x=0, y=0, z=0, occupancy=" ", tempFactor=" ", charge=" "):
        self.rectype = rectype
        self.id = id
        self.name = name
        self.altLoc = altLoc
        self.resName = resName
        self.chainID = chainID
        self.resSeq = resSeq
        self.iCode = iCode
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.tempFactor = tempFactor
        self.charge = charge


class Bond(object):

    """Bond Class"""

    def __init__(self, atom1=0, atom2=0):
        self.atom1 = atom1
        self.atom2 = atom2


class Residue(object):

    """Residue Class"""

    def __init__(self, id=-1, name="", atoms=None, chain=""):
        self.id = id
        self.name = name
        if atoms == None:
            self.atoms = []
        else:
            self.atoms = atoms
        self.chain = chain


class Chain(object):

    """Chain Class"""

    def __init__(self, id=-1, name="", residues=None):
        self.id = id
        self.name = name
        if residues == None:
            self.residues = []
        else:
            self.residues = residues


class Molecule(object):

    """Molecule Class"""

    def __init__(self, id=0, name="", atoms=None, bonds=None, residues=None, chains=None):
        self.id = id
        self.name = name

        if atoms == None:
            self.atoms = {}
        else:
            self.atoms = atoms

        if bonds == None:
            self.bonds = []
        else:
            self.bonds = bonds

        if residues == None:
            self.residues = {}
        else:
            self.residues = residues

        if chains == None:
            self.chains = []
        else:
            self.chains = chains

    def residue_total(self):
        return len(self.residues)

    def atom_total(self):
        return len(self.atoms)

    def bond_total(self):
        return len(self.bonds)

    def chain_total(self):
        return len(self.chains)


class pyPDB(object):

    """PDB Class"""

    def __init__(self, filename):
        self.filename = filename
        self.molecule = None
        self.selectedAtoms = []
        self.reduced = []
        self._readFile()
        self.verbose = False

    def _readFile(self):
        m = Molecule(self.filename)
        m.name = os.path.splitext(self.filename)[0].lower()

        f = open(self.filename, 'r').read().replace('\r\n', '\n')

        l = 0
        temp_chain = []
        chain_no = 1
        for line in f.splitlines():
            l += 1

            if (line[0:4] == 'ATOM' or line[0:6] == 'HETATM'):
                # get atom information
                atom = self._readAtom(line)
                # add atom to molecule atoms
                m.atoms[atom.id] = atom
                if atom.resSeq not in m.residues.keys():
                    # new residue
                    r = Residue()
                    r.id = atom.resSeq
                    r.name = atom.resName
                    r.atoms = [atom.id]
                    chain_name = line[21:22]
                    r.chain = chain_name
                    m.residues[r.id] = r
                    temp_chain.append(r)

                else:
                    # new atom to residue
                    m.residues[atom.resSeq].atoms.append(atom)

            if line[0:6] == 'CONECT':
                bonds_in_line = self._readBonds(line)
                for bond in bonds_in_line:
                    m.bonds.append(bond)

            if 'TER' in line:
                c = Chain()
                c.name = line[21:22]
                c.residues = temp_chain
                c.id = chain_no
                m.chains.append(c)
                temp_chain = []
                chain_no = chain_no + 1

        if m.bond_total() == 0:
            print 'Warning: No CONECT info, so no bond analysis.'

        if 'TER' not in f:
            print 'Warning: No TER statement, so no chains are built.\n'

        self.molecule = m

    def _readAtom(self, line):
        a = Atom()
        a.rectype = line[0:6]  # ATOM or HETATM
        iid = line[7:11].strip()
        a.id = int(iid)
        a.name = line[12:14].strip()
        a.altLoc = line[16]
        a.resName = line[17:20]
        a.chainID = line[21]
        a.resSeq = int(line[22:26])
        a.iCode = line[26]
        a.x = float(line[30:37])
        a.y = float(line[38:45])
        a.z = float(line[46:53])
        a.occupancy = line[54:59]
        a.tempFactor = line[60:65]
        a.charge = line[78:89]
        return a

    def _readBonds(self, line):
        fields = line.split()
        bonds = []
        n = 2
        while n < len(fields):
            bond = Bond()
            bond.atom1 = int(fields[1])
            bond.atom2 = int(fields[n])
            bonds.append(bond)
            n += 1

        return bonds

    def distanceBetweenAtoms(self, atomid1, atomid2):
        import numpy
        atom1 = self.molecule.atoms[atomid1]
        atom2 = self.molecule.atoms[atomid2]

        a = numpy.array((atom1.x, atom1.y, atom1.z))
        b = numpy.array((atom2.x, atom2.y, atom2.z))
        dist = numpy.linalg.norm(a - b)

        return int(dist * 100) / 100.00

    def atomsWithinDistanceOfAtom(self, atomid, distance):

        referenceAtom = self.molecule.atoms[atomid]

        atomsWithinDistance = []
        atomDistances = []
        self.selectedAtoms = []
        for key in self.molecule.atoms:
            if self.distanceBetweenAtoms(atomid, self.molecule.atoms[key].id) <= distance:
                if self.molecule.atoms[key].id != atomid:
                    atomsWithinDistance.append(self.molecule.atoms[key])
                    d = self.distanceBetweenAtoms(
                        atomid, self.molecule.atoms[key].id)
                    atomDistances.append(d)
                    self.selectedAtoms.append(self.molecule.atoms[key])

        return (atomsWithinDistance, atomDistances)

    def toJSON(self):
        ret = '{ \n'
        ret += '\t "atom_total": {0},\n'.format(self.molecule.atom_total())
        ret += '\t "residue_total": {0},\n'.format(self.molecule.residue_total())
        ret += '\t "bond_total": {0}'.format(self.molecule.bond_total())
        ret += '\n}'
        return ret

    def distanceMap(self):
        n1 = 0
        dist_map = []
        for atom in self.molecule.atoms:
            atom1 = self.molecule.atoms[atom]
            temp_distances = []

            for a2 in self.molecule.atoms:
                atom2 = self.molecule.atoms[a2]
                temp_distances.append(
                    self.distanceBetweenAtoms(atom1.id, atom2.id))

            dist_map.append(temp_distances)

        return dist_map

    def plotDistanceMap(self, save=False, directory='', close=False):
        import numpy
        import matplotlib.pyplot as plt

        m = self.distanceMap()
        matrix = numpy.matrix(m)

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect('equal')
        plt.title('Distance Map')
        extent = self.molecule.atom_total() + 0.5
        plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.hot,
                   extent=(0.5, extent, 0.5, extent))
        plt.colorbar()

        if save == True:
            plt.savefig('{}distance_map.pdf'.format(directory))

        if close == True:
            plt.close()
        else:
            plt.show()

    def selectAtom(self, atomid):
        atom = self.molecule.atoms[atomid]

        alreadySelected = False
        for atoms in self.selectedAtoms:
            if atomid == atoms.id:
                alreadySelected = True

        if alreadySelected == False:
            self.selectedAtoms.append(atom)

        return self  # enables chaining

    def selectAtoms(self, atomids=[]):
        for atom in atomids:
            alreadySelected = False
            for atoms in self.selectedAtoms:
                if atom == atoms.id:
                    alreadySelected = True

            if alreadySelected == False:
                self.selectedAtoms.append(self.molecule.atoms[atom])

    def reduce(self):
        for atom in self.molecule.atoms:
            if 'H' not in self.molecule.atoms[atom].name:
                self.reduced.append(self.molecule.atoms[atom])

        return self.reduced

    def listResiduesFromAtoms(self, atoms):
        residues = []
        for atom in atoms:
            if atom not in residues:
                residues.append(self.molecule.residues[atom.residue_id])

        temp_residue_list = []
        for residue in residues:
            if residue.id not in temp_residue_list:
                temp_residue_list.append(residue.id)

        return temp_residue_list

    def toAmberMask(self, key='residues'):
        ret = ''
        i = 1

        if(key == 'residues'):
            for residue in self.listResiduesFromAtoms(self.selectedAtoms):
                if i == len(self.listResiduesFromAtoms(self.selectedAtoms)):
                    comma = ''
                else:
                    comma = ','

                ret += '{0}{1}'.format(residue, comma)
                i = i + 1
            return ret

        elif(key == 'atoms'):
            for atom in self.selectedAtoms:
                if i == len(self.selectedAtoms):
                    comma = ''
                else:
                    comma = ','

                ret += '{0}{1}'.format(atom.id, comma)
                i = i + 1
            return ret

    def removeSelection(self):
        self.selectedAtoms = []
        return self

    def writePDB(self):
        if not self.selectedAtoms:
            atomsToWrite = self.molecule.atoms
        else:
            atomsToWrite = []
            for atom in self.selectedAtoms:
                atomsToWrite.append(atom.id)

        for atom in sorted(atomsToWrite, key=lambda k: k):
            a = self.molecule.atoms[atom]
            print self._get_atom_line(a)

    def _get_atom_line(self, a):

        # COLUMNS        DATA  TYPE    FIELD        DEFINITION
        # ---------------------------------------------------------------------
        #  0 -  5        Record name   "ATOM  "
        #  6 - 10        Integer       serial       Atom  serial number.
        # 12 - 15        Atom          name         Atom name.
        # 16             Character     altLoc       Alternate location indicator.
        # 17 - 19        Residue name  resName      Residue name.
        # 21             Character     chainID      Chain identifier.
        # 22 - 25        Integer       resSeq       Residue sequence number.
        # 26             AChar         iCode        Code for insertion of residues.
        # 30 - 37        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        # 38 - 45        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        # 46 - 53        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        # 54 - 59        Real(6.2)     occupancy    Occupancy.
        # 60 - 65        Real(6.2)     tempFactor   Temperature  factor.
        # 76 - 77        LString(2)    element      Element symbol, right-justified.
        # 78 - 89        LString(2)    charge       Charge  on the atom.

        #          1         2         3         4         5         6         7         8
        # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
        # MODEL        1
        # ATOM      1  N   LEU A  25      80.669  55.349  53.905  1.00 39.12           N
        # ATOM      1  N   LEU A   2      80.660  55.340  53.900  1.0

        args = (a.rectype, a.id, a.name, a.altLoc, a.resName,
                a.chainID, a.resSeq, a.iCode,
                a.x, a.y, a.z, a.occupancy)

        return "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s" % args

    def translateCoordinates(self, translateVector):
        if not self.selectedAtoms:
            atoms_list = self.molecule.atoms
        else:
            atoms_list = []
            for atom in self.selectedAtoms:
                atoms_list.append(atom.id)

        for atom in atoms_list:
            a = self.molecule.atoms[atom]
            a1 = numpy.array([a.x, a.y, a.z])
            a2 = numpy.array(translateVector)
            s = numpy.add(a1, a2)
            self.molecule.atoms[atom].x = s[0]
            self.molecule.atoms[atom].y = s[1]
            self.molecule.atoms[atom].z = s[2]

        return self

    def reorderResidues(self):
        counter = 1
        _log(self.verbose, 'Reordering Residues:', colour="red", bold=True)
        for i in self.molecule.residues:
            res = self.molecule.residues[i]
            _log(self.verbose, '{}-{} ----> {}'.format(res.name, res.id, counter))
            res.id = counter
            counter = counter + 1
        return self


    def describeResidues(self):
        description = "{}\n----------------------------------\n".format(
            os.path.basename(self.molecule.name))
        for i in self.molecule.residues:
            description += "Residue ID: {:3g} -- Residue Name: {} -- Chain ID: {}\n".format(
                self.molecule.residues[i].id, self.molecule.residues[i].name, self.molecule.residues[i].chain)
        return description


def mergePDBs(pdbs, output):
    # TODO: Deal with two residues having same ID but are on different chains
    with open(output, 'w') as outfile:
        for fname in pdbs:
            with open(fname) as infile:
                for line in infile:
                    if line[0:4] == 'ATOM':
                        outfile.write(line)
            outfile.write("TER\n")

    return output

def _log(verbose=True, message="", colour=None, background=None, bold=False, underline=False, inverted=False, run=False):

    if verbose:

        colours = {
            'black':    '90',
            'red':      '91',
            'green':    '92',
            'yellow':   '93',
            'blue':     '94',
            'magenta':  '95',
            'cyan':     '96',
            'white':    '97'
        }

        backgrounds = {
            'default':  '49',
            'black':    '100',
            'red':      '101',
            'green':    '102',
            'yellow':   '103',
            'blue':     '104',
            'magenta':  '105',
            'cyan':     '106',
            'white':    '107'
        }

        if bold:                   message = '\033[1m' + message + '\033[21m'
        if underline:              message = '\033[4m' + message + '\033[24m'
        if background is not None: message = '\033[' + backgrounds[background] + 'm' + message + '\033[49m'
        if colour is not None:     message = '\033[' + colours[colour] + 'm' + message + '\033[0m'
        if inverted:               message = '\033[7m' + message + '\033[27m'

        if run:
            print message,
        else:
            print message

    return

def unZip(archive, uncompressed):
    import gzip
    f = gzip.open(archive, 'r')
    g = open(uncompressed, 'w')
    g.writelines(f.readlines())
    f.close()
    g.close()
    os.remove(archive)
    return True

def downloadPDB(pdbCode, output=""):
    import urllib

    pdb = "{pdbid}.pdb.gz".format(pdbid=pdbCode)
    url = "http://www.rcsb.org/pdb/files/{pdb}".format(pdb=pdb)

    urllib.urlretrieve(url, pdb)

    if output == "":
        output_path = "{pdbid}.pdb".format(pdbid=pdbCode)
    else:
        output_path = output

    unZip(pdb, output_path)

    return output_path


if __name__ == '__main__':

# load pdb
    #p = pyPDB('pdbs/gly.pdb')

# if you need to download the pdb, you can load it straight away
    #p2 = pyPDB(downloadPDB('1P47', 'pdbs/1P47.pdb'))

# we can merge two pdb files:
    p3 = pyPDB('pdbs/multiple.pdbqt')
    p3.verbose = True

# after merging, we probably need to reorder the residues:
#    p3.reorderResidues()

# and also describe the residues:
    # p3.describeResidues()

# translate the coordinates (or selection)
# p.translateCoordinates([10,5,1])
# for atom in p.molecule.atoms:
#    a = p.molecule.atoms[atom]
#    print "[{0:g}, {1:g}, {2:g}]".format(a.x, a.y, a.z)

# select one atom
#   p.selectAtom(4)

# select multiple atoms individually (this continues after the previous one)
#   p.selectAtom(5).selectAtom(6)

# or select multiple atoms all in one go
#   p.selectAtoms([4, 5, 6])

# the 'p' pyPDB instance now has a selectedAtoms attribute that is iterable:
#   for atom in p.selectedAtoms:
#       print '{}{}'.format(atom.name, atom.id)

# calculate a distance map
#   print p.distanceMap()

# and also plot it
#   p.plotDistanceMap(save=False, close=True)

# calculate the distance between two atoms
#   print p.distanceBetweenAtoms(8, 9)

# calculate atoms within a given distance of another atom
#   print p.atomsWithinDistanceOfAtom(10, 3)

# you can iterate over something like the above such as:
#   atomsWithinDistance = p.atomsWithinDistanceOfAtom(10, 3)
#   i = 0
#   for x in atomsWithinDistance[0]:
#       print 'Atom {}{} is within {} of {}{}: {}'.format(x.name, x.id, 3,
#           p.molecule.atoms[10].name, 10, atomsWithinDistance[1][i])
#       i += 1

# or even make an amber mask:
#   print p.toAmberMask('atoms')

# output a description of 'p' as json
#   print p.toJSON()

# reduce a pdb:
#   p.reduce()

# ...which can be iterated over:
#   for atom in p.reduce():
#       print '{}{}'.format(atom.name, atom.id)

# the selection (or all atoms if no selection) can be written to a pdb file
#   p.writePDB()

# the selection can be removed using
#   p.removeSelection()
