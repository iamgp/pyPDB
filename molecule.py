class Atom(object):
    """Atom Class"""
    def __init__(self, id=-1, element="", coords=None, residue_id=-1, residue_name=""):
        self.id = id
        self.element = element
        self.residue_id = residue_id
        self.residue_name = residue_name
        if coords == None:
            self.coords = [0,0,0]
        else:
            self.coords = coords

class Bond(object):
    """Bond Class"""
    def __init__(self, atom1=0, atom2=0):
        self.atom1 = atom1
        self.atom2 = atom2

class Residue(object):
    """Residue Class"""
    def __init__(self, id=-1, name="", atoms=None):
        self.id = id
        self.name = name
        if atoms == None:
            self.atoms = []
        else:
            self.atoms = atoms

class Molecule(object):
    """Molecule Class"""
    def __init__(self, id=0, name="", atoms=None, bonds=None, residues=None):
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

    def residue_total(self):
        return len(self.residues)

    def atom_total(self):
        return len(self.atoms)

    def bond_total(self):
        return len(self.bonds)
