# чтобы разделить на блоки
import numpy
import math

block_tokens = ('@<TRIPOS>ATOM', '@<TRIPOS>BOND', '@<TRIPOS>SUBSTRUCTURE')

atom_tokens = ('atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id', 'subst_name', 'charge')
bond_tokens = ('bond_id', 'origin_atom_id', 'target_atom_id', 'bond_type')
substructure_tokens = ('subst_id', 'subst_name', 'root_atom', 'subst_type', 'dict_type', 'chain')

block_subtokens = (atom_tokens, bond_tokens, substructure_tokens)


class MolParser(object):
    def __init__(self, file_name_):
        self.file_name = file_name_
        self.text = self.read_file(file_name_)
        self.blocks = self.read_by_blocks()
        self.Mol2 = self.read_by_tokens()

    def read_file(self, file_name_):
        with open(file_name_, 'r') as f:
            return f.read()

    def read_by_blocks(self):
        blocks = []
        for i in range(len(block_tokens)):
            blocks.append(self.text.split(block_tokens[i])[1])
        for i in range(len(block_tokens) - 1):
            blocks[i] = blocks[i].split(block_tokens[i + 1])[0]
        return blocks

    def read_by_tokens(self):
        Mol2 = {}
        for b_k, b_v in enumerate(self.blocks):
            block = []
            lines = b_v.split('\n')
            lines = (list(filter(lambda x: x != '', lines)))
            for j in lines:
                elements = j.split('\t')
                if len(elements) == 1: elements = j.split(' ')
                attributes = {}
                for k, v in enumerate(block_subtokens[b_k]):
                    attributes[v] = elements[k]
                # line__ = {int([*attributes.values()][0]): attributes}
                block.append(attributes)
            Mol2[block_tokens[b_k]] = block
        return Mol2

    @property
    def atoms(self):
        return self.Mol2['@<TRIPOS>ATOM']

    @property
    def bonds(self):
        return self.Mol2['@<TRIPOS>BOND']

    @property
    def substruct(self):
        return self.Mol2['@<TRIPOS>SUBSTRUCTURE']


class Coord(object):
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    @property
    def X(self):
        return self.x

    @property
    def Y(self):
        return self.y

    @property
    def Z(self):
        return self.z


class Energy(object):
    def __init__(self, mol2_):
        self.Mol2 = mol2_
        self.Angles = {}
        self.sum_ = 0

    def calcBondLength(self):
        for i in range(len(self.Mol2.bonds)):
            origin_atom_id = int(self.Mol2.bonds[i]['origin_atom_id'])
            target_atom_id = int(self.Mol2.bonds[i]['target_atom_id'])
            dot_1 = Coord(self.Mol2.atoms[origin_atom_id - 1]['x'], self.Mol2.atoms[origin_atom_id - 1]['y'],
                          self.Mol2.atoms[origin_atom_id - 1]['z'])
            dot_2 = Coord(self.Mol2.atoms[target_atom_id - 1]['x'], self.Mol2.atoms[target_atom_id - 1]['y'],
                          self.Mol2.atoms[target_atom_id - 1]['z'])
            self.Mol2.bonds[i]['bond_length'] = self.Distance(dot_1, dot_2)

    def ChangeNames(self):
        with open('TriposForceField(MMFF94).txt', 'r') as f:
            table = f.read().lower()

        lines = table.split('\n')
        names_map = []
        for i in lines:
            names_map.append(tuple(i.split('=')))
        for i in range(len(self.Mol2.atoms)):
            atom_type = self.Mol2.atoms[i]['atom_type'].lower()
            if atom_type.find(' ') != -1:
                atom_type.remove(' ')
            for j in names_map:
                if j[0] == atom_type:
                    self.Mol2.atoms[i].update({'atom_type': j[1]})

    def CalcANgles(self):

        for i in range(len(self.Mol2.bonds)):
            origin_atom_id = int(self.Mol2.bonds[i]['origin_atom_id'])
            target_atom_id = int(self.Mol2.bonds[i]['target_atom_id'])
            new_target_atom_id = -1
            for j in range(len(self.Mol2.bonds)):
                if j != i:
                    if int(self.Mol2.bonds[j]['origin_atom_id']) == target_atom_id:
                        new_target_atom_id = int(self.Mol2.bonds[j]['target_atom_id'])
                        sum_cur = self.CalcAngle_(origin_atom_id, target_atom_id, new_target_atom_id)
                        self.sum_ += sum_cur
            if new_target_atom_id == -1:
                for j in range(len(self.Mol2.bonds)):
                    if j != i:
                        if int(self.Mol2.bonds[j]['target_atom_id']) == target_atom_id:
                            new_target_atom_id = int(self.Mol2.bonds[j]['origin_atom_id'])
                            sum_cur = self.CalcAngle_(origin_atom_id, target_atom_id, new_target_atom_id)
                            self.sum_ += sum_cur

    def CalcAngle_(self, origin_atom_id, target_atom_id, new_target_atom_id):
        dot_1 = Coord(self.Mol2.atoms[origin_atom_id - 1]['x'], self.Mol2.atoms[origin_atom_id - 1]['y'],
                      self.Mol2.atoms[origin_atom_id - 1]['z'])
        dot_2 = Coord(self.Mol2.atoms[target_atom_id - 1]['x'], self.Mol2.atoms[target_atom_id - 1]['y'],
                      self.Mol2.atoms[target_atom_id - 1]['z'])
        dot_3 = Coord(self.Mol2.atoms[new_target_atom_id - 1]['x'], self.Mol2.atoms[new_target_atom_id - 1]['y'],
                      self.Mol2.atoms[new_target_atom_id - 1]['z'])
        return self.CalcAngle(dot_1, dot_2, dot_3)

    def Distance(self, dot_1, dot_2):
        return math.sqrt(
            math.pow(dot_2.x - dot_1.x, 2) + math.pow(dot_2.y - dot_1.y, 2) + math.pow(dot_2.z - dot_1.z, 2))

    def CalcAngle(self, dot_1, dot_2, dot_3):
        ba = [dot_1.x - dot_2.x, dot_1.y - dot_2.y, dot_1.z - dot_2.z]
        bc = [dot_1.x - dot_3.x, dot_1.y - dot_3.y, dot_1.z - dot_3.z]
        cosine_angle = numpy.dot(ba, bc) / (numpy.linalg.norm(ba) * numpy.linalg.norm(bc))
        angle = numpy.arccos(cosine_angle)
        pAngle = numpy.degrees(angle)
        return pAngle


def main():
    Mol2 = MolParser('mol.mol2')

    E = Energy(Mol2)
    E.CalcANgles()
    E.ChangeNames()
    print(E.sum_)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

