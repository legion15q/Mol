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


class CreateSmiles(object):
    def __init__(self, mol2_):
        self.Mol2 = mol2_
        self.smiles_ = ''
        self.start_bond_id = self.find_start_atom()
        self.create_smiles(self.start_bond_id)
        print(self.smiles_)


    def find_start_atom(self):
        start = None
        for i in self.Mol2.bonds:
            bond_id = i['bond_id']
            origin = i['origin_atom_id']
            found_start = True
            for j in self.Mol2.bonds:
                target = j['target_atom_id']
                if origin == target:
                    found_start = False
                    continue
            if found_start is True:
                return int(bond_id)


    def calc_count_of_output_bonds(self, origin):
        k = 0
        for j in self.Mol2.bonds:
            cur_origin = j['origin_atom_id']
            if cur_origin == origin:
                k += 1
        return k

    def create_smiles(self, bond_id):
        if bond_id is None:
            print(self.smiles_)
            exit(1)
        start_origin_id = self.Mol2.bonds[bond_id - 1]['origin_atom_id']
        if bond_id == self.start_bond_id:
            start_origin_name = self.Mol2.atoms[int(start_origin_id) - 1]['atom_name']
            self.smiles_ += start_origin_name

        k = self.calc_count_of_output_bonds(start_origin_id)
        while k != 0:
            if k > 1:
                orig_ = self.Mol2.bonds[bond_id - 1]['target_atom_id']
                new_bond_id = self.get_next_bond_id_with_origin(bond_id, orig_)
                if self.Mol2.bonds[bond_id - 1]['bond_type'] == '2':
                    self.smiles_ += '='
                self.smiles_ += '('
                self.smiles_ += self.Mol2.atoms[int(orig_) - 1]['atom_name']
                if new_bond_id is not None:
                    self.create_smiles(new_bond_id)
                else:
                    self.smiles_ += ')'
                k -= 1
            else:
                new_bond_id = self.get_next_bond_id_with_origin(bond_id, start_origin_id)
                if new_bond_id is not None:
                    self.create_smiles(new_bond_id)
                else:
                    k -= 1
                    if self.Mol2.bonds[int(self.start_bond_id) - 1]['origin_atom_id'] != self.Mol2.bonds[int(bond_id) - 1]['origin_atom_id']:
                        self.smiles_ += ')'
                    break
                k -= 1


    def get_next_bond_id_with_origin(self, bond_id, orig_) -> int:
        for i in range(int(bond_id), len(self.Mol2.bonds)):
            new_origin = self.Mol2.bonds[i]['origin_atom_id']
            if orig_ == new_origin:
                return int(self.Mol2.bonds[i]['bond_id'])




def main():
    Mol2 = MolParser('mol.mol2')
    CreateSmiles(Mol2)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
