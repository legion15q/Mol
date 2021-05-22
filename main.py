# чтобы разделить на блоки
import numpy

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
                line_ = {}
                for k, v in enumerate(block_subtokens[b_k]):
                    line_[v] = elements[k]
                block.append(line_)
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


def main():
    obj = MolParser('mol.mol2')
    print(obj.atoms)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
