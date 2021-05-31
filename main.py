import numpy
from scipy.optimize import minimize, brute
import math
from energy import bs_energy, bs_energy_K_r_mean, bs_energy_req_mean, bending_energy, bending_energy_K_q_mean, \
    bending_energy_qeq_mean, vdw_params, vdw_params_ek_mean, vdw_params_Rj_mean

# чтобы разделить на блоки
block_tokens = ('@<TRIPOS>ATOM', '@<TRIPOS>BOND', '@<TRIPOS>SUBSTRUCTURE')

atom_tokens = ('atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id', 'subst_name', 'charge')
bond_tokens = ('bond_id', 'origin_atom_id', 'target_atom_id', 'bond_type')
substructure_tokens = ('subst_id', 'subst_name', 'root_atom', 'subst_type', 'dict_type', 'chain')

block_subtokens = (atom_tokens, bond_tokens, substructure_tokens)

eps_0 = 8.85418781762e-12
k_clmb = 1 / (4 * math.pi * eps_0)
k_clmb = 8.9875517873681764 * (10 ** 9)


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


def Distance(dot_1, dot_2):
    return math.sqrt(
        math.pow(dot_2.x - dot_1.x, 2) + math.pow(dot_2.y - dot_1.y, 2) + math.pow(dot_2.z - dot_1.z, 2))


def CalcAngle(dot_1, dot_2, dot_3):
    a = [dot_1.x, dot_1.y, dot_1.z]
    b = [dot_2.x, dot_2.y, dot_2.z]
    c = [dot_3.x, dot_3.y, dot_3.z]
    ba = [aa - bb for aa, bb in zip(a, b)]
    bc = [cc - bb for cc, bb in zip(c, b)]

    nba = math.sqrt(sum((x ** 2.0 for x in ba)))
    ba = [x / nba for x in ba]

    nbc = math.sqrt(sum((x ** 2.0 for x in bc)))
    bc = [x / nbc for x in bc]

    scale = sum((aa * bb for aa, bb in zip(ba, bc)))

    angle = math.acos(scale)
    return math.degrees(angle)


def foo():
    foo.counter += 1


foo.counter = 0


class Energy(object):
    def __init__(self, mol2_):
        self.Mol2 = mol2_
        self.ChangeNames()
        self.unic_chains = []
        self.lst_of_unic_chains = []
        self.map_unic_chains_to_angles = []
        self.map_unic_chains_to_r = []
        self.map_unic_long_chains = []

    def ChangeNames(self):
        with open('TriposForceField(MMFF94).txt', 'r') as f:
            table = f.read()

        lines = table.split('\n')
        names_map = []
        for i in lines:
            names_map.append(tuple(i.split('=')))
        for i in range(len(self.Mol2.atoms)):
            atom_type = self.Mol2.atoms[i]['atom_type']
            if atom_type.find(' ') != -1:
                ''.join(atom_type.split(' '))
            for j in names_map:
                if j[0] == atom_type:
                    self.Mol2.atoms[i].update({'atom_type': j[1]})

    def calc_all_long_chains(self):
        self.lst_of_unic_chains = []
        bonds_tuples_lst = []
        self.map_unic_long_chains = []
        for i in range(len(self.Mol2.bonds)):
            origin_atom_id = int(self.Mol2.bonds[i]['origin_atom_id'])
            target_atom_id = int(self.Mol2.bonds[i]['target_atom_id'])
            temp_tuple = (origin_atom_id, target_atom_id)
            bonds_tuples_lst.append(temp_tuple)
        #   for k in range(len(bonds_tuples_lst)):
        k = 0
        temp_lst = [0]
        unic_chains = self.unic_chains
        while len(temp_lst) != 0:
            temp_lst = []
            for i in range(len(unic_chains)):
                chain_set = set(unic_chains[i])
                for j in range(len(bonds_tuples_lst)):
                    bond_set = set(bonds_tuples_lst[j])
                    if len(chain_set.intersection(bond_set)) == 1:
                        if len({unic_chains[i][0]}.intersection(bond_set)) == 1:
                            temp_lst.append(
                                tuple(bond_set - {unic_chains[i][0]}) + unic_chains[i])
                        if len({unic_chains[i][-1]}.intersection(bond_set)) == 1:
                            temp_lst.append(
                                unic_chains[i] + tuple(bond_set - {unic_chains[i][-1]}))
            for i in temp_lst:
                for j in temp_lst:
                    if temp_lst.count(j) > 1:
                        temp_lst.remove(j)
                    elif i != j:
                        if set(j) == set(i):
                            temp_lst.remove(j)
            self.map_unic_long_chains.append(temp_lst)
            unic_chains = self.map_unic_long_chains[k]
            k += 1

        ends_of_triplets = []
        ends_of_duplets = []
        for i in self.map_unic_chains_to_r:
            ends_of_duplets.append((list(i.keys())[0][0], list(i.keys())[0][-1]))

        for i in self.unic_chains:
            ends_of_triplets.append((i[0], i[-1]))

        for i in self.map_unic_long_chains:
            for j in i:
                self.lst_of_unic_chains.append(j)

        for i in self.lst_of_unic_chains:
            for j in self.lst_of_unic_chains:
                if ({i[0], i[-1]} == {j[0], j[-1]}) and (len(i) <= len(j)) and (i != j):
                    self.lst_of_unic_chains.remove(j)

        for j in self.lst_of_unic_chains:
            for i in self.lst_of_unic_chains:
                for k in ends_of_triplets:
                    if {i[0], i[-1]} == set(k):
                        self.lst_of_unic_chains.remove(i)
                for k in ends_of_duplets:
                    if {i[0], i[-1]} == set(k):
                        self.lst_of_unic_chains.remove(i)

    def calc_vdw_coloubm_energy_energy(self, x=None):

        self.calc_all_long_chains()
        vdw_coloubm_energy = 0
        ind = 0
        for j in self.lst_of_unic_chains:

            dot_1 = Coord(self.Mol2.atoms[j[0] - 1]['x'], self.Mol2.atoms[j[0] - 1]['y'],
                          self.Mol2.atoms[j[0] - 1]['z'])
            dot_2 = Coord(self.Mol2.atoms[j[-1] - 1]['x'], self.Mol2.atoms[j[-1] - 1]['y'],
                          self.Mol2.atoms[j[-1] - 1]['z'])
            name1 = self.Mol2.atoms[j[0] - 1]['atom_type']
            name2 = self.Mol2.atoms[j[-1] - 1]['atom_type']

            R_i, e_i = vdw_params_Rj_mean, vdw_params_ek_mean
            if vdw_params.get(tuple([name1])) is not None:
                coef = vdw_params.get(tuple([name1]))
                R_i, e_i = coef['R*j'], coef['e_k']
            R_j, e_j = vdw_params_Rj_mean, vdw_params_ek_mean
            if vdw_params.get(tuple([name2])) is not None:
                coef = vdw_params.get(tuple([name2]))
                R_j, e_j = coef['R*j'], coef['e_k']
            R_ij = R_i + R_j
            e_ij = math.sqrt(e_i * e_j)
            A_ij = e_ij * (R_ij ** 12)
            B_ij = 2 * e_ij * (R_ij ** 6)
            if x is None:
                R = Distance(dot_1, dot_2)
            else:
                R = x[ind]
                ind += 1

            if len(j) == 4:
                vdw_coloubm_energy += 0.5 * A_ij / (R ** 12) - 0.5 * B_ij / (R ** 6) + \
                                      0.83 * k_clmb * (10 ** -9) * float(self.Mol2.atoms[j[0] - 1]['charge']) * float(
                    self.Mol2.atoms[j[-1] - 1][
                        'charge']) / R
            else:
                vdw_coloubm_energy += A_ij / (R ** 12) - B_ij / (R ** 6) + \
                                      k_clmb * (10 ** -9) * \
                                      float(self.Mol2.atoms[j[0] - 1]['charge']) * float(self.Mol2.atoms[j[-1] - 1][
                                                                                             'charge']) / R

        return vdw_coloubm_energy

    def calcBondEnergy(self, x=None):
        self.map_unic_chains_to_r = []
        bse = 0
        ind = 0
        N = len(self.Mol2.atoms)
        if x is not None:
            with open('trajectory.txt', 'a') as f:
                f.write(str(N) + '\n')
                f.write(str(foo.counter) + '\n')
                ind = 0
                for m in range(N):
                    atom_id = ind / 3
                    f.write(str(self.Mol2.atoms[int(atom_id)]['atom_name']) + ' ' + ' '.join(
                        [str(x[ind + m]) for m in range(0, 3)]) + '\n')
                    ind += 3
            f.close()
            foo()

        for i in range(len(self.Mol2.bonds)):
            origin_atom_id = int(self.Mol2.bonds[i]['origin_atom_id'])
            target_atom_id = int(self.Mol2.bonds[i]['target_atom_id'])
            dot_1 = Coord(self.Mol2.atoms[origin_atom_id - 1]['x'], self.Mol2.atoms[origin_atom_id - 1]['y'],
                          self.Mol2.atoms[origin_atom_id - 1]['z'])
            dot_2 = Coord(self.Mol2.atoms[target_atom_id - 1]['x'], self.Mol2.atoms[target_atom_id - 1]['y'],
                          self.Mol2.atoms[target_atom_id - 1]['z'])
            if x is None:
                r = Distance(dot_1, dot_2)
            else:
                or_atom_id = (origin_atom_id - 1) * 3
                tg_atom_id = (target_atom_id - 1) * 3
                r = Distance(Coord(x[or_atom_id], x[or_atom_id + 1], x[or_atom_id + 2]),
                             Coord(x[tg_atom_id], x[tg_atom_id + 1], x[tg_atom_id + 2]))
                ind += 1

            self.Mol2.bonds[i]['bond_length'] = r
            name1 = self.Mol2.atoms[origin_atom_id - 1]['atom_type']
            name2 = self.Mol2.atoms[target_atom_id - 1]['atom_type']
            K_r, req = bs_energy_K_r_mean, bs_energy_req_mean
            if bs_energy.get(tuple((name1, name2))) is not None:
                coef = bs_energy.get(tuple((name1, name2)))
                K_r, req = coef['K_r'], coef['req']

            bse += K_r * (r - req) ** 2
            self.map_unic_chains_to_r.append({(origin_atom_id, target_atom_id): r})

        return bse / 2

    def CalcAnglesEnergy(self, x=None):
        self.unic_chains = []
        self.map_unic_chains_to_angles = []
        # не забыть углы из таблицы перевести в градусы
        for i in range(len(self.Mol2.bonds)):
            origin_atom_id = int(self.Mol2.bonds[i]['origin_atom_id'])
            target_atom_id = int(self.Mol2.bonds[i]['target_atom_id'])
            for j in range(len(self.Mol2.bonds)):
                if j != i:
                    if int(self.Mol2.bonds[j]['origin_atom_id']) == target_atom_id:
                        new_target_atom_id = int(self.Mol2.bonds[j]['target_atom_id'])
                        self.unic_chains.append((origin_atom_id, target_atom_id, new_target_atom_id))

                    if int(self.Mol2.bonds[j]['origin_atom_id']) == origin_atom_id:
                        new_target_atom_id = int(self.Mol2.bonds[j]['target_atom_id'])
                        self.unic_chains.append((new_target_atom_id, origin_atom_id, target_atom_id))

            for j in range(len(self.Mol2.bonds)):
                if j != i:
                    if int(self.Mol2.bonds[j]['target_atom_id']) == target_atom_id:
                        new_target_atom_id = int(self.Mol2.bonds[j]['origin_atom_id'])
                        self.unic_chains.append((origin_atom_id, target_atom_id, new_target_atom_id))

        for i in self.unic_chains:
            for j in self.unic_chains:
                if self.unic_chains.count(j) > 1:
                    self.unic_chains.remove(j)
                elif i != j:
                    if set(j) == set(i):
                        self.unic_chains.remove(j)
        ind = 0
        for i in self.unic_chains:
            if x is None:
                angle = self.CalcAngle_(i[0], i[1], i[2])
            else:
                angle = x[ind]
                ind += 1
            self.map_unic_chains_to_angles.append({i: angle})

        # print(self.map_unic_chains_to_angles)

        be = 0
        for i in self.map_unic_chains_to_angles:
            name1 = self.Mol2.atoms[list(i.keys())[0][0] - 1]['atom_type']
            name2 = self.Mol2.atoms[list(i.keys())[0][1] - 1]['atom_type']
            name3 = self.Mol2.atoms[list(i.keys())[0][2] - 1]['atom_type']
            K_q, qeq = bending_energy_K_q_mean, bending_energy_qeq_mean
            if bending_energy.get(tuple((name1, name2, name3))) is not None:
                coef = bending_energy.get(tuple((name1, name2, name3)))
                K_q, qeq = coef['K_q'], coef['qeq']
            # print(list(i.keys())[0], '(', name1, name2, name3, ')', K_q, qeq)

            be += K_q * (list(i.values())[0] - qeq) ** 2
        return be / 2

    def CalcAngle_(self, origin_atom_id, target_atom_id, new_target_atom_id):
        dot_1 = Coord(self.Mol2.atoms[origin_atom_id - 1]['x'], self.Mol2.atoms[origin_atom_id - 1]['y'],
                      self.Mol2.atoms[origin_atom_id - 1]['z'])
        dot_2 = Coord(self.Mol2.atoms[target_atom_id - 1]['x'], self.Mol2.atoms[target_atom_id - 1]['y'],
                      self.Mol2.atoms[target_atom_id - 1]['z'])
        dot_3 = Coord(self.Mol2.atoms[new_target_atom_id - 1]['x'], self.Mol2.atoms[new_target_atom_id - 1]['y'],
                      self.Mol2.atoms[new_target_atom_id - 1]['z'])
        return CalcAngle(dot_1, dot_2, dot_3)


def constraint(x):
    return x


def main():
    Mol2 = MolParser('azz.mol2')
    E = Energy(Mol2)
    print(E.calcBondEnergy())
    print(E.CalcAnglesEnergy())
    print(E.calc_vdw_coloubm_energy_energy())
    n = len(E.Mol2.atoms) * 3
    # res1 = minimize(lambda x: CalcSumOfBondEnergy(x), x0=numpy.array([0] * n), method='SLSQP', bounds=[(0, None)] * n)
    # print(res1.x)
    m = len(E.map_unic_chains_to_angles)
    # res2 = minimize(lambda x: CalcSumOfAngleEnergy(x), x0=numpy.array([0] * m), method='SLSQP', bounds=[(0, None)] * m)
    # print(res2.x)
    s = len(E.lst_of_unic_chains)
    # res3 = minimize(lambda x: calc_sum_vdw_coloubm_energy_energy(x), x0=numpy.array([0.1] * s), method='SLSQP', bounds=[(0, None)] * s)
    # print(res3.x)

    x0_ = []
    for i in range(len(E.Mol2.atoms)):
        x_ = E.Mol2.atoms[i]['x']
        y_ = E.Mol2.atoms[i]['y']
        z_ = E.Mol2.atoms[i]['z']
        x0_.append(x_)
        x0_.append(y_)
        x0_.append(z_)
    res = minimize(
        lambda x:
        numpy.sum((
            E.calcBondEnergy(x[:n]),
            E.CalcAnglesEnergy(x[n:(n + m)]),
            E.calc_vdw_coloubm_energy_energy(x[(n + m):])
        )),
        x0=numpy.array(x0_ + [10] * (m + s)),

        method='SLSQP',

        options={
            'maxiter': 25,
            'ftol': 1e-2,
        },

        bounds=[(0., None)] * n + [(0., 180.)] * m + [(0., None)] * s,
        constraints={"fun": lambda x:
        numpy.sum((
            E.calcBondEnergy(x[:n]),
            E.CalcAnglesEnergy(x[n:(n + m)]),
            E.calc_vdw_coloubm_energy_energy(x[(n + m):])
        )), "type": "ineq"}

    )
    print(res)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
