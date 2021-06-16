import math
import numpy as np
numbers = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
simbols = ['+', '-', '=', '[', ']', '(', ')']


class SmilesData(object):
    def __init__(self, file_name_):
        self.file_name = file_name_
        self.smiles_ = self.read_file(file_name_).split()
        self.unic_chains = self.get_unic_chains(self.smiles_)
        self.docs = self.calc_smiles_vector(self.smiles_)

     #   self.calc_smiles_vector(['C=CN=C1'])
     #   print(self.unic_chains)

    def read_file(self, file_name_):
        with open(file_name_, 'r') as f:
            return f.read()

    def get_unic_chains(self, lst_) -> set:
     #   global unic_chain
        unic_chains = set()
        for i in lst_:
            for j in range(len(i)):
                if numbers.count(i[j]) > 0:
                    if numbers.count(i[j - 1]) > 0:
                        unic_chain += i[j]
                        # for m in range(len(unic_atom)):
                        #     if simbols.count(unic_atom[m]) > 0:
                        #         self.unic_chains.add(unic_atom[pos:m])
                        #         pos = m
                        unic_chains.add(unic_chain)
                    else:
                        unic_chain = ''
                        k = j - 1
                        while k != -1:
                            temp = i[k]
                            if numbers.count(temp) == 0:
                                unic_chain += i[k]
                            else:
                                break
                            k -= 1
                        unic_chain = unic_chain[::-1]
                        unic_chain += i[j]
                        pos = 0
                        # for m in range(len(unic_atom)):
                        #     if simbols.count(unic_atom[m]) > 0:
                        #         self.unic_chains.add(unic_atom[pos:m])
                        #         pos = m
                        unic_chains.add(unic_chain)
        return unic_chains

    def calc_smiles_vector(self, smiles_):
        matr = []
        for i in smiles_:
            vec = []
            item = self.get_unic_chains([i])
            for j in self.unic_chains:
                if j in item:
                    vec.append(1)
                else:
                    vec.append(0)
            matr.append(vec)
        return matr
        # print(matr)


class Ranking(object):
    def __init__(self, file_name, request):
        self.a = SmilesData(file_name)
        self.f = self.a.calc_smiles_vector(request)[0]
        self.p = {}
        self.ranking_docs()
        self.get_result()



    def ranking_docs(self):
        for i in range(len(self.a.docs)):
            self.p[ self.cos_sim(self.f, self.a.docs[i]) ]= self.a.smiles_[i]

    def get_result(self):
        x = sorted(self.p.items(), key=lambda x: x[0], reverse=True)
        for i in x:
            print(i)

    def cos_sim(self, vec_1: list, vec_2: list):
        vec_1 = np.array(vec_1)
        vec_2 = np.array(vec_2)
        vec_1_norm = vec_1/math.sqrt(sum(vec_1**2))
        vec_2_norm = vec_2/math.sqrt(sum(vec_2**2))
        return np.dot(vec_1_norm, vec_2_norm)


def main():
    a = Ranking('GDB13_Subset-ABCDEFGH.smi', ['O1N=NC2'])


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
