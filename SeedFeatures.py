from collections import Counter
from InteractionRichPresentation import *

class SeedFeatures(object):

    def __init__(self, seed):
        self.seed = seed
        self.seed.replace_T_U()
        if self.is_canonical() :
            self.canonic = "Canonic"
        elif self.is_simple_non_canonical() :
            self.canonic = "Simple Non-Canonic"
        else:
            self.canonic = "None"

    def is_canonical (self) :
        #exact W-C pairing of 2-7 or 3-8 nts of mirna
        mirna = self.seed.mir_inter
        mrna = self.seed.mrna_inter

        pairs2_7 = self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c']
        pairs3_8 = self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_c']
        return pairs2_7==6 or pairs3_8==6

    def is_simple_non_canonical (self):
        #pairing at positions 2-7 or 3-8, allowing G-U pairs and up to one bulged or mismatched nucleotide
        def evaluate_simple_non_canonical_condition (dict):
            return (dict['count_mismatch'] <= 1) and ((dict['count_c'] + dict['count_w']) >= 5)

        mirna = self.seed.mir_inter
        mrna = self.seed.mrna_inter
        pairs2_7_dict = self.seed_complementary(mirna[-7:-1], mrna[-7:-1])
        pairs2_7_is_simple_non_canonical = evaluate_simple_non_canonical_condition (pairs2_7_dict)
        pairs3_8_dict = self.seed_complementary(mirna[-8:-2], mrna[-8:-2])
        pairs3_8_is_simple_non_canonical = evaluate_simple_non_canonical_condition(pairs3_8_dict)

        return pairs2_7_is_simple_non_canonical or pairs3_8_is_simple_non_canonical



    def extract_seed_features(self):
        self.seed_match_type()



    def seed_complementary(self, seq1, seq2):
        count_c = 0
        count_w = 0
        count_mismatch = 0
        c = ['AU', 'UA', 'GC', 'CG']
        w = ['GU', 'UG']

        for i in range(len(seq1)):
            ss = seq1[i] + seq2[i]
            if ss in c:
                count_c += 1
            elif ss in w:
                count_w += 1
            else:
                count_mismatch+=1
        result = {'count_c': count_c, 'count_w': count_w, 'count_mismatch' : count_mismatch}
        return result

    def seed_match_type(self):  # 26
        c4 = ['AU', 'UA', 'GC', 'CG']
        w2 = ['GU', 'UG']
        smt_dic = {'Seed_match_8mer': 0,
                   'Seed_match_8merA1': 0,
                   'Seed_match_7mer1': 0,
                   'Seed_match_7mer2': 0,
                   'Seed_match_7merA1': 0,
                   'Seed_match_6mer1': 0,
                   'Seed_match_6mer2': 0,
                   'Seed_match_6mer3': 0,
                   'Seed_match_6mer1GU1': 0,
                   'Seed_match_6mer2GU1': 0,
                   'Seed_match_6mer3GU1': 0,
                   'Seed_match_6mer1GU2': 0,
                   'Seed_match_6mer2GU2': 0,
                   'Seed_match_6mer3GU2': 0,
                   'Seed_match_6mer1GU3': 0,
                   'Seed_match_6mer2GU3': 0,
                   'Seed_match_6mer3GU3': 0,
                   'Seed_match_6mer1GU4': 0,
                   'Seed_match_6mer2GU4': 0,
                   'Seed_match_6mer3GU4': 0,
                   'Seed_match_6mer1GU5': 0,
                   'Seed_match_6mer2GU5': 0,
                   'Seed_match_6mer3GU5': 0,
                   'Seed_match_6mer1GU6': 0,
                   'Seed_match_6mer2GU6': 0,
                   'Seed_match_6mer3GU6': 0}
        mirna = self.seed.mir_inter
        mrna = self.seed.mrna_inter

        # # Seed_match_8mer
        if self.seed_complementary(mirna[-8:], mrna[-8:])['count_c'] == 8:
            smt_dic['Seed_match_8mer'] = 1

        # # Seed_match_8merA1
        if self.seed_complementary(mirna[-8:-1], mrna[-8:-1])['count_c'] == 7 and mirna[-1] == 'A' and mirna[-1] + mrna[
            -1] not in c4:
            smt_dic['Seed_match_8merA1'] = 1
        if self.seed_complementary(mirna[-8:-1], mrna[-8:-1])['count_c'] == 7 and mrna[-1] == 'A' and mirna[-1] + mrna[
            -1] not in c4:
            smt_dic['Seed_match_8merA1'] = 1

        # # Seed_match_7mer1
        if self.seed_complementary(mirna[-7:], mrna[-7:])['count_c'] == 7:
            smt_dic['Seed_match_7mer1'] = 1

        # # Seed_match_7mer2
        if self.seed_complementary(mirna[-8:-1], mrna[-8:-1])['count_c'] == 7:
            smt_dic['Seed_match_7mer2'] = 1

        # # Seed_match_7merA1
        if self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c'] == 6 and mirna[-1] == 'A' and mirna[-1] + mrna[
            -1] not in c4:
            smt_dic['Seed_match_7merA1'] = 1
        if self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c'] == 6 and mrna[-1] == 'A' and mirna[-1] + mrna[
            -1] not in c4:
            smt_dic['Seed_match_7merA1'] = 1

        # # Seed_match_6mer1, Seed_match_6mer2, Seed_match_6mer3
        if self.seed_complementary(mirna[-6:], mrna[-6:])['count_c'] == 6:
            smt_dic['Seed_match_6mer1'] = 1
        if self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c'] == 6:
            smt_dic['Seed_match_6mer2'] = 1
        if self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_c'] == 6:
            smt_dic['Seed_match_6mer3'] = 1

        # # Seed_match_6mer1GU1,2,3,4,5,6
        if self.seed_complementary(mirna[-6:], mrna[-6:])['count_c'] == 5 and self.seed_complementary(mirna[-6:], mrna[-6:])[
            'count_w'] == 1:
            smt_dic['Seed_match_6mer1GU1'] = 1
        if self.seed_complementary(mirna[-6:], mrna[-6:])['count_c'] == 4 and self.seed_complementary(mirna[-6:], mrna[-6:])[
            'count_w'] == 2:
            smt_dic['Seed_match_6mer1GU2'] = 1
        if self.seed_complementary(mirna[-6:], mrna[-6:])['count_c'] == 3 and self.seed_complementary(mirna[-6:], mrna[-6:])[
            'count_w'] == 3:
            smt_dic['Seed_match_6mer1GU3'] = 1
        if self.seed_complementary(mirna[-6:], mrna[-6:])['count_c'] == 2 and self.seed_complementary(mirna[-6:], mrna[-6:])[
            'count_w'] == 4:
            smt_dic['Seed_match_6mer1GU4'] = 1
        if self.seed_complementary(mirna[-6:], mrna[-6:])['count_c'] == 1 and self.seed_complementary(mirna[-6:], mrna[-6:])[
            'count_w'] == 5:
            smt_dic['Seed_match_6mer1GU5'] = 1
        if self.seed_complementary(mirna[-6:], mrna[-6:])['count_c'] == 0 and self.seed_complementary(mirna[-6:], mrna[-6:])[
            'count_w'] == 6:
            smt_dic['Seed_match_6mer1GU6'] = 1

        # # Seed_match_6mer2GU1,2,3,4,5,6
        if self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c'] == 5 and \
                self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_w'] == 1:
            smt_dic['Seed_match_6mer2GU1'] = 1
        if self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c'] == 4 and \
                self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_w'] == 2:
            smt_dic['Seed_match_6mer2GU2'] = 1
        if self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c'] == 3 and \
                self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_w'] == 3:
            smt_dic['Seed_match_6mer2GU3'] = 1
        if self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c'] == 2 and \
                self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_w'] == 4:
            smt_dic['Seed_match_6mer2GU4'] = 1
        if self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c'] == 1 and \
                self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_w'] == 5:
            smt_dic['Seed_match_6mer2GU5'] = 1
        if self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_c'] == 0 and \
                self.seed_complementary(mirna[-7:-1], mrna[-7:-1])['count_w'] == 6:
            smt_dic['Seed_match_6mer2GU6'] = 1

        # # Seed_match_6mer3GU1,2,3,4,5,6
        if self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_c'] == 5 and \
                self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_w'] == 1:
            smt_dic['Seed_match_6mer3GU1'] = 1
        if self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_c'] == 4 and \
                self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_w'] == 2:
            smt_dic['Seed_match_6mer3GU2'] = 1
        if self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_c'] == 3 and \
                self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_w'] == 3:
            smt_dic['Seed_match_6mer3GU3'] = 1
        if self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_c'] == 2 and \
                self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_w'] == 4:
            smt_dic['Seed_match_6mer3GU4'] = 1
        if self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_c'] == 1 and \
                self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_w'] == 5:
            smt_dic['Seed_match_6mer3GU5'] = 1
        if self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_c'] == 0 and \
                self.seed_complementary(mirna[-8:-2], mrna[-8:-2])['count_w'] == 6:
            smt_dic['Seed_match_6mer3GU6'] = 1


        self.smt_dic = smt_dic
        self.seed_type = [seed_type for seed_type, onehot in smt_dic.items() if onehot == 1]

    def tostring(self):
        classstr = ""
        classstr = classstr + "is canonic : {}\n".format(self.canonic)
        classstr = classstr + "seed type :  {}\n".format(self.seed_type)

        return classstr


    def __str__(self):
        return self.tostring()
