from collections import Counter
from InteractionRichPresentation import *

class MatchingFeatures(object):

    def __init__(self, seed):
        self.seed = seed
        self.seed.replace_T_U()
        if self.is_canonical() :
            self.canonic = "Canonic"
        elif self.is_simple_non_canonical() :
            self.canonic = "Simple Non-Canonic"
        else:
            self.canonic = "None"



    def extract_matching_features(self):
        self.seed_match_type()

    def miRNA_match_position(mir, mr_site):  # 20
        mirna = mir.upper().replace('T', 'U')
        mrna = mr_site.upper().replace('T', 'U')

        AU = ['AU', 'UA']
        GC = ['GC', 'CG']
        GU = ['GU', 'UG']

        if len(mirna) < 20:
            mirna += '-' * 20
            mrna += '-' * 20
        mmp_dic = {}
        for i in range(21)[1:]:
            key = 'miRNA_match_position' + str(i + 1)
            pair = mirna[-i] + mrna[-i]
            if pair in AU:
                mmp_dic[key] = 2
            elif pair in GC:
                mmp_dic[key] = 1
            elif pair in GU:
                mmp_dic[key] = 3
            elif '-' in pair:
                mmp_dic[key] = 5
            else:
                mmp_dic[key] = 4
        return mmp_dic

    def miRNA_pairing_count(mir, mr_site):  # 6*3=18
        mirna = mir.upper().replace('T', 'U')
        mrna = mr_site.upper().replace('T', 'U')

        AU = ['AU', 'UA']
        GC = ['GC', 'CG']
        GU = ['GU', 'UG']
        MM = ['AA', 'AG', 'AC', 'UU', 'UC', 'GA', 'GG', 'CA', 'CU', 'CC']

        mpc_dic = {'Seed_GC': 0,
                   'Seed_AU': 0,
                   'Seed_GU': 0,
                   'Seed_mismatch': 0,
                   'Seed_bulge': 0,
                   'Seed_bulge_nt': 0,
                   'Total_GC': 0,
                   'Total_AU': 0,
                   'Total_GU': 0,
                   'Total_mismatch': 0,
                   'Total_bulge': 0,
                   'Total_bulge_nt': 0,
                   'X3p_GC': 0,
                   'X3p_AU': 0,
                   'X3p_GU': 0,
                   'X3p_mismatch': 0,
                   'X3p_bulge': 0,
                   'X3p_bulge_nt': 0}
        for i in range(len(mirna) + 1)[1:]:
            pair = mirna[-i] + mrna[-i]
            if pair in AU:
                mpc_dic['Total_AU'] += 1
                if -9 < i < -1:
                    mpc_dic['Seed_AU'] += 1
                if i <= -9:
                    mpc_dic['X3p_AU'] += 1
            elif pair in GC:
                mpc_dic['Total_GC'] += 1
                if -9 < i < -1:
                    mpc_dic['Seed_GC'] += 1
                if i <= -9:
                    mpc_dic['X3p_GC'] += 1
            elif pair in GU:
                mpc_dic['Total_GU'] += 1
                if -9 < i < -1:
                    mpc_dic['Seed_GU'] += 1
                if i <= -9:
                    mpc_dic['X3p_GU'] += 1
            elif pair in MM:
                mpc_dic['Total_mismatch'] += 1
                if -9 < i < -1:
                    mpc_dic['Seed_mismatch'] += 1
                if i <= -9:
                    mpc_dic['X3p_mismatch'] += 1
            elif '-' in pair:
                mpc_dic['Total_bulge_nt'] += 1
                if -9 < i < -1:
                    mpc_dic['Seed_bulge_nt'] += 1
                if i <= -9:
                    mpc_dic['X3p_bulge_nt'] += 1
        mirna = 'A' + mirna
        for i in range(len(mirna) + 1)[1:]:
            if mirna[-i] == '-' and mirna[-i - 1] != '-':
                mpc_dic['Total_bulge'] += 1
                if -9 < i < -1:
                    mpc_dic['Seed_bulge'] += 1
                if i <= -9:
                    mpc_dic['X3p_bulge'] += 1
        return mpc_dic

    def tostring(self):
        classstr = ""
        classstr = classstr + "is canonic : {}\n".format(self.canonic)
        classstr = classstr + "seed type :  {}\n".format(self.seed_type)

        return classstr


    def __str__(self):
        return self.tostring()
