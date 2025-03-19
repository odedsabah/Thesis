import pandas as pd
import os
import sys


class Metaphlane4_abundance:

    def __init__(self,Path_MGX_MP):
        self.Path_MGX_MP_CD = Path_MGX_MP
        self.species_abundance = None
        self.Taxa_Abundance = None
        self.species_main = None

    def get_data(self, Project_name):
        '''
        This function reads the information from a metaphlan file and stores it in two lists (bacteria_name, percent_list)
        '''
        self.species_abundance = dict()
        samples = [x for x in os.listdir(self.Path_MGX_MP_CD) if x.endswith(".txt")]
        for sample in samples:
            sample_name = sample.split('.')[0]
            for k in self.species_abundance:
                self.species_abundance[k] = 0
            # print(sample_name)
            with open(f'{self.Path_MGX_MP_CD}/{sample}') as s:
                f = s.readlines()
            for line in f:
                #k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Coprococcus|s__Coprococcus_comes 2|1239|186801|186802|186803|33042|410072        0.00264
                if line.startswith("k__Bacteria"):
                    try:
                        self.species_abundance[line.strip().split("\t")[0].split("|")[-1]] = line.split()[2]
                    except:
                        print(f'this values: {line.split()[2]} is not abundance ')
            self.Taxa_Abundance = pd.DataFrame((self.species_abundance.items()), columns=['Abundance',sample_name]).set_index('Abundance')

            self.species_main = pd.concat([self.species_main,self.Taxa_Abundance], axis=1, join="outer")
        self.species_main.to_csv(f"~/metaanalysis/datasets/species_abundance/{Project_name}_species_abundance.csv")
        # self.species_main.to_csv(f"~/{Project_name}_species_abundance.csv")


def main():
    if len(sys.argv) != 3:
        quit("\nUsage: " + sys.argv[0] + "<Path_MGX_MP4> <Project_name> \n\n")

    Path_MGX_MP = sys.argv[1]
    Project_name = sys.argv[2]

    M_A = Metaphlane4_abundance(Path_MGX_MP)
    M_A.get_data(Project_name)

if __name__ == '__main__':
    main()


