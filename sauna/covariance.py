import sandy
import os
import pandas as pd
import numpy as np
import math
import multiprocessing
import pathlib 

from .auxiliary import *

"""List of commonly used covariances in terms of the MT numbers:

- MT1   - total (n,t)
- MT2   - elastic scattering (n,n)
- MT4   - inelastic scattering (n,n')
- MT16  - inelastic scattering (n,2n)
- MT18  - fission (n,f)
- MT102 - radioactive capture (n,gamma)
- MT103 - neutron-proton (n,p)
- MT104 - neutron-deutron (n,H2)
- MT105 - neutron-triton (n,H3)
- MT106 - neutron-helium-3 (n,He3)
- MT107 - neutron-alpha (n,alpha)
- MT251 - average scattering cosine (<mu>)
- MT452 - total neutron multiplicity (nu_t)
- MT455 - total neutron multiplicity (nu_d)
- MT1018 - fission spectrum (chi) - not a part of the standard
           MT number list of the ENDF-6 format, taken for
           simplicity as it is done in the SCALE code system

"""

MTS = {1, 2, 4, 16, 18, 102, 103, 104, 105, 106, 107, 251, 452, 455, 456}


class Covariances():
    """A collection of covariances to interact with.
    To properly create a Covariances instance it is necessary
    to set the path to a folder.

    Attributes
    ----------
    library : str, optional
        Library of the covariances. If not specified, the assumed library
        is ENDF/B-VII.1, and it is processed accordingly with possible 
        errors prone to some libraries.
    group_structure : numpy.ndarray
        Group structure used in the elements (Covariance instances) of 
        Covariances to process the data to.
    covariances : list
        List of the Covariance instances.
    zams : numpy.ndarray
        List of zams present in the Covariances instance.

    """

    def __init__(self, library=''):
        self.library = library
        self._group_structure = []
        self._covariances = []

    def __repr__(self):
        return (f"{self.__class__.__name__}({self.library!r}, {(len(self.group_structure)-1)!r})")

    @property
    def library(self):
        return self._library
    
    @library.setter
    def library(self, library):
        self._library = library

    @property
    def group_structure(self):
        return self._group_structure
    
    @group_structure.setter
    def group_structure(self, group_structure):
        group_number = len(group_structure)
        if group_number >=2 &  group_number <= 1501:
            self._group_structure = np.array(group_structure)
        else:
            raise ValueError(f'The number of group must be between 1 and 1500, but {group_number-1} is provided')

    @property
    def covariances(self):
        return self._covariances
    
    @covariances.setter
    def covariances(self, covariances):
        self._covariances = covariances

    def append(self, covariance):
        self._covariances.append(covariance)

    def extend(self, covariances):
        self._covariances.extend(covariances)

    def from_endf(self, file):
        """Generate covariance matrices from an ENDF-6 file.

        Parameters
        ----------
        file : str
            The relative path to an ENDF-6 file to process, e.g.,
            '../NuclearData/BROND-3.1/n_1325_13-Al-27.dat'.

        Return
        ------
        covs_in_file: list
            List of Covariance instances in an ENDF-6 file.

        """

        # Read and write the log file
        text_file = open(f'{self.library}.txt', "a+")

        # Get the text from an ENDF-6 file
        path = pathlib.Path(file)
        text = get_text(path)
        covs_in_file = []
        try:
            # Get an ENDF-6 tape after proccessing
            tape = sandy.Endf6.from_text(text)
            mat = tape.mat[0]
            zam = get_zam(tape)


            print("-------------------------------")
            print("Processing", path)
            print("-------------------------------")

            #These conditions are specified otherwise NJOY does not allow processing some files
            errorr = process_file(tape, self.library, group_structure=self.group_structure)

            try:
                cov33 = errorr["errorr33"].get_cov()
                local_mf = cov33.data
                mts = local_mf.index.levels[1]

                for mt_row in set(mts) & MTS:
                    for mt_col in set(mts) & MTS:

                        # Check whether it is a zero matrix
                        local_df = local_mf.loc[(mat, mt_row), (mat,mt_col)]
                        if local_df.eq(0).all(axis=None) == False:

                            # Avoid double accounting
                            if mt_row <= mt_col:
                                covariance = Covariance()
                                covariance.library = self.library
                                covariance.group_structure = self.group_structure
                                covariance.mf = 33
                                covariance.mat = mat
                                covariance.zam_1 = zam
                                covariance.zam_2 = zam
                                covariance.reaction_1 = mt_row
                                covariance.reaction_2 = mt_col
                                covariance.dataframe = local_df
                                covs_in_file.append(covariance)   

                if (self.library == "ENDF-B-VI.8") & (zam == 922380):
                    local_df = local_mf.loc[(mat, 51), (mat,51)]
                    if local_df.eq(0).all(axis=None) == False:

                        covariance = Covariance()
                        covariance.library = self.library
                        covariance.group_structure = self.group_structure
                        covariance.mf = 33
                        covariance.mat = mat
                        covariance.zam_1 = zam
                        covariance.zam_2 = zam
                        covariance.reaction_1 = 4
                        covariance.reaction_2 = 4
                        covariance.dataframe = local_df
                        covs_in_file.append(covariance)   

            except KeyError:
                text_file.write(f'{path} has got a KeyError exception. (check MF33 presence)\n')
            except TypeError as inst:
                text_file.write(f'{path} has got a TypeError exception. (MF33, {inst})\n')

            try:
                cov31 = errorr["errorr31"].get_cov()
                local_mf = cov31.data
                mts = local_mf.index.levels[1]
                for mt_row in set(mts) & MTS:
                    for mt_col in set(mts) & MTS:

                        # Check whether it is a zero matrix
                        local_df = local_mf.loc[(mat, mt_row), (mat,mt_col)]
                        if local_df.eq(0).all(axis=None) == False:

                            # Avoid double accounting
                            if mt_row <= mt_col:
                                covariance = Covariance()
                                covariance.library = self.library
                                covariance.group_structure = self.group_structure
                                covariance.mf = 31
                                covariance.mat = mat
                                covariance.zam_1 = zam
                                covariance.zam_2 = zam
                                covariance.reaction_1 = mt_row
                                covariance.reaction_2 = mt_col
                                covariance.dataframe = local_df
                                covs_in_file.append(covariance)                       
            except KeyError:
                text_file.write(f'{path} has got a KeyError exception. (check MF31 presence)\n')
            except TypeError as inst:
                text_file.write(f'{path} has got a TypeError exception. (MF31, {inst})\n')

            try:
                cov35 = errorr["errorr35"].get_cov()
                local_mf = cov35.data
                covariance = Covariance()
                covariance.library = self.library
                covariance.group_structure = self.group_structure
                covariance.mf = 35
                covariance.mat = mat
                covariance.zam_1 = zam
                covariance.zam_2 = zam
                covariance.reaction_1 = 1018
                covariance.reaction_2 = 1018
                covariance.dataframe = local_mf.loc[(mat, 18), (mat, 18)]  
                covs_in_file.append(covariance) 
            except KeyError:
                text_file.write(f'{path} has got a KeyError exception. (check MF35 presence)\n')
            except TypeError as inst:
                text_file.write(f'{path} has got a TypeError exception. (MF35, {inst})\n')    

            try:
                cov34 = errorr["errorr34"].get_cov()
                local_mf = cov34.data
                covariance = Covariance()
                covariance.library = self.library
                covariance.group_structure = self.group_structure
                covariance.mf = 34
                covariance.mat = mat
                covariance.zam_1 = zam
                covariance.zam_2 = zam
                covariance.reaction_1 = 251
                covariance.reaction_2 = 251
                covariance.dataframe = local_mf.loc[(mat, 251), (mat,251)]  
                covs_in_file.append(covariance) 

            except KeyError:
                text_file.write(f'{path} has got a KeyError exception. (check MF34 presence)\n')
            except TypeError as inst:
                text_file.write(f'{path} has got a TypeError exception. (MF34, {inst})\n')

        except ValueError as inst:
            text_file.write(f'{path} has got a ValueError exception. ({inst})\n')
        except pd.errors.ParserError as inst:
            text_file.write(f'{path} has got a pandas.errors.ParserError issue. Processing issue. ({inst})\n')
        except IndexError as inst:
            text_file.write(f'{path} has got an IndexError issue. Processing issue. ({inst})\n')

        text_file.close()   

        return covs_in_file

    def from_endfs(self, folder, extension='.dat', parallel = True):
        """Generate covariance matrices from ENDF-6 files in a folder.

        Parameters
        ----------
        folder : str
            The relative path to a folder with ENDF-6 files to process, e.g.,
            '../NuclearData/BROND-3.1/'.
        extension : str, optional
            This defines the extension of the ENDF-6 files in the folder for
            processing, e.g., '.dat', '.tendl'.
        parallel : bool, optional
            This defines whether multithreading applied.

        """

        # Get all the paths to files inside a fodler
        path_folder = pathlib.Path(folder).glob(f'*{extension}')
        paths = [p for p in path_folder if p.is_file()]

        # Get over the whole folder for each file
        if parallel == True:
            with multiprocessing.Pool(processes = multiprocessing.cpu_count()) as pool:
                results = pool.map(self.from_endf, paths)
                for res in results:
                    self.extend(res)
        else:
            for p in paths:
                res = self.from_endf(p)
                self.extend(res)
        
        print('Total number of covariances:', len(self.covariances))
        print("-------------------------------")
        print(f'The {self.library} library has been processed.')
        print("-------------------------------")

    def from_abbn(self, file):
        """Generate covariance matrices from an ABBN file.
        It assumes each file contains data about only one
        nuclide or element.

        Parameters
        ----------
        file : str
            The relative path to an ABBN file to be imported, e.g.,
            '../NuclearData/ABBN/n_1325_13-Al-27.dat'.

        """

        """Dictionary of the nuclides present in the ABBN
        covariance data and the corresponding ZAM number.

        """

        ABBN_NUCLIDES = {'H'   : [10010, 10020, 10030],
                        'LI-7' : [30070],
                        'BE'   : [40100],
                        'B-10' : [50100],
                        'C'    : [60000, 60120, 60130, 60140],
                        'N'    : [70140, 70150],
                        'O'    : [80160, 80170, 80180],
                        'NA'   : [110230],
                        'MG'   : [120240, 120250, 120260],
                        'AL'   : [130270],
                        'SI'   : [140280, 140290, 140300],
                        'CR'   : [240500, 240520, 240530],
                        'MN'   : [250550],                 
                        'FE'   : [260540, 250560, 260570],
                        'NI'   : [280580, 280600],       
                        'CU'   : [290630,290650],
                        'ZN'   : [300640, 300660, 300670, 300680, 300700],
                        'ZR'   : [400900, 400910, 400920, 400930, 400940, 400950, 400960],
                        'MO'   : [420920, 420940, 420950, 420960, 420970, 420980, 421000],
                        'W'    : [741800, 741820, 741830, 741840, 741860],
                        'PB'   : [822040, 822060, 822070, 822080],
                        'BI'   : [832090],
                        'TH28' : [902280],
                        'TH30' : [902300],
                        'TH32' : [902320],
                        'PA31' : [912310],
                        'PA33' : [912330],
                        'U233' : [922330],
                        'U234' : [922340],
                        'U235' : [922350],
                        'U236' : [922360],
                        'U238' : [922380],
                        'NP37' : [932370],
                        'PU38' : [942380],
                        'PU39' : [942390],
                        'PU40' : [942400],
                        'PU41' : [942410],
                        'PU42' : [942420],
                        'AM41' : [952410],
                        'AM43' : [952430],
                        'CM42' : [962420],
                        'CM43' : [962430],
                        'CM44' : [962440],
                        'CM45' : [962450],
                        'CM46' : [962460],
                        'CM48' : [962480]}

        """Dictionary of the reactions in the ABBN covariance data
        and the corresponding MT number.

        """

        ABBN_REACTIONS = {    "1" :   1,
                              "2" :   2,
                              "4" :   4,
                             "18" :  18,
                            "101" : 101,
                            "102" : 102,
                            "452" : 452}

        nuclides_1 = []
        nuclides_2 = []
        raw_covs = []
        group_numbers = []
        reactions_1 = []
        reactions_2 = []
        # Additional lines besides covariance values
        ADDITIONAL_LINES = 10

        lines = open(file, "r").readlines()
        length = len(lines)
        if lines[1].strip() == '*':
            l = 1
        else:
            l = 0
        while l < length:
            group_number = int(lines[l + 1][43:46].strip())
            reaction = lines[l + 4][14:17].strip()

            # How many groups are present
            given_number = int(lines[l + 6].split()[-1])

            # Get uncertainty vector and correlations
            uncertainty_vector = np.zeros(group_number)
            cor_matrix = np.zeros((group_number, group_number))
            start_group = group_number - given_number
            for row in range(given_number):
                variance = float(lines[l + 8 + row].split()[1])
                uncertainty_vector[start_group + row] = np.sqrt(variance)

                # Populate the correlation matrix
                for col in range(0, row+1):
                    if col == row:
                        cor_matrix[start_group + row, start_group + col] = int(lines[l + 8 + row].split()[col + 2]) / 100
                    else:
                        cor_matrix[start_group + row, start_group + col] = int(lines[l + 8 + row].split()[col + 2]) / 100
                        cor_matrix[start_group + col, start_group + row] = cor_matrix[start_group + row, start_group + col]

            # Get a base for forming Covariance instances
            nuclide_1 = ABBN_NUCLIDES[lines[l + 4].split()[1][4:]]
            nuclide_2 = ABBN_NUCLIDES[lines[l + 4].split()[1][4:]]
            cov_matrix = cor_matrix * np.outer(uncertainty_vector, uncertainty_vector)
            for nuclide in nuclide_1:
                raw_covs.append(cov_matrix)
                nuclides_1.append(nuclide)
                nuclides_2.append(nuclide)
                group_numbers.append(group_number)
                reactions_1.append(ABBN_REACTIONS[reaction])
                reactions_2.append(ABBN_REACTIONS[reaction])

            # move l to next cov
            l += ADDITIONAL_LINES + given_number

        # Remove redundunt covariances
        if (101 in reactions_1) & (102 in reactions_1):
            index = reactions_1.index(101)
            raw_covs.pop(index)
            nuclides_1.pop(index)
            nuclides_2.pop(index)
            group_numbers.pop(index)
            reactions_1.pop(index)
            reactions_2.pop(index)
        elif (101 in reactions_1) & (102 not in reactions_1):
            index = reactions_1.index(101)
            reactions_1[index] = 102
            reactions_2[index] = 102
     
        for rc in range(len(raw_covs)):
            covariance = Covariance()
            covariance.library = self.library
            covariance.zam_1 = nuclides_1[rc]
            covariance.zam_2 = nuclides_2[rc]
            covariance.group_structure = self.group_structure
            covariance.reaction_1 = reactions_1[rc]
            covariance.reaction_2 = reactions_2[rc]
            covariance.dataframe = pd.DataFrame(raw_covs[rc])
            # Set the MF number based upon the MT number
            if covariance.reaction_1 == 251:
                covariance.mf = 34
            elif (covariance.reaction_1 == 452) | (covariance.reaction_1 == 455) | (covariance.reaction_1 == 456):
                covariance.mf = 31
            elif (covariance.reaction_1 == 1018):
                covariance.mf = 35
            else:
                covariance.mf = 33
            self.covariances.append(covariance) 

    def from_abbns(self, folder, extension='.TAB'):
        """Import covariance matrices from ABBN files in a folder
        to the Covariances instances.

        Parameters
        ----------
        folder : str
            Relative path to a folder with the abbn covariance files, e.g.,
            '../NuclearData/ABBN_COV'.

        """
        
        path_folder = pathlib.Path(folder).glob(f'*{extension}')
        paths = [p for p in path_folder if p.is_file()]
        
        print('Importing covariance data from the ABBN files.')
        for path in paths:
            try:
                self.from_abbn(path)
            except:
                print(f'An issue was encountered while importing {path}. Check the file.')
        
        print('The covariances have been imported successfully.')

    def from_excel(self, file):
        """Import a covariance matrix from an .xlsx file 
        to the Covariances instance.

        Parameters
        ----------
        file : str
            Path to a .xlsx 

        """

        fname = file.name.split('.')[0]
        zam_1, zam_2, first_mt, second_mt = fname.split('-')

        covariance = Covariance()
        covariance.library = self.library
        covariance.group_structure = self.group_structure
        covariance.zam_1 = int(zam_1)
        covariance.zam_2 = int(zam_2)
        covariance.reaction_1 = int(first_mt)
        covariance.reaction_2 = int(second_mt)
        covariance.dataframe = pd.read_excel(file).fillna(0)

        # Set the MF number based upon the MT number
        if first_mt == 251:
            covariance.mf = 34
        elif (first_mt == 452) | (first_mt == 455) | (first_mt == 456):
            covariance.mf = 31
        elif (first_mt == 1018):
            covariance.mf = 35
        else:
            covariance.mf = 33

        self.covariances.append(covariance) 
             
    def from_excels(self, folder):
        """Import covariance matrices from .xlsx files in a folder
        to the Covariances instances.

        Parameters
        ----------
        folder : str
            Relative path to a folder with .xlsx files to processn, e.g.,
            '../NuclearData/BROND-3.1/'.

        """
        
        covtype = '.xlsx'
        path_folder = pathlib.Path(folder).glob(f'*{covtype}')
        paths = [p for p in path_folder if p.is_file()]
        
        print('Importing covariance data from Excels.')
        for path in paths:
            self.from_excel(path)
        
        print('The covariances have been imported successfully.')

    def from_commara(self, file):
        """Import covariance matrices from a COMMARA file to
        the Covariances instances.

        Parameters
        ----------
        file : str
            Relative path to a COMMARA file to import covariances, e.g.,
            '../NuclearData/COMMARA.cov'.

        """

        """Dictionary of the nuclides present in the COMMARA
        covariance data and the corresponding ZAM number.

        """

        COMMARA_NUCLIDES = {'H1'      :  10010,
                            'H2'      :  10020,
                            'He4'     :  20040,
                            'Li6'     :  30060,
                            'Li7'     :  30070,
                            'Be9'     :  40090,
                            'B10'     :  50100,
                            'B11'     :  50110,
                            'CARBON'  :  60000,
                            'N15'     :  70150, 
                            'O16'     :  80160,
                            'F19'     :  90190,
                            'Na23'    : 110230,
                            'Mg24'    : 120240,
                            'Mg25'    : 120250,
                            'Mg26'    : 120260,
                            'Al27'    : 130270,
                            'Si28'    : 140280,
                            'Si29'    : 140290,
                            'Si30'    : 140300,
                            'Cr50'    : 240500,
                            'Cr52'    : 240520,
                            'Cr53'    : 240530,
                            'Mn55'    : 250550,
                            'Fe54'    : 260540,
                            'Fe56'    : 250560,
                            'Fe57'    : 260570,
                            'Ni58'    : 280580,
                            'Ni60'    : 280600,
                            'Zr90'    : 400900,
                            'Zr91'    : 400910,
                            'Zr92'    : 400920,
                            'Zr93'    : 400930,
                            'Zr94'    : 400940,
                            'Zr95'    : 400950,
                            'Zr96'    : 400960,
                            'Nb95'    : 410950,
                            'Mo92'    : 420920,
                            'Mo94'    : 420940,
                            'Mo95'    : 420950,
                            'Mo96'    : 420960,
                            'Mo97'    : 420970,
                            'Mo98'    : 420980,
                            'Mo100'   : 421000,
                            'Tc99'    : 430990,
                            'Ru101'   : 441010,
                            'Ru102'   : 441020,
                            'Ru103'   : 441030,
                            'Ru104'   : 441040,
                            'Ru106'   : 441060,
                            'Rh103'   : 451030,
                            'Pd105'   : 461050,
                            'Pd106'   : 461060,
                            'Pd107'   : 461070,
                            'Pd108'   : 461080,
                            'Ag109'   : 471090,
                            'I127'    : 531270,
                            'I129'    : 531290,
                            'Xe131'   : 541310,
                            'Xe132'   : 541320,
                            'Xe134'   : 541340,
                            'Cs133'   : 551330,
                            'Cs135'   : 551350,
                            'La139'   : 571390,
                            'Ce141'   : 581410,
                            'Pr141'   : 591410,
                            'Nd143'   : 601430,
                            'Nd145'   : 601450,
                            'Nd146'   : 601460,
                            'Nd148'   : 601480,
                            'Pm147'   : 611470,
                            'Sm149'   : 621490,
                            'Sm151'   : 621510,
                            'Sm152'   : 621520,
                            'Eu153'   : 631530,
                            'Eu155'   : 631550,
                            'Gd155'   : 641550,
                            'Gd156'   : 641560,
                            'Gd157'   : 641570,
                            'Gd158'   : 641580,
                            'Gd160'   : 641600,
                            'Er166'   : 681660,
                            'Er167'   : 681670,
                            'Er168'   : 681680,
                            'Er170'   : 681700,
                            'Pb204'   : 822040,
                            'Pb206'   : 822060,
                            'Pb207'   : 822070,
                            'Pb208'   : 822080,
                            'Bi209'   : 832090,
                            'Th232'   : 902320,
                            'U233'    : 922330,
                            'U234'    : 922340,
                            'U235'    : 922350,
                            'U236'    : 922360,
                            'U238'    : 922380,
                            'Np237'   : 932370,
                            'Pu238'   : 942380,
                            'Pu239'   : 942390,
                            'Pu240'   : 942400,
                            'Pu241'   : 942410,
                            'Pu242'   : 942420,
                            'Am241'   : 952410,
                            'Am242m'  : 952421,
                            'Am243'   : 952430,
                            'Cm242'   : 962420,
                            'Cm243'   : 962430,
                            'Cm244'   : 962440,
                            'Cm245'   : 962450,
                            'Cm246'   : 962460}

        """Dictionary of the reactions in the COMMARA covariance data
        and the corresponding MT number

        """

        COMMARA_REACTIONS = {
                            "ELASTIC"   :     2,
                            "INELASTIC" :     4,
                            "NxN"       :    16,
                            "FISSION"   :    18,
                            "CAPTURE"   :   102,
                            "P1ELASTIC" :   251,
                            "NU"        :   452,
                            "CHI"       :  1018}
     
        lines = open(file, "r").readlines()

        nuclides_1 = []
        nuclides_2 = []
        raw_covs = []
        group_numbers = []
        reactions_1 = []
        reactions_2 = []

        for l in range(len(lines)):
            # if this is True, it containes numbers
            line = lines[l]
            if line.upper().isupper():
                
                if l != 0:
                    group_numbers.append(-1+int(np.sqrt(1+len(numbers))))
                    raw_covs.append(numbers)
                    
                numbers = []
                splitted_line = line.split()
                reactions_1.append(splitted_line[0])
                nuclides_1.append(splitted_line[1])
                reactions_2.append(splitted_line[2])
                nuclides_2.append(splitted_line[3])

            else:

                if '-' in line:
                    line = line.replace('-', ' -')
                splitted_values = [float(x) for x in line.split()]
                numbers.extend(splitted_values)

        # Append the last one
        group_numbers.append(-1+int(np.sqrt(1+len(numbers))))
        raw_covs.append(numbers)

        for rc in range(len(raw_covs)):
            covariance = Covariance()
            covariance.library = self.library
            covariance.zam_1 = int(COMMARA_NUCLIDES[nuclides_1[rc]])
            covariance.zam_2 = int(COMMARA_NUCLIDES[nuclides_2[rc]])
            covariance.reaction_1 = int(COMMARA_REACTIONS[reactions_1[rc]])
            covariance.reaction_2 = int(COMMARA_REACTIONS[reactions_2[rc]])
            covariance.group_structure = self.group_structure
            
            # inverse the position of the values
            diag_1 = raw_covs[rc][:group_numbers[rc]][::-1]
            diag_2 = raw_covs[rc][group_numbers[rc]:2*group_numbers[rc]][::-1]
            corr_matrix = np.array(raw_covs[rc][2*group_numbers[rc]:])[::-1].reshape(group_numbers[rc],group_numbers[rc])

            covariance.dataframe = pd.DataFrame(np.diag(diag_1)@corr_matrix@np.diag(diag_2))

            # Set the MF number based upon the MT number
            if covariance.reaction_1 == 251:
                covariance.mf = 34
            elif (covariance.reaction_1 == 452) | (covariance.reaction_1 == 455) | (covariance.reaction_1 == 456):
                covariance.mf = 31
            elif (covariance.reaction_1 == 1018):
                covariance.mf = 35
            else:
                covariance.mf = 33

            self.covariances.append(covariance) 

    def from_coverx(self, ampxcovconverter, file):
        """Import covariance matrices from a COVERX (the format
        of AMPX of the SCALE code system) file to the
        Covariances instances.

        Parameters
        ----------
        ampxcovconverter : str
            Path to the AmpxCOVConverter executable, which translates
            a COVERX file from the binary format to ASCII, the file
            name is not included e.g.
            'home/example/SCALE-6.2.4-Source/build/install/bin/AmpxCOVConverter'.
        file : str
            Path to a COVERX file

        Notes
        -----
            A large number of groups takes considerably more time for
            the AmpxCOVConverter work. The 56-group structure is recommended to use.
        """

        HEAD_LINES = 2
        MAX_GROUP_COLUMNS = 25
        CORR_HEAD = 3

        print('Converting AMPX binaries via AmpxCOVConverter ...')
        os.system(f'{ampxcovconverter}  -i {file} -f new -p yes')
        print('Converting AMPX ASCII data ...')

        lines = open(f'{file}.toc', "r").readlines()

        def scaleid_to_zam(scale_id):
            """Translate a SCALE ID to the corresponding
            ZAM value.

            Parameters
            ----------
            scale_id : int
                SCALE ID of a nuclide

            Return
            ------
            zam : int
                ZAM value from a SCALE ID

            Notes
            -----
                The SCALE IDs contain more types (e.g. S(alpha, beta)),
                which are not considered in this implementation. The 
                IDs are provided https://scale-manual.ornl.gov/XSLib.html .

            """  
                
            str_id = str(scale_id)
            length = len(str_id)

            if (length == 4) | (length == 5) | (length == 6):
                zam = scale_id * 10
            elif length == 7:
                first_digit = int(str_id[0])
                zam = (scale_id - first_digit*1e6)*10 + first_digit
            else:
                print(f'Import of the SCALE ID = {scale_id} is not supported.')

            return int(zam)

        # Get first lines ids
        first_line_ids = []
        for i in range(len(lines)):
            # Each covariance matrix description starts with ' Material 1'
            if lines[i].startswith(' Material 1'):
                first_line_ids.append(i)

                # Get materuals and reactions for the current covariance
                param_line = lines[i].split(',')
                zam_1 = scaleid_to_zam(int(param_line[0].split()[2]))
                first_mt = int(param_line[0].split()[5])
                zam_2 = scaleid_to_zam(int(param_line[1].split()[2]))
                second_mt = int(param_line[1].split()[5])
                
                if first_mt > second_mt:
                    first_mt , second_mt = second_mt, first_mt
                    zam_1, zam_2 = zam_2, zam_1

                # Get the number of groups
                group_number = 0
                while True:
                    if '*** correlation matrix ***' in lines[group_number + i]:
                        group_number -= HEAD_LINES
                        break 
                    group_number += 1
                
                # Get the group structure and two uncertainty vectors
                uncertainty_vector_1 = [] # col
                uncertainty_vector_2 = [] # row
                groups = []
                for e in range(HEAD_LINES + i, HEAD_LINES + group_number + i):
                    line = lines[e]

                    uncertainty_vector_1.append(float(line.split()[5]))
                    uncertainty_vector_2.append(float(line.split()[6]))

                    if e == (HEAD_LINES + i): groups.append(float(line.split()[1]))
                    groups.append(float(line.split()[2]))

                # Get number of part covariance matrix divided by
                number_of_parts = math.ceil(group_number / MAX_GROUP_COLUMNS)
                part_structure = []
                columns_left = group_number - (number_of_parts - 1)*MAX_GROUP_COLUMNS
                for p in range(number_of_parts):
                    if p == number_of_parts - 1:
                        part_structure.append(columns_left)
                    else:
                        part_structure.append(MAX_GROUP_COLUMNS)

                # Get correlations
                first_part = HEAD_LINES + group_number + HEAD_LINES + CORR_HEAD + i
                matrix = []
                for p in range(number_of_parts):
                    matrix_part = []

                    # Append rows of part p
                    for g in range(first_part + (CORR_HEAD + group_number)*p, first_part + group_number + (CORR_HEAD + group_number)*p):
                        row = []
                        splitted_line = lines[g].split()
                        for i in range(1, part_structure[p] + 1):
                            row.append(float(splitted_line[i]))
                        matrix_part.append(row)
                    
                    # Append data of values of part p
                    matrix_part = np.array(matrix_part)
                    if p == 0:
                        matrix = matrix_part
                    else:
                        matrix = np.append(matrix, matrix_part, axis = 1)

                # Convert the correlation matrix into the covariance matrix
                matrix = np.diag(uncertainty_vector_1)@matrix@np.diag(uncertainty_vector_2) / 1e3
                


                covariance = Covariance()
                covariance.library = self.library
                covariance.zam_1 = zam_1
                covariance.zam_2 = zam_2
                covariance.reaction_1 = first_mt
                covariance.reaction_2 = second_mt
                covariance.group_structure = groups[::-1]
                covariance.dataframe = pd.DataFrame(np.rot90(matrix, 2))

                # Set the MF number based upon the MT number
                if first_mt == 251:
                    covariance.mf = 34
                elif (first_mt == 452) | (first_mt == 455) | (first_mt == 456):
                    covariance.mf = 31
                elif (first_mt == 1018):
                    covariance.mf = 35
                else:
                    covariance.mf = 33

                self.covariances.append(covariance) 

        os.remove(f'{file}.new')
        os.remove(f'{file}.new.toc')
        os.remove(f'{file}.toc')

    def check_eigenvalues(self, eps=1e-5):
        """Check the mathermatical corectness of the 
        covariance matrices whether it is positive semidefinite.
        The method provides an output as a number of incorrect covariances
        and the data itself.

        """          

        number_of_positive = 0
        number_of_nonpositive = 0
        number_of_symmetric = 0
        number_of_nonsymmetric = 0
        number_of_complex = 0
        incorrect_ratios = []
        text_file = open(f'{self.library}-eigenvalues.txt', "a+")

        for cov in self.covariances:         
            if cov.reaction_1 == cov.reaction_2:
                number_of_symmetric += 1
                matrix = cov.dataframe.to_numpy()
                eigenvalues = np.linalg.eigvals(matrix)
                is_positive = np.all(eigenvalues >= 0)
                is_symmetric = np.all(np.abs(matrix-matrix.T) < 1e-8)
                

                if not is_symmetric:
                    number_of_nonsymmetric += 1

                if is_positive:
                    number_of_positive += 1
                    text_file.write(f'The matrix  {cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2} is positive semi-definite. The eigenvalues - \n[{eigenvalues}] \n\n')
                
                elif np.any(np.iscomplex(eigenvalues)):
                    number_of_complex += 1
                    first_complex = get_complex(eigenvalues)
                    ratio = np.abs(eigenvalues[first_complex])/np.max(eigenvalues).real
                    incorrect_ratios.append(ratio)
                    text_file.write(f'The matrix  {cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2} contains complex values. \
                                    \nThe first complelx eigenvalue at g={first_complex+1}. \
                                    \nThe ratio is {ratio}. \
                                    \nThe eigenvalues - \n[{eigenvalues}] \n\n')                    
                else:
                    number_of_nonpositive += 1
                    first_negative = get_negative(eigenvalues)
                    ratio = np.abs(eigenvalues[first_negative])/np.max(eigenvalues)
                    incorrect_ratios.append(ratio)
                    text_file.write(f'The matrix {cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2} is not positive semi-definite. \
                                    \nThe first negative eigenvalue at g={first_negative+1}. \
                                    \nThe ratio is {ratio}. \
                                    \nThe eigenvalues - \n[{eigenvalues}] \n\n') 

                

        text_file.write(f'The number of non-symmetric is {number_of_nonsymmetric} among {number_of_symmetric} checked matrices, which must be symmetric')
        text_file.write(f'The number of positive semi-definite matrices is {number_of_positive} of {number_of_symmetric} symmetric matrices among a total number of {len(self.covariances)} matrices')
        text_file.write(f'The number of matrices containing complex eigenvalues is {number_of_complex} of {number_of_symmetric} symmetric matrices among a total number of {len(self.covariances)} matrices')
        text_file.write(f'The number of non-positive and not complex matrices is {number_of_nonpositive} of {number_of_symmetric} symmetric matrices among a total number of {len(self.covariances)} matrices')
        text_file.write(f'The number of matrices with incorrect eigenvalues assuming tolerance={eps} is {len([i for i in incorrect_ratios if i > eps])} of {len(incorrect_ratios)} incorrect matrices among a total number of {len(self.covariances)} matrices')
        text_file.close()

        print('-------------------------------------------')
        print(f'The number of non-symmetric is {number_of_nonsymmetric} among {number_of_symmetric} checked matrices, which must be symmetric')
        print(f'The number of positive semi-definite matrices is {number_of_positive} of {number_of_symmetric} symmetric matrices among a total number of {len(self.covariances)} matrices')
        print(f'The number of matrices containing complex eigenvalues is {number_of_complex} of {number_of_symmetric} symmetric matrices among a total number of {len(self.covariances)} matrices')
        print(f'The number of matrices with incorrect eigenvalues assuming tolerance={eps} is {len([i for i in incorrect_ratios if i > eps])} of {len(incorrect_ratios)} incorrect matrices among a total number of {len(self.covariances)} matrices')
        
        print('-------------------------------------------')

    def check_corrs(self, eps = 1e-5):
        """Check whether the correlation matrices are correct or
        not, i.e. the values are between -1 and 1. The method
        provides an output as a number of incorrect correlation
        matrices and the data itself.

        Parameters
        ----------
        eps : float
            Allowed deviation from -1 and 1 to be ignored during
            a checking. That is, the matrices are assumed incorrect
            if at least one element has a value of <-1-eps or >1+eps .

        """

        number_of_symmetric = 0
        number_of_corr = 0
        text_file = open(f'{self.library}-correlations.txt', "a+")

        for cov in self.covariances:     
            if cov.reaction_1 == cov.reaction_2:
                array_of_covs = cov.dataframe.to_numpy()
                number_of_symmetric += 1
                corr = cov_to_corr(array_of_covs)

                if np.any(corr > 1+eps) | np.any(corr < -1-eps):
                    number_of_corr += 1
                    text_file.write(f'The correlation matrix {cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2} is not correct.  \n{corr}\n')
                    
        text_file.write(f'The number of incorrect correlation matrices is {number_of_corr} of {number_of_symmetric} symmetric matrices among a total number of {len(self.covariances)} matrices')
        text_file.close()

        print('-------------------------------------------')
        print(f'The number of incorrect correlation matrices is {number_of_corr} of {number_of_symmetric} symmetric matrices among a total number of {len(self.covariances)} matrices')
        print('-------------------------------------------')

    def export_corrs(self, save_to='./correlations/', fix_corr = False, eps = 1e-5):
        """Auxiliary method to export the correlation matrices based upon
        the covariance matrices.

        Parameters
        ----------
        save_to: str
            Relative path to save the processed correlations, e.g., './BROND-3.1'.

        """  

        # Create the directory if it does not exist
        os.makedirs(save_to, exist_ok=True)

        for cov in self.covariances:
            if cov.reaction_1 == cov.reaction_2:
                array_of_covs = cov.dataframe.to_numpy()
                corr = cov_to_corr(array_of_covs)

                if fix_corr:
                    corr = np.where(corr >  1+eps,  1, corr)                 
                    corr = np.where(corr < -1-eps, -1, corr)

            temp_df = cov.dataframe.copy()
            temp_df[:] = corr

            if os.path.exists(f'./{save_to}/{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2}.xlsx') == False:   
                temp_df.to_excel(f'./{save_to}/{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2}.xlsx', index=False)
            else:
                print(f'{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2}.xlsx exists')

        print('-------------------------------------------------')
        print('The correlations have been exported successfully.')
        print('-------------------------------------------------')

    def limit_covs(self):
        """Limit variances and covariances to 100% to avoid
        unphysical or overestimated values. It assumes that there is
        no change in correlations while changing the variances
        and covariances.

        """       

        number_of_matrices = 0
        text_file = open(f'{self.library}-limited.txt', "a+")

        for cov in self.covariances:
            if cov.reaction_1 == cov.reaction_2:
                # Get the covs for given zam if reaction_1 is not equal to reaction_2
                zam_covs = [i for i in self.get_by_zam(cov.zam_1) if i.reaction_1 != i.reaction_2]

                # Get the indices where the values are incorrect
                array_of_covs = cov.dataframe.to_numpy()
                diag = np.sqrt(np.diag(array_of_covs))
                incorrect_indices = np.where(diag > 1)[0]

                # Fix data for incorrect indices
                if len(incorrect_indices) > 0:
                    number_of_matrices += 1
                    text_file.write(f'{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2} has been fixed. \n')

                    for i in incorrect_indices:
                        # Firstly, fix covariances when the reactions are not the same
                        for zam_cov in zam_covs:
                            if zam_cov.reaction_1 == cov.reaction_1:
                                zam_cov.dataframe.values[i,:] /= diag[i]  
                            elif zam_cov.reaction_2 == cov.reaction_1:
                                zam_cov.dataframe.values[:,i] /= diag[i]  

                        # Secondly fix covariances when the reactions are the same
                        array_of_covs[i,:] /= diag[i]
                        array_of_covs[:,i] /= diag[i]
                        cov.dataframe.values[:,:] = array_of_covs

        text_file.close()

        print('-------------------------------------------')
        print(f'The number of matrices with over 100% values is {number_of_matrices} of a total number of {len(self.covariances)} matrices.')
        print('-------------------------------------------')

    def get_by_zam(self, zam):
        """Get a list of Covariance instances by a ZAM value
        from a Covariances instance.

        Parameters
        ----------
        zam : int
            ZAM value (ZAM = Z*10000 + A*10 + M)

        Returns
        -------
        numpy.ndarray
            List of Covariances based upon first ZAM.

        Notes
        -----
        Cross-material correlations are accounted only
        when the first ZAM has the correlation.
        The method is usually used when it iterated over all the
        ZAMs, therefore, cross-material correlations, corresponding
        to the second ZAM, are naturally accounted.
        """

        array = np.array([cov for cov in self.covariances if cov.zam_1 == zam])

        return array

    def get_by_reaction(self, mt):
        """Get a list of Covariance instances by a MT value
        from a Covariances instance.

        Parameters
        ----------
        mt : int
            MT value

        Returns
        -------
        numpy.ndarray
            List of Covariances based upon first MT.

        Notes
        -----
        Cross-reaction correlations are accounted only
        when the first reaction has the correlation.

        """
        array = np.array([cov for cov in self.covariances if cov.reaction_1 == mt])

        return array

    def get_by_params(self, zam_1, zam_2, reaction_1, reaction_2):
        """Get a Covariance instance by the ZAM, Reaction 1 MT
        number, and Reaction 2 MT numbers from a Covariances instance.

        Parameters
        ----------
        zam_1 : int
            ZAM value (ZAM = Z*10000 + A*10 + M)
        zam_2 : int
            ZAM value (ZAM = Z*10000 + A*10 + M)
        reaction_1 : int
            MT number for the first reactio
        reaction_2 : int
            MT number for the second reaction
            
        Returns
        -------
        Covariance

        """
        
        return next(cov for cov in self.covariances if (cov.zam_1 == zam_1) & (cov.zam_2 == zam_2) & (cov.reaction_1 == reaction_1) & (cov.reaction_2 == reaction_2))

    def to_excels(self, save_to='./covariances/', fix_corr = False, eps = 1e-5):
        """Export covariances in the Coviariances instance
        to Excel files.

        Parameters
        ----------
        save_to: str, optional
            Relative path to save the processed covariances, e.g., './BROND-3.1'.
        fix_corr: bool, optional
            Fix those matrices which correlations are mathematically incorrect.
            It sets the values over 1+err to 1 and the values less than -1-err
            to zero. err is asuumed equal to 1e-5 since the values in ENDF-6 files
            are limited to 1e-6.

        """  

        # Create the directory if it does not exist
        os.makedirs(save_to, exist_ok=True)

        
        for cov in self.covariances:
            if fix_corr:
                temp_df = cov.dataframe.copy()
                if (cov.reaction_1 == cov.reaction_2) & (cov.zam_1 == cov.zam_2): 
                    temp_df[:] = fix_corrs(cov, eps)

                if os.path.exists(f'./{save_to}/{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2}.xlsx') == False:   
                    temp_df.to_excel(f'./{save_to}/{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2}.xlsx', index=False)
                else:
                    print(f'{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2}.xlsx exists')
            else:
                if os.path.exists(f'./{save_to}/{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2}.xlsx') == False:   
                    cov.dataframe.to_excel(f'./{save_to}/{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2}.xlsx', index=False)
                else:
                    print(f'{cov.zam_1}-{cov.zam_2}-{cov.reaction_1}-{cov.reaction_2}.xlsx exists')

        print('------------------------------------------------')
        print('The covariances have been exported successfully.')
        print('------------------------------------------------')

class Covariance():
    """Contains a covariance matrix and corresponding information.

    Attributes
    ----------
    library : str, optional
        Name of the library
    mat : int, optional
        MAT number according to the corresponding ENDF-6 file.
    zam_1 : int
        ZAM number for the corresponding nuclide acorrding to
        the following formula: ZAM = Z*10000 + A*10 + M .
    zam_2 : int
        ZAM number 2 for the corresponding nuclide acorrding to
        the following formula: ZAM = Z*10000 + A*10 + M.
    mf : int, optional
        MF number according to the corresponding ENDF-6 file.
    reaction_1 : int
        MT number for the first reaction
    reaction_2 : int
        MT number for the second reaction
    group_structure : numpy.ndarray
        Energy grid for the covariance matrix. Contains the number of
        elements equal to G, the number of energy groups, + 1.
    dataframe : pandas.DataFrame
        Dataframe of the covariance matrix.
    
    """

    def __init__(self):
        self._library = None
        self._mat = None
        self._zam_1 = None
        self._zam_2 = None
        self._file = None
        self._reaction_1 = None
        self._reaction_2 = None
        self._group_structure = None
        self._dataframe = None

    def __repr__(self):
        return (f"{self.__class__.__name__}({self.zam_1!r}-{self.zam_2!r}, {self.reaction_1!r}-{self.reaction_2!r}, {(len(self.group_structure)-1)!r})")

    @property
    def library(self):
        return self._library
    
    @library.setter
    def library(self, library):
        self._library = library

    @property
    def mat(self):
        return self._mat
    
    @mat.setter
    def mat(self, mat):
        self._mat = mat

    @property
    def zam_1(self):
        return self._zam_1

    @property
    def zam_2(self):
        return self._zam_2

    @zam_1.setter
    def zam_1(self, zam_1):
        self._zam_1 = zam_1

    @zam_2.setter
    def zam_2(self, zam_2):
        self._zam_2 = zam_2

    @property
    def mf(self):
        return self._mf
    
    @mf.setter
    def mf(self, mf):
        self._mf = mf

    @property
    def reaction_1(self):
        return self._reaction_1
    
    @reaction_1.setter
    def reaction_1(self, reaction_1):
        self._reaction_1 = reaction_1
    
    @property
    def reaction_2(self):
        return self._reaction_2
    
    @reaction_2.setter
    def reaction_2(self, reaction_2):
        self._reaction_2 = reaction_2
    
    @property
    def group_structure(self):
        return self._group_structure
    
    @group_structure.setter
    def group_structure(self, group_structure):
        self._group_structure = group_structure

    @property
    def dataframe(self):
        return self._dataframe
    
    @dataframe.setter
    def dataframe(self, dataframe):
        self._dataframe = dataframe    
