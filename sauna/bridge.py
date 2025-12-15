from dataclasses import dataclass
import serpentTools
import math
import numpy as np
import pathlib

@dataclass
class SensitivityResult:
    functional: str
    zam: int
    reaction: int
    group_structure: list
    sensitivity_vector: list
    uncertainty_vector: list
    sensitivity: float
    uncertainty: float

class Serpent():
    """A static class responsible for the routines to get sensitivities
    from Serpent.

    """

    """Dictionary of commonly used functionals for sensitivity
    and uncertainty analysis in Serpent:

    - keff   - Eigenvalue
    - beff   - Effective delayed neutron fraction
    - leff   - Effective prompt neutron lifetime (not Lambda-eff,
               but effective prompt neutron generation time:
               l-eff = k-eff*Lambda-eff)
    - lambda - Effective precursor decay constant

    """

    FUNCTIONALS = {'keff'         : 'Eigenvalue',
                   'leff'         : 'l-eff',
                   'beff'         : 'beta-eff',
                   'beffGroup1'   : 'beta-eff-1',
                   'beffGroup2'   : 'beta-eff-2',
                   'beffGroup3'   : 'beta-eff-3',
                   'beffGroup4'   : 'beta-eff-4',
                   'beffGroup5'   : 'beta-eff-5',
                   'beffGroup6'   : 'beta-eff-6',
                   'beffGroup7'   : 'beta-eff-7',
                   'beffGroup8'   : 'beta-eff-8',
                   'lambda'       : 'lambda-eff',
                   'lambdaGroup1' : 'lambda-eff-1',
                   'lambdaGroup2' : 'lambda-eff-2',
                   'lambdaGroup3' : 'lambda-eff-3',
                   'lambdaGroup4' : 'lambda-eff-4',
                   'lambdaGroup5' : 'lambda-eff-5',
                   'lambdaGroup6' : 'lambda-eff-6',
                   'lambdaGroup7' : 'lambda-eff-7',
                   'lambdaGroup8' : 'lambda-eff-8'}

    """Dictionary of reactions constained as strings in
    Serpent results.
    
    """

    REACTION_TO_MT = {
        'total xs'     : 1,
        'ela leg mom 1': 251,
        'nubar total'  : 452,
        'nubar delayed': 455,
        'nubar prompt' : 456,
        'chi total'    : 1018,
        'chi delayed'  : 1455,
        'chi prompt'   : 1456}


    def __new__(cls):
        # This class is not intended to be instantiated.
        raise TypeError('A static class cannot be instantiated.')

    @staticmethod
    def _get_mt(reaction):
        """Get an MT number from Serpent perturbation.

        Parameters
        ----------
        pert : str
            Serpent perturbation name

        Returns
        -------
        int
            MT number corresponding to the given perturbation
        """   

        if reaction in Serpent.REACTION_TO_MT:
            return Serpent.REACTION_TO_MT[reaction]
        
        else:
            return int(reaction.split()[1])

    @staticmethod 
    def read(file):
        """A serpentTools reader wrapper method to get data from a sensitivity
        Serpent output file.

        Parameters
        ----------
        s: str
            Path to a _sens0.m file
        
        Returns
        -------
        List of SensitivityResult instances

        """ 
        
        path = pathlib.Path(file)
        s = serpentTools.read(path)
        serpent_functionals = list(s.sensitivities.keys())
        zams = s.zais
        zam_length = len(zams)
        group_structure = [e*1e6 for e in s.energies[:]]
        print('The Serpent sensitivity file contains the following functionals:', serpent_functionals)
        sensitivities = []
        for f in serpent_functionals:
            for zai in np.arange(zam_length):
                zam = list(zams.items())[zai][0]
                for pert in s.perts.keys():
                    try:
                        sens = s.sensitivities[f]
                        slice = sens[s.materials['total'], s.zais[zam], s.perts[pert]]
                        value, std = slice[:, 0], abs(slice[:, 1] * slice[:, 0])

                        int_s = s.energyIntegratedSens[f]
                        int_slice = int_s[s.materials['total'],  s.zais[zam],  s.perts[pert]]
                        int_value, int_std = int_slice[0], np.abs(int_slice[1] * int_slice[0])

                        if pert == 'chi total':
                            int_value = np.sum(value)
                        
                        # This addresses some cases where group values are nan
                        if np.isnan(int_value):
                            value = [0 if np.isnan(group_value) else group_value for group_value in value]
                            std   = [0 if np.isnan(group_std)   else group_std   for group_std   in std ]
                            int_value = np.sum(value)
                            int_std = np.sqrt(np.sum([group_std**2 for group_std in std]))

                        sensitivity = SensitivityResult(functional         = Serpent.FUNCTIONALS[f],
                                                  zam                = zam,
                                                  reaction           = Serpent._get_mt(pert),
                                                  group_structure    = group_structure,
                                                  sensitivity_vector = value,
                                                  uncertainty_vector = std,
                                                  sensitivity        = int_value,
                                                  uncertainty        = int_std)
                        sensitivities.append(sensitivity)
                        
                    except:
                        # In case get_mt is unable to get correct value
                        continue

        return sensitivities
    
class SCALE():
    """"A static class responsible for the routines to get sensitivities
    from SCALE.
    
    """

    # The number of lines in the HEAD of an SDF
    SDF_HEAD_NUMBER = 5

    def __new__(cls):
        # This class is not intended to be instantiated.
        raise TypeError('A static class cannot be instantiated.')

    @staticmethod
    def _scaleid_to_zam(scale_id):
        """Convert a SCALE ID to the corresponding
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
            IDs are provided https://scale-manual.ornl.gov/XSLib.html

        """  
        
        str_id = str(scale_id)
        length = len(str_id)

        if length in (4, 5, 6):
            zam = scale_id * 10
        elif length == 7:
            first_digit = int(str_id[0])
            zam = (scale_id - first_digit*1e6)*10 + first_digit
        else:
            print(f'Import of the SCALE ID = {scale_id} is not supported.')

        return int(zam)
    
    @staticmethod 
    def read(file, type, functional):
        """Read the SCALE output.

        Parameters
        ----------
        file: str
            Path to a .sdf file
        type : str
            Type of the SDF. type='A' if the results of a deterministic
            calculation are imported. type='B' if the resuls of a 
            stochastic calculation are imported.
            'B' is set by default, 'A' is not implemented yet.
        functional : str
            Since each .sdf file provides data for one functional,
            this argument should be passed.

        Returns
        -------
        list
            List of SensitivityResult instances

        """ 
        path = pathlib.Path(file)
        
        lines = open(path, "r").readlines()
        
        # Values for a generic .sdf file 
        group_number = int(lines[1].split()[0])
        energy_lines_number = math.ceil((group_number + 1)/5)
        group_lines_number = math.ceil(group_number/5)
        lines_btw_values = 4 + group_lines_number*2 
        profile_number = int(lines[2].split()[0])
        program = lines[-11].split()[1]

        # Get group structure
        groups = []
        for l in range(SCALE.SDF_HEAD_NUMBER, SCALE.SDF_HEAD_NUMBER + energy_lines_number):
            line = lines[l]
            groups.extend([float(x) for x in line.split()])
        
        group_structure = groups[::-1]

        print('Number of profiles:', profile_number)
        print('Number of groups:', group_number)
        print('Program:', program)

        if (type != 'B') & ((program == 'kenovi') | (program == 'kenova')):
            print(f"A stochastic code {program} was used. Setting type = 'B'.")
        elif (type != 'A') & (type != 'B'):
            raise ValueError("Please, set type='A' if the SDF is from a deterministic solution of SCALE, \
                            e.g. TSUNAMI-2D, or set type='B' if the SDF is \
                            from a stochastic solution of SCALE, e.g. TSUNAMI-3D")
        
        sensitivities = []
        for p in np.arange(profile_number):
            
            profile_parms = lines[SCALE.SDF_HEAD_NUMBER + energy_lines_number + lines_btw_values*p].split()
            reaction = int(profile_parms[3])
            
            # Id of the volume in the SCALE input
            region_id = lines[SCALE.SDF_HEAD_NUMBER + energy_lines_number + lines_btw_values*p + 1].split()[0]

            # Set this to '0' to extract only sensitivities of the whole system
            if region_id == '0':
                zam = SCALE._scaleid_to_zam(int(profile_parms[2]))
                
                integral_sensitivity = lines[SCALE.SDF_HEAD_NUMBER + energy_lines_number + lines_btw_values*p + 3].split()
                int_value, int_std = float(integral_sensitivity[0]), float(integral_sensitivity[1])
                
                value, std = [], []
                for i in range(group_lines_number):
                    value_line = lines[SCALE.SDF_HEAD_NUMBER + energy_lines_number + lines_btw_values*p + 4 + i]
                    std_line = lines[SCALE.SDF_HEAD_NUMBER + energy_lines_number + lines_btw_values*p + 4 + group_lines_number + i]
                    value.extend([float(x) for x in value_line.split()])
                    std.extend([float(x) for x in std_line.split()])
                value, std = value[::-1], std[::-1]

                sensitivity = SensitivityResult(functional = functional,
                                          zam = zam,
                                          reaction = reaction,
                                          group_structure = group_structure,
                                          sensitivity_vector = value,
                                          uncertainty_vector = std,
                                          sensitivity = int_value,
                                          uncertainty = int_std)
                
                sensitivities.append(sensitivity)

        return sensitivities
    