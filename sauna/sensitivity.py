import numpy as np
import pandas as pd
import os

from  .bridge import Serpent, SCALE


"""Dictionary of used reactions in terms of the MT numbers in
current implementation:

- MT1    - total (n,t)
- MT2    - elastic scattering (n,n)
- MT4    - inelastic scattering (n,n')
- MT16   - inelastic scattering (n,xn)
- MT18   - fission (n,f)
- MT102  - radioactive capture (n,γ)
- MT103  - neutron-proton (n,p)
- MT104  - neutron-deutron (n,D)
- MT105  - neutron-triton (n,T)
- MT106  - neutron-helium-3 (n,He3)
- MT107  - neutron-alpha (n,α)
- MT251  - average scattering cosine (<μ>)
- MT452  - total neutron multiplicity (ν)
- MT455  - delayed neutron multiplicity (ν-delayed)
- MT456  - prompt neutron multiplicity (ν-prompt)
- MT1002 - thermal elastic scattering laws (S(α,β), it is not a part of the
           established MT numbers) 
- MT1004 - thermal inelastic scattering laws (S(α,β), it is not a part of the
           established MT numbers) 
- MT1018 - fission spectrum (χ) - not a part of the standard
           MT number list of the ENDF-6 format, taken for
           simplicity (the value of 1018 is taken from SCALE)
- MT1455 - delayed fission spectrum (χ-delayed) - not a part of the standard
           MT number list of the ENDF-6 format, taken for
           simplicity
- MT1456 - prompt fission spectrum (χ-prompt) - not a part of the standard
           MT number list of the ENDF-6 format, taken for
           simplicity
"""

REACTIONS = {1   : "(n,t)"      ,
             2   : "(n,n)"      ,
             4   : "(n,n')"     ,
             16  : "(n,2n)"     ,
             18  : "(n,f)"      ,
             102 : "(n,γ)"      ,
             103 : "(n,p)"      ,
             104 : "(n,D)"      ,
             105 : "(n,T)"      ,
             106 : "(n,He3)"    ,
             107 : "(n,α)"      ,
             251 : "(<μ>)"      ,
             452 : "(ν)"        ,
             455 : "(ν-d)",
             456 : "(ν-p)" ,
             1002: "(S(α,β))-el"  ,
             1004: "(S(α,β))-inel"  ,
             1018: "(χ)"        ,
             1455: "(χ-d)",
             1456: "(χ-p)"}

class Sensitivities():
    """A collection of Sensitivity instances to interact with.

    Attributes
    ----------
    group_structure : numpy.ndarray
        Energy grid for the sensitivity vector, it contains the number of
        elements equal to G, the number of energy groups, + 1.
    sensitivities: list
        List of the Sensitivity instances
    functionals : list
        List of the functionals included into the sensetivities list.
    zams : list
        List of the ZAM values includeed into the sensetivities list.
        
    """

    def __init__(self):
        self._group_structure = []
        self._sensitivities = []
        self._functionals = []
        self._zams = []

    def __repr__(self):
        return (f"{self.__class__.__name__}({len(self.functionals)!r}, {len(self.zams)!r}, {len(self.sensitivities)!r}, {(len(self.group_structure)-1)!r})")

    @property
    def group_structure(self):
        return self._group_structure

    @group_structure.setter
    def group_structure(self, group_structure):
        group_number = len(group_structure)
        if group_number >=2 &  group_number <= 1501:
            self._group_structure = group_structure
        else:
            raise ValueError(f'The number of group must be between 1 and 1500 due to the NJOY limitations, \
                             but {group_number-1} is provided')

    @property
    def sensitivities(self):
        return self._sensitivities
    
    @sensitivities.setter
    def sensitivities(self, sensitivities):
        self._sensitivities = sensitivities

    @property
    def functionals(self):
        return self._functionals

    @property
    def zams(self):
        return self._zams
    
    def append(self, sensitivity):
        self._sensitivities.append(sensitivity)
        if sensitivity.functional not in self.functionals:
            self._functionals.append(sensitivity.functional)  
        if sensitivity.zam not in self.zams:
            self._zams.append(sensitivity.zam) 
        
    def extend(self, sensitivities):
        self._sensitivities.extend(sensitivities)
        for sens in sensitivities:
            if sens.functional not in self.functionals:
                self._functionals.append(sens.functional)  
            if sens.zam not in self.zams:
                self._zams.append(sens.zam) 

    @classmethod
    def reactivity_difference(cls, sensitivities_nom, sensitivities_pert, keff_nom, keff_pert, dkeff_nom = 0, dkeff_pert = 0, functional = 'Reactivity-difference'):
        """Calculate senitivity of a reactivity difference to quantify 
        the influence of the inputs on the reactivity effects.
        
        Parameters
        ----------
        sensitivities_nom : Sensitivities
            Sensitivities instance for nominal state sensitivities
        sensitivities_pert : Sensitivities
            Sensitivities instance for perturbed state sensitivities
        keff_nom : float
            Eigenvalue in the nominal state
        keff_pert : float
            Eigenvalue in the perturbed state
        dkeff_nom : float, optional
            Statistical eigenvalue uncertainty in the nominal state
        dkeff_pert : float, optional
            Statistical eigenvalue uncertainty in the perturbed state
        functional : str, optional
            Name of the functional to be used for the new Sensitivities
            instance

        Returns
        -------
        Sensitivities
            Reactivity-difference sensitivities

        """   

        reactivity_sensitivities = Sensitivities()
        reactivity_sensitivities.group_structure = sensitivities_nom.group_structure
        
        # Number of elements for filling incorrect values
        out = np.ones(len(reactivity_sensitivities.group_structure) - 1)

        # Append for all the senses present in the nominal data
        nom_zams  = sensitivities_nom.zams
        for zam in nom_zams:
            senses_nom  = sensitivities_nom.get_by_zam(zam)
        
            for sens_nom in senses_nom:

                sens_pert = sensitivities_pert.get_by_params('Eigenvalue', zam, sens_nom.reaction)
                
                # Create a Sensitivity instance to put the reactivity sensitivity data
                reactivity_sensitivity = Sensitivity()
                reactivity_sensitivity.group_structure = sens_nom.group_structure
                reactivity_sensitivity.functional      = functional
                reactivity_sensitivity.zam             = sens_nom.zam
                reactivity_sensitivity.reaction        = sens_nom.reaction

                # Assign a sensitivity vector and its uncertainty for reactivity difference
                reactivity_sensitivity.sensitivity_vector = (sens_pert.sensitivity_vector*keff_nom - sens_nom.sensitivity_vector*keff_pert)/abs(keff_nom-keff_pert)
                
                numerator_1 = (sens_pert.sensitivity_vector*dkeff_nom)**2 + (sens_pert.uncertainty_vector*keff_nom)**2
                numerator_2 = (sens_nom.sensitivity_vector*dkeff_pert)**2 + (sens_nom.uncertainty_vector*keff_pert)**2
                numerator   =  numerator_1 + numerator_2
                denominator = dkeff_nom**2  + dkeff_pert**2
                first_term  = np.divide(numerator,(sens_pert.sensitivity_vector*keff_nom - sens_nom.sensitivity_vector*keff_pert)**2, out = out, 
                                        where=(sens_pert.sensitivity_vector*keff_nom - sens_nom.sensitivity_vector*keff_pert)!=0)
                second_term = denominator/(keff_nom-keff_pert)**2

                reactivity_sensitivity.uncertainty_vector = abs(reactivity_sensitivity.sensitivity_vector) * np.sqrt(first_term + second_term)

                # Assign an integrated sensitivity and its uncertainty for reactivity difference
                if (sens_pert.sensitivity*keff_nom - sens_nom.sensitivity*keff_pert) != 0:
                    reactivity_sensitivity.sensitivity = (sens_pert.sensitivity*keff_nom - sens_nom.sensitivity*keff_pert)/abs(keff_nom-keff_pert)
                    numerator_1 = (sens_pert.sensitivity*dkeff_nom)**2 + (sens_pert.uncertainty*keff_nom)**2
                    numerator_2 = (sens_nom.sensitivity*dkeff_pert)**2 + (sens_nom.uncertainty*keff_pert)**2
                    numerator   = numerator_1 + numerator_2
                    denominator = dkeff_nom**2  + dkeff_pert**2
                    first_term  = numerator/(sens_pert.sensitivity*keff_nom - sens_nom.sensitivity*keff_pert)**2
                    second_term = denominator/(keff_nom-keff_pert)**2

                    reactivity_sensitivity.uncertainty = abs(reactivity_sensitivity.sensitivity) * np.sqrt(first_term + second_term)

                else:
                    reactivity_sensitivity.sensitivity = 0
                    reactivity_sensitivity.uncertainty = 0

                reactivity_sensitivities.append(reactivity_sensitivity)

        # Append all the senses non-present in the nominal data but present in perturbed data
        pert_zams = list(set(sensitivities_pert.zams) - set(sensitivities_nom.zams))
        for zam in pert_zams:

            senses_pert  = sensitivities_pert.get_by_zam(zam)

            for sens_pert in senses_pert:
                # Create a Sensitivity instance to put the reactivity sensitivity data
                reactivity_sensitivity = Sensitivity()
                reactivity_sensitivity.group_structure = sens_pert.group_structure
                reactivity_sensitivity.functional      = functional
                reactivity_sensitivity.zam             = sens_pert.zam
                reactivity_sensitivity.reaction        = sens_pert.reaction

                # Assign a sensitivity vector and its uncertainty for reactivity difference
                reactivity_sensitivity.sensitivity_vector = sens_pert.sensitivity_vector*keff_nom / abs(keff_nom-keff_pert)

                numerator   =  (sens_pert.sensitivity_vector*dkeff_nom)**2 + (sens_pert.uncertainty_vector*keff_nom)**2
                denominator = dkeff_nom**2  + dkeff_pert**2
                first_term  = np.divide(numerator,(sens_pert.sensitivity_vector*keff_nom)**2, out = out, 
                                        where=(sens_pert.sensitivity_vector*keff_nom) != 0)
                second_term = denominator/(keff_nom-keff_pert)**2

                reactivity_sensitivity.uncertainty_vector = abs(reactivity_sensitivity.sensitivity_vector) * np.sqrt(first_term + second_term)

                # Assign an integrated sensitivity and its uncertainty for reactivity difference
                if sens_pert.sensitivity != 0:
                    reactivity_sensitivity.sensitivity = sens_pert.sensitivity*keff_nom / abs(keff_nom-keff_pert)

                    numerator   =  (sens_pert.sensitivity*dkeff_nom)**2 + (sens_pert.uncertainty*keff_nom)**2
                    denominator = dkeff_nom**2  + dkeff_pert**2
                    first_term  = numerator / (sens_pert.sensitivity*keff_nom)**2
                    second_term = denominator / (keff_nom-keff_pert)**2

                    reactivity_sensitivity.uncertainty = abs(reactivity_sensitivity.sensitivity) * np.sqrt(first_term + second_term)
                    
                else:
                    reactivity_sensitivity.sensitivity = 0
                    reactivity_sensitivity.uncertainty = 0

                reactivity_sensitivities.append(reactivity_sensitivity)

        return reactivity_sensitivities

    @classmethod
    def beta_eff(cls, sensitivities_nom, sensitivities_prompt, keff_nom, keff_prompt,  functional = 'beta-eff'):
        """Calculate senitivity of an effective delayed neutron fraction via
        the nominal eigenvalue and prompt eigenvalue
        
        Parameters
        ----------
        sensitivities_num : Sensitivities
            Sensitivities of eigenvalue
        sensitivities_Lambda : Sensitivities
            Sensitivities of prompt eigenvalue
        keff_nom : float
            Nominal eigenvalue
        keff_prompt : float
            Prompt eigenvalue
        functional : str, optional
            Name of the functional to be used for the new Sensitivities
            instance

        Returns
        -------
        Sensitivities
            Ratio sensitivities
        
        """        
        sensitivities_b = Sensitivities()
        sensitivities_b.group_structure = sensitivities_nom.group_structure
        denom_name = sensitivities_prompt.functionals[0]

        for zam in sensitivities_nom.zams:
            zam_senses = sensitivities_prompt.get_by_zam(zam)

            for nom_sens in zam_senses:
                prompt_sens = sensitivities_prompt.get_by_params(denom_name, zam, nom_sens.reaction)
                    
                sensitivity_b = Sensitivity()
                sensitivity_b.group_structure = nom_sens.group_structure
                sensitivity_b.functional      = functional
                sensitivity_b.zam             = nom_sens.zam
                sensitivity_b.reaction        = nom_sens.reaction
                sensitivity_b.sensitivity_vector = keff_prompt/(keff_nom-keff_prompt)*(nom_sens.sensitivity_vector - prompt_sens.sensitivity_vector)
                sensitivity_b.sensitivity = keff_prompt/(keff_nom-keff_prompt)*(nom_sens.sensitivity - prompt_sens.sensitivity)
                    
                sensitivity_b.uncertainty_vector = [0]*(len(sensitivities_b.group_structure) - 1)
                sensitivity_b.uncertainty = 0

                sensitivities_b.append(sensitivity_b)

        return sensitivities_b
    
    @classmethod
    def ratio(cls, sensitivities_num, sensitivities_denom, functional = 'Ratio'):
        """Calculate senitivity of a ratio of two arbitrary sensitivities.
        For instnace, it can be used to calculate the effective neutron generation
        time via the effective life-time and the eigenvalue
        
        Parameters
        ----------
        sensitivities_num : Sensitivities
            Sensitivities of a functional in the numerator
        sensitivities_Lambda : Sensitivities
            Sensitivities of a functional in the denomerator
        functional : str, optional
            Name of the functional to be used for the new Sensitivities
            instance

        Returns
        -------
        Sensitivities
            Ratio sensitivities
        
        """        
        sensitivities_r = Sensitivities()
        sensitivities_r.group_structure = sensitivities_num.group_structure
        denom_name = sensitivities_denom.functionals[0]

        for zam in sensitivities_num.zams:
            zam_senses = sensitivities_denom.get_by_zam(zam)

            for num_sens in zam_senses:
                denom_sens = sensitivities_denom.get_by_params(denom_name, zam, num_sens.reaction)

                sensitivity_r = Sensitivity()
                sensitivity_r.group_structure = num_sens.group_structure
                sensitivity_r.functional      = functional
                sensitivity_r.zam             = num_sens.zam
                sensitivity_r.reaction        = num_sens.reaction
                sensitivity_r.sensitivity_vector = num_sens.sensitivity_vector - denom_sens.sensitivity_vector
                sensitivity_r.sensitivity = num_sens.sensitivity - denom_sens.sensitivity
                    
                sensitivity_r.uncertainty_vector = np.sqrt(num_sens.uncertainty_vector**2 + denom_sens.uncertainty_vector**2)
                sensitivity_r.uncertainty = np.sqrt(num_sens.uncertainty**2 + denom_sens.uncertainty**2)

                sensitivities_r.append(sensitivity_r)

        return sensitivities_r

    @classmethod
    def promt_decay(cls, sensitivities_prompt, sensitivities_l, k_prompt, functional = 'alpha'):
        """Calculate senitivity of a prompt decay constant also known as 
        Rossi-alpha or prompt alpha eigenvalue via the sensitivities of the
        prompt eigenvalue and effective neutron life-time.
        
        Parameters
        ----------
        sensitivities_prompt : Sensitivities
            Sensitivities of the prompt eigenvalue, i.e. with no
            delayed neutrons in the simulation
        sensitivities_Lambda : Sensitivities
            Sensitivities of the effective life-time (l-eff)
        k_prompt : float
            Prompt eigenvalue
        functional : str, optional
            Name of the functional to be used for the new Sensitivities
            instance

        Returns
        -------
        Sensitivities
            Prompt alpha sensitivities
        
        """        

        sensitivities_a = Sensitivities()
        sensitivities_a.group_structure = sensitivities_prompt.group_structure

        for zam in sensitivities_prompt.zams:
            zam_senses = sensitivities_prompt.get_by_zam(zam)

            for prompt_sens in zam_senses:
                l_sens = sensitivities_l.get_by_params('l-eff', zam, prompt_sens.reaction)
                    
                sensitivity_a = Sensitivity()
                sensitivity_a.group_structure = prompt_sens.group_structure
                sensitivity_a.functional      = functional
                sensitivity_a.zam             = prompt_sens.zam
                sensitivity_a.reaction        = prompt_sens.reaction
                sensitivity_a.sensitivity_vector = k_prompt/(k_prompt-1)*prompt_sens.sensitivity_vector - l_sens.sensitivity_vector
                sensitivity_a.sensitivity = k_prompt/(k_prompt-1)*prompt_sens.sensitivity - l_sens.sensitivity
                    
                sensitivity_a.uncertainty_vector = [0]*(len(sensitivities_a.group_structure) - 1)
                sensitivity_a.uncertainty = 0

                sensitivities_a.append(sensitivity_a)
        
        return sensitivities_a

    @classmethod
    def breeding_ratio(cls, sensitivities_gamma, sensitivities_fission, R_gamma, R_fission, discard_scattering = False, functional = 'Breeding Ratio'):
        """Calculate senitivity of a breeding ratio via the sensitivities
        of the fissile-gamma/fertile-gamma ratio and the fissile-fission/
        fertile-gamma ratio. The method is intented to avoid the limitation
        of the modern stochastic codes to tally only single reaction rate in
        both numerator and denominator.
        
        Parameters
        ----------
        sensitivities_gamma : Sensitivities
            Sensitivities instance for the fissile-gamma/fertile-gamma ratio,
            e.g. the ratio of Pu239(n,\gamma) and U238(n,\gamma) reaction
            rates. Make sure it contains only one functional
        sensitivities_fission : Sensitivities
            Sensitivities instance for the fissile-fission/fertile-gamma ratio,
            e.g. the ratio of Pu239(n,f) and U238(n,\gamma) reaction rates.
            Make sure it contains only one functional
        R_gamma : float
            Nominal value of the fissile-gamma/fertile-gamma ratio
        R_fission : float
            Nominal value of the fissile-fission/fertile-gamma ratio
        discard_scattering : bool, optional
            Define whether to include or not the scattering sensitivities
            since they might have unreasonably high values in some cases
        functional : str, optional
            Name of the functional to be used for the new Sensitivities
            instance

        Returns
        -------
        Sensitivities
            Breeding ratio sensitivities
        
        """

        br_sensitivities = Sensitivities()
        br_sensitivities.group_structure = sensitivities_gamma.group_structure
        fission_name = sensitivities_fission.functionals[0]

        for zam in sensitivities_gamma.zams:
            zam_senses = sensitivities_gamma.get_by_zam(zam)

            for gamma_sens in zam_senses:
                if (discard_scattering == False) | ((gamma_sens.reaction != 2) & (gamma_sens.reaction != 4)):
                    fis_sens = sensitivities_fission.get_by_params(fission_name, zam, gamma_sens.reaction)
                    
                    br_sensitivity = Sensitivity()
                    br_sensitivity.group_structure = gamma_sens.group_structure
                    br_sensitivity.functional      = functional
                    br_sensitivity.zam             = gamma_sens.zam
                    br_sensitivity.reaction        = gamma_sens.reaction
                    br_sensitivity.sensitivity_vector = -1/(1+R_fission/R_gamma)*gamma_sens.sensitivity_vector - 1/(1+R_gamma/R_fission)*fis_sens.sensitivity_vector
                    br_sensitivity.sensitivity = -1/(1+R_fission/R_gamma)*gamma_sens.sensitivity - 1/(1+R_gamma/R_fission)*fis_sens.sensitivity
                    
                    br_sensitivity.uncertainty_vector = [0]*(len(br_sensitivities.group_structure) - 1)
                    br_sensitivity.uncertainty = 0

                    br_sensitivities.append(br_sensitivity)

        return br_sensitivities

    def from_serpent(self, file):
        """Import sensitivity data of all the reactions from a Serpent output 
        file (_sens0.m) to the Sensitivities instance.

        Parameters
        ----------
        file: str
            Path to the Serpent output, e.g., './MOX3600_sens0.m'.

        """  

        sensitivities = Serpent.read(file)
        for sens in sensitivities:
            if sens.reaction in REACTIONS.keys():
                sensitivity = Sensitivity()
                sensitivity.functional =    sens.functional
                sensitivity.group_structure = sens.group_structure
                sensitivity.zam = sens.zam
                sensitivity.reaction = sens.reaction
                sensitivity.sensitivity_vector = sens.sensitivity_vector
                sensitivity.uncertainty_vector = sens.uncertainty_vector
                sensitivity.sensitivity = sens.sensitivity
                sensitivity.uncertainty = sens.uncertainty
                self.append(sensitivity)

    def from_scale(self, file, type='B', functional = 'Eigenvalue'):
        """Import sensitivities from an SDF (Sensitivity Data
        File) file of SCALE to the Sensitivities instance.

        Parameters
        ----------
        file : str
            Relative path to an SDF file to get covariances, e.g.,
            '../SCALE/MET1000/MET1000.sdf'.  
        type : str, optional
            Type of the SDF. type='A' if the results of a deterministic
            calculation are imported. type='B' if the resuls of a 
            stochastic calculation are imported.
            'B' is set by default, 'A' is not implemented yet.
        functional : str, optional
            Since each .sdf file provides data for one functional,
            this argument should be passed. The default value is 
            'Eigenvalue'.

        Notes
        -----
            The description of SDF is provided in 
            https://scale-manual.ornl.gov/tsunami-ip-appAB.html .

        """  
              
        sensitivities = SCALE.read(file, type, functional)

        for sens in sensitivities:
            if sens.reaction in REACTIONS.keys():
                sensitivity = Sensitivity()
                sensitivity.functional =    sens.functional
                sensitivity.group_structure = sens.group_structure
                sensitivity.zam = sens.zam
                sensitivity.reaction = sens.reaction
                sensitivity.sensitivity_vector = sens.sensitivity_vector
                sensitivity.uncertainty_vector = sens.uncertainty_vector
                sensitivity.sensitivity = sens.sensitivity
                sensitivity.uncertainty = sens.uncertainty
                self.append(sensitivity)


        print('The data have been imported successfully.')

    def get_by_functional(self, functional):
        """Get a list of Sensitivity instances by a functional 
        from a Sensetivities instance.

        Parameters
        ----------
        functional : str
            Functional name           
        Returns
        -------
        numpy.ndarray
            Sensitivity instances with the same functionals in an numpy.ndarray .

        """

        return np.array([sensitivity for sensitivity in self.sensitivities if sensitivity.functional == functional])

    def get_by_zam(self, zam):
        """Get a list of Sensitivity instances by a ZAM value
        from a Sensetivities instance.

        Parameters
        ----------
        zam : int
            ZAM value (ZAM = Z*10000 + A*10 + M)

        Returns
        -------
        list
            Sensitivity instances with the same ZAM in an numpy.ndarray .

        """

        return [sensitivity for sensitivity in self.sensitivities if sensitivity.zam == zam]

    def get_by_reaction(self, reaction):
        """Get a list of Sensitivity instances by an MT number
        from a Sensetivities instance.

        Parameters
        ----------
        reaction : int
            MT number

        Returns
        -------
        list
            Sensitivity instances with the same MT number in an numpy.ndarray .

        """

        return [sensitivity for sensitivity in self.sensitivities if sensitivity.reaction == reaction]

    def get_by_params(self, functional, zam, reaction):
        """Get a list of Sensitivity instances by a functional,
        ZAM, reaction from a Sensetivities instance.

        Parameters
        ----------
        functional : str
            Functional name
        zam : int
            ZAM value (ZAM = Z*10000 + A*10 + M)
        reaction : int
            MT number        
            
        Returns
        -------
        Sensitivity
            Sensitivity instance with requested functional, zam, and reaction.

        """
        
        try:
            return next(sensitivity for sensitivity in self.sensitivities if (sensitivity.functional == functional) & (sensitivity.zam == zam) & (sensitivity.reaction == reaction))
        
        except:
            sensitivity = Sensitivity()
            sensitivity.group_structure = self.group_structure
            sensitivity.functional = functional
            sensitivity.zam = zam
            sensitivity.reaction = reaction
            sensitivity.sensitivity = 0.
            sensitivity.uncertainty = 0.
            sensitivity.sensitivity_vector = np.zeros(len(sensitivity.group_structure) - 1)
            sensitivity.uncertainty_vector = np.zeros(len(sensitivity.group_structure) - 1)

            return sensitivity

    def to_dataframe(self, functional, sort=True):
        """Export sensitivity data of all the relative sensitivities 
        from a Sensitivities instance to a dataframe. It generates diffirent
        sheet for given functional.

        Parameters
        ----------
        name: functional
            Name of the functional
        sort: bool, optional
            Sort by absolute sensitivities. The default value is
            True.       

        Returns
        -------
        sensitivity_df : pandas.DataFrame
        
        """ 

        sensitivity_df = pd.DataFrame(columns = ['Nuclide', 'MT', 'Reaction', 'Sensitivity', 'Statistical uncertainty'])
        temp_list = [sens for sens in self.sensitivities if sens.functional == functional]
        for sens in temp_list:
            new_row =   {
                'Nuclide':  f'{sens.zam}',
                'MT': f'MT{sens.reaction}',
                'Reaction': REACTIONS[sens.reaction],
                'Sensitivity': sens.sensitivity,
                'Statistical uncertainty': sens.uncertainty            
            }
            sensitivity_df.loc[len(sensitivity_df)] = new_row
        
        if sort:
            sensitivity_df = sensitivity_df.reindex(sensitivity_df['Sensitivity'].abs().sort_values(ascending=False).index).reset_index(drop=True)  

        return sensitivity_df

    def to_excel(self, name='Sensitivity.xlsx', sort = True):
        """Export sensitivity data of all the integral relative sensitivities 
        from a Sensitivities instance. Generates diffirent sheet for each 
        functional.

        Parameters
        ----------
        name: str
            Path where to save the sensitivities, e.g., 'MOX3600_sens.xlsx'.
            The default value is 'Sensitivity.xlsx'.
        sort: bool, optional
            Sort by absolute sensitivities. The default value is 'True'.
        """  
   
        for functional in self.functionals:
            sensitivity_df = self.to_dataframe(functional, sort)

            if os.path.exists(name):
                with pd.ExcelWriter(name, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:  
                    sensitivity_df.to_excel(writer, sheet_name = functional)
            else:
                with pd.ExcelWriter(name, engine='openpyxl', mode='w') as writer:
                    sensitivity_df.to_excel(writer, sheet_name = functional)

class Sensitivity():
    """Contains a sensitivity vector and corresponding information for a functional.

    Attributes
    ----------
    functional : str
        Name of the functional (output parameter), e.g. eigenvalue, beta-eff, etc.
    zam : int
        ZAM number for the corresponding nuclide acorrding to
        the following formula: ZAM = Z*10000 + A*10 + M. 
    reaction : int
        MT number for a reaction from the MTS list.
    sensitivty : float
        Integral relative sensitivity (can be interpreted as a
        sum of each group senseitivty).
    uncertainty : float
        Statistical uncertainty of the sensitivity attribute.
    group_structure : numpy.ndarray
        Energy grid for the sensitivity vector. Contains the number of
        elements equal to G, the number of energy groups, + 1.
    sensitivity_vector : np.ndarray
        Sensitivity vector, i.e. array of relative sensetivities for each
        group.
    uncertainties : numpy.ndarray
        Array of statistical uncertainties for each value in the list of the 
        sensitivity_vector attribute.

    """

    def __init__(self):
        self._functional = None
        self._zam = None
        self._reaction = None
        self._sensitivity = None
        self._uncertainty = None
        self._group_structure = []
        self._sensitivity_vector = []
        self._uncertainties = None

    def __repr__(self):
        return (f"{self.__class__.__name__}({self.functional!r}, {self.zam!r}, {self.reaction!r}, {(len(self.group_structure)-1)!r})")

    @property
    def functional(self):
        return self._functional
    
    @functional.setter
    def functional(self, functional):
        self._functional = functional

    @property
    def zam(self):
        return self._zam

    @zam.setter
    def zam(self, zam):
        self._zam = zam

    @property
    def reaction(self):
        return self._reaction

    @reaction.setter
    def reaction(self, reaction):
        self._reaction = reaction

    @property
    def sensitivity(self):
        return self._sensitivity

    @sensitivity.setter
    def sensitivity(self, sensitivity):
        self._sensitivity = sensitivity        

    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, uncertainty):
        self._uncertainty = uncertainty    

    @property
    def group_structure(self):
        return self._group_structure

    @group_structure.setter
    def group_structure(self, group_structure):
        group_number = len(group_structure)
        if group_number >=2 &  group_number <= 1501:
            self._group_structure = group_structure
        else:
            raise ValueError(f'The number of group must be between 1 and 1500 due to NJOY limitations, \
                             but {group_number-1} is provided')

    @property
    def sensitivity_vector(self):
        return self._sensitivity_vector

    @sensitivity_vector.setter
    def sensitivity_vector(self, sensitivity_vector):
        self._sensitivity_vector = np.array(sensitivity_vector) 

    @property
    def uncertainty_vector(self):
        return self._uncertainty_vector

    @uncertainty_vector.setter
    def uncertainty_vector(self, uncertainty_vector):
        self._uncertainty_vector = np.array(uncertainty_vector) 

    def from_serpent(self, file):
        """Export sensitivity data of a single reaction from a Serpent output 
        file (_sens0.m) to the Sensitivity instance.

        Parameters
        ----------
        file: str
            Path to the Serpent output, e.g., './MOX3600_sens0.m'.

        Notes
        -----
        This method is not intended for multiple calls for creating
        different Sensitivity instances, because each call reads the whole
        Serpent output and processes it. A similar method of Sensitivies
        should be used instead to populate the instances in a single
        call.

        """  

        sensitivities = Serpent.read(file)

        return next(sensitivity for sensitivity in sensitivities if (sensitivity.functional == self.functional) & (sensitivity.zam == self.zam) & (sensitivity.reaction == self.reaction))

    def from_scale(self, file):
        """Export sensitivity data of a single reaction from a SCALE output 
        file (.sdf) to the Sensitivity instance.

        Parameters
        ----------
        file: str
            Relative path to the SCALE output, e.g., './MET1000.sdf'.
       
        Notes
        -----
        This method is not intended for multiple calls for creating
        different Sensitivity instances, because each call reads the whole
        SCALE output and processes it. A similar method of Sensitivies
        should be used instead to populate the instances in a single
        call.

        """  

        sensitivities = Serpent.read(file)

        return next(sensitivity for sensitivity in sensitivities if (sensitivity.functional == self.functional) & (sensitivity.zam == self.zam) & (sensitivity.reaction == self.reaction))
