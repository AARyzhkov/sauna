import numpy as np
import pandas as pd
import os
import scipy.optimize
import copy


class Analysis():
    """A static class responsible for the routines, related
    to sensitivity and uncertainty analysis. This class
    governs the intercations among sensitivities and
    uncertainties.

    """

    def __new__(cls):
        raise TypeError('A static class cannot be instantiated.')

    @classmethod
    def cov_to_unc(cls, covariance):
        """Calculate uncertainty from covariance
        with the preservation of the sign.
        
        Parameters
        ----------
        covariance: float
            Covariance between two values
        
        Returns
        -------
        float
            Uncertainty including sign of covariance

        """      

        if covariance >= 0:
            uncertainty = np.sqrt(covariance)
        else:
            uncertainty = -np.sqrt(abs(covariance))
        
        return uncertainty

    @classmethod
    def unc_to_cov(cls, uncertainty):
        """Calculate covariance from uncertainty
        with the preservation of the sign.
        
        Parameters
        ----------
        uncertainty: float
            Uncertainty between two values
        
        Returns
        -------
        float
            Covariance including sign of uncertainty

        """     
        
        if uncertainty >= 0:
            covariance = uncertainty**2
        else:
            covariance = -abs(uncertainty)**2  

        return covariance

    @classmethod
    def get_std(cls, vec_1, std_1, cov_matrix, vec_2, std_2):
        """Calculate statistical uncertainty of uncertainty based upon the
        sensitivities and their statisticl uncertainty.W

        Parameters
        ----------
        vec_1: numpy.ndarray
            The first sensitivity vector
        std_1: numpy.ndarray
            Statistical uncertainty of the first sensitivity vector
        cov_matrix: numpay.ndarray
            Covariance matrix 
        vec_2: numpy.ndarray
            The second sensitivity vector
        std_2: numpy.ndarray
            Statistical uncertainty of the second sensitivity vector

        Returns
        -------
        float
            Statistical uncertainty of the functional uncertainty

        """   

        if np.all(vec_1 ==  vec_2):
            numerator   = np.square(vec_1) @ np.square(cov_matrix) @ np.square(std_1)
            denominator = vec_1 @ cov_matrix @ vec_1
        else:
            numerator = np.square(std_1) @ np.square(cov_matrix) @ np.square(vec_2) + np.square(vec_1) @ np.square(cov_matrix) @ np.square(std_2)
            denominator = 4 * vec_1 @ cov_matrix @ vec_2

        std = np.sqrt(abs(numerator / denominator))

        return std     

    @classmethod
    def get_individual_std(cls, sensitivity_1, std_1, covariance, sensitivity_2, std_2):
        """Calculate statistical uncertainty of uncertainty based upon the
        sensitivities and their statisticl uncertainty.

        Parameters
        ----------
        sensitivity_1: float
            The first sensitivity value
        std_1: float
            Statistical uncertainty of the first sensitivity value
        cov_matrix: float
            Covariance value 
        sensitivity_2: float
            The second sensitivity value
        std_2: float
            Statistical uncertainty of the second sensitivity value

        Returns
        -------
        float
            Statistical uncertainty of the functional uncertainty

        """   

        if sensitivity_1 * covariance * sensitivity_2 != 0:
            if (sensitivity_1 ==  sensitivity_2):
                std = np.sqrt(np.abs(covariance)) * std_1
            else:
                numerator = np.square(std_1*covariance*sensitivity_2) + np.square(sensitivity_1*covariance*std_2)
                denomerator = abs(4*sensitivity_1*covariance*sensitivity_2)
                std = np.sqrt(abs(numerator / denomerator))
        else: 
            std = 0

        return std    

    @classmethod
    def sandwich(cls, sensitivity_1, covariance, sensitivity_2):
        """Propagate a functional uncertainty via the Sandwich rule.
         
        Parameters
        ----------
        sensitivity_1 : Sensitivity
            First Sensitivity instance
        covariance : Covariance
            Covariance instance
        sensitivity_2 : Sensitivity
            Second Sensitivity instance
    
        Return
        ------
        float, float
            Return a tuple of the functional uncertainty and its statistical uncertainty.

        """

        # Get sensitivity vectors and its uncertainties 
        vec_1, std_1 = sensitivity_1.sensitivity_vector, sensitivity_1.uncertainty_vector
        vec_2, std_2 = sensitivity_2.sensitivity_vector, sensitivity_2.uncertainty_vector

        cov_matrix = covariance.dataframe.to_numpy()

        covariance = vec_1 @ cov_matrix @ vec_2

        if covariance != 0: 
            std = cls.get_std(vec_1, std_1, cov_matrix, vec_2, std_2)
        else:
            std = 0

        # Artificially set std no more than uncertainty
        uncertainty = cls.cov_to_unc(covariance)
        if std > abs(uncertainty):
            std = abs(uncertainty)

        return uncertainty, std

    @classmethod
    def compare(cls, sensitivities_1, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'E', covariances = None, reactions = None, save_to = None):
        """Calculate similarity indices
        
        Parameters
        ----------
        sensitivities_1 : Sensitivities
            Sensitivities instance for the model (application)
        sensitivities_2 : Sensitivities
            Sensitivities instance for another model (experiment) that
            a comparision is conducted with
        functional_1 : str, optional
            First functional name. The defualt values is 'Eigenvalue'.
        functional_2 : str, optional 
            Second functional name. The defualt values is 'Eigenvalue'.
        type: str, optional
            Type of the similarity index to compare. It supports the
            {'E', c_k', 'G'} indices. . The defualt values is 'E'.
        covariances: Covariances, optional
            Covariances instance to use for similarity assessment,
            it must be specified if only the c_k type of similarity
            is calculated, otherwise it is not used.
        reactions: list, optional
            List of MT numbers to take into account for similarity assessment.
            None is default and accounts all the reactions. This argument dictates to
            take only certain MT numbers. It is currently applicable only for G.
        save_to : str, optional
            Path to save the breakdown for c_k. The default value is 'c_k.xlsx'

        Returns
        -------
        float
            Similarity index according to the set type

        Notes
        -----
        Inherently provides uncertainties for each system

        """   

        # Check if covariances are provided
        if (type == 'c_k') & (covariances == None):
            raise ValueError(f"Covariances for the c_k method must be provided.")

        if type == "E":

            numerator = 0
            denominator_a = 0
            denominator_e = 0

            for sensitivity_a in sensitivities_1.get_by_functional(functional_1):
                sensitivity_e = sensitivities_2.get_by_params(functional_2, sensitivity_a.zam, sensitivity_a.reaction)
                numerator += sensitivity_a.sensitivity_vector * sensitivity_e.sensitivity_vector
                denominator_a += sensitivity_a.sensitivity_vector * sensitivity_a.sensitivity_vector
            for sensitivity_e in sensitivities_2.get_by_functional(functional_2):
                denominator_e += sensitivity_e.sensitivity_vector * sensitivity_e.sensitivity_vector

            return np.sum(numerator) / np.sqrt(np.sum(denominator_a) * np.sum(denominator_e))   
        
        elif type == "G":   
            numerator = 0
            denominator = 0
            if reactions == None:
                for sensitivity_1 in sensitivities_1.get_by_functional(functional_1):           
                
                    vec_a = sensitivity_1.sensitivity_vector
                    vec_e = sensitivities_2.get_by_params(functional_2, sensitivity_1.zam, sensitivity_1.reaction).sensitivity_vector
                    
                    for g in range(len(vec_a)):

                        denominator += np.abs(vec_a[g])

                        if (np.sign(vec_a[g]) != np.sign(vec_e[g])):
                            numerator += np.abs(vec_a[g])
                        elif (np.abs(vec_a[g]) >= np.abs(vec_e[g])):
                            numerator += np.abs(vec_a[g] - vec_e[g])
                    
            else:
                sensitivities = [sensitivity for sensitivity in sensitivities_1.sensitivities if (sensitivity.functional == functional_1)  & (sensitivity.reaction in reactions)]
                for sensitivity_1 in sensitivities:
                    vec_a = sensitivity_1.sensitivity_vector
                    vec_e = sensitivities_2.get_by_params(functional_2, sensitivity_1.zam, sensitivity_1.reaction).sensitivity_vector                    
                    
                    for g in range(len(vec_a)):

                        denominator += np.abs(vec_a[g])

                        if (np.sign(vec_a[g]) != np.sign(vec_e[g])):
                            numerator += np.abs(vec_a[g])
                        elif (np.abs(vec_a[g]) >= np.abs(vec_e[g])):
                            numerator += np.abs(vec_a[g] - vec_e[g])                    
                    
            index = 1 - numerator/denominator 
        
            return index 
        elif type == "c_k":

            # Create an uncertainty dataframe to populated
            c_k_df = pd.DataFrame(columns = ['Nuclide 1', 'Reaction 1', 'Nuclide 2', 'Reaction 2', 'Individual c_k'])

            numerator = 0
            denominator_a = 0
            denominator_e = 0
            
            # Get the list of all zams containing in both the application and 
            # experimental Sensitivities instances
            zams = sensitivities_1.zams + list(set(sensitivities_2.zams) - set(sensitivities_1.zams))
            for zam in zams:
                covs = covariances.get_by_zam(zam)

                for cov in covs:
                    zam_1 = cov.zam_1
                    zam_2 = cov.zam_2
                    first_mt  = cov.reaction_1
                    second_mt = cov.reaction_2

                    if (first_mt != 1) & (second_mt != 1) & (first_mt != 455) & (first_mt != 456) & (second_mt != 455) & (second_mt != 456):
                        sensitivity_a1 = sensitivities_1.get_by_params(functional_1, zam_1, first_mt)
                        sensitivity_a2 = sensitivities_1.get_by_params(functional_1, zam_2, second_mt)

                        sensitivity_e1 = sensitivities_2.get_by_params(functional_2, zam_1, first_mt)
                        sensitivity_e2 = sensitivities_2.get_by_params(functional_2, zam_2, second_mt)
                        
                        uncertainty_ae, _ = cls.sandwich(sensitivity_a1, cov, sensitivity_e2)
                        cov_ae = cls.unc_to_cov(uncertainty_ae)
                        
                        numerator += cov_ae
                        individual_contribution = cov_ae

                        uncertainty_a, _  = cls.sandwich(sensitivity_a1, cov, sensitivity_a2)
                        cov_a = cls.unc_to_cov(uncertainty_a)
                        denominator_a += cov_a

                        uncertainty_e, _  = cls.sandwich(sensitivity_e1, cov, sensitivity_e2)
                        cov_e = cls.unc_to_cov(uncertainty_e)
                        denominator_e += cov_e

                        if (first_mt != second_mt) | (zam_1 != zam_2):       
                           uncertainty_ea, _ = cls.sandwich(sensitivity_e1, cov, sensitivity_a2)
                           cov_ea = cls.unc_to_cov(uncertainty_ea)
                           numerator += cov_ea
                           individual_contribution += cov_ea
                           denominator_a += cov_a
                           denominator_e += cov_e
                        
                        # Create a new row to append it to the dataframe
                        new_row =   {'Nuclide 1'      : f'{sensitivity_a1.zam}',
                                     'Reaction 1'     : f'MT{sensitivity_a1.reaction}',
                                     'Nuclide 2'      : f'{sensitivity_a2.zam}',
                                     'Reaction 2'     : f'MT{sensitivity_a2.reaction}',
                                     'Individual c_k' : individual_contribution}   
                        
                        c_k_df.loc[len(c_k_df)] = new_row

            denominator = np.sqrt(denominator_a*denominator_e)            
            index = numerator/denominator

            # Add the total uncertainty row 
            total_c_k = np.sum(c_k_df['Individual c_k'])

            total_row = {'Nuclide 1'      : 'total',
                         'Reaction 1'     : 'total',
                         'Nuclide 2'      : 'total',
                         'Reaction 2'     : 'total',
                         'Individual c_k' : total_c_k,}   
            c_k_df.loc[len(c_k_df)] = total_row

            # Sort the uncertainties by its absolute values
            c_k_df = c_k_df.reindex(c_k_df['Individual c_k'].abs().sort_values(ascending=False).index).reset_index(drop=True)    
            c_k_df['Individual c_k'] = c_k_df['Individual c_k'].div(denominator)
            
            # Export to Excel
            if save_to != None:
                name = save_to
                if os.path.exists(name):
                    with pd.ExcelWriter(name, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:  
                        c_k_df.to_excel(writer)
                else:
                    with pd.ExcelWriter(name, engine='openpyxl', mode='w') as writer:
                        c_k_df.to_excel(writer) 

            
            print(f"{functional_1} uncertainty of the first  model is {np.sqrt(denominator_a)*100:2f}%.")
            print(f"{functional_2} uncertainty of the second model is {np.sqrt(denominator_e)*100:2f}%.")
            print(f"The ck index is equal to {index}.")
            if denominator_a < denominator_e:
                print(f"The possible reduction is approximately equal to {np.sqrt(denominator_a*(1-index**2))*100:2f}%.")
            else:
                print(f"The application uncertainty is less than or equal to experimental.")
            return index
        elif type == "c_r":
            raise NotImplementedError(f"The set method of comparison '{type}' is not implemented yet.")
        elif type == "R":
            raise NotImplementedError(f"The completeness parameter '{type}' is not implemented yet.")
        elif type == "Z":
            raise NotImplementedError(f"The penalty assessment '{type}' is not implemented yet.")
        else:
            raise NotImplementedError(f"The set method of comparison '{type}' is not implemented.")
  
    @classmethod
    def get_breakdown(cls, sensitivities, covariances, save_to = None, by_total = False):
        """Propagate the uncertainties of the functionals in the Sensitivities
        instance via the Sandwich rule and produce the dataframes with the
        uncertainty breakdown by uncertainty sources.
         
        Parameters
        ----------
        sensitivities : Sensitivities
            Sensitivities instance
        covariances : Covariances
            Covariances instance
        save_to : str, optioanl
            Path to save the breakdown
        by_total: bool, optional
            Whether to compute the uncertainty by total cross section. 
            The default values is False.

        Return
        ------
        dict
            Return a dictionary in the form {str : pandas.DataFrame}
            where the key is the functional, and item is the uncertainty breakdown. 

        Notes
        -----
        The implementation currently does not account the inherent decay constant 
        uncertainty for lambda-eff — only impact from the data in MF31-35.

        """

        # Create a dictionary of dataframes for different functionals
        dataframes = {}

        # Populate the uncertainty dataframe for each nuclide and reaction
        # Get necessary functionals from sensetivities
        for functional in sensitivities.functionals:

            # Create an uncertainty dataframe to populated
            uncertainty_df = pd.DataFrame(columns = ['Nuclide 1', 'Reaction 1', 'Nuclide 2', 'Reaction 2', 'Uncertainty [%]', 'Statistical Uncertainty [%]'])

            # Get necessary zams from sensetivities
            for zam in sensitivities.zams:

                # Get covs based upon the zams
                covs = covariances.get_by_zam(zam)

                # Get sensitivities for each cov and propagate the uncertainties
                for cov in covs:
                    zam_1 = cov.zam_1
                    zam_2 = cov.zam_2
                    first_mt  = cov.reaction_1
                    second_mt = cov.reaction_2
                    if (first_mt != 1) & (second_mt != 1) & (by_total == False):
                        if (first_mt != 455) & (first_mt != 456) & (second_mt != 455) & (second_mt != 456):

                            # Get Sensitivity instances for the corresponding covariance matrix
                            sensitivity_1 = sensitivities.get_by_params(functional, zam_1, first_mt)
                            sensitivity_2 = sensitivities.get_by_params(functional, zam_2, second_mt)

                            # Get the uncertainty intrduced by given reactions to a functional
                            uncertainty, std = cls.sandwich(sensitivity_1, cov, sensitivity_2)

                            # Account that cross-reaction covariance must be doubled 
                            if (first_mt != second_mt) | (zam_1 != zam_2):
                                uncertainty, std = uncertainty*np.sqrt(2), std*np.sqrt(2)

                            # Create a new row to append it to the dataframe
                            new_row =   {'Nuclide 1'      : f'{sensitivity_1.zam}',
                                         'Reaction 1'     : f'MT{sensitivity_1.reaction}',
                                         'Nuclide 2'      : f'{sensitivity_2.zam}',
                                         'Reaction 2'     : f'MT{sensitivity_2.reaction}',
                                         'Uncertainty [%]': uncertainty * 100,
                                         'Statistical Uncertainty [%]' : std * 100}   

                            # Populate the uncertainty dataframe for each reaction
                            uncertainty_df.loc[len(uncertainty_df)] = new_row

                        # Account delayed neutron fraction if present for beta-eff and lambda-eff
                        # with no prompt data since they can be accounted whether via both prompt or total
                        # nu-fission, and the second one is taken here
                        elif (('beta-eff' in functional) | ('lambda-eff' in functional)) & (first_mt == 455 | second_mt == 455) & (first_mt != 456 | second_mt != 456):
                            # Get Sensitivity instances for the corresponding covariance matrix
                            sensitivity_1 = sensitivities.get_by_params(functional, zam_1, first_mt)
                            sensitivity_2 = sensitivities.get_by_params(functional, zam_2, second_mt)

                            # Get the uncertainty intrduced by given reactions to a functional
                            uncertainty, std = cls.sandwich(sensitivity_1, cov, sensitivity_2)

                            # Account that cross-reaction covariance must be doubled 
                            if (first_mt != second_mt) | (zam_1 != zam_2):
                                uncertainty, std = uncertainty*np.sqrt(2), std*np.sqrt(2)

                            # Create a new row to append it to the dataframe
                            new_row =   {'Nuclide 1'      : f'{sensitivity_1.zam}',
                                         'Reaction 1'     : f'MT{sensitivity_1.reaction}',
                                         'Nuclide 2'      : f'{sensitivity_2.zam}',
                                         'Reaction 2'     : f'MT{sensitivity_2.reaction}',
                                         'Uncertainty [%]': uncertainty * 100,
                                         'Statistical Uncertainty [%]' : std * 100}   

                            # Populate the uncertainty dataframe for each reaction
                            uncertainty_df.loc[len(uncertainty_df)] = new_row

                    elif (((first_mt == 1) & (second_mt == 1)) | ((first_mt == 1018) & (second_mt == 1018)) | ((first_mt == 251) & (second_mt == 251))| ((first_mt == 452) & (second_mt == 452))) & (by_total == True):
                        sensitivity_1 = sensitivities.get_by_params(functional, zam_1, first_mt)
                        sensitivity_2 = sensitivities.get_by_params(functional, zam_2, second_mt)
                        uncertainty, std = cls.sandwich(sensitivity_1, cov, sensitivity_2)

                        # Create a new row to append it to the dataframe
                        new_row =   {'Nuclide 1'      : f'{sensitivity_1.zam}',
                                     'Reaction 1'     : f'MT{sensitivity_1.reaction}',
                                     'Nuclide 2'      : f'{sensitivity_2.zam}',
                                     'Reaction 2'     : f'MT{sensitivity_2.reaction}',
                                     'Uncertainty [%]': uncertainty * 100,
                                     'Statistical Uncertainty [%]' : std * 100}   

                        # Populate the uncertainty dataframe for each reaction
                        uncertainty_df.loc[len(uncertainty_df)] = new_row

            # Add the total uncertainty row 
            total_uncertainty = np.sqrt(np.sum(i*i if i >=0 else -i*i for i in uncertainty_df['Uncertainty [%]']))
            stat_uncertainty  = np.sqrt(np.sum(uncertainty_df['Statistical Uncertainty [%]'][i]**2 * uncertainty_df['Uncertainty [%]'][i]**2 / total_uncertainty**2 for i in range(len(uncertainty_df['Statistical Uncertainty [%]']))))

            total_row = {'Nuclide 1'   : 'total',
                         'Reaction 1'  : 'total',
                         'Nuclide 2'   : 'total',
                         'Reaction 2'  : 'total',
                         'Uncertainty [%]' : total_uncertainty,
                         'Statistical Uncertainty [%]' : stat_uncertainty}   
            uncertainty_df.loc[len(uncertainty_df)] = total_row

            # Sort the uncertainties by its absolute values
            uncertainty_df = uncertainty_df.reindex(uncertainty_df['Uncertainty [%]'].abs().sort_values(ascending=False).index).reset_index(drop=True)       
            dataframes[functional] = uncertainty_df
        
            # Export to Excel
            if save_to != None:
                name = save_to
                if os.path.exists(name):
                    with pd.ExcelWriter(name, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:  
                        uncertainty_df.to_excel(writer, sheet_name = functional)
                else:
                    with pd.ExcelWriter(name, engine='openpyxl', mode='w') as writer:
                        uncertainty_df.to_excel(writer, sheet_name = functional) 

        return dataframes

    @classmethod
    def get_detailed(cls, sensitivities, covariances, save_to = None, by_total = False):
        """Get detailed (groupwise) uncertainty breakdown. 
         
        Parameters
        ----------
        sensitivities : Sensitivities
            Sensitivities instance
        covariances : Covariances
            Covariances instance
        save_to : str, optional
            Path to save the breakdown
        by_total: bool, optional
            Whether to compute the uncertainty by total cross section. 
            The default values is False.

        Return
        ------
        dict
            Return a dictionary in the form {str : pandas.DataFrame}
            where the key is the functional, and item is the uncertainty breakdown.

        Notes
        -----
        It works slow due to many lines when impoorting to Excel, which is the 
        bottleneck here.
        The method currently does not account the inherent decay constant 
        uncertainty for lambda-eff — only impact from the data in MF31-35.

        """
        
        # Create a dictionary of dataframes for different functionals
        dataframes = {}
        group_number = len(sensitivities.group_structure) - 1
        indices = np.arange(group_number)
        sqrt_two = np.sqrt(2)
        def get_rows(functional, zam_1, zam_2, first_mt, second_mt):
            """An auxiliary method to populate the dataframe 

            Parameters
            ----------
            functional : str
                Functional name
            zam_1 : int
                ZAM value
            zam_2 : int
                ZAM value    
            first_mt : int
                First MT number
            second_mt : int
                Second MT number

            Return
            ------
            dict
                Return a dictionary to append them to the dataframe

            """
            # Get Sensitivity instances for the corresponding covariance matrix
            sensitivity_1 = sensitivities.get_by_params(functional, zam_1, first_mt)
            sensitivity_2 = sensitivities.get_by_params(functional, zam_2, second_mt)
            is_mts_equal = (first_mt == second_mt)
            vec_1 = sensitivity_1.sensitivity_vector
            vec_2 = sensitivity_2.sensitivity_vector
            std_1 = sensitivity_1.uncertainty_vector
            std_2 = sensitivity_2.uncertainty_vector
            cov_matrix = cov.dataframe.to_numpy()

            temp_rows = []
            for i in indices:
                for j in indices:
                    uncertainty = cls.cov_to_unc(vec_1[i] * cov_matrix[i][j] * vec_2[j])
                    std = cls.get_individual_std(vec_1[i], std_1[i], cov_matrix[i][j], vec_2[j], std_2[j])

                    # Artificially set std not more than uncertainty
                    if std > abs(uncertainty):
                        std = abs(uncertainty)

                    # Account that cross-reaction covariance must be doubled 
                    if not is_mts_equal:
                         uncertainty, std = uncertainty * sqrt_two, std * sqrt_two

                    # Create a new row to append it to the dataframe
                    new_row =   {'Nuclide 1'      : f'{zam_1}',
                                 'Reaction 1'     : f'MT{first_mt}',
                                 'Group 1'        : group_number - i,
                                 'Nuclide 2'      : f'{zam_2}',
                                 'Reaction 2'     : f'MT{second_mt}',
                                 'Group 2'        : group_number - j,
                                 'Uncertainty [%]': uncertainty * 100,
                                 'Statistical Uncertainty [%]' : std * 100}
                    temp_rows.append(new_row)
            return  temp_rows

        # Populate the uncertainty dataframe for each nuclide and reaction
        # Get necessary functionals from sensetivities
        if by_total == False:
            for functional in sensitivities.functionals:
                rows = []
                if ('beta-eff' == functional) | ('lambda-eff' == functional):
                    for zam in sensitivities.zams:
                        covs = covariances.get_by_zam(zam)
                        for cov in covs:
                            zam_1 = cov.zam_1
                            zam_2 = cov.zam_2
                            first_mt  = cov.reaction_1
                            second_mt = cov.reaction_2

                            if (first_mt != 1) & (second_mt != 1):
                                is_nu_t = (first_mt != 455) & (first_mt != 456) & (second_mt != 455) & (second_mt != 456)
                                if is_nu_t: 
                                    rows.extend(get_rows(functional, zam_1, zam_2, first_mt, second_mt)) 
                                elif (first_mt == 455 | second_mt == 455) & (first_mt != 456 | second_mt != 456): 
                                    rows.extend(get_rows(functional, zam_1, zam_2, first_mt, second_mt))   
                else:
                    for zam in sensitivities.zams:
                        covs = covariances.get_by_zam(zam)
                        for cov in covs:
                            zam_1 = cov.zam_1
                            zam_2 = cov.zam_2
                            first_mt  = cov.reaction_1
                            second_mt = cov.reaction_2

                            if (first_mt != 1) & (second_mt != 1):
                                is_nu_t = (first_mt != 455) & (first_mt != 456) & (second_mt != 455) & (second_mt != 456)
                                if is_nu_t: rows.extend(get_rows(functional, zam_1, zam_2, first_mt, second_mt))  
                        # Create and populated the final dataframe
                uncertainty_df = pd.DataFrame(rows, columns = ['Nuclide 1', 'Reaction 1', 'Group 1', 'Nuclide 2', 'Reaction 2', 'Group 2', 'Uncertainty [%]', 'Statistical Uncertainty [%]'])                             
                
                # Add the total uncertainty row 
                total_uncertainty = np.sqrt(np.sum(i*i if i >=0 else -i*i for i in uncertainty_df['Uncertainty [%]']))
                stat_uncertainty  = np.sqrt(np.sum(uncertainty_df['Statistical Uncertainty [%]'][i]**2 * uncertainty_df['Uncertainty [%]'][i]**2 / total_uncertainty**2 for i in range(len(uncertainty_df['Statistical Uncertainty [%]']))))

                total_row = {'Nuclide 1'   : 'total',
                            'Reaction 1'  : 'total',
                            'Group 1'     : 'total',
                            'Nuclide 2'   : 'total',
                            'Reaction 2'  : 'total',
                            'Group 2'     : 'total',
                            'Uncertainty [%]' : total_uncertainty,
                            'Statistical Uncertainty [%]' : stat_uncertainty} 
                
                uncertainty_df.loc[len(uncertainty_df)] = total_row
                
                # Sort the uncertainties by its absolute values
                uncertainty_df = uncertainty_df.reindex(uncertainty_df['Uncertainty [%]'].abs().sort_values(ascending=False).index).reset_index(drop=True)       
                dataframes[functional] = uncertainty_df
        else:
            for functional in sensitivities.functionals:
                rows = []
                if ('beta-eff' == functional) | ('lambda-eff' == functional):
                    for zam in sensitivities.zams:
                        covs = covariances.get_by_zam(zam)
                        for cov in covs:
                            zam_1 = cov.zam_1
                            zam_2 = cov.zam_2
                            first_mt  = cov.reaction_1
                            second_mt = cov.reaction_2

                            if ((first_mt == 1) & (second_mt == 1)) | ((first_mt == 1018) & (second_mt == 1018)) | ((first_mt == 251) & (second_mt == 251)) | ((first_mt == 452) & (second_mt == 452)):
                                rows.extend(get_rows(functional, zam_1, zam_2, first_mt, second_mt))  
                else:
                    for zam in sensitivities.zams:
                        covs = covariances.get_by_zam(zam)
                        for cov in covs:
                            zam_1 = cov.zam_1
                            zam_2 = cov.zam_2
                            first_mt  = cov.reaction_1
                            second_mt = cov.reaction_2

                            if ((first_mt == 1) & (second_mt == 1)) | ((first_mt == 1018) & (second_mt == 1018)) | ((first_mt == 251) & (second_mt == 251)) | ((first_mt == 452) & (second_mt == 452)):
                                rows.extend(get_rows(functional, zam_1, zam_2, first_mt, second_mt))  
                # Create and populated the final dataframe
                uncertainty_df = pd.DataFrame(rows, columns = ['Nuclide 1', 'Reaction 1', 'Group 1', 'Nuclide 2', 'Reaction 2', 'Group 2', 'Uncertainty [%]', 'Statistical Uncertainty [%]'])                             
                
                # Add the total uncertainty row 
                total_uncertainty = np.sqrt(np.sum(i*i if i >=0 else -i*i for i in uncertainty_df['Uncertainty [%]']))
                stat_uncertainty  = np.sqrt(np.sum(uncertainty_df['Statistical Uncertainty [%]'][i]**2 * uncertainty_df['Uncertainty [%]'][i]**2 / total_uncertainty**2 for i in range(len(uncertainty_df['Statistical Uncertainty [%]']))))

                total_row = {'Nuclide 1'   : 'total',
                            'Reaction 1'  : 'total',
                            'Group 1'     : 'total',
                            'Nuclide 2'   : 'total',
                            'Reaction 2'  : 'total',
                            'Group 2'     : 'total',
                            'Uncertainty [%]' : total_uncertainty,
                            'Statistical Uncertainty [%]' : stat_uncertainty} 
                
                uncertainty_df.loc[len(uncertainty_df)] = total_row

                # # Sort the uncertainties by its absolute values
                uncertainty_df = uncertainty_df.reindex(uncertainty_df['Uncertainty [%]'].abs().sort_values(ascending=False).index).reset_index(drop=True)       
                dataframes[functional] = uncertainty_df

        # Export to Excel
        if save_to != None:
            with pd.ExcelWriter(save_to, engine='openpyxl', mode='w') as writer:
                for sheet, df in dataframes.items():
                    print('Importing the dataframe for: {sheet}...')
                    df.to_excel(writer, sheet_name = sheet) 

        return dataframes

    @classmethod
    def get_concentration_uncertainty(cls, sensitivities, uncertainty, targets, background_zams, fraction, fraction_type = 'ao', uncertainty_type = 'normal'):
        """Propagate the uncertainties of the functionals in the Sensitivities
        instance based upon the idea that macroscopic XS are strictly defined
        as a product of the nuclide concentration and the corresponding total
        microscopic XS, and it does not matter whether it is the sensitivity
        to the total XS or the concentration.
         
        Parameters
        ----------
        sensitivities : Sensitivities
            Sensitivities instance
        uncertainty : float
            Uncertainty in the concentration: 68% of the confidence interval
            for the 'normal' distribution or 1/2 of the range of the 
            'uniform' and 'triangular' distributions.
        targets: numpy.ndarray of ints
            ZAM values for the nuclide, which uncertainty
            is provided in the uncertainty parameter.
        background_zams: numpy.ndarray
            Nuclides to be considered in the calculation and exclude the target nuclide. 
            Make sure that correct nuclides (ZAMs) for the region, where the influence
            of uncertainties in concentrations is assessed. For example, ZAMs
            for cladding must not be provided if the Pu239 fraction uncecrtainty
            influence in the HM fuel composition is of interest. 
        fraction: float
            Fraction of the target nuclide among the nuclides parameter.
        fraction_type: string, optional
            Type of the fraction which the uncertainty is calculated for.
            It support only the atomic fraction 'ao', which is the same for
            concentrations, and the weight fraction 'wo'.  The default
            value is 'ao'.
        uncertainty_type: string, optional
            Type of the uncertainty distribution. It is assumed as
            one of the following {'normal', 'uniform', 'triangular'}.
            The default value is 'normal'.

        Return
        ------
        dict
            Return a dictionary in the form {str : float}
            where the key is the functional, and item is the uncertainty
            based upon constrained sensitivity.
  
        """

        nuclides = background_zams

        # Make sure that target values are not included in the nuclides list
        for target in targets:
            if target in nuclides:
                nuclides.remove(target)
                print(f'{target} has been removed from the background_zams array.')

        concentration_uncertainties = {}

        for functional in sensitivities.functionals:
            
            unconstrained_sensitivity = 0
            for target in targets:
                unconstrained_sensitivity += sensitivities.get_by_params(functional, target, 1).sensitivity

            constrainer = 0
            for nuclide in nuclides:
                constrainer += sensitivities.get_by_params(functional, nuclide, 1).sensitivity

            # Get constrained sensitivities:
            if fraction_type == 'ao':
                constrained_sensitivity = unconstrained_sensitivity - fraction/(1-fraction)*constrainer
                concentration_uncertainty = np.abs(constrained_sensitivity * uncertainty)
                concentration_uncertainties[functional] = concentration_uncertainty
            elif fraction_type == 'wo':
                constrained_sensitivity = (unconstrained_sensitivity - constrainer) / (1 - fraction)
                concentration_uncertainty = np.abs(constrained_sensitivity*uncertainty)
                concentration_uncertainties[functional] = concentration_uncertainty
            else:
                raise NotImplementedError(f"The fraction type '{fraction_type} is not implemented', only 'ao' and 'wo' are supported.")
            
        if uncertainty_type == 'normal':
            pass
        elif uncertainty_type == 'uniform':
            for functional_uncertainty in concentration_uncertainties.keys():
                concentration_uncertainties[functional_uncertainty] /= (2*np.sqrt(3))
        elif uncertainty_type == 'triangular':
            for functional_uncertainty in concentration_uncertainties.keys():
                concentration_uncertainties[functional_uncertainty] /= np.sqrt(6)
        else:
            raise NotImplementedError(f"The distribution type '{uncertainty_type}' is not implemented, only 'normal', 'uniform', and 'triangular' are supported.")
            
        return concentration_uncertainties

    @classmethod
    def get_density_uncertainty(cls, sensitivities, uncertainty, targets, uncertainty_type = 'normal'):
        """Propagate the uncertainties of the functionals in the Sensitivities
        instance based upon the idea that macroscopic XS are strictly defined
        as a product of the nuclide concentration and the corresponding total
        microscopic XS, and it does not matter whether it is the sensitivity
        to the total XS or the concentration.
         
        Parameters
        ----------
        sensitivities : Sensitivities
            Sensitivities instance
        uncertainty : float
            Uncertainty in the concentration: 68% of the confidence interval
            for the 'normal' distribution or 1/2 of the range of the 
            'uniform' or 'triangular' distributions
        targets: numpy.ndarray of ints
            ZAM values for the nuclide, which uncertainty
            is provided in the uncertainty parameter. Make sure the
            target nuclides are located in the volume of interest
        uncertainty_type: string, optional
            Type of the uncertainty distribution. It is assumed as
            one of the following {'normal', 'uniform', 'triangular'}.
            The default value is 'normal'.

        Return
        ------
        dict
            Return a dictionary in the form {str : float}
            where the key is the functional, and item is the uncertainty
            based upon total sensitivity.
        
        """

        uncertainties = {}

        for functional in sensitivities.functionals:
            sensitivity = 0
            for target in targets:
                sensitivity += sensitivities.get_by_params(functional, target, 1).sensitivity
            
            uncertainty = np.abs(sensitivity * uncertainty)
            uncertainties[functional] = uncertainty

        if uncertainty_type == 'normal':
            pass
        elif uncertainty_type == 'uniform':
            for functional_uncertainty in uncertainties.keys():
                uncertainties[functional_uncertainty] /= (2*np.sqrt(3))
        elif uncertainty_type == 'triangular':
            for functional_uncertainty in uncertainties.keys():
                uncertainties[functional_uncertainty] /= np.sqrt(6)
        else:
            raise NotImplementedError(f"The distribution type '{uncertainty_type}' is not implemented, only 'normal', 'uniform', and 'triangular' are supported.")
            
        return uncertainties

    @classmethod
    def tars(cls, model_sensitivities, covariances, tars, number_of_reactions=10, lower_boundary=0.005, method='trust-constr', cost_type='A', energy_costs = None, maxiter=1000, tol = 1e-5):
        """Calculates target accuracy requirements based upon
        a Sensitivities instance, an Uncertainties instance
        and TAR values.
         
        Parameters
        ----------
        model_sensitivities : list
            List of Sensitivities instances
        covariances : Covariances
            Covariances instance
        tars: dict
            Dictionary of target accuracies in the form of
            {functional : TAR, ...}, e.g. {'Eigenvalue' : 0.005}.
        number_of_reactions : int, optional
            Number of top symmetric covariances to take into account in the
            optimization process. The other uncertainties are assumed
            to be constant. The actual number of the going to be tweaked
            covariances is more since it the parameter does not include the
            cross-reaction covariances. Especially, the actual number can be
            significantly increased if two or more functionals considered
            with different uncertainty sources (up to *len(functionals)).
            Default value is 10.
        lower_boundary : float, optional
            Lower boundary for uncertainty. Each uncertainty cannnot
            be decreased lower than this value. The current
            implementation assumes that if the value is already
            lower than the set value, the value will not be changed.
            The method does not handle lower boundary equal zero, to
            avoid that, it is recommended to provide the lower boundary
            reasonably higher than zero (i.e. is not going to be
            reached), e.g. 1e-4. The default value is 0.005.
        method : string, optional
            Method to be applied by the optimizer. Only 'trust-constr'
            is currently implemented as a more reliable approach:
            https://github.com/scipy/scipy/issues/9640#issuecomment-451918238 .
            The default value is 'trust-constr'.
        cost_type : string, optional
            Parameter identifying cost parameters {'A', 'B', 'C'}. The 
            coefficients are defined according to the NEA/OECD WPEC/SG26
            report. The cost for other parameters (not mentioned in the
            report) is set to maximum, i.e. the (n,n') cost is set to (\chi)).
            The default value is 'A', assuming there is no difference among
            reactions.
        energy_cost : list, optional
            List of energy weights identifying the cost at each energy group that
            is multipled with cost_type to get the total cost in the TAR excercise.
            The default value is None, meaning there is no energy dependence of the
            cost function.
        maxiter: float, optional
            Number of maximum iterations to be passed to the SciPy optimizer. The default
            value is 1000.
        tol: float, optional
            Tolerance value to be passed to the SciPy optimizer. The default value is
            1e-5.
            
        Return
        ------
        Covariances
            Return a Covariances instance with updated data corresponding to
            TARs.

        """   

        # The code is structured as follows:
        # 1 Create a list of uncertainties and other params for the optimization problem
        # 2 Create auxiliary functions to run the problem
        # 3 Solve the problem using the auxiliry functions
        # 4 Move data from the results to the new covariances

        if lower_boundary < 0:
            raise ValueError('Make the lower_boundary parameter more than or equal to zero.')
        

        # Variables for iteration
        group_number = len(covariances.group_structure) - 1
        base_indices = np.arange(number_of_reactions)
        
        # Prepare individual sensitivities and uncertainties to be used
        # for the optimization
        zam_senses = []
        uncertainties = []
        zams = []
        reactions = []
        zam_covs = []
        zam_mt = []
        all_reactions = set()
        dfs = []
        costs = []
        
        if energy_costs == None:
            energy_costs.extend([1]*group_number)

        # Cost Function
        def define_cost(zam, reaction, cost_type='A'):
            if (zam in (922350, 922380, 942390)) & (reaction in (18, 102, 452)):
                if cost_type == 'A':
                    cost_coefficient = 1
                elif cost_type == 'B':
                    cost_coefficient = 1
                elif cost_type == 'C':
                    cost_coefficient = 1
            elif (zam >= 900000) & (reaction in (18, 102, 452)):
                if cost_type == 'A':
                    cost_coefficient = 1
                elif cost_type == 'B':
                    cost_coefficient = 2
                elif cost_type == 'C':
                    cost_coefficient = 2
            elif (zam < 900000) & (reaction == 102):
                if cost_type == 'A':
                    cost_coefficient = 1
                elif cost_type == 'B':
                    cost_coefficient = 1
                elif cost_type == 'C':
                    cost_coefficient = 1
            elif reaction == 2:
                if cost_type == 'A':
                    cost_coefficient = 1
                elif cost_type == 'B':
                    cost_coefficient = 1
                elif cost_type == 'C':
                    cost_coefficient = 1
            elif reaction == 4:
                if cost_type == 'A':
                    cost_coefficient = 1
                elif cost_type == 'B':
                    cost_coefficient = 3
                elif cost_type == 'C':
                    cost_coefficient = 10
            else:
                if cost_type == 'A':
                    cost_coefficient = 1
                elif cost_type == 'B':
                    cost_coefficient = 3
                elif cost_type == 'C':
                    cost_coefficient = 10

            return cost_coefficient

        # Get dataframe for the uncertainty sources
        # Get MTs present in the dataframe to avoid accounting others
        # Remove cross-correlations, cross-correations are fixed during next step
        for sensitivities in model_sensitivities:
            model_dfs = {}
            for functional in tars:
                df = cls.get_breakdown(sensitivities, covariances)[functional]
                all_reactions.update(df['Reaction 1'].drop_duplicates().sort_values().str.replace('MT', '').drop(df[df['Nuclide 1']=='total'].index, inplace=False).astype(int))
                df = df.drop(df[(df['Reaction 1'] != df['Reaction 2']) | (df['Nuclide 1'] != df['Nuclide 2'])].index, inplace=False)[1:(number_of_reactions+1)].reset_index(drop=True)
                model_dfs[functional] = df
            dfs.append(model_dfs)

        # Populate list of uncertainties_to_minimize based upon the number of reactions
        for s, sensitivities in enumerate(model_sensitivities):
            for functional in tars:
                for row in base_indices:
                    mt  = int(dfs[s][functional]['Reaction 1'][row][2:])
                    zam = int(dfs[s][functional]['Nuclide 1'][row])

                    if zam not in zams:
                        zam_covs.extend(covariances.get_by_zam(zam))
                        
                    if (zam, mt) not in zam_mt:
                        zam_mt.append((zam, mt))
                        zams.append(zam)
                        reactions.append(mt)
                        costs.extend([define_cost(zam, mt, cost_type)]*group_number)


                        cov = covariances.get_by_params(zam, zam, mt, mt)

                        # Get the uncertainties from the covariance matrix
                        diag = np.diag(cov.dataframe.to_numpy())
                        uncertainties.extend(np.sqrt(diag))
        
        # Redefine energy costs 
        energy_costs = energy_costs * int(len(uncertainties)/group_number)
        costs = np.multiply(costs,energy_costs)

        for sensitivities in model_sensitivities:
            temp_senses = []
            for functional in tars:
                for zam in set(zams):
                    for reaction in all_reactions:
                        temp_senses.append(sensitivities.get_by_params(functional, zam, reaction))
            zam_senses.append(temp_senses)

        # Populate list of covs which are going to be tweaked
        zam_covs = [cov for cov in zam_covs if cov.reaction_1 in all_reactions]


        # Update indices for all the functionals present
        indices = np.arange(len(zam_mt))

        # Prepare covariance list which is going to be accessed during the
        # optimization
        temp_covs = []
        for row in indices:
            mt  = reactions[row]
            zam = zams[row]
            temp_covs.append([cov for cov in zam_covs if (((cov.zam_1 == zam) & (cov.reaction_1 == mt)) | ((cov.zam_2 == zam) & (cov.reaction_2 == mt)))])

        # Constraints for uncertainties: (unc_i)_min <= unc_i <= (unc_i)_0
        lower_boundaries = [unc if unc < lower_boundary else lower_boundary for unc in uncertainties]
        upper_boundaries = uncertainties.copy()

         # Sandwich formula functions to be constrained to TAR
        def sandwich_constraint(uncertainties_to_minimize, functional, model_index, sensitivities):
            """An auxiliary function to propagate uncertainties
            from given uncertainties, sensitivities are taken from the
            parental method - tars(cls,...).

            Parameters
            ----------
            uncertainties_to_minimize : numpy.ndarray
                values to propagate the uncertainities
            
            Return
            ------
            float
                Return the functional uncertainty 

            """  
            sum = 0

            for row in indices:
            
                mt  = reactions[row]
                zam = zams[row]
    
                base_diag = uncertainties[int(row*group_number):int((row+1)*group_number)]
                new_diag = uncertainties_to_minimize[int(row*group_number):int((row+1)*group_number)]
    
                for cov in temp_covs[row]:
                    symm_coef = 1
                    if (cov.zam_1 == zam) & (cov.reaction_1 == mt):
                        sens_1 = next(sens for sens in zam_senses[model_index] if (sens.functional == functional) & (sens.zam == zam) & (sens.reaction == mt))
                        first_base = base_diag
                        first_new = new_diag
                    else:
                        sens_1 = sensitivities.get_by_params(functional, cov.zam_1, cov.reaction_1)
                        first_base = 1
                        first_new = 1
                        symm_coef = 2

                    if (cov.zam_2 == zam) & (cov.reaction_2 == mt):
                        sens_2 = next(sens for sens in zam_senses[model_index] if (sens.functional == functional) & (sens.zam == zam) & (sens.reaction == mt))
                        second_base = base_diag
                        second_new = new_diag        
                    else:
                        sens_2 = sensitivities.get_by_params(functional, cov.zam_2, cov.reaction_2)
                        second_base = 1
                        second_new = 1
                        symm_coef = 2

                    matrix     = cov.dataframe.to_numpy()
                    base_outer = np.outer(first_base, second_base)
                    cor        = np.divide(matrix, base_outer, out=np.zeros_like(matrix), where= base_outer != 0 )
                    new_outer  = np.outer(first_new, second_new)
                    new_cov_mat  = cor * new_outer   
                    sandwich = sens_1.sensitivity_vector @ new_cov_mat @ sens_2.sensitivity_vector
                    sum +=  symm_coef * sandwich

            uncertainty = sum
            print(f'Model {model_index} {functional}:', np.sqrt(uncertainty).real)
            return uncertainty

        # Non-linear constraints for trust-constr
        tar_constraints = []
        for s, sensitivities in enumerate(model_sensitivities):
            for functional in tars:
                tar_constraint = scipy.optimize.NonlinearConstraint(lambda x, functional=functional, model_index=s, sensitivities=sensitivities: sandwich_constraint(x,functional,model_index,sensitivities) , 0, tars[functional]**2,
                                                            jac = 'cs', hess = scipy.optimize.SR1())
                tar_constraints.append(tar_constraint)
        
        # Weight functions to be minimized
        def weight_function(uncertainties_to_minimize):
            """An auxiliary function to calculates the sum of
            the weighted inverse variances weight_i/uncertainties_i**2.

            Parameters
            ----------
            uncertainties_to_minimize : numpy.ndarray
                values to propagate the uncertainities

            Return
            ------
            float
                Return the sum 

            """  

            # The function assumes that the the initial uncertainties
            # do not depend on uncertainties_to_minimize and returns
            # always zero for them. If the values are set at zero only
            # if the uncertainties_to_minimize, the optimizer yield 
            # inaccurate results
            return np.sum([0 if unc == 0 else lambda_i/target**2 for unc, target, lambda_i in zip(uncertainties, uncertainties_to_minimize, costs)])

        def weight_jacobian(uncertainties_to_minimize):
            """An auxiliary function to calculates the Jacobian of
            the sum of the weighted inverse variances, which is
            -2*weight_i/uncertainties_i**3.

            Parameters
            ----------
            uncertainties_to_minimize : numpy.ndarray
                values to propagate the uncertainities

            Return
            ------
            numpy.ndarray
                Return the Jacobian


            """  

            # The function assumes that the the initial uncertainties
            # do not depend on uncertainties_to_minimize and returns
            # always zero for them. If the values are set at zero only
            # if the uncertainties_to_minimize, the optimizer yield 
            # inaccurate results
          
            return [0 if unc == 0 else -2*lambda_i/target**3 for unc, target, lambda_i in zip(uncertainties, uncertainties_to_minimize, costs)]

        if method == 'trust-constr':
            
            bounds = scipy.optimize.Bounds(lower_boundaries, upper_boundaries)

            res = scipy.optimize.minimize(weight_function, uncertainties, method = 'trust-constr',
                                                       jac = weight_jacobian, hess = scipy.optimize.SR1(),
                                                       constraints = tar_constraints,
                                                       options = {'disp': True,  'maxiter': maxiter},
                                                       tol =  tol,
                                                       bounds = bounds)     
        else:
            raise NotImplementedError(f'The set method {method} is not implemented.')

        optimized_uncertainties = [0 if unc == 0 else target for unc, target in zip(uncertainties, res.x)]
        print(optimized_uncertainties)

        # Completely copies the main covariances to avoid rewriting them
        target_covariances = copy.deepcopy(covariances)
        
        # Get covs with the zams of interest
        zam_covs = []
        [zam_covs.extend(target_covariances.get_by_zam(zam)) for zam in set(zams)]

        changed_covariances = []
        for row in indices:
            mt  = reactions[row]
            zam = zams[row]
            changed_covariances.append((zam,mt))

            temp_covs = [cov for cov in zam_covs if ((cov.zam_1 == zam) & (cov.reaction_1 == mt)) | ((cov.zam_2 == zam) & (cov.reaction_2 == mt))]
            base_diag = uncertainties[int(row*group_number):int((row+1)*group_number)]
            new_diag = optimized_uncertainties[int(row*group_number):int((row+1)*group_number)]
            
            for cov in temp_covs:
                if (cov.zam_1 == zam) & (cov.reaction_1 == mt):
                    first_new = new_diag
                    first_base = base_diag
                else:
                    first_new = 1
                    first_base = 1

                if (cov.zam_2 == zam) & (cov.reaction_2 == mt):
                    second_new = new_diag
                    second_base = base_diag
                else:
                    second_new = 1
                    second_base = 1

                matrix     = cov.dataframe.to_numpy()
                base_outer = np.outer(first_base, second_base)
                cor        = np.divide(matrix, base_outer, out=np.zeros_like(matrix), where=base_outer != 0)
                new_outer  = np.outer(first_new, second_new)
                cov.dataframe[:] = cor * new_outer   

        base_cost = weight_function(uncertainties)
        final_cost = weight_function(optimized_uncertainties)
        cost = final_cost - base_cost

        print('The total number of symmmetric covariances tweaked:', len(zam_mt))
        print('The tweaked covariances:', changed_covariances)
        print('The base cost function:', base_cost)
        print('The cost function:', final_cost)
        print('The net cost:', cost)

        return target_covariances


















