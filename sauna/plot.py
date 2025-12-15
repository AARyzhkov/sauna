import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import shutil

usetex = True if shutil.which('latex') else False
if usetex:
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['text.usetex']  = True

from sauna.auxiliary import cov_to_corr

# The plot parameters are hardcoded as Elsevier recommends:
# https://www.elsevier.com/about/policies-and-standards/author/artwork-and-media-instructions/artwork-sizing
# The expected font size should be about 7 pt (6 pt sub/supersctips), but in some cases it may reach 10 pt
# Single column width is 9 cm (3.54 in.), and two column width is 19 cm (7.48 in.)
# 1000 dpi is set for every plot
# Some parameters are subjective since the requirements are not declared to them

LINEWIDTH = 0.5
XMIN = 1E-5
XMAX = 2E+7
LINESTYLES = ['-', '--','-.', ':']
plt.rcParams['axes.linewidth']    = 0.2
plt.rcParams['grid.linewidth']    = 0.2
plt.rcParams['xtick.major.width'] = 0.2
plt.rcParams['xtick.minor.width'] = 0.2
plt.rcParams['ytick.major.width'] = 0.2
plt.rcParams['ytick.minor.width'] = 0.2
plt.rcParams['font.size']       = 7
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['axes.labelsize']  = 7
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['axes.autolimit_mode'] = 'round_numbers'

"""Dictionary of commonly used reactions in terms of the MT numbers
to translate them to the LaTeX-like formalism :

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
- MT1002 - thermal elastic scattering laws (S(α,β)), it is not a part of the
           established MT numbers) 
- MT1004 - thermal inelastic scattering laws (S(α,β)), it is not a part of the
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

LATEX_REACTIONS = {1   : "(n,t)"      ,
                   2   : "(n,n)"      ,
                   4   : "(n,n')"     ,
                   16  : "(n,2n)"     ,
                   18  : "(n,f)"      ,
                   102 : "(n,\gamma)" ,
                   103 : "(n,p)"      ,
                   104 : "(n,D)"      ,
                   105 : "(n,T)"      ,
                   106 : "(n,^{3}He)"    ,
                   107 : "(n,\\alpha)",
                   251 : "(\langle \mu \\rangle)",
                   452 : "(\\nu)"            ,
                   455 : "(\\nu_{d})"    ,
                   456 : "(\\nu_{p})"     ,
                   1002: "(S_{el}(\\alpha,\\beta))" ,
                   1004: "(S_{in}(\\alpha,\\beta))" ,
                   1018: "(\chi)"            ,
                   1455: "(\chi_{d})"    ,
                   1456: "(\chi_{p})"}

ATOMIC_NUMBERS = { 1: 'H' ,
                   2: 'He',
                   3: 'Li',
                   4: 'Be',
                   5: 'B' ,
                   6: 'C' ,
                   7: 'N' ,
                   8: 'O' ,
                   9: 'F' ,
                  10: 'Ne',
                  11: 'Na',
                  12: 'Mg',
                  13: 'Al',
                  14: 'Si',
                  15: 'P' ,
                  16: 'S' ,   
                  17: 'Cl', 
                  18: 'Ar',
                  19: 'K' ,
                  20: 'Ca',
                  21: 'Sc',
                  22: 'Ti',
                  23: 'V' ,
                  24: 'Cr',
                  25: 'Mn',
                  26: 'Fe',
                  27: 'Co',  
                  28: 'Ni',  
                  29: 'Cu',  
                  30: 'Zn',
                  31: 'Ga',
                  32: 'Ge',
                  33: 'As',  
                  34: 'Se',  
                  35: 'Br', 
                  36: 'Kr', 
                  37: 'Rb',
                  38: 'Sr',
                  39: 'Y' ,   
                  40: 'Zr',  
                  41: 'Nb',  
                  42: 'Mo',
                  43: 'Tc',
                  44: 'Ru',
                  45: 'Rh',  
                  46: 'Pd',  
                  47: 'Ag',  
                  48: 'Cd', 
                  49: 'In',
                  50: 'Sn',
                  51: 'Sb',  
                  52: 'Te',  
                  53: 'I' ,   
                  54: 'Xe', 
                  55: 'Cs',
                  56: 'Ba',
                  57: 'La',  
                  58: 'Ce',  
                  59: 'Pr',  
                  60: 'Nd',
                  61: 'Pm',
                  62: 'Sm',
                  63: 'Eu',  
                  64: 'Gd',  
                  65: 'Tb',  
                  66: 'Dy',
                  67: 'Ho',
                  68: 'Er',
                  69: 'Tm',  
                  70: 'Yb',  
                  71: 'Lu',  
                  72: 'Hf',
                  73: 'Ta',
                  74: 'W' ,
                  75: 'Re',  
                  76: 'Os',  
                  77: 'Ir',  
                  78: 'Pt', 
                  79: 'Au',
                  80: 'Hg',
                  81: 'Tl',  
                  82: 'Pb',  
                  83: 'Bi', 
                  84: 'Po', 
                  85: 'At',
                  86: 'Rn',
                  87: 'Fr', 
                  88: 'Ra',
                  89: 'Ac', 
                  90: 'Th',
                  91: 'Pa',
                  92: 'U' ,
                  93: 'Np',  
                  94: 'Pu',  
                  95: 'Am',  
                  96: 'Cm', 
                  97: 'Bk',
                  98: 'Cf',
                  99: 'Es', 
                 100: 'Fm', 
                 101: 'Md', 
                 102: 'No',
                 103: 'Lr',
                 104: 'Rf', 
                 105: 'Db', 
                 106: 'Sg', 
                 107: 'Bh',
                 108: 'Hs',
                 109: 'Mt', 
                 110: 'Ds', 
                 111: 'Rg', 
                 112: 'Cn',
                 113: 'Nh',
                 114: 'Fl', 
                 115: 'Mc', 
                 116: 'Lv', 
                 117: 'Ts',
                 118: 'Og' }

class Plot():
    """"A static class responsible for plotting-related
    routines.
    
    """

    def __new__(cls):
        raise TypeError('A static class cannot be instantiated')
    
    @staticmethod
    def plot_sensitivity(sensitivities, name='sensitivity_profile', format = 'svg', show_integral = True, show_uncertainty  = True, normalization_type = 0, annotations = None, width=3.54):
        """Plot sensitivity profile
        
        Parameters
        ----------
        sensitivities : List
            List of the Sensitivity instances
        name : str, optional
            Name of the plot to save. The default value is
            'sensitivity_profile'.
        format : str, optional 
            Format of the plot to save. The default value is
            'svg'.
        show_integral : bool, optional
            Whether show the integral value in the legend 
            or not. The default value is True.
        show_uncertainty : bool, optional
            Whether show the statistical uncertainty. The default
            value is True.
        normalization_type : int, optional
            Types for normalizing the sensitivity profiles:
            0 - no normalization
            1 - per unit energy (can be recommended when different group structures are considered)
            2 - per unit lethargy
            The default value is 0.
        annotations : list, optional
            List of additional comments to the legened on the plot. The number
            of elements have to correspond to the sensitivities attribue length.
        width : float, optional
            The width of the plot in inches. The default value is 3.53.

        Notes
        -----
        Sensitivities with different energy structure is not 
        recommended to be compared since its each group 
        can be represented as a sum of a number of group
        with smaller sensitivities (if an extreme case is
        considered, the sensitivity at a distinct
        point of energy is equal to zero). Normalizing per 
        unit energy provides a better plot for a comparision
        than one standard or another one per unit lethargy.
        
        """      

        fig, ax = plt.subplots(figsize=(width, 2.655))
        s = 0

        # Annotation check
        a = 0
        if annotations != None:
            if len(annotations) != len(sensitivities):
                raise ValueError('The numbers of annotations and sensitivities are not the same.')        

        for sensitivity in sensitivities:
            profile = sensitivity.sensitivity_vector  
            std = sensitivity.uncertainty_vector 
            energies = sensitivity.group_structure
            reaction = sensitivity.reaction
            int_sensitivity = sensitivity.sensitivity
            int_uncertainty = sensitivity.uncertainty * 1e4

            zam = sensitivity.zam
            charge = round(zam / 1e4)
            mass_number = str(round((zam - charge*1e4)/10))
            if int(str(zam)[-1]): mass_number.append('m')
            if mass_number == '0': mass_number = '\mathrm{nat}'
            symbol = ATOMIC_NUMBERS[charge]

            # Normalize the plotting data
            if   normalization_type == 0:
                normaliziation = 1
            elif normalization_type == 1:
                normaliziation = [energies[i+1]-energies[i] for i in range(len(energies)-1)]
            elif normalization_type == 2:
                normaliziation = [energies[i+1]/energies[i] for i in range(len(energies)-1)]
            else:
                raise ValueError(f'The type of normalization was set to {normalization_type}, while 1, 2, or 3 are allowed.')

            reaction = LATEX_REACTIONS[reaction]

            label = f'$^{{{mass_number}}}${symbol}${reaction}$'

            if show_integral == True:
                if np.abs(int_sensitivity) > 10**(-3):
                    label += f'\n Integral: {int_sensitivity:.4f}({int_uncertainty:.0f})'
                else:
                    sens_exp = f'{int_sensitivity:.2e}'.split('e')
                    unc_exp  = f'{(int_uncertainty/10**float(sens_exp[1])/1e2):.0f}'.split('e')
                    if '-0' in sens_exp[1]:
                        sens_exp[1] = sens_exp[1].replace('-0', '-')
                    label += f'\n Integral: {sens_exp[0]}({unc_exp[0]})$\cdot$10$^{{{sens_exp[1]}}}$'

            # Add annotations if passed
            if annotations != None:
                label += f'{annotations[a]}'
                a += 1 

            s += 1 
            ax.stairs(profile/normaliziation, energies, label=label, alpha=0.7, linewidth = LINEWIDTH, linestyle=LINESTYLES[s%len(LINESTYLES)-1])
            if show_uncertainty == True:
                ax.fill_between(energies[1:], (profile+std)/normaliziation, (profile-std)/normaliziation, step='pre', alpha=0.2)
        
        ax.grid("both")
        ax.set_axisbelow(True)
        ax.set_xscale("log")
        ax.set_xlim(XMIN, XMAX)

        # Ticker settings
        ax.set_xticks([10**i for i in range(-5, 8, 1)])
        ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(['' if i%2==0 else f'$10^{{{i}}}$' for i  in range(-5, 8, 1)]))
        ax.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=np.arange(2, 10)*0.1, numticks=130))
        ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

        ax.legend(shadow = False, fontsize = 'medium', fancybox = True)
        ax.set_xlabel("Energy [eV]", fontsize='medium')
        ax.autoscale_view()
        
        # Setting the Y-axis label
        if   normalization_type == 0:
            ax.set_ylabel("Sensitivity", fontsize='medium')
        elif normalization_type == 1:
            ax.set_ylabel("Sensitivity per unit energy [eV$^{-1}$]", fontsize='medium')
        elif normalization_type == 2:
            ax.set_ylabel("Sensitivity per unit lethargy", fontsize='medium')
        else:
            raise ValueError(f'The type of normalization was set to {normalization_type}, while 1, 2, or 3 are allowed.')
        
        fig.savefig(f"{name}.{format}", format = format, dpi = 1000, bbox_inches='tight')

    @staticmethod
    def plot_uncertainty(covariances, name='uncertainty_profile', format = 'svg', annotations = None, width=3.54):
        """Plot the diagonal values of a covariance
        matrix.
        
        Parameters
        ----------
        covariances : List
            List of the Covariances instances
        name : str, optional
            Name of the plot to save. The default value is
            'uncertainty_profile'.
        format : str, optional 
            Format of the plot to save. The default value is
            'svg'.
        annotations : List
            List of values to add at the end of the legend if
            the plots have to be differentiated when the same
            uncertainties from different sources are used.
        width : str, optional
            The width of the plot in inches. The default value is 3.53.

        """      

 
        fig, ax = plt.subplots(figsize=(width, 2.655))
        s = 0
        use_ylog = False
        contain_cross = False

        # Annotation check
        a = 0
        if annotations != None:
            if len(annotations) != len(covariances):
                raise ValueError('The numbers of annotations and covariances are not the same.')
          

        for covariance in covariances:
            
            energies = covariance.group_structure
            reaction_1 = covariance.reaction_1
            reaction_2 = covariance.reaction_2

            zam_1 = covariance.zam_1
            zam_2 = covariance.zam_2
            charge_1 = round(zam_1 / 1e4)
            charge_2 = round(zam_2 / 1e4)
            mass_number_1 = str(round((zam_1 - charge_1*1e4)/10))
            mass_number_2 = str(round((zam_2 - charge_2*1e4)/10))
            if int(str(zam_1)[-1]): mass_number_1.append('m')
            if mass_number_1 == '0': mass_number_1 = '\mathrm{nat}'
            symbol_1 = ATOMIC_NUMBERS[charge_1]
            if int(str(zam_2)[-1]): mass_number_2.append('m')
            if mass_number_2 == '0': mass_number_1 = '\mathrm{nat}'
            symbol_2 = ATOMIC_NUMBERS[charge_2]

            reaction_1 = LATEX_REACTIONS[reaction_1]
            reaction_2 = LATEX_REACTIONS[reaction_2]
            
            label = f'$^{{{mass_number_1}}}${symbol_1}${reaction_1}$'
            if (zam_1 == zam_2) & (reaction_1 == reaction_2):
                profile = np.sqrt(np.diag(covariance.dataframe)) * 1e2
            else:
                profile = np.array([np.sqrt(i)*1e2 if i >= 0 else -np.sqrt(-i)*1e2 for i in np.diag(covariance.dataframe)])
                label += f'$-$$^{{{mass_number_2}}}${symbol_2}${reaction_2}$'
                contain_cross = True

            # Make log scale if the uncertainties too high
            if np.any(profile > 200):
                use_ylog = True

            # Add annotations if present
            if annotations != None:
                label += f'{annotations[a]}'
                a += 1 

            s += 1
            ax.stairs(profile, energies, label=label, alpha=0.7, linewidth = LINEWIDTH, linestyle=LINESTYLES[s%len(LINESTYLES)-1])

        ax.grid("both")
        ax.set_axisbelow(True)
        ax.set_xscale("log")
        ax.set_xlim(XMIN, XMAX)

        if use_ylog & (contain_cross == False):
            ax.set_yscale("log")

        # Ticker settings for the log scale
        ax.set_xticks([10**i for i in range(-5, 8, 1)])
        ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(['' if i%2==0 else f'$10^{{{i}}}$' for i  in range(-5, 8, 1)]))
        ax.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=np.arange(2, 10)*0.1, numticks=130))
        ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

        ax.legend(shadow = False, fontsize = 'medium', fancybox = True)
        ax.set_xlabel("Energy [eV]", fontsize='medium')
        ax.set_ylabel("Uncertainty [$\%$]", fontsize='medium')
        ax.autoscale_view()

        fig.savefig(f"{name}.{format}", format = format, dpi = 1000, bbox_inches='tight')

    @staticmethod
    def plot_cov(covariance,  name='cov_matrix', format = 'svg', tick_step = 1, type = 'matplotlib'):
        """Plot a covariance matrix from a Covariance
        instance.
        
        Parameters
        ----------
        covariance : Cobariance
            Covariances instance
        name : str
            Name of the plot
        format : str
            Format of the plot
        tick_step : int
            Number of ticks between numbers
        type : str
            Type of plot: {'matplotlib', 'seaborn'}

        """      

        fig, ax = plt.subplots(figsize=(3.54, 2.655))

        energies = covariance.group_structure
        group_number = len(energies) - 1
        reaction_1 = covariance.reaction_1
        reaction_2 = covariance.reaction_2
        covmatrix = np.rot90(covariance.dataframe, 2)

        zam_1 = covariance.zam_1
        zam_2 = covariance.zam_2
        charge_1 = round(zam_1 / 1e4)
        charge_2 = round(zam_2 / 1e4)
        mass_number_1 = str(round((zam_1 - charge_2*1e4)/10))
        mass_number_2 = str(round((zam_1 - charge_2*1e4)/10))
        if int(str(zam_1)[-1]): mass_number_1.append('m')
        if mass_number_1 == '0': mass_number_1 = '\mathrm{nat}'
        symbol_1 = ATOMIC_NUMBERS[charge_1]
        if int(str(zam_2)[-1]): mass_number_2.append('m')
        if mass_number_2 == '0': mass_number_2 = '\mathrm{nat}'
        symbol_2 = ATOMIC_NUMBERS[charge_2]

        reaction_1 = LATEX_REACTIONS[reaction_1]
        reaction_2 = LATEX_REACTIONS[reaction_2]

        label_1 = f'$^{{{mass_number_1}}}${symbol_1}${reaction_1}$'
        label_2 = f'$^{{{mass_number_2}}}${symbol_2}${reaction_2}$'

        if type == 'seaborn':
            ticklabels = [i     for i in range(1, group_number+1, 1)][::-1]
            xticks     = [i+0.5 for i in range(0, group_number, tick_step)]
            yticks     = [i-0.5 for i in range(tick_step, group_number+tick_step, tick_step)]

            cmap = sns.diverging_palette(240, 10, as_cmap=True)
            sns.heatmap(covmatrix, cmap=cmap, center=0, square=True, linewidths=.5,
                        xticklabels=ticklabels, yticklabels=ticklabels)

            ax.set_xticks(xticks)
            ax.set_yticks(yticks)
            ax.set_xlabel(f'{label_2}', fontsize='medium')
            ax.set_ylabel(f'{label_1}', fontsize='medium')

        elif type == 'matplotlib':
            xticks     = [i+1 for i in range(0, group_number, tick_step)]
            yticks     = [i+1 for i in range(int(tick_step/2), group_number, tick_step)]
            plot = ax.pcolormesh(np.arange(0.5, group_number+1, 1), np.arange(0.5, group_number+1, 1), covmatrix, cmap='RdBu_r', edgecolors='white', linewidths=0.1, rasterized = False, norm=mpl.colors.CenteredNorm())
            ax.set_xticks(xticks)
            ax.invert_xaxis()
            ax.set_yticks(yticks)
            ax.set_xlabel(f'{label_2}', fontsize='medium')
            ax.set_ylabel(f'{label_1}', fontsize='medium')
            ax.set(frame_on=False)
            ax.set_aspect(1) 

            cb = fig.colorbar(plot)
            cb.outline.set_visible(False)

            # Align tick labels at
            tick_labels = cb.ax.get_yticklabels()
            [label.set_horizontalalignment('right') for label in tick_labels]
            cb.ax.tick_params(axis='y', pad = 30)

        else:
            print(f'There is no type "{type}" implemented. Set "seaborn" or "matplotlib"')

        plt.subplots_adjust(bottom=0.2)
        fig.savefig(f"{name}.{format}", format = format, dpi=1000)

    @staticmethod
    def plot_cor(covariance, covariance_1 = None, covariance_2 = None, name='corr_matrix', format = 'svg'):
        """Plot a correlation matrix from a Covariance
        instance as an NJOY's COVR-like plot.
        
        Parameters
        ----------
        covariance : Covariance
            Covariance instance to be plotted
        covariance_1 : Covariance, optional
            Covariance instance for reaction_1 of covariance to be passed
            if a channel-channhel correlation is plotted.
        covariance_2 : Covariance, optional
            Covariance instance for reaction_2 of covariance to be passed
            if a channel-channhel correlation is plotted.
        name : str, optional
            Name of the plot. The default value is 'corr_matrix'.
        format : str, optional
            Format of the plot to save. The default value is 'svg'.

        """      

        axes = []
        fig = plt.figure(figsize=(3.54, 3.9648))
        gs  = plt.GridSpec(5,4, height_ratios=[0.2, 0.05, 0.9, 0.1, 0.05], width_ratios=[0.01, 0.85, 0.05, 0.2])
        axes.append(fig.add_subplot(gs[2,1])) # heatmap
        axes.append(fig.add_subplot(gs[2,3])) # uncertainty (right)
        axes.append(fig.add_subplot(gs[0,1])) # uncertainty (top)
        axes.append(fig.add_subplot(gs[4,1])) # colormap
        
        energies = covariance.group_structure
        group_number = len(energies) - 1
        zam_1 = covariance.zam_1
        zam_2 = covariance.zam_2
        reaction_1 = covariance.reaction_1
        reaction_2 = covariance.reaction_2
        covmatrix = covariance.dataframe
        
        if (reaction_1 == reaction_2) & (zam_1 == zam_2):
            diag_1 = np.sqrt(np.diag(covmatrix))
            diag_2 = diag_1
            cormatrix = np.rot90(cov_to_corr(covmatrix), 2)
        else:
            if (zam_1 != covariance_1.zam_1) | (zam_2 != covariance_1.zam_2) | (zam_1 != covariance_2.zam_1) | (zam_2 != covariance_2.zam_2):
                raise ValueError('The covariances are for incorrect nuclides. Check covariance_1 and covariance_2.')
            if (covariance_1 == None) | (covariance_2 == None):
                raise ValueError('Non-symmetric covariances are provided, provide the appropriate symmetric matrices.')
            elif (covariance_1.reaction_1 != reaction_1) | (covariance_1.reaction_2 != reaction_1):   
                raise ValueError('The first symmetric matrix is incorrect.')
            elif (covariance_2.reaction_1 != reaction_2) | (covariance_2.reaction_2 != reaction_2):
                raise ValueError('The second symmetric matrix is incorrect.')
            diag_1 = np.sqrt(np.diag(covariance_1.dataframe))
            diag_2 = np.sqrt(np.diag(covariance_2.dataframe))
            cormatrix = np.rot90(covmatrix/np.outer(diag_1,diag_2), 2)

        
        zam_1 = covariance.zam_1
        zam_2 = covariance.zam_2
        charge_1 = round(zam_1 / 1e4)
        charge_2 = round(zam_2 / 1e4)
        mass_number_1 = str(round((zam_1 - charge_2*1e4)/10))
        mass_number_2 = str(round((zam_1 - charge_2*1e4)/10))
        if int(str(zam_1)[-1]): mass_number_1.append('m')
        if mass_number_1 == '0': mass_number_1 = '\mathrm{nat}'
        symbol_1 = ATOMIC_NUMBERS[charge_1]
        if int(str(zam_2)[-1]): mass_number_2.append('m')
        if mass_number_2 == '0': mass_number_2 = '\mathrm{nat}'
        symbol_2 = ATOMIC_NUMBERS[charge_2]

        reaction_1 = LATEX_REACTIONS[reaction_1]
        reaction_2 = LATEX_REACTIONS[reaction_2]

        label_1 = f'$^{{{mass_number_1}}}${symbol_1}${reaction_1}$'
        label_2 = f'$^{{{mass_number_2}}}${symbol_2}${reaction_2}$'

        profile_1 = diag_1 * 1e2
        profile_2 = diag_2 * 1e2

        # Correlation plot
        xticks = [10**i for i in range(-5, 8, 1)]
        yticks = xticks
        tick_labels = ['' if i%2==0 else f'$10^{{{i}}}$' for i  in range(-5, 8, 1)]
        plot = axes[0].pcolormesh(energies[::-1], energies[::-1], cormatrix, cmap='RdBu_r', edgecolors='white', linewidths=0.1, rasterized = False, norm=mpl.colors.CenteredNorm())
        axes[0].set_xlim(XMIN, XMAX)
        axes[0].set_ylim(XMIN, XMAX)

        axes[0].set_xscale("log")
        axes[0].set_yscale("log")

        axes[0].set_xticks(xticks, labels=tick_labels)
        axes[0].set_yticks(yticks, labels=tick_labels)
    
        axes[0].invert_yaxis()
        
        # Uncertainty plot (right)
        axes[1].hist(x=energies[:-1], bins=energies, weights=profile_1, linestyle='solid', alpha=0.8, histtype='step', align='mid', orientation='horizontal', linewidth = LINEWIDTH/2)
        axes[1].set_ylim(XMIN, XMAX)

        axes[1].set_ylabel("Energy [eV]", fontsize='medium', rotation=-90, labelpad=10.0)
        axes[1].set_yscale("log")
        axes[1].set_yticks(yticks, labels=tick_labels)
        axes[1].invert_yaxis()


        # Uncertainty plot (top)
        axes[2].hist(x=energies[:-1], bins=energies, weights=profile_2, linestyle='solid', alpha=0.8, histtype='step', align='mid', orientation='vertical', linewidth = LINEWIDTH/2)
        axes[2].set_xlim(XMIN, XMAX)
        axes[2].set_xlabel("Energy [eV]", fontsize='medium', labelpad=2.0)
        axes[2].set_xscale("log") 
        axes[2].set_xticks(yticks, labels=tick_labels)

        axes[0].set_xlabel(f'{label_2}', fontsize='medium')
        axes[0].set_ylabel(f'{label_1}', fontsize='medium')
        axes[0].set_xticks(xticks)
        axes[0].set_yticks(yticks)
        axes[0].set(frame_on=False)
        axes[0].tick_params(length=2)
        axes[0].set_aspect(1)

        cb = mpl.colorbar.Colorbar(ax = axes[3], mappable = plot, orientation='horizontal', ticks=mpl.ticker.MaxNLocator(5), ticklocation = 'bottom')
        cb.outline.set_visible(False)

        axes[1].xaxis.set_major_locator(mpl.ticker.MaxNLocator(3))
        axes[1].tick_params(length=2, rotation=-90)
        axes[1].set_xlabel("Uncertainty [$\%$]", fontsize='medium')
        axes[2].yaxis.set_major_locator(mpl.ticker.MaxNLocator(3))
        axes[2].set_ylabel("Uncertainty [$\%$]", fontsize='medium')
        axes[2].tick_params(length=2) 

        fig.savefig(f"{name}.{format}", format=format, dpi = 1000)

# TODO: Add plot uncertainty + sensitivity plot
# TODO: Add model-corr plot