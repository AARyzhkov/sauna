---
icon: person-chalkboard
---

# Defining Sensitivities

## Importing sensitivities

The main class for working with sensitivities here is $$\texttt{Sensitivities}$$. The class is intended to handle $$\texttt{Sensitivity}$$ instances. A $$\texttt{Sensitivities}$$ instance can be populated manually by creating $$\texttt{Sensitivity}$$ instances and defining their parameters though the approach is rather inefficient. Consequently, the most conventional way to create a proper $$\texttt{Sensitivities}$$ instance is to import the data from neutron transport simulation results such as Serpent or SCALE.

```python
serpent_sensitivities = sauna.Sensitivities()
serpent_sensitivities.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
serpent_sensitivities.from_serpent('../../models/MET1000_FC_56_sens0.m')
```

```python
scale_sensitivities = sauna.Sensitivities()
scale_sensitivities.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
scale_sensitivities.from_scale('../../models/MET1000_FC_NFP.2025.02.28T10.09.44.sdf', 'B', 'Eigenvalue')
```

{% hint style="warning" %}
The `from_scale()` method does not mean it is allowed only for SCALE results: it works with `.sdf` files, which can be produced by other tools such as MCNP. Although SCALE is a code system, and the sensitivities are computed via the KENO-VI transport code in the frame of the TSUNAMI-3D module of SCALE, it was decided to set for generalizing that for the people who are unfamiliar with the SCALE structure but are aware of the code system.
{% endhint %}

One has to provide the group structure to make sure the other functions of SAUNA work properly. It can be accessed from the `sauna.GROUP_STRUCTURES`  dictionary or provided as a list of energy boundaries from the lowest to the largest G+1 values. The predefined group structures are provided below:

* [ABBN-28/299](https://inis.iaea.org/collection/NCLCollectionStore/_Public/28/075/28075558.pdf)
* [ECCO-33](http://serpent.vtt.fi/mediawiki/index.php/ECCO_33-group_structure)
* [WPEC-SG46](https://www.epj-conferences.org/10.1051/epjconf/202328414012)
* [SCALE-44](https://www-nds.iaea.org/publications/indc/indc-czr-0001.pdf)
* [SCALE-56/252](https://info.ornl.gov/sites/publications/Files/Pub213297.pdf)
* [UFLIB-70](https://info.ornl.gov/sites/publications/Files/Pub213297.pdf)

## Accessing sensitivities

In some cases, it is necessary to get a $$\texttt{Sensitivity}$$ instance with certain parameters. For this reason, some methods are provided. Three main parameters are of interest for asking sensitivities: the functional name, ZAM, and MT. The available functional can be known by accessing the `functional` attribute of a $$\texttt{Sensitivities}$$ instance, and the same goes for ZAMs.

```
serpent_sensitivities.functionals
serpent_sensitivities.zams
```

ZAM values are standard and defined through a >=5-digit number. The values are used to define a nuclide. The first digits are taken by a number of the atomic number. Next, three digits is taken by the mass number though they are zero if the nuclide represents a natural composition such as 60000 for $$^{\text{nat}}\text{C}$$. The last digit is taken by the number of the metastable state, which is usually zero though in some cases it is not, for instance, 952421 for $$^{242\text{m}}\text{Am}$$. The MT numbers are governed by [the ENDF-6 format](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf). In addition, all the MT numbers can also be found on [the NEA website](https://www.oecd-nea.org/dbdata/data/manual-endf/endf102_MT.pdf). Only a part of the numbers is relevant and provided below:

```
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
- MT1002 - thermal elastic scattering laws (S(α,β))
- MT1004 - thermal inelastic scattering laws (S(α,β))
- MT1018 - fission spectrum (χ) - the value of 1018 is taken from SCALE
- MT1455 - delayed fission spectrum (χ-delayed)
- MT1456 - prompt fission spectrum (χ-prompt)
```

{% hint style="warning" %}
MT=1004, 1018, 1455, and 1456 are not a part of the standard MT number list of the ENDF-6 format and are taken for simplicity. The standard MT numbers are limited by 999.
{% endhint %}

Here are the examples of accessing $$\texttt{Sensitivity}$$ instances.

```python
serpent_sensitivities.get_by_functional('Eigenvalue') # returns a list of Sensitivity instances
serpent_sensitivities.get_by_zam(922380) # returns a list of Sensitivity instances
scale_sensitivities.get_by_reaction(4) # returns a list of Sensitivity instances
scale_sensitivities.get_by_parameters('Eigenvalue', 922380, 4) # returns Sensitivity instance
```

## Exporting sensitivities

For an analysis, the data can be imported whether to a dataframe leveraging the pandas package capabilities or imported to an Excel file using one of the following two lines, respectively.

```python
sensitivitiy_df = scale_sensitivities.to_dataframe('Eigenvalue', sort=True)
scale_sensitivities.to_excel(name='scale_sensitivities.xlsx', sort=True)
```
