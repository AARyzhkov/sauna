---
icon: person-chalkboard
---

# Defining Covariances

Defining covariances is conducted in a similar to sensitivities manner: the covariances are stored in a $$\texttt{Covariances}$$ instance as a list of $$\texttt{Covariance}$$ instances. The conventional way to get covariances is to import them from different sources.  Currently, five formats are supported: ENDF-6, Excel, COMMARA, ABBN, and COVERX.

## ENDF-6

Firstly, the ENDF-6 format stores the covariance data in a pointwise format requiring processing to get groupwise covariances. NJOY is used for this purpose, the interaction between SAUNA and NJOY relies on the SANDY package interface, which is done in the following way.

```python
covariances = sauna.Covariances()
covariances.library = 'ENDF-B-VII.1'
covariances.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
covariances.from_endfs('./path/to/folder/with/endf/files', extension='.dat', parallel=True)
```

## Excel

The library is given explicitly to properly handle unique features of different libraries that have to be accounted for processing. To avoid processing each time, it is suggested to export the covariances into the Excel format for using them later. It can be done via the following line.

```python
covariances.to_excels('./path/where/to/save/covariances/')
```

To get the data from the Excel files the $$\texttt{from_excels()}$$ method is used. For instance, to get the provided with SAUNA data, the following line can be used with the path to the folder with the Excel files.

```python
covariances.from_excels('../../covariances/ENDF-B-VII.1-56/')
```

## COMMARA

The COMMARA data are stored in one file, consequently, it only requires path to the file.

```python
covariances.from_commara('../../covariances/COMMARA.cov')
```

## ABBN

The ABBN format stores the data in a different file for each nuclide, therefore the path to a folder has to be provided.

```
covariances.from_abbns('./path/to/folder/with/ABBN/covariances')
```

## AMPX

Finally, the latest supported format is COVERX of AMPX. COVERX is not supported directly: it relies on the AMPX module, TOC, which the path has to be provided to.

```python
ampx_path = 'path/to/SCALE-6.2.4-Source/build/install/bin/AmpxCOVConverter'
coverx_path = 'scale.rev08.56groupcov7.1'
covariances.from_coverx(ampx_path, coverx_path)
```

{% hint style="warning" %}
The process takes a substantial amount of time to import the data. It can be recommended to export the covariances via `to_excels()` for a later usage.
{% endhint %}
