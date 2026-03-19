from __future__ import annotations

import pathlib
import numpy as np
from numpy.typing import NDArray


def get_text(file: str) -> str:
    """Extract the text from an ENDF-6 file and fix it

    Parameters
    ----------
    file : str
        The relative path to an ENDF-6 file to process, e.g.,
        '../NuclearData/BROND-3.1/n_5040_50-Sn-117.dat'
    
    Returns
    -------
    str
        The fixed text of an ENDF-6 file 

    """
    
    # Read the ENDF-6 file
    path = pathlib.Path(file)
    with open(path, 'r+') as file:
        text = file.read()

    # Remove errors from the text of an ENDF-6 file
    text = fix_endf(text)

    return text


def fix_endf(text: str) -> str:
    """Fix ENDF-6 files to allow pandas to parse them correctly.
    This does not interact with the original ENDF-6 files. They
    are fixed internally.

    Parameters
    ----------
    text : str
        The text from the ENDF-6 file
    
    Returns
    -------
    str
        The fixed an ENDF-6 file as string

    """

    # These two fix some exceptions with the BROND-3.1 ENDF files
    if '+ ' in text:
        text = text.replace('+ ', '+0')
    if '- ' in text:
        text = text.replace('- ', '-0')
    
    # This fixes some exceptions with the BROND-3.1 and CENDL-3.2 ENDF files
    if '?' in text:
        text = text.replace('?', ' ')
    
    return text


def get_zam(tape: "sandy.Endf6") -> int:  # type: ignore[name-defined]
    """Get ZAM (ZA * 10 + META) of a nuclide from tape
     
    Parameters
    ----------
    tape: sandy.Endf6
        An ENDF-6 tape before processing

    Returns
    -------
    int
        ZAM of the nuclide

    Notes
    -----
    This method is used to properly name files

    """
    # Get the MAT number
    mat: int = tape.mat[0]

    # Get data from the tape
    info: dict = tape.read_section(mat, 1, 451)
    meta: int = info["LISO"]
    za: int = int(info["ZA"])
    zam: int = za * 10 + meta

    return zam


def process_file(
    tape: "sandy.Endf6",  # type: ignore[name-defined]
    library: str,
    group_structure: list[float],
) -> "sandy.endf.errorr":  # type: ignore[name-defined]
    """Process an ENDF-6 tape through NJOY/ERRORR for a given library and group structure.
     
    Parameters
    ----------
    tape: sandy.Endf6
        An ENDF-6 tape to process
    library : str
        The argument defines how the method interacts with peculiarities of
        different libraries. The special exceptions are provided for
        'ENDF/B-VIII.0', 'JEFF-4T2.2', 'JEFF-4T3', and 'BROND-3.1'. There are
        no exceptions for other state-of-the-art libraries. 
    group_structure : list of float
        Group structure to process an ENDF-6 file in eV

    Returns
    -------
    sandy.endf.errorr
       Processed from an ENDF-6 file tape

    Notes
    -----
    Not provided

    """   

    zam: int = get_zam(tape)
    # These conditions are specified otherwise NJOY does not allow processing the following files
    if library == 'ENDF/B-VIII.0':
        if zam in [922380, 922350]:                               
            errorr = tape.get_errorr(temperature=293.6, err=0.001, mubar=False,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)
        elif zam in [641570, 501240, 501170]:
            errorr = tape.get_errorr(temperature=293.6, err=0.001, nubar=False, chi=False,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)                                
        else:
            errorr = tape.get_errorr(temperature=293.6, err=0.001, 
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)
    elif library == 'JEFF-4T2.2':
        if zam in [922390, 922400]:
            errorr = tape.get_errorr(temperature=293.6, err=0.001, nubar=False,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)
        else:
            errorr = tape.get_errorr(temperature=293.6, err=0.001, 
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)
    elif library == 'JEFF-4T3':
        if zam == 922380:
            errorr = tape.get_errorr(temperature=293.6, err=0.001, mubar=False, chi=False,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)
        elif zam in [922390, 922400]:
            errorr = tape.get_errorr(temperature=293.6, err=0.001, nubar=False,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)
        else:
            errorr = tape.get_errorr(temperature=293.6, err=0.001,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False) 
    elif library == 'BROND-3.1':
        if zam in [80180, 150310, 250550, 420940, 420950, 420960, 420970, 420980, 501170, 501240, 641570, 741800, 741830]:
            errorr = tape.get_errorr(temperature=293.6, err=0.001, nubar=False, chi=False,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)                
        elif zam in [922380, 932390]:
            errorr = tape.get_errorr(temperature=293.6, err=0.001, chi=False,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)
        elif zam == 942390:
            errorr = tape.get_errorr(temperature=293.6, err=0.001, mubar=False, chi=False,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False)
        else:
            errorr = tape.get_errorr(temperature=293.6, err=0.001,
                                     minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                     errorr_kws=dict(ek=group_structure),verbose=False) 
    elif library in ['ENDF/B-VII.1','JENDL-5','JEFF-3.3','CENDL-3.2']:
        errorr = tape.get_errorr(temperature=293.6, err=0.001, 
                                 minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                 errorr_kws=dict(ek=group_structure),verbose=False)
    else:
        errorr = tape.get_errorr(temperature=293.6, err=0.001, 
                                 minimal_processing=True, groupr_kws=dict(ek=group_structure),
                                 errorr_kws=dict(ek=group_structure),verbose=False)

    return errorr


def cov_to_corr(array_of_covariances: NDArray[np.floating]) -> NDArray[np.floating]:
    """Generate a correlation matrix from a symmetric covariance matrix.

    Parameters
    ----------
    array_of_covariances : numpy.ndarray
        Covariance matrix

    Return
    ------
    numpy.ndarray
        Correlation matrix derived from the covariance matrix

    """
    diag: NDArray[np.floating] = np.sqrt(np.diag(array_of_covariances))
    outer: NDArray[np.floating] = np.outer(diag, diag)
    correlations: NDArray[np.floating] = np.divide(
        array_of_covariances, outer,
        out=np.zeros_like(array_of_covariances),
        where=outer != 0,
    )
    return correlations


def corr_to_cov(
    array_of_correlations: NDArray[np.floating],
    variances: NDArray[np.floating],
) -> NDArray[np.floating]:
    """Generate a covariance matrix from a symmetric correlation matrix.

    Parameters
    ----------
    array_of_correlations : numpy.ndarray
        Correlation matrix
    variances : numpy.ndarray
        Variances for the corresponding correlations

    Return
    ------
    numpy.ndarray
        Covariance matrix derived from the correlation matrix

    """
    diag: NDArray[np.floating] = np.sqrt(variances)
    outer_diag: NDArray[np.floating] = np.outer(diag, diag)
    covariances: NDArray[np.floating] = array_of_correlations * outer_diag
    return covariances


def fix_corrs(cov: "Covariance", eps: float = 1e-5) -> NDArray[np.floating]:  # type: ignore[name-defined]
    """Make sure that a symmetric covariance matrix does
    not have non-mathematical correlations.

    Parameters
    ----------
    cov : Covariance
        A symmetric Covariance instance (reaction_1 must equal reaction_2)
    eps : float, optional
        Tolerance for clipping correlations outside [-1, 1]. Default is 1e-5.

    Returns
    -------
    numpy.ndarray
        Fixed symmetric covariance matrix as a numpy array

    Raises
    ------
    ValueError
        If a non-symmetric (cross-reaction) covariance is provided.
    
    """      

    if cov.reaction_1 == cov.reaction_2:
        array_of_covs: NDArray[np.floating] = cov.dataframe.to_numpy()

        corr: NDArray[np.floating] = cov_to_corr(array_of_covs)

        corr = np.where(corr > 1+eps, 1, corr)                 
        corr = np.where(corr < -1-eps, -1, corr)

        cov = corr_to_cov(corr, np.diag(cov.dataframe.to_numpy()))

        return cov
    else:
        raise ValueError(f'Non-symmetric matrix is provided.')


def get_negative(list: list[float]) -> int | None:
    """Return the index of the first negative value in a list.

    Parameters
    ----------
    list : list of float
        Input list to search for a negative value

    Returns
    -------
    int or None
        Index of the first negative element, or None if none found

    """
    for i in range(len(list)):
        if list[i] < 0:
            return i


def get_complex(list: list[complex]) -> int | None:
    """Return the index of the first complex value in a list.

    Parameters
    ----------
    list : list of complex
        Input list to search for a complex value

    Returns
    -------
    int or None
        Index of the first complex element, or None if none found

    """
    for i in range(len(list)):
        if np.iscomplex(list[i]):
            return i