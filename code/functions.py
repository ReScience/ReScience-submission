#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------
# Copyright (c) 2016, Rafael Neto Henriques
#
# Contributors : Rafael Neto Henriques (rafaelnh21@gmail.com)
# -------------------------------------------------------------------------
# References:
#
# Hoy, A.R., Koay, C.G., Kecskemeti, S.R., Alexander, A.L., 2014
# Optimization of a free water elimination two-compartmental model
# for diffusion tensor imaging. NeuroImage 103, 323-333.
# doi: 10.1016/j.neuroimage.2014.09.053
# -------------------------------------------------------------------------


import numpy as np

from dipy.reconst.dti import (design_matrix, decompose_tensor,
                              from_lower_triangular, lower_triangular)
from dipy.reconst.vec_val_sum import vec_val_vect
from dipy.core.ndindex import ndindex

import scipy.optimize as opt

# -------------------------------------------------------------------------
# Weigthed linear least squares fit procedure
# -------------------------------------------------------------------------


def wls_iter(design_matrix, sig, S0, Diso=3e-3, mdreg=1.5e-3,
             min_signal=1.0e-6, piterations=3):
    """ Applies weighted linear least squares fit of the water free elimination
    model to single voxel signals.

    Parameters
    ----------
    design_matrix : array (g, 7)
        Design matrix holding the covariants used to solve for the regression
        coefficients.
    sig : array (g, )
        Diffusion-weighted signal for a single voxel data.
    S0 : float
        Non diffusion weighted signal (i.e. signal for b-value=0).
    Diso : float, optional
        Value of the free water isotropic diffusion. Default is set to 3e-3
        $mm^{2}.s^{-1}$. Please ajust this value if you are assuming different
        units of diffusion.
    mdreg : float, optimal
        Tissue compartment mean diffusivity regularization threshold.
        If tissue's mean diffusivity is almost near the free water diffusion
        value, the diffusion signal is assumed to be only free water diffusion
        (i.e. volume fraction will be set to 1 and tissue's diffusion
        parameters are set to zero). Default md_reg was set to
        1.5e-3 $mm^{2}.s^{-1}$ according to [1]_.
    min_signal : float
        The minimum signal value. Needs to be a strictly positive
        number. Default: minimal signal in the data provided to `fit`.
    piterations : inter, optional
        Number of iterations used to refine the precision of f. Default is set
        to 3 corresponding to a precision of 0.01.

    Returns
    -------
    All parameters estimated from the free water tensor model.
    Parameters are ordered as follows:
        1) Three diffusion tensor's eigenvalues
        2) Three lines of the eigenvector matrix each containing the
           first, second and third coordinates of the eigenvector
        3) The volume fraction of the free water compartment

    References
    ----------
    .. [1] Henriques, R.N., Rokem, A., Garyfallidis, E., St-Jean, S., Peterson,
           E.T., Correia, M.M., 2017. Re: Optimization of a free water
           elimination two-compartmental model for diffusion tensor imaging.
           ReScience
    """
    W = design_matrix

    # Define weights
    S2 = np.diag(sig**2)

    # Defining matrix to solve fwDTI wls solution
    WTS2 = np.dot(W.T, S2)
    inv_WT_S2_W = np.linalg.pinv(np.dot(WTS2, W))
    invWTS2W_WTS2 = np.dot(inv_WT_S2_W, WTS2)

    # Process voxel if it has significant signal from tissue
    if np.mean(sig) > min_signal and S0 > min_signal:
        # General free-water signal contribution
        fwsig = np.exp(np.dot(design_matrix,
                              np.array([Diso, 0, Diso, 0, 0, Diso, 0])))

        df = 1  # initialize precision
        flow = 0  # lower f evaluated
        fhig = 1  # higher f evaluated
        ns = 9  # initial number of samples per iteration
        for p in range(piterations):
            df = df * 0.1
            fs = np.linspace(flow+df, fhig-df, num=ns)  # sampling f
            SFW = np.array([fwsig, ]*ns)  # repeat contributions for all values
            FS, SI = np.meshgrid(fs, sig)
            SA = SI - FS*S0*SFW.T
            # SA < 0 means that the signal components from the free water
            # component is larger than the total fiber. This cases are present
            # for inapropriate large volume fractions (given the current S0
            # value estimated). To overcome this issue negative SA are replaced
            # by data's min positive signal.
            SA[SA <= 0] = min_signal
            y = np.log(SA / (1-FS))
            all_new_params = np.dot(invWTS2W_WTS2, y)
            # Select params for lower F2
            SIpred = (1-FS)*np.exp(np.dot(W, all_new_params)) + FS*S0*SFW.T
            F2 = np.sum(np.square(SI - SIpred), axis=0)
            Mind = np.argmin(F2)
            params = all_new_params[:, Mind]
            f = fs[Mind]  # Updated f
            flow = f - df  # refining precision
            fhig = f + df
            ns = 19

        if mdreg is None:
            evals, evecs = decompose_tensor(from_lower_triangular(params))
            fw_params = np.concatenate((evals, evecs[0], evecs[1], evecs[2],
                                        np.array([f])), axis=0)
        else:
            # MD regularization - if tissue's md is larger than mdreg,
            # the voxel will be classified as containing only free water
            md = (params[0] + params[2] + params[5]) / 3
            if md > mdreg:
                fw_params = np.zeros(13)
                fw_params[12] = 1.0
            else:
                evals, evecs = decompose_tensor(from_lower_triangular(params))
                fw_params = np.concatenate((evals, evecs[0], evecs[1],
                                            evecs[2], np.array([f])), axis=0)
    else:
        fw_params = np.zeros(13)

    return fw_params


def wls_fit_tensor(gtab, data, Diso=3e-3, mask=None, min_signal=1.0e-6,
                   piterations=3, mdreg=1.5e-3):
    r""" Computes weighted least squares (WLS) fit to calculate self-diffusion
    tensor using a linear regression model [1]_.

    Parameters
    ----------
    gtab : a GradientTable class instance
        The gradient table containing diffusion acquisition parameters.
    data : ndarray ([X, Y, Z, ...], g)
        Data or response variables holding the data. Note that the last
        dimension should contain the data. It makes no copies of data.
    Diso : float, optional
        Value of the free water isotropic diffusion. Default is set to 3e-3
        $mm^{2}.s^{-1}$. Please ajust this value if you are assuming different
        units of diffusion.
    mask : array, optional
        A boolean array used to mark the coordinates in the data that should
        be analyzed that has the shape data.shape[:-1]
    min_signal : float
        The minimum signal value. Needs to be a strictly positive
        number. Default: 1.0e-6.
    piterations : inter, optional
        Number of iterations used to refine the precision of f. Default is set
        to 3 corresponding to a precision of 0.01.
    mdreg : float, optimal
        Tissue compartment mean diffusivity regularization threshold.
        If tissue's mean diffusivity is almost near the free water diffusion
        value, the diffusion signal is assumed to be only free water diffusion
        (i.e. volume fraction will be set to 1 and tissue's diffusion
        parameters are set to zero). Default md_reg was set to
        1.5e-3 $mm^{2}.s^{-1}$ according to [1]_.

    Returns
    -------
    fw_params : ndarray (x, y, z, 13)
        Matrix containing in the last dimension the free water model parameters
        in the following order:
            1) Three diffusion tensor's eigenvalues
            2) Three lines of the eigenvector matrix each containing the
               first, second and third coordinates of the eigenvector
            3) The volume fraction of the free water compartment.

    References
    ----------
    .. [1] Henriques, R.N., Rokem, A., Garyfallidis, E., St-Jean, S., Peterson,
           E.T., Correia, M.M., 2017. Re: Optimization of a free water
           elimination two-compartmental model for diffusion tensor imaging.
           ReScience
    """
    fw_params = np.zeros(data.shape[:-1] + (13,))
    W = design_matrix(gtab)

    # Prepare mask
    if mask is None:
        mask = np.ones(data.shape[:-1], dtype=bool)
    else:
        if mask.shape != data.shape[:-1]:
            raise ValueError("Mask is not the same shape as data.")
        mask = np.array(mask, dtype=bool, copy=False)

    # Prepare S0
    S0 = np.mean(data[..., gtab.b0s_mask], axis=-1)

    index = ndindex(mask.shape)
    for v in index:
        if mask[v]:
            params = wls_iter(W, data[v], S0[v], min_signal=min_signal,
                              Diso=3e-3, piterations=piterations, mdreg=mdreg)
            fw_params[v] = params

    return fw_params

# -------------------------------------------------------------------------
# non-linear least squares fit procedure
# -------------------------------------------------------------------------


def _nls_err_func(tensor_elements, design_matrix, data, Diso=3e-3,
                  cholesky=False, f_transform=False):
    """ Error function for the non-linear least-squares fit of the tensor water
    elimination model.

    Parameters
    ----------
    tensor_elements : array (8, )
        The six independent elements of the diffusion tensor followed by
        -log(S0) and the volume fraction f of the water elimination
        compartment. Note that if cholesky is set to true, tensor elements are
        assumed to be written as Cholesky's decomposition elements. If
        f_transform is true, volume fraction f has to be converted to
        ft = arcsin(2*f - 1) + pi/2
    design_matrix : array
        The design matrix
    data : array
        The voxel signal in all gradient directions
    Diso : float, optional
        Value of the free water isotropic diffusion. Default is set to 3e-3
        $mm^{2}.s^{-1}$. Please ajust this value if you are assuming different
        units of diffusion.
    cholesky : bool, optional
        If true, the diffusion tensor elements were decomposed using cholesky
        decomposition. See fwdti.nls_fit_tensor
        Default: False
    f_transform : bool, optional
        If true, the water volume fraction was converted to
        ft = arcsin(2*f - 1) + pi/2, insuring f estimates between 0 and 1.
        See fwdti.nls_fit_tensor
        Default: True
    """
    tensor = np.copy(tensor_elements)
    if cholesky:
        tensor[:6] = cholesky_to_lower_triangular(tensor[:6])

    if f_transform:
        f = 0.5 * (1 + np.sin(tensor[7] - np.pi/2))
    else:
        f = tensor[7]

    # This is the predicted signal given the params:
    y = (1-f) * np.exp(np.dot(design_matrix, tensor[:7])) + \
        f * np.exp(np.dot(design_matrix,
                          np.array([Diso, 0, Diso, 0, 0, Diso, tensor[6]])))

    # Compute the residuals
    return data - y


def _nls_jacobian_func(tensor_elements, design_matrix, data, Diso=3e-3,
                       cholesky=False, f_transform=False):
    """The Jacobian is the first derivative of the least squares error
    function.

    Parameters
    ----------
    tensor_elements : array (8, )
        The six independent elements of the diffusion tensor followed by
        -log(S0) and the volume fraction f of the water elimination
        compartment. Note that if f_transform is true, volume fraction f is
        converted to ft = arcsin(2*f - 1) + pi/2
    design_matrix : array
        The design matrix
    Diso : float, optional
        Value of the free water isotropic diffusion. Default is set to 3e-3
        $mm^{2}.s^{-1}$. Please ajust this value if you are assuming different
        units of diffusion.
    f_transform : bool, optional
        If true, the water volume fraction was converted to
        ft = arcsin(2*f - 1) + pi/2, insuring f estimates between 0 and 1.
        See fwdti.nls_fit_tensor
        Default: True
    """
    tensor = np.copy(tensor_elements)
    if f_transform:
        f = 0.5 * (1 + np.sin(tensor[7] - np.pi/2))
    else:
        f = tensor[7]

    t = np.exp(np.dot(design_matrix, tensor[:7]))
    s = np.exp(np.dot(design_matrix,
                      np.array([Diso, 0, Diso, 0, 0, Diso, tensor[6]])))
    T = (f-1.0) * t[:, None] * design_matrix
    S = np.zeros(design_matrix.shape)
    S[:, 6] = f * s

    if f_transform:
        df = (t-s) * (0.5*np.cos(tensor[7]-np.pi/2))
    else:
        df = (t-s)
    return np.concatenate((T - S, df[:, None]), axis=1)


def nls_iter(design_matrix, sig, S0, Diso=3e-3, mdreg=2.7e-3,
             min_signal=1.0e-6, cholesky=False, f_transform=True,
             jac=True):
    """ Applies non linear least squares fit of the water free elimination
    model to single voxel signals.

    Parameters
    ----------
    design_matrix : array (g, 7)
        Design matrix holding the covariants used to solve for the regression
        coefficients.
    sig : array (g, )
        Diffusion-weighted signal for a single voxel data.
    S0 : float
        Non diffusion weighted signal (i.e. signal for b-value=0).
    Diso : float, optional
        Value of the free water isotropic diffusion. Default is set to 3e-3
        $mm^{2}.s^{-1}$. Please ajust this value if you are assuming different
        units of diffusion.
    min_signal : float
        The minimum signal value. Needs to be a strictly positive
        number.
    cholesky : bool, optional
        If true it uses cholesky decomposition to insure that diffusion tensor
        is positive define.
        Default: False
    f_transform : bool, optional
        If true, the water volume fractions is converted during the convergence
        procedure to ft = arcsin(2*f - 1) + pi/2, insuring f estimates between
        0 and 1.
        Default: True
    jac : bool
        Use the Jacobian? Default: False

    Returns
    -------
    All parameters estimated from the free water tensor model.
    Parameters are ordered as follows:
        1) Three diffusion tensor's eigenvalues
        2) Three lines of the eigenvector matrix each containing the
           first, second and third coordinates of the eigenvector
        3) The volume fraction of the free water compartment.
    """
    # Initial guess
    params = wls_iter(design_matrix, sig, S0, min_signal=min_signal, Diso=Diso)

    # Process voxel if it has significant signal from tissue
    if np.mean(sig) > min_signal and S0 > min_signal:
        # converting evals and evecs to diffusion tensor elements
        evals = params[:3]
        evecs = params[3:12].reshape((3, 3))
        dt = lower_triangular(vec_val_vect(evecs, evals))

        # Cholesky decomposition if requested
        if cholesky:
            dt = lower_triangular_to_cholesky(dt)

        # f transformation if requested
        if f_transform:
            f = np.arcsin(2*params[12] - 1) + np.pi/2
        else:
            f = params[12]

        # Use the Levenberg-Marquardt algorithm wrapped in opt.leastsq
        start_params = np.concatenate((dt, [-np.log(S0), f]), axis=0)
        if jac:
            this_tensor, status = opt.leastsq(_nls_err_func, start_params[:8],
                                              args=(design_matrix, sig, Diso,
                                                    cholesky, f_transform),
                                              Dfun=_nls_jacobian_func)
        else:
            this_tensor, status = opt.leastsq(_nls_err_func, start_params[:8],
                                              args=(design_matrix, sig, Diso,
                                                    cholesky, f_transform))

        # Invert the cholesky decomposition if this was requested
        if cholesky:
            this_tensor[:6] = cholesky_to_lower_triangular(this_tensor[:6])

        # Invert f transformation if this was requested
        if f_transform:
            this_tensor[7] = 0.5 * (1 + np.sin(this_tensor[7] - np.pi/2))

        # The parameters are the evals and the evecs:
        evals, evecs = decompose_tensor(from_lower_triangular(this_tensor[:6]))
        params = np.concatenate((evals, evecs[0], evecs[1], evecs[2],
                                 np.array([this_tensor[7]])), axis=0)
    return params


def nls_fit_tensor(gtab, data, mask=None, Diso=3e-3,
                   min_signal=1.0e-6, f_transform=True, cholesky=False,
                   jac=True):
    """
    Fit the water elimination tensor model using the non-linear least-squares.

    Parameters
    ----------
    gtab : a GradientTable class instance
        The gradient table containing diffusion acquisition parameters.
    data : ndarray ([X, Y, Z, ...], g)
        Data or response variables holding the data. Note that the last
        dimension should contain the data. It makes no copies of data.
    mask : array, optional
        A boolean array used to mark the coordinates in the data that should
        be analyzed that has the shape data.shape[:-1]
    Diso : float, optional
        Value of the free water isotropic diffusion. Default is set to 3e-3
        $mm^{2}.s^{-1}$. Please ajust this value if you are assuming different
        units of diffusion.
    min_signal : float
        The minimum signal value. Needs to be a strictly positive
        number. Default: 1.0e-6.
    f_transform : bool, optional
        If true, the water volume fractions is converted during the convergence
        procedure to ft = arcsin(2*f - 1) + pi/2, insuring f estimates between
        0 and 1.
        Default: True
    cholesky : bool, optional
        If true it uses cholesky decomposition to insure that diffusion tensor
        is positive define.
        Default: False
    jac : bool
        Use the Jacobian? Default: False

    Returns
    -------
    fw_params : ndarray (x, y, z, 13)
        Matrix containing in the last dimension the free water model parameters
        in the following order:
            1) Three diffusion tensor's eigenvalues
            2) Three lines of the eigenvector matrix each containing the
               first, second and third coordinates of the eigenvector
            3) The volume fraction of the free water compartment

    References
    ----------
    .. [1] Henriques, R.N., Rokem, A., Garyfallidis, E., St-Jean, S., Peterson,
           E.T., Correia, M.M., 2017. Re: Optimization of a free water
           elimination two-compartmental model for diffusion tensor imaging.
           ReScience
    """
    # Analyse compatible input cases
    if jac is True and cholesky is True:
        raise ValueError("Cholesky decomposition is not compatible with jac.")

    fw_params = np.zeros(data.shape[:-1] + (13,))
    W = design_matrix(gtab)

    # Prepare mask
    if mask is None:
        mask = np.ones(data.shape[:-1], dtype=bool)
    else:
        if mask.shape != data.shape[:-1]:
            raise ValueError("Mask is not the same shape as data.")
        mask = np.array(mask, dtype=bool, copy=False)

    # Prepare S0
    S0 = np.mean(data[..., gtab.b0s_mask], axis=-1)

    # Loop data fitting through all voxels
    index = ndindex(mask.shape)
    for v in index:
        if mask[v]:
            params = nls_iter(W, data[v], S0[v], Diso=Diso,
                              min_signal=min_signal,
                              f_transform=f_transform,
                              cholesky=cholesky, jac=jac)
            fw_params[v] = params

    return fw_params


def lower_triangular_to_cholesky(tensor_elements):
    """ Perfoms Cholesky decompostion of the diffusion tensor

    Parameters
    ----------
    tensor_elements : array (6,)
        Array containing the six elements of diffusion tensor's lower
        triangular.
    Returns
    -------
    cholesky_elements : array (6,)
        Array containing the six Cholesky's decomposition elements
        (R0, R1, R2, R3, R4, R5) [1]_.
    References
    ----------
    .. [1] Koay, C.G., Carew, J.D., Alexander, A.L., Basser, P.J.,
           Meyerand, M.E., 2006. Investigation of anomalous estimates of
           tensor-derived quantities in diffusion tensor imaging. Magnetic
           Resonance in Medicine, 55(4), 930-936. doi:10.1002/mrm.20832
    """
    R0 = np.sqrt(tensor_elements[0])
    R3 = tensor_elements[1] / R0
    R1 = np.sqrt(tensor_elements[2] - R3**2)
    R5 = tensor_elements[3] / R0
    R4 = (tensor_elements[4] - R3*R5) / R1
    R2 = np.sqrt(tensor_elements[5] - R4**2 - R5**2)

    return np.array([R0, R1, R2, R3, R4, R5])


def cholesky_to_lower_triangular(R):
    """ Convert Cholesky decompostion elements to the diffusion tensor elements

    Parameters
    ----------
    R : array (6,)
        Array containing the six Cholesky's decomposition elements
        (R0, R1, R2, R3, R4, R5) [1]_.

    Returns
    -------
    tensor_elements : array (6,)
        Array containing the six elements of diffusion tensor's lower
        triangular.

    References
    ----------
    .. [1] Koay, C.G., Carew, J.D., Alexander, A.L., Basser, P.J.,
           Meyerand, M.E., 2006. Investigation of anomalous estimates of
           tensor-derived quantities in diffusion tensor imaging. Magnetic
           Resonance in Medicine, 55(4), 930-936. doi:10.1002/mrm.20832
    """
    Dxx = R[0]**2
    Dxy = R[0]*R[3]
    Dyy = R[1]**2 + R[3]**2
    Dxz = R[0]*R[5]
    Dyz = R[1]*R[4] + R[3]*R[5]
    Dzz = R[2]**2 + R[4]**2 + R[5]**2
    return np.array([Dxx, Dxy, Dyy, Dxz, Dyz, Dzz])


# -------------------------------------------------------------------------
# Supplementary function
# -------------------------------------------------------------------------


def nls_iter_bounds(design_matrix, sig, S0, Diso=3e-3,
                    min_signal=1.0e-6, bounds=None, jac=True):
    """ Applies non-linear least-squares fit with constraints of the water free
    elimination model to single voxel signals.

    Parameters
    ----------
    design_matrix : array (g, 7)
        Design matrix holding the covariants used to solve for the regression
        coefficients.
    sig : array (g, )
        Diffusion-weighted signal for a single voxel data.
    S0 : float
        Non diffusion weighted signal (i.e. signal for b-value=0).
    Diso : float, optional
        Value of the free water isotropic diffusion. Default is set to 3e-3
        $mm^{2}.s^{-1}$. Please ajust this value if you are assuming different
        units of diffusion.
    min_signal : float
        The minimum signal value. Needs to be a strictly positive
        number.
    bounds : 2-tuple of arrays with 14 elements, optional
        Lower and upper bounds on fwdti model variables and the log of
        non-diffusion signal S0. Use np.inf with an appropriate sign to
        disable bounds on all or some variables. When bounds is set to None
        the following default variable bounds is used:
            ([0., -Diso, 0., -Diso, -Diso, 0., 0., np.exp(-10.)],
             [Diso, Diso, Diso, Diso, Diso, Diso, 1., np.exp(10.)])
    jac : bool
        Use the Jacobian? Default: False

    Returns
    -------
    All parameters estimated from the free water tensor model.
    Parameters are ordered as follows:
        1) Three diffusion tensor's eigenvalues
        2) Three lines of the eigenvector matrix each containing the
           first, second and third coordinates of the eigenvector
        3) The volume fraction of the free water compartment.

    References
    ----------
    .. [1] Henriques, R.N., Rokem, A., Garyfallidis, E., St-Jean, S., Peterson,
           E.T., Correia, M.M., 2017. Re: Optimization of a free water
           elimination two-compartmental model for diffusion tensor imaging.
           ReScience
    """
    # Initial guess
    params = wls_iter(design_matrix, sig, S0,
                      min_signal=min_signal, Diso=Diso)

    # Set bounds
    if bounds is None:
        bounds = ([0., -Diso, 0., -Diso, -Diso, 0., -10., 0],
                  [Diso, Diso, Diso, Diso, Diso, Diso, 10., 1])
    else:
        # In the helper subfunctions it was easier to have log(S0) first than
        # the water volume. Therefore, we have to reorder the boundaries if
        # specified by the user
        S0low = np.log(bounds[0][7])
        S0hig = np.log(bounds[1][7])
        bounds[0][7] = bounds[0][6]
        bounds[1][7] = bounds[1][6]
        bounds[0][6] = S0low
        bounds[1][6] = S0hig

    # Process voxel if it has significant signal from tissue
    if np.mean(sig) > min_signal and S0 > min_signal:
        # converting evals and evecs to diffusion tensor elements
        evals = params[:3]
        evecs = params[3:12].reshape((3, 3))
        dt = lower_triangular(vec_val_vect(evecs, evals))
        f = params[12]

        # Use the Levenberg-Marquardt algorithm wrapped in opt.leastsq
        start_params = np.concatenate((dt, [-np.log(S0), f]), axis=0)
        lb = np.array(bounds[0])
        ub = np.array(bounds[1])
        start_params[start_params < lb] = lb[start_params < lb]
        start_params[start_params > ub] = ub[start_params > ub]
        if jac:
            out = opt.least_squares(_nls_err_func, start_params[:8],
                                    args=(design_matrix, sig,
                                          Diso, False, False),
                                    jac=_nls_jacobian_func,
                                    bounds=bounds)
        else:
            out = opt.least_squares(_nls_err_func, start_params[:8],
                                    args=(design_matrix, sig,
                                          Diso, False, False),
                                    bounds=bounds)
        this_tensor = out.x

        # The parameters are the evals and the evecs:
        evals, evecs = decompose_tensor(from_lower_triangular(this_tensor[:6]))
        params = np.concatenate((evals, evecs[0], evecs[1], evecs[2],
                                 np.array([this_tensor[7]])), axis=0)
    return params


def nls_fit_tensor_bounds(gtab, data, mask=None, Diso=3e-3,
                          min_signal=1.0e-6, bounds=None, jac=True):
    """
    Fit the water elimination tensor model using the non-linear least-squares
    with constraints

    gtab : a GradientTable class instance
        The gradient table containing diffusion acquisition parameters.
    data : ndarray ([X, Y, Z, ...], g)
        Data or response variables holding the data. Note that the last
        dimension should contain the data. It makes no copies of data.
    mask : array, optional
        A boolean array used to mark the coordinates in the data that should
        be analyzed that has the shape data.shape[:-1]
    Diso : float, optional
        Value of the free water isotropic diffusion. Default is set to 3e-3
        $mm^{2}.s^{-1}$. Please ajust this value if you are assuming different
        units of diffusion.
    min_signal : float
        The minimum signal value. Needs to be a strictly positive
        number. Default: 1.0e-6.
    bounds : 2-tuple of arrays with 14 elements, optional
        Lower and upper bounds on fwdti model variables and the log of
        non-diffusion signal S0. Use np.inf with an appropriate sign to
        disable bounds on all or some variables. When bounds is set to None
        the following default variable bounds is used:
            ([0., -Diso, 0., -Diso, -Diso, 0., 0., np.exp(-10.)],
             [Diso, Diso, Diso, Diso, Diso, Diso, 1., np.exp(10.)])
    jac : bool
        Use the Jacobian? Default: False

    Returns
    -------
    fw_params : ndarray (x, y, z, 13)
        Matrix containing in the last dimension the free water model parameters
        in the following order:
            1) Three diffusion tensor's eigenvalues
            2) Three lines of the eigenvector matrix each containing the
               first, second and third coordinates of the eigenvector
            3) The volume fraction of the free water compartment

    References
    ----------
    .. [1] Henriques, R.N., Rokem, A., Garyfallidis, E., St-Jean, S., Peterson,
           E.T., Correia, M.M., 2017. Re: Optimization of a free water
           elimination two-compartmental model for diffusion tensor imaging.
           ReScience
    """
    fw_params = np.zeros(data.shape[:-1] + (13,))
    W = design_matrix(gtab)

    # Prepare mask
    if mask is None:
        mask = np.ones(data.shape[:-1], dtype=bool)
    else:
        if mask.shape != data.shape[:-1]:
            raise ValueError("Mask is not the same shape as data.")
        mask = np.array(mask, dtype=bool, copy=False)

    # Prepare S0
    S0 = np.mean(data[:, :, :, gtab.b0s_mask], axis=-1)

    index = ndindex(mask.shape)
    for v in index:
        if mask[v]:
            params = nls_iter_bounds(W, data[v], S0[v], Diso=Diso,
                                     min_signal=min_signal, bounds=bounds,
                                     jac=jac)
            fw_params[v] = params

    return fw_params
