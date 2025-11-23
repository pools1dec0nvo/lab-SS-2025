#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inicialização para o trabalho de laboratório de Sinais e Sistemas

Created on Nov 2022, updated Nov. 2024

@author: pmqa
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.signal import freqz

# Unit step function
u = lambda t: (t >= 0).astype(float)

# Unit ramp function
s = lambda t: (t >= 0) * t

# Approximate unit step function
uD = lambda t, D: (1 / D) * (s(t + D / 2) - s(t - D / 2))

# Numeric derivative of the unit step
duD = lambda t, D, h: (uD(t + h / 2, D) - uD(t - h / 2, D)) / h

# Approximate delta function
deltaD = lambda t, D: duD(t, D, 1e-6)

def seqsin(fs, *args):
    """
    Generate a sequence of sinusoids with specified frequencies and durations.
    
    Parameters:
    fs   - Sampling frequency (Hz)
    *args - Variable argument list with pairs of (frequency, duration) for each sinusoid.
            A frequency of zero corresponds to a period of silence.
    
    Returns:
    x    - Array containing the sequence of sinusoids
    """
    # Check if the number of arguments is even
    if len(args) % 2 != 0:
        raise ValueError("seqsin: Number of (frequency, duration) arguments must be even")
    
    # Duration in samples for the initial and final tapers
    n_in = round(0.001 * fs)  # Initial taper duration
    n_out = round(0.01 * fs)  # Final taper duration
    
    # Create fade-in and fade-out tapers
    fade_in = 0.5 * (1 + np.cos(np.pi / n_in * np.arange(-n_in, 0)))
    fade_out = 0.5 * (1 + np.cos(np.pi / n_out * np.arange(-1, -n_out - 1, -1)))
    
    # Initialize the output signal
    x = np.array([])
    
    # Loop through frequency-duration pairs
    for i in range(0, len(args), 2):
        freq = args[i]
        duration = args[i + 1]
        nsamples = round(duration * fs)
        
        # Generate the sinusoid or silence
        if freq == 0:
            sinusoid = np.zeros(nsamples)
        else:
            sinusoid = np.sin(2 * np.pi * freq / fs * np.arange(nsamples))
        
        # Apply fade-in and fade-out tapers
        if nsamples > n_in + n_out:
            sinusoid[:n_in] *= fade_in
            sinusoid[-n_out:] *= fade_out
            
        # Append to the output signal
        x = np.concatenate((x, sinusoid))
    
    return x

def rir(fs, mic, n, r, rm, src):
    """
    Room Impulse Response using the mirror image method.
    
    Parameters:
    fs  - sample rate
    mic - list or array [x, y, z] coordinates of the microphone
    n   - virtual source count factor (program will account for (2*n+1)^3 sources)
    r   - reflection coefficient for the walls, typically -1 < r < 1
    rm  - list or array [x, y, z] dimensions of the room
    src - list or array [x, y, z] coordinates of the sound source
    
    Returns:
    h   - impulse response vector
    """
    # Indices for sequence
    nn = np.arange(-n, n + 1)
    
    # Calculate rms and srcs components
    rms = nn + 0.5 - 0.5 * np.float_power(-1, nn)
    srcs = np.float_power(-1, nn)
    
    # Calculate x, y, z positions (Equation 2, 3, 4)
    xi = srcs * src[0] + rms * rm[0] - mic[0]
    yj = srcs * src[1] + rms * rm[1] - mic[1]
    zk = srcs * src[2] + rms * rm[2] - mic[2]

    # Create 3D matrices from x, y, z positions (meshgrid equivalent)
    i, j, k = np.meshgrid(xi, yj, zk, indexing='ij')
    
    # Calculate distances (Equation 5)
    d = np.sqrt(i**2 + j**2 + k**2)
    
    # Convert distance to time samples (Equation 6, adjusted)
    time = np.round(fs * d / 343).astype(int) + 1
    
    # Compute reflection coefficients (Equation 9)
    e, f, g = np.meshgrid(nn, nn, nn, indexing='ij')
    c = r ** (np.abs(e) + np.abs(f) + np.abs(g))
    
    # Calculate e as the equivalent to Equation 10
    e = c / d
    
    # Flatten and convert to sparse matrix (Equivalent to Equation 11)
    h_sparse = csr_matrix((e.ravel(), (time.ravel(), np.zeros_like(time.ravel()))))
    
    # Convert to dense array and normalize
    h = h_sparse.toarray().flatten()
    h = h / np.max(np.abs(h))
    
    return h

def cftransform(x, i0=0, n=512, fs=1.0, *args):
    """
    Evaluate the continuous-time Fourier transform of a sampled continuous-time signal.

    Parameters:
    x (ndarray): Vector of input samples
    i0 (int): Index of sample corresponding to time t = 0
    n (int): Number of evaluation points
    fs (float): Sample rate
    *args: Optional arguments for freqz ('whole')

    Returns:
    X (ndarray): (Sampled) continuous-time Fourier transform
    f (ndarray): Frequency axis (Hz)
    """

    # Compute the frequency response of x
    f, X = freqz(x, a=1, worN=n, fs=fs, *args)
    X = (X / fs) * np.exp(1j * 2 * np.pi * f * i0 / fs)

    # Convert frequencies to range [-fs/2, fs/2]
    f = (f + fs/2) % fs - fs/2

    return X, f

from system2 import *

from system3 import *

print("\n\nSinais e Sistemas - trabalho de laboratório: inicialização concluída.\n")
