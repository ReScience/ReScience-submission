from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import uncertainpy as un


onset_threshold = 0.55          # Fraction of the normalized voltage
end_threshold = -0.1            # Relative to the onset_threshold
burst_threshold = 60            # ms
min_spike_amplitude = 10        # mV


def duration(time, voltage):
    """
    Calculate the duration of a series of events in a voltage trace.

    Parameters
    ----------
    time : array_like
        Time array, in ms.
    voltage : array_like
        Voltage array from a neuron.

    Returns
    -------
    duration : numpy.array
        A numpy array containing the duration of each spike (in ms) in `voltage`.

    Notes
    -----
    Normalizes the voltage before the spikes are calculated.
    """

    # Find spikes in the normalized voltage trace
    spikes = un.features.Spikes(time,
                                voltage,
                                threshold=onset_threshold,
                                end_threshold=end_threshold,
                                trim=False,
                                normalize=True,
                                min_amplitude=min_spike_amplitude)

    # Calculate the duration of each spike
    duration = []
    for spike in spikes:
        duration.append(spike.time[-1] - spike.time[0])

    return np.array(duration)



def burstiness(event_durations, burst_threshold=burst_threshold):
    """
    Calculate the burstiness from a series of event durations by finding the
    fraction of events that is above `burst_threshold`.

    Parameters
    ----------
    event_durations : array
        A numpy array containing the duration of each spike in `voltage`.
    burst_threshold : float, optional
        Default is 70 ms.

    Returns
    -------
    burstiness_factor : float, None
        The fraction of events that is above `burst_threshold`. Returns None
        if there are no events.
    """
    if not event_durations.size:
        return None

    # Find all events greater than the burst_threshold
    burst_occurrences = event_durations > burst_threshold

    # Calculate the fraction of these events
    burstiness_factor = np.sum(burst_occurrences)/len(burst_occurrences)

    return burstiness_factor



def burstiness_factor(time, spikes, info):
    """
    Calculate the burstiness factor from a uncertainpy.Spikes object.
    Compatible with uncertainpy.SpikingFeatures.

    Parameters
    ----------
    time : {None, numpy.nan, array_like}
        Time values of the model. If no time values it is None or numpy.nan.
    spikes : uncertainpy.Spikes
        Spikes found in the model result.
    info : dictionary
        Not used, but is required as input

    Returns
    -------
    time : None
    burstiness_factor : float, None
        The fraction of events that is above `burst_threshold`. Returns None
        if there are no events.
    """
    event_durations = []
    for spike in spikes:
        event_durations.append(spike.time[-1] - spike.time[0])

    event_durations = np.array(event_durations)

    burstiness_factor = burstiness(event_durations, burst_threshold=burst_threshold)

    return None, burstiness_factor

