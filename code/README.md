### Code repository

This folder contains four Python scripts:
  *  **netPD.py:** Main parameters of simulation is defined in this code. Use this script to run the simulation.
  *  **netParams.py:** Networks structure parameters.
  *  **netModels.py:** Function responsible to make connections between different populations of
                       neurons and to connect stimulation (background or DC currents).
  *  **neuronModels.py:** Neuron model equations and parameters.

Running the scripts:
--------------------

To run the model for the spontaneous activity in figure 2, type from the console:

``` 
python netPD.py
```

To run the experiment about dependence of spontaneous acitivity on input (figure 5), change the "bg_type" variable in line 41 or "stim" variable in line 45 from file netPD.py and run the simulation with the same command.

Platform Information:
---------------------

Platform: linux

gcc (GCC) 5.4.0

Python: 2.7.5 

brian: 2.0.1

NumPy: 1.11.1

SciPy: 0.19.0

Matplotlib: 2.0.0

Pandas:  0.19.2
