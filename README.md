# Getting Started with the Human Sinoatrial Node Model 

<p>Welcome to the Human Sinoatrial Node Model repository! This project implements the mathematical model of the human sinus node action potential based on the work by Fabbri et al. [<a href="#ref1"><strong>1</strong></a>]. This guide will help you get started with using the model.</p>

##

<p>The membrane potential equation for the human SAN cell model is given by:</p>

<p><img src="https://raw.githubusercontent.com/CellularSyntax/cellularsyntax.github.io/main/fabbri_san_membrane_voltage_equation.PNG" style=""></p>

<p>The currents in this equation, such as the <strong><em>I<sub>f</sub></em></strong> (Funny current), <strong><em>I<sub>CaL</sub></em></strong> (L-type calcium current), <strong><em>I<sub>CaT</sub></em></strong> (T-type calcium current), <strong><em>I<sub>Kr</sub></em></strong> (Rapid delayed rectifier potassium current), <strong><em>I<sub>Ks</sub></em></strong> (Slow delayed rectifier potassium current), <strong><em>I<sub>K,ACh</sub></em></strong> (Potassium current modulated by acetylcholine), <strong><em>I<sub>to</sub></em></strong> (Transient outward potassium current), <strong><em>I<sub>Na</sub></em></strong> (Sodium current), <strong><em>I<sub>NaK</sub></em></strong> (Sodium-potassium pump current), <strong><em>I<sub>NaCa</sub></em></strong> (Sodium-calcium exchanger current), and <strong><em>I<sub>Kur</sub></em></strong> (Ultra-rapid delayed rectifier potassium current), collectively contribute to the intricate electrical behavior of the sinoatrial node. They are crucial for its proper functioning in generating and conducting cardiac impulses.</p>

<p>Here, \( \frac{dV}{dt} \) is the rate of change of the membrane potential with respect to time, \( C \) is the membrane capacitance, and the terms inside the parentheses represent the various ionic currents that play a vital role in the depolarization and repolarization phases of the action potential.</p>

<p>This equation elegantly encapsulates the interactions and dynamics of the SAN cell membrane, highlighting the role of each ion channel in generating the characteristic rhythmic electrical activity of the SAN</p>


## Prerequisites 

<p>To use this model, you'll need:</p>
<ul>
  <li>Python (version 2.6 or later)</li>
  <li>NumPy library</li>
  <li>JSON library</li>
  <li>SciPy library</li>
  <li>Matplotlib library</li>
</ul>

## Installation 

<ol>
  <li>Clone this repository to your local machine:</li>
</ol>
<pre><code>$ git clone https://github.com/CellularSyntax/HumanSAN-Fabbri2017.git
$ cd HumanSAN-Fabbri2017
</code></pre>
<ol start="2">
  <li>Install the required libraries:</li>
</ol>
<pre><code>$ pip install numpy scipy matplotlib</code></pre>

## How to Use 

<ol>
  <li>Open the <code>config/config.json</code> file to view and modify the model parameters and initial conditions.</li>
  <li>Explore the <code>tests/test_san_model.json</code> file for an example of how to use the <code>SinoatrialNode</code> class to simulate the model and obtain results.</li>
  <li>The main implementation of the model is in the <code>model/SinoatrialNode.py</code> file. The class <code>SinoatrialNode</code> contains methods to initialize the model, update its conditions, calculate constants, and more.</li>
</ol>

### Minimal Example

To get started with the Sinoatrial Node Model, you can use the following minimal example to simulate the model and visualize the results using Python:

```python
from model.SinoAtrialNode import SinoAtrialNode
import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from helper_functions import parse_model_parameters

# Load the JSON file containing model parameters
constants, initial_conditions, constant_desc, init_cond_desc = parse_model_parameters("../config/config.json")

# Set the simulation duration and create the SinoAtrialNode object
sim_dur = 2
san = SinoAtrialNode(constant_descriptions=constant_desc,
                     state_descriptions=init_cond_desc,
                     constants=constants,
                     initial_conditions=initial_conditions)

# Solve the model using solve_ivp
sol = solve_ivp(san.calculate_derivatives, [0, sim_dur], list(san.y), method='BDF', rtol=1e-6,
                t_eval=np.arange(0, sim_dur, 1e-4), vectorized=False)

# Visualize the results
plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label="Transmembrane Potential")
plt.xlabel("Time (ms)")
plt.ylabel("Voltage (mV)")
plt.title("Sinoatrial Node Model Simulation")
plt.show()
```
<p>This code will produce the following output. </p>
<img src="https://raw.githubusercontent.com/CellularSyntax/cellularsyntax.github.io/main/simulation_plot.png" style="">

## Model Parameters and References 

<p>The values for constants and initial conditions in this JSON structure were obtained from the CellML model available <a href="https://models.cellml.org/e/568/HumanSAN_Fabbri_Fantini_Wilders_Severi_2017.cellml/view"><b>here</b></a>. These values were derived from the original publication:</p>
<ol>
  <li id="ref1">Fabbri, Alan, et al. "Computational analysis of the human sinus node action potential: model development and effects of mutations." <i>The Journal of Physiology</i> 595.7 (2017): 2365-2396. (DOI: <a href="https://doi.org/10.1113/JP273259">10.1113/JP273259</a>)</li>
</ol>

## Blog Article 

<p>For a detailed explanation of how this code was implemented and insights into modeling the transmembrane potential, please refer to the corresponding <a href="https://cellularsyntax.github.io/2023/08/12/modeling-the-human-sinoatrial-node.html"><b>blog article</b></a>.</p>

## Citation 

<p>If you use this model in your research or projects, please consider citing the original work by Fabbri et al. [<a href="#ref1"><strong>1</strong></a>] and acknowledging this repository.</p>

## References 
<ol>
  <li id="ref1">Fabbri, Alan, et al. "Computational analysis of the human sinus node action potential: model development and effects of mutations." <i>The Journal of Physiology</i> 595.7 (2017): 2365-2396. (DOI: <a href="https://doi.org/10.1113/JP273259">10.1113/JP273259</a>)</li>
</ol>
