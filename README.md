# Getting Started with the Sinoatrial Node Model 

<p>Welcome to the Sinoatrial Node Model repository! This project implements the mathematical model of the human sinus node action potential based on the work by Fabbri et al. [<a href="#ref1"><strong>1</strong></a>]. This guide will help you get started with using the model.</p>

## Prerequisites 

<p>To use this model, you'll need:</p>
<ul>
  <li>Python (version X.X or later)</li>
  <li>NumPy library</li>
  <li>JSON library</li>
</ul>

## Installation 

<ol>
  <li>Clone this repository to your local machine:</li>
</ol>
<pre><code>git clone https://github.com/yourusername/your-repo-name.git
cd your-repo-name
</code></pre>
<ol start="2">
  <li>Install the required libraries:</li>
</ol>
<pre><code>pip install numpy
</code></pre>

## How to Use 

<ol>
  <li>Open the <code>config/config.json</code> file to view and modify the model parameters and initial conditions.</li>
  <li>Explore the <code>tests/test_san_model.json</code> file for an example of how to use the <code>SinoatrialNode</code> class to simulate the model and obtain results.</li>
  <li>The main implementation of the model is in the <code>model/SinoatrialNode.py</code> file. The class <code>SinoatrialNode</code> contains methods to initialize the model, update its conditions, calculate constants, and more.</li>
</ol>

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
