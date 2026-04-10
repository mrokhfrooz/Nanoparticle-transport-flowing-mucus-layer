# Nanoparticle-transport-flowing-mucus-layer
Computational modeling of nanoparticle transport in cilia-driven mucus flows, focusing on diffusivity, flow structure, and healthy vs. pathological conditions.

This repository is provided to support the submitted manuscript entitled:

**"Nanoparticle transport in a flowing mucus layer"**, submitted to *Nanoscale*.

---

## Cilia-driven flows submodule

The subfolder **`Cilia-driven-flows-streamlines`** computes the flow fields based on Section **2.3 (Flow kinematics)** of the manuscript.

The core simulation code  
`cilia_driven_flow_Choudhury_JFM2023.cpp`  
computes the flow field using the model described in:

> Choudhury, Anjishnu, et al.  
> "On the role of viscoelasticity in mucociliary clearance: a hydrodynamic continuum approach."  
> *Journal of Fluid Mechanics*, 971 (2023): A33.

---

## Governing equations

The implementation solves Equations (4.1a–b), (4.2), (4.3a), (A1a–d), and (A3a–d) from the original paper to obtain the velocity field:

\[
\mathbf{U} = (u, v)
\]

where:
- \(u\): velocity component in the x-direction  
- \(v\): velocity component in the y-direction  

---

## Numerical setup

- The computational domain is discretized on a uniform square grid.
- The resulting velocity field can be exported in `.csv` or Excel format for post-processing and visualization.
- Example output files include:  
  `square-6um-box-beta=0.csv`

---

## Visualization

The following Python scripts are provided for post-processing:

- `velocity_field_streamline.py` → generates streamline visualizations  
- `velocity_field_vectors.py` → plots velocity vector fields  

---

## Output

- Flow fields are printed during execution.
- Data can be saved for external visualization and analysis.

---

## Contact

For questions, please contact:

- Mohammad Reza Rokhforouz: mrokhforouz92@gmail.com  
