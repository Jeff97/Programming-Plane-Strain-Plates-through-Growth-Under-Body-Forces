Here, we provide some supplementary files for the manuscript[1]. 

- Movie 1 shows the simulated deformation process of different cases. To download the movie, please click the item, then click "View raw". 

- Each folder contains the ABAQUS input files and the corresponding UMAT subroutines

[1] **Shape-Programming Framework for Hyperelastic Plane-Strain Plates through Differential Growth in the Presence of Body Forces**, Zhanfeng Li, Xiaohu Yao, Mokarram Hossain, Jiong Wang*



# Simulation related files (needs expertise)

Each folder also contains necessary files to run the simulation using ABAQUS. To run the simulation, please make sure 

- Visual Studio and OneAPI are linked to ABAQUS, such that it can be used for subroutine development. [Here](https://www.researchgate.net/publication/349991987_Linking_ABAQUS_20192020_and_Intel_oneAPI_Base_Toolkit_FORTRAN_Compiler) is a linking guide from the internet. 

- Submit the job through ABAQUS COMMAND window


# Theory and numerical verification

This study aims to establish a framework for controlling the shape of hyperelastic plane-strain plates under body forces by prescribing non-uniform growth fields. 

We adopt a nonlinear plate theory of growth and then derive the asymptotic analytical solutions that reveal the relationship between the applied body forces, the required growth functions, and the geometric properties of the desired target shapes. 

![Theory](https://github.com/Jeff97/Programming-Plane-Strain-Plates-through-Growth-Under-Body-Forces/blob/main/Theory.jpg)

Building on these analytical results, we propose a shape-programming framework that involves calculating growth fields and designing appropriate loading paths. The framework is examined by 3D finite element simulations for representative examples, showing good quantitative agreement. 

![SimulationResults](https://github.com/Jeff97/Programming-Plane-Strain-Plates-through-Growth-Under-Body-Forces/blob/main/SimulationResults.jpg)

# Experimental results

Experiments based on silicone plates further validate the feasibility of the proposed shape programming approach. 

To achieve this, we employ an inverse design methodology. We first define the target shape for a plate clamped at one end, which is a flat plane in this case. From this target shape, we then calculate the required growth field based on the derived analytical formulae. This field is engineered to counteracts the bending deformation due to gravitational force. The balance of mechanical equilibrium ensures the plate settles into its intended shape. 

![ExperimentMethod](https://github.com/Jeff97/Programming-Plane-Strain-Plates-through-Growth-Under-Body-Forces/blob/main/ExperimentMethod.jpg)

To validate our theoretical framework, we designed six experimental settings. These cases explore the influence of some key factors on the shape-programming accuracy and effectiveness. We investigate two target shapes: a flat configuration and a downward-bending semicircle. Two different silicone elastomers, Ecoflex 00-30 and Dragon Skin 20, are selected to assess the effect of material properties. Furthermore, we fabricated samples with two different thicknesses, 10 mm and 5 mm, to investigate the role of structural stiffness.

![ExperimentResult1](https://github.com/Jeff97/Programming-Plane-Strain-Plates-through-Growth-Under-Body-Forces/blob/main/ExperimentResult1.jpg)

![ExperimentResult2](https://github.com/Jeff97/Programming-Plane-Strain-Plates-through-Growth-Under-Body-Forces/blob/main/ExperimentResult2.jpg)

Our results demonstrate that body force is not merely an external load but also a design variable in morphing structures. 