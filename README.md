# galaxyMock

**galaxyMock** is a toolset to create a self-consistent model of a galaxy's stellar distribution, incorporating multiple stellar components, and to transform it into a realistic observed mock galaxy.

## Update History

**Last updated on:** November 5th, 2024  
**Author:** Guillaume Thomas \(*guillaume.thomas.astro@gmail.com*\)

---

## Requirements

- Python 3 is required.
- The installation script will automatically install any required Python packages.

---

## Installation

To install **galaxyMock**, simply run the following command in your terminal:


> ./install.sh


---

## Repository Structure

This repository contains the following files and directories:

- **make_model.py**: Generates a realistic stellar distribution model for a galaxy.
- **make_mock.py**: Projects the stellar model onto the sky and creates an observed mock galaxy, resembling real observational data.
- **iso/**: Contains isochrones from the Padova stellar library. Additional isochrone sets can be added to this folder as needed.
- **out/**: Contains the generated models and mocks.
- **examples/**: Contains example files that demonstrate how to use the package:
  - **mock_sculptor_GaiaDR3.fits**: Example of an observed-like mock of the Sculptor dSph model in Gaia DR3 based on the model of *Arroyo-Polonio et al. 2025*.
  - **params.ini**: Configuration file that is used to set parameters for generating the model of the galaxy.
  - **params_proj.ini**: Configuration file that is used to set parameters for making the observed-like mock galaxy.

---

## Example Usage

Here’s how to use the tools in **galaxyMock**:

### 1. **Generate a Galaxy Model**

You can use the **make_model.py** script to generate a realistic stellar distribution model of a galaxy. This will require modifying the **exemples/params.ini** configuration file to set the parameters for the galaxy model, such as stellar components, mass, age, etc.

**Example command to generate the model:**
> python make_model.py exemples/params.ini


This command reads the **exemples/params.ini** file and creates a stellar model similar to the Sculptor model of *Arroyo-Polonio et al. 2025*, which will be saved in a file called **out/model.fits**.

### 2. **Create a Mock Observation**

After generating a stellar model, you can use the **make_mock.py** script to project the model onto the sky and generate a mock observation. This will simulate real-world data, such as what might be observed by Gaia.

**Example command to generate the mock:**
> python make_mock.py out/model.fits exemples/params_proj.ini


This will take the galaxy model **out/model.fits** and project it into an observed mock, saved as **mock.fits**. This mock data will have the characteristics of a real-world observation, such as Gaia DR3.

### 3. **Example Files**

- **examples/mock_sculptor_GaiaDR3.fits**: This is the mock observation corresponding to the Sculptor galaxy model, generated to simulate Gaia DR3 observations.

You can modify the parameters in the **params_proj.ini** file and re-run the scripts to generate different models and mocks as needed.

---



