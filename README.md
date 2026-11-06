# biomass
[![License: MIT](https://shields.io)](https://opensource.org/licenses/MIT)
![Python](https://img.shields.io/badge/python-3.14-blue.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23992410.svg)](https://doi.org/10.5281/zenodo.23992410)

This package provides unit aware calculations of biomass for tree species and tall shrubs found in the 
Southwestern United States. The package loads data, searches for columns, defines units, and then calculates 
species level biomass using allometric equations. The package focuses on allometric equations that can 
provide both total biomass, as well as foliar biomass and available canopy fuel. 

## Acknowledgements
This project was supported by a grant from ESTCP (RC23-8363): Advancing Wildfire, Ecosystem Science and 
Technology Capacity Within an Innovation Landscapes Network

## Installation
1. Clone the repository
```bash
git clone https://github.com/SWFireInnovation/biomass.git
```
1. create and activate the conda environment
```bash
conda create --file enivronment.yaml
conda activate calc_biomass
```

The package is now ready to use.

## Usage Example
The labs folder contains many example Jupyter Labs used during the development of the project. Files labeled 
`*_methods.ipynb` contain the final methods used with this package for specific projects. They are good 
examples of how to use the package with some basic QC.

Note that this package will return no data for species without an equation. If a species is dominant in an 
area, but has no equation, it is recommended that an equation for another species be substituted. See 
Adjusting species_equ.yaml below.

### Fill out data_info.yaml
After installation, the first step is to change the file extension of `data_info.template`, to a `.yaml` file. 
Next, fill out the fields with the correct file paths and column names for your project:

```yaml
data_path: c:\Users\me\my_tree_data.xlsx
units:
  input:
    My DBH column: cm
    My Height column: m
    My CBH column: m
  out:
    weight: kg
```

### Adjusting species_equ.yaml
There is not an equation for every species. For example, there are allometric equations for 1 pinyon species. 
It could be reasonably argued that the biomass allometry of other pinyon species are similar enough to 
substitute the 1 pinyon species with a published equation. Each user must decide for themselves what 
substitutions are appropriate for their uses.

If a dominant species in your data has no equation, it will return no biomass. It is recommended that an 
equation from another species be substituted.

To change substitutions or equations used, alter `species_equ.yaml`. Here is an example where a Border Pinyon 
is substituting an equation for a 2-needle pinyon:

```yaml
PIDI3:
  name: Border Pinyon
  equ: PJ_Grier
  spp_eq: PIED
```

Note that equations are grouped into python classes by the paper that they are taken from. So, in this 
example, `PJ_Grier` is a python class containing all the equations from a paper published on Pinyon/Juniper 
species where Grier is the last name of the first author. A list of available python classes (i.e. equation 
groups) is listed below, or can be found in the file `calc_biomass/allometry.py`

### Python Terminal Use
First load the data:

```python
from calc_biomass import manage_data

data = manage_data.LoadData(data_path_yaml='../data_info.yaml')
data.load_data(data.data_path)
```

Next, calculate the desired biomass component for every tree in your dataset with an assigned equation:

```python
wt = manage_data.CalcBiomass(data)
# what equations should be used for these species?
wt.add_equ_column()

wt.calc_biomass('total')

wt.calc_biomass('foliage')

wt.calc_biomass('avl_canfuel')
```

## Available Allometric Equation Sets
Allometric equations are all in the sub-module `calc_biomass/allometry.py`. Each paper has its own class, 
and most or all species published in that paper have been programmed in that class. Class names are labeled by 
the types of trees included and the last name of the first author. For example `OakWoodland_Cojnacky` contains 
equations for oak woodland species in a paper published by Chojnacky. 

A current list of species and equation sets includes:

* OakWoodland_Chojnacky
  * *Quercus gambelii*
  * *Quercus emoryii*
  * *Quercus* (Uniformly used for *Q. arizonica, Q. hypoleucoides*, and *Q. oblongifolia* )
* PJ_Grier
  * *Pinus edulis*
  * *Juniperous monosperma*
* BCtimber_Standish
  * *Abies grandis*
  * *Abies lasiocarpa*
  * *Betula papyrifera*
  * *Larix occidentalis*
  * *Picea engelmannii*
  * *Pinus contorta* 
  * *Pinus monticola*
  * *Pinus ponderosa*
  * *Populus tremuloides*
  * *Populus trichocarpa*
  * *Pseudotsuga menziesii* (interior)
  * *Thuja plicata*
  * *Tsuga heterophyla*
  * *Tsuga mertensia*
* Mesquite_McClaren
  * *Prosopis veluntina* 
* MEXpines_VarasLarreta
  * *Juniperus deppeana*
  * *Pinus engelmannii*
  * *Pinus leiophylla*
  * *Pinus strobiformus*

## License
[MIT License](https://github.com/SWFireInnovation/biomass/blob/main/LICENSE)