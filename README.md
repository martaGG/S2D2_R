# S2D2: Small-scale Significant DBSCAN Detection

R code for the detection of small pristine structure with 3-sigma significance associated to position/velocity of stars (or young stellar objects) in starForming Regions. Developed under the SFM (StarFormMapper) EU project.

S2D2 uses DBSCAN detect the smallest significant structure in a spatial/spatial-kinematic space. The procedure proposes, in structured regions, calculation of the epsilon and Nmin parameters (as described in González et al., 2020) for DBSCAN to retrieve the smallest structures in the region with a user defined level of significance. If the region is not structured, or the user wants it, eps and Nmin can be supplied, and Dbscan will be performed with these parameters, without, however, guaranteeing any level of significance.

Usage example:

```
> source('main.R')
```

## Input

csv file with a table containing the input parameters. Must be provided in the adress hardcoded in main.R

Example file content:
```
filename, dim, coord, eps, Nmin, Qlim, signif
Taurus_RaDec.csv,2,Ra Dec,,, 0.7,99.85
```

### Input parameters

#### filename
Path for file with coordinates of stars/objects in your region. Uploaded 2D example file with Taurus coordinates.

V1: csv file with the coordinates: See options and expected values in the coord parameter section.

#### dim
Dimension of the space of search.
- V1: Integer, only=2.
- Future: limited options ('2D','3D','2+2D', '3+2D'...)

#### coord
Coordinate frame of the input, depending on the dimension. Cannot be left empty.
- 2D
  - 'Ra Dec' (without apostrophes) for great circle distance, assuming right ascension and declination in degrees.
  - any other string will calculate the euclidean distance

#### eps
Scale parameter supplied to DBSCAN, associated with the size of the structures to search.

- Default: empty field to search for the smallest scale in structured regions.
- If a float eps is supplied, Nmin must also be supplied.


#### Nmin
Number of points supplied to DBSCAN, associated with the density of a neighbourhood with radius eps.

- Default:empty field to calculate the Nmin guaranteeing  the required significance (see parameter signif).

- If an integer Nmin is supplied, eps must also be supplied.

#### Qlim
limit of Q parameter for considering a region structured, and calculate automatically the eps and Nmin values according to the procedure, as described in González et al. 2020.

We note that the classical limit Q parameter for structured regions is 0.8, and that a conservative limit of 0.7 avoids the possibility of analyising a region withous structure (as described in Gonzalez et al. 2020)

#### signif
If float, user supplied value (in percentage) required for structure retrieval
If empty or invalid, default strict value, 99.85 (>3 sigma).
## Output
- Console outputs some values of variables and status 
- Ascii file named as the input file with the extension .out containing the coordinates of each star in the region (in the  user supplied values) and an additional column with an integer representing the number of substructure assigned. Indexes follow the R convention,so the value 0 represents noise stars (those not assigned to any cluster).

- pdf file named as the input file with the extension .pdf with a plot of the region where:
  - grey stars are noise.
  - stars in significant structures are overplotted in colour. Each nest will be plotted in a different colour taken from a viridis colour table with as many different shades as NESTs, so the specific colours will depend on the amount of structures retrieved. 

## Requirements
R with libraries:
	- fpc
	- astrolibR
	- stats
	- viridis
	
## Description

### 2D
We refer to Gonzalez et al. 2020 and references therein for a complete description of the procedure.

### Structured regions
We will consider that a starForming region is structured when the Q parameter (Cartwright & Whitworth, 2003) is lower than the user supplied limit. 
In that case, and if the user has not provided eps and Nmin (default) we use our procedure to calculate them and obtain the smallest scale significant structure in the region.

#### eps calculation: Small-scale
We calculate the length scale for DBSCAN (epsilon) using the One point correlation function, OPCF (Joncour et al. 2017) , comparing the first nearest neighbour distance distribution of the sample with the first nearest neighbour distance distribution of a homogeneous random distribution (CSR, or complete spatial randomness) with intensity equal to the local density derived from the mean of the 6th neighbour distribution of the sample.

#### Nmin calculation: Significance
We iteratively calculate the significance of a structure of that scale and a fixed number of points k until we reach the significance value provided by the reader. The significance of a structure of size eps and k points, as described in Joncour et al 2018, is given by the the probability of having k-1 nearest neighbours in an eps neighbourhood under a homogeneous random distribution with intensity rho.

#### DBSCAN detection

We run scikit.learn’s DBSCAN for the previously calculated eps and Nmin. 


### Unstructured regions

If the Q parameter is larger than Qlim, and the user has not supplied an eps and Nmin, we display an error message (we cannot guarantee that the region is structured), and suggest the user to try the procedure providing eps and Nmin.

#### DBSCAN detection
We run fpc package DBSCAN for the user defined eps and Nmin.

## Acknowledging this
Please cite Gonzalez et al 2020 if you use this code. 

## References

To be completed
- Cartwright & Withworth, 2003
- González et al 2020.
- Joncour et al. 2017
- Joncour et al. 2018


