### GCA
Processing functions for Airborne Laser Scanning (ALS) and the Global Canopy Atlas (GCA).

All code is released under a GPLv3 license (https://www.gnu.org/licenses/). All publicly available products are released under a CC BY 4.0 license (https://creativecommons.org/licenses/by/4.0/).

### METADATA
Always use this code with a metadata file. The metadata file needs to be provided in csv format and end in _metadata.csv to be read by the pipeline. A typical name would be "GCA_metadata.csv" or "_GCA_metadata.csv". Details on how to fill out the metadata file can be found in the documentation (cf. DOCUMENTATION AND REFERENCE), and example can be found online on Zenodo (cf. WORKED EXAMPLES). Key information to be always provided: coordinate reference system (acq_crs) and acquisition dates (acq_mindt, acq_maxdt). Note that the pipeline may not run if these are incorrect/empty.

### SOFTWARE 
The pipeline is dependent on LAStools (https://rapidlasso.de). For most use cases, this requires a license (cf. https://rapidlasso.de/pricing/). Always use the most recent version of the pipeline and download a recent version of LAStools before processing. When systematic processing errors occur, it's useful to try both the most recent and a slightly older version of LAStools to rule out issues due to a recent update.

### DOCUMENTATION AND REFERENCE
Documentation can be found in INFO files (cf. ALS_processing folder) and in this publication: 

Fischer, F. J., Jackson, T., Vincent, G., & Jucker, T. (2024). Robust characterisation of forest structure from airborne laser scanning—A systematic assessment and sample workflow for ecologists. Methods in Ecology and Evolution, 15, 1873–1888. https://doi.org/10.1111/2041-210X.14416 

### WORKED EXAMPLES

A recent worked example with v.1.0.2 of the pipeline is available on Zenodo. It uses scans from multiple years of Dutch ALS campaigns (AHN2-5): https://zenodo.org/records/14722001

Data and code for the original paper in Methods in Ecology and Evolution can be found here: https://zenodo.org/records/10878070. Note that the pipeline has been updated since then, with improved error handling. It can also now be run on Linux, so using newer versions in this repository is strongly recommended.
