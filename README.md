# GetSubgenome_M3
These scripts used for dividing subgenomes

The mash-distance (M-distance) can be read in the acticle:

Ondov, B.D., Treangen, T.J., Melsted, P. et al. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biol 17, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x

The corresponding software <mash> used in the Div_m3 can be installed in the site https://github.com/marbl/mash

You can use the old version (shell language version) => Div_M3_pip.sh or

You can use the updated version (python language version) => GetSubgenome_M3.py

I recommend using the Python version, as it is highly convenient and extremely fast — subgenome partitioning can be completed within minutes. 

Below is a detailed description of GetSubgenome_M3.py.

Mash Chromosome Comparison and Visualization Tool

Description
This tool performs comparative analysis of chromosomes across genomes using Mash distances, with special support for polyploid genomes (diploid, triploid, tetraploid, etc.). It provides both computational analysis and visualization capabilities through a command-line interface.

Features
Supports any ploidy level (configurable number of subgenomes)

Complete pipeline from FASTA files to publication-ready visualizations

Modular design allowing separate execution of calculation and visualization steps

Colorblind-friendly visualizations with dynamic color schemes

Parallel processing for improved performance

Installation
Prerequisites
Software Dependencies:

seqkit (v2.0+ recommended)

Mash (v2.0+ recommended)

Python 3.8+

Python Packages:

pip install pandas matplotlib seaborn scipy numpy
Setup
Clone or download the script

Make the script executable:

bash
chmod +x mash_analysis.py
Usage
Command Structure
./GetSubgenome_M3.py <command> [options] [arguments]
Available Commands
all - Run complete pipeline (calculation + visualization)

./GetSubgenome_M3.py all [options] <ref_abb> <ref_fasta> <qry_abb> <qry_fasta>
calculate - Perform calculations only

./GetSubgenome_M3.py calculate [options] <ref_abb> <ref_fasta> <qry_abb> <qry_fasta>
plot - Generate visualization from pre-calculated results

./GetSubgenome_M3.py plot [options] <data_file>
Common Options
Option	Description	Default
-s, --subgenomes	Number of subgenomes to consider	2
-t, --threads	Number of parallel threads	4
-o, --output-pdf	Custom output PDF filename	Auto-generated
Examples
Basic diploid analysis:

bash
./GetSubgenome_M3.py all HSAP human.fasta OSAT rice.fasta
Tetraploid analysis with 8 threads:

bash
./GetSubgenome_M3.py all -s 4 -t 8 ATHA arabidopsis.fasta BRRA brassica.fasta
Calculate only (no visualization):

bash
./GetSubgenome_M3.py calculate HSAP human.fasta MMUS mouse.fasta
Visualize pre-calculated results:

bash
./GetSubgenome_M3.py plot HSAP_MMUS_mashDistance.filter.Gadd -o comparison.pdf
Output Files
Calculation Results:

[prefix]_mashDistance - Raw Mash distance results

[prefix]_mashDistance.filter.Gadd - Processed results with subgenome assignments

Visualization:

[prefix]_mashDistance.filter.Gadd.pdf - Publication-ready PDF visualization

Advanced Configuration
Handling Complex Genomes
For polyploid genomes, specify the expected number of subgenomes:

bash
# For hexaploid  analysis:
./GetSubgenome_M3.py all -s 6 REF REF.fasta QUERY QUERY.fasta
Customizing Visualization
The script automatically:

Generates distinct colors for each subgenome

Uses different markers (circles, squares, triangles, etc.)

Adjusts figure size based on ploidy level

Staggers labels to avoid overlap

Troubleshooting
Dependency Errors:

Ensure seqkit and mash are in your PATH

Verify Python package versions with pip freeze

Memory Issues:

For large genomes, reduce sketch size with -s parameter in build_mash_db()

Visualization Problems:

For >6 subgenomes, consider modifying generate_color_scheme() for more distinct colors

License
[ MIT License]

Citation
If you use this tool in your research, please cite:
xxx
