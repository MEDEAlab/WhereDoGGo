### Description
```
                           ▄              ▄
                          ▌▒█           ▄▀▒▌
                          ▌▒▒█        ▄▀▒▒▒▐
                         ▐▄▀▒▒▀▀▀▀▄▄▄▀▒▒▒▒▒▐
                       ▄▄▀▒░▒▒▒▒▒▒▒▒▒█▒▒▄█▒▐
                     ▄▀▒▒▒░░░▒▒▒░░░▒▒▒▀██▀▒▌
                    ▐▒▒▒▄▄▒▒▒▒░░░▒▒▒▒▒▒▒▀▄▒▒▌
                    ▌░░▌█▀▒▒▒▒▒▄▀█▄▒▒▒▒▒▒▒█▒▐
                   ▐░░░▒▒▒▒▒▒▒▒▌██▀▒▒░░░▒▒▒▀▄▌
                   ▌░▒▄██▄▒▒▒▒▒▒▒▒▒░░░░░░▒▒▒▒▌
                  ▌▒▀▐▄█▄█▌▄░▀▒▒░░░░░░░░░░▒▒▒▐
                  ▐▒▒▐▀▐▀▒░▄▄▒▄▒▒▒▒▒▒░▒░▒░▒▒▒▒▌
                  ▐▒▒▒▀▀▄▄▒▒▒▄▒▒▒▒▒▒▒▒░▒░▒░▒▒▐
                   ▌▒▒▒▒▒▒▀▀▀▒▒▒▒▒▒░▒░▒░▒░▒▒▒▌
                   ▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒░▒░▒░▒▒▄▒▒▐
                    ▀▄▒▒▒▒▒▒▒▒▒▒▒░▒░▒░▒▄▒▒▒▒▌
                      ▀▄▒▒▒▒▒▒▒▒▒▒▄▄▄▀▒▒▒▒▄▀
                        ▀▄▄▄▄▄▄▀▀▀▒▒▒▒▒▄▄▀
                           ▒▒▒▒▒▒▒▒▒▒▀▀
__        ___                   ____         ____  ____      ___
\ \      / / |__   ___ _ __ ___|  _ \  ___  / ___|/ ___| ___|__ \
 \ \ /\ / /| '_ \ / _ \ '__/ _ \ | | |/ _ \| |  _| |  _ / _ \ / /
  \ V  V / | | | |  __/ | |  __/ |_| | (_) | |_| | |_| | (_) |_|
   \_/\_/  |_| |_|\___|_|  \___|____/ \___/ \____|\____|\___/(_)
```

**WhereDoGGo?** (Where Does the Genome Go?) is a user-friendly pipeline for prokaryotic phylogenomics. Its main function is to perform accurate placement of genomes in the taxon of interest all the way to domain level (Bacteria or Archaea), while addressing sources of error in (deep) phylogenies.It includes all necessary steps for any phylogenetic analysis, from picking and downloading genomes, to concatenating protein marker sets, and running a multitude of phylogenetic analyses accounting for various types of errors. WhereDoGGo? comprises of four main modules/steps, plus a number of auxiliary scripts. Below is a short description for each module:

1. **doggo_fetch**: `doggo_fetch.py` will download genomes and create a local genome database based on GTDB taxonomy, genome stats, and user input.

   - **Options**:
     ```
     -i, --input (required): INPUT must be a tab-delimited file as per GTDB\'s metadata files for Bacteria or Archaea.
     -lvl, --level (required): LEVEL must be the highest taxonomic level for genomes to be selected as presented in GTDB taxonomy strings e.g., p__Asgardarchaeota (i.e., get genomes from within LEVEL).
     -res, --resolution (required): RESOLUTION must be a taxonomic level lower than LEVEL. For each RESOLUTION in LEVEL, NUMBER genomes are picked (or as many as available). Must be p (phylum), c (class), o (order), f (family), g (genus), or all..
     -n, --number (required): NUMBER must be the number of genomes to be picked per RESOLUTION. Must be a positive integer or all. If RESOLUTION is all, NUMBER must also be all.
     -ig, --ignore (optional): IGNORE must be a text file containing genome assembly accessions (one per line, versionless e.g., GCA_011362025) that will not be picked.
     ```

   - **Example usage**:
     ```
     python doggo_fetch.py -i metadata.tsv -lvl p__Asgardarchaeota -res f -n 2 -ig ignore_list.txt
     #This will fetch the two highest-quality genomes for each family in p__Asgardarchaeota. It will also try to maximize the taxonomic range, by not repeating a genus (the level below family), unless one genome from every genus has already been picked.
     ```

   **Note**: The metadata input file must be parsed, i.e., it must contain one genome per assembly accession number. Parsed metadata files (only the representative genome of each species) for Bacteria and Archaea are included in the metadata/ directory. Ignore lists for Bacteria and Archaea are included in the ignore/ directory.

2. **doggo_herd**: `doggo_herd.py` will create a database for a set of local genomes (contig files).

   - **Options**:
     ```
     -loc, --location  (required): LOCATION must be a directory containing all genome files. The trailing slash will be added, if not included.
     -ext, --extension  (required): EXTENSION must be the file extension of the genome files. The leading dot will be added, if not included.
     -prj, --project (optional): PROJECT must be a project code consisting of alphanumeric characters and/or underscores only. This code will be used as the filename stem of the generated files. If not provided, it will default to five random alphanumeric characters.
     -code25, --code25 (optional): CODE25 must be provided if the local genomes are SR1 or Gracilibacteria (c__JAEDAM01). All ORF predictions will use Genetic Code 25.
     ```

   - **Example usage**:
     ```
     doggo_herd.py -loc /genome_files -ext .fna -prj project_1
     ```

3. **doggo_sniff**: `doggo_sniff.py` will run homology searches, alignments, trimming, concatenating (only markers with sequences > 50% taxa), and all the intermediate parsing steps.

   - **Options**:
     ```
     -db, --databases (required): DATABASES must be one or more multi-FASTA files (e.g., from the output of doggo_fetch or doggo_herd) against which the hmm searches will be run.
     -hmm, --hmm (required): HMM must be a directory containing the .hmm files (HMM profiles) for the markers.
     -con, --concatenation (optional): CONCATENATION must be the filename stem of the concatenation output containing alphanumeric characters and/or underscores only. If not provided, it will default to five random alphanumeric characters.
     -cut, --cutoffs (optional): CUTOFFS must be a tab-delimited file with two columns, the marker HMM profile filename and its domain bitscore cutoff as it would be input in HMMER. The HMM profile filename must only contain alphanumeric characters and/or underscores and use the .hmm extension. The domain bitscore must be a number, with or without decimals. For any files not included or if this argument is not provided, default domain bitscore cutoff is 30.
     -f, --fuse (optional): FUSE will run an additional step to fuse fragmented adjacent sequences (up to three fragments, four or more will be ignored), based on their accessions. WARNING: This option is still experimental, use it with caution and manually check your final alignments and concatenation.
     ```

   - **Example usage**:
     ```
     python doggo_sniff.py -db databases.faa -hmm /hmm_folder -con proj_1 -cut cutoffs.txt -f
     ```

     **Note**: A variation of doggo_sniff.py, called doggo_sniff_fftns.py has been included. It is identical except that it runs multiple alignments with the faster but less accurate FFT-NS-2 algorithm of MAFFT instead of E-INS-i.
     **Note 2**: HMM profiles and their manuall manually determined cutoffs are included in the hmm/ and cutoffs/ directories respectively, for the GTDB marker sets (bacteria120, archaea53) and 16 contiguous ribosomal proteins (for Bacteria and Archaea).

4. **doggo_zoomies**: `doggo_zoomies.py` will run all the different phylogenetic analyses in IQ-TREE. (TODO: Add a section on the errors that the different analyses address.)

   - **Options**:
     ```
     -i, --input (required): INPUT must be the input FASTA file.
     -MFP, --MFP (optional): MFP will run IQ-TREE with MODELFINDER to select the model (matrices: LG, WAG, JTT; frequencies: FU,F,FO).
     -C60, --C60 (optional):C60 will run IQ-TREE with the C60 mixture model under the matrix picked by MFP (and with the MFP phylogeny as guide tree under the PMSF approximation) or LG (if the MFP log and phylogeny not available) and 10 FreeRate categories i.e., LG (or WAG or JTT )+C60+R10.
     -SR4, --SR4 (optional): SR4 will recode data to the 4-state SR alphabet and run IQ-TREE as with the MFP option but with the GTR4 matrix.
     -SR4C60, --SR4C60 (optional): SR4C60 will recode data to the 4-state SR alphabet and run IQ-TREE with the SR4C60 model. If available, the SR4 phylogeny will be used as guide tree.
     -D6, --D6 (optional): DD6 will recode data to the 6-state Dayhoff alphabet and run IQ-TREE as with the MFP option but with the GTR6 matrix.
     -nex, --nex (optional, unless the SR4C60 option or DESAT option with the SR4C60 argument is selected): NEX must be the input Nexus file for running mixture model analyses on recoded datasets.
     -GHOST, --GHOST (optional): GHOST will run IQ-TREE with the GHOST heterotachy model with the maximum number of mixture categories calculated based on the dataset's number of sequences and positions.
     -desat, --desat (optional): DESAT will create a set of desaturated datasets at 5% incements (from 5% to 90% of positions retained) for each of the arguments given and run IQ-TREE with the appropriate options. For any argument given, the corresponding option must also be picked or have been run previously to calculate site-specific rates for the desaturation. DESAT arguments can be MFP, C60, SR4, or SR4C60.
     -AU, --AU (optional): AU will take a user-defined set of leaves and, if they are monophyletic, place them on all possible positions on the tree, optimize branch lengths, and run the Approximately Unbiased test on them versus one or more of the phylogenies run previously. The latter are specified with the same names as their respective options. AU arguments can be MFP, C60, SR4, SR4C60, or D6.
     -AUclade, --AUclade (optional, unless the AU option is provided): AUclade must be a text file with the clade (one leaf name per line) whose position on the phylogeny will be tested.
     -leaves, --leaves (required): LEAVES must be a text file with complete leaf names that will be used to converting them in the phylogenies e.g., an .assembliesnames file. All instances of tab-delimitation will be converted to spaces.

     ```

   - **Example usage**:
     ```
     python doggo_zoomies.py -i concatenation.faa -MFP -C60 -SR4 -SR4C60 -nex nexus_file.nex -desat -AU -AUclade clade.txt -leaves leaf_names.assembliesnames
     ```

   **Note**: All scripts used by the four modules can be run individually as well. Check the individual scripts for usage information.

5. **Auxiliary scripts**: WhereDoGGo? includes a number of auxiliary scripts to facilitate the analyses.
    pickgenomes_dry.py: This script will do a dry run (i.e., no downloads) of checking which genomes would be downloaded by a doggo_fetch run. Use it to determine the parameters of your doggo_fetch run.
    ```
    Usage: python pickgenomes_dry.py <input_tsv> <tax_level> <tax_resolution> <number> <ignore_list>
    <input_tsv> must be a tab-delimited file as per GTDB\'s metadata files for Bacteria or Archaea. (required)
    <tax_level> must be the highest taxonomic level for genomes to be selected as presented in GTDB taxonomy strings e.g. p__Asgardarchaeota (i.e., get genomes from within <tax_level>). (required)
    <tax_resolution> must be a taxonomic level lower than <tax_level>. For each <tax_resolution> in <tax_level>, the script picks <number> genomes (or as many as available). It must be p (phylum), c (class), o (order), f (family), g (genus), or all. "all" picks all genomes for <tax_level>. (required)
    <number> must be the number of genomes to be picked per <tax_resolution> or "all". If <tax_resolution> is "all", <number> must also be "all". (required)
    <ignore_list> must be a text file containing genome assembly accessions (once per line, versionless, e.g. GCA_011362025) that should be not included in the pickgenomes process (optional)
    ```
    subsampledb.sh: This script will create a local database that is a subset of a prexisting database. Use it to avoid having to rerun doggo_fetch if, for example, you've already downloaded all genomes from Bacteria and/or Archaea.
    ```
    Usage: subsampledb.sh <assemblies> <inputdb>
    <assemblies> must be a text file of versionless assemblies that will be subsampled from inputdb (1/line). (required)
    <inputdb> must be a local database in FASTA format. (required)
    ```
    parsematadata.py: This script will parse a file containing GTDB metadata and output only the lines containing the representative genome for each species. Use it for the Methods section of publications using WhereDoGGo?.
    ```
    Usage: python parsemetadata.py <input_file> <output_file>
    <input_file> must be tab-delimited GTDB metadata. (required)
    <output_file> must be the name of the output file that will contain the representative genome of each species cluster. (required)
    ```
    checkversions.py: This script will print the version of Python, the libraries, and other programs that WhereDoGGo? will call.
    ```
    Usage: python checkversions.py
    ```

### Dependencies

We highly recommend creating a fresh conda environment to install WhereDoGGo? dependencies. This will ensure that no package incompatibilities will arise. Here is a list with all WhereDoGGo? dependencies and links for downloading them (we suggest conda and/or pip):
IMPORTANT NOTE: Currently WhereDoGGo? is up-to-date to work with Python 3.12 and Biopython 1.83 and there will be issues if you try running it with older versions.

1. [Biopython (v1.83 or higher)]
 ```
 https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython
 Cock, P. J. A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25, 1422–1423 (2009).
 ```
2. [ETE3]
 ```
 https://anaconda.org/etetoolkit/ete3 or https://pypi.org/project/ete3/
 Huerta-Cepas, J., Serra, F. & Bork, P. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. Mol Biol Evol 33, 1635–1638 (2016).
 ```
3. [NumPy]
 ```
 https://numpy.org/install/
 Harris, C. R. et al. Array programming with NumPy. Nature 585, 357–362 (2020).
 ```
4. [pandas]
 ```
 https://github.com/pandas-dev/pandas
 McKinney, W. Data structures for statistical computing in Python. in SciPy vol. 445 51–56 (2010).
 ```
5. [ncbi-datasets-cli]
 ```
 https://github.com/ncbi/datasets
 (use github link to cite)
 ```
6. [Pyrodigal]
 ```
 https://github.com/althonos/pyrodigal
 Larralde, M. Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. JOSS 7, 4296 (2022).
 ```
7. [HMMER]
 ```
 https://anaconda.org/bioconda/hmmer
 Eddy, S. R. Accelerated Profile HMM Searches. PLoS Comput Biol 7, e1002195 (2011).
 ```
8. [Pullseq]
 ```
 https://anaconda.org/bioconda/pullseq
 https://github.com/bcthomas/pullseq
 ```
9. [MAFFT]
 ```
 https://anaconda.org/bioconda/mafft
 #WhereDoGGo? uses E-INS-i and FFT-NS-2.
 ```
10. [BMGE]
 ```
 https://anaconda.org/bioconda/bmge
 Criscuolo, A. & Gribaldo, S. BMGE (Block Mapping and Gathering with Entropy): a new software for selection of phylogenetic informative regions from multiple sequence alignments. BMC Evol Biol 10, 210 (2010).
 ```
11. [IQ-TREE 2 (v2.2 or higher)]
 ```
 https://anaconda.org/bioconda/iqtree
 Minh, B. Q. et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Molecular Biology and Evolution 37, 1530–1534 (2020).
 Kalyaanamoorthy, S., Minh, B. Q., Wong, T. K. F., Von Haeseler, A. & Jermiin, L. S. ModelFinder: fast model selection for accurate phylogenetic estimates. Nat Methods 14, 587–589 (2017).
 Hoang, D. T., Chernomor, O., Von Haeseler, A., Minh, B. Q. & Vinh, L. S. UFBoot2: Improving the Ultrafast Bootstrap Approximation. Molecular Biology and Evolution 35, 518–522 (2018).
 Guindon, S. et al. New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0. Systematic Biology 59, 307–321 (2010).
 Wang, H.-C., Minh, B. Q., Susko, E. & Roger, A. J. Modeling site heterogeneity with posterior mean site frequency profiles accelerates accurate phylogenomic estimation. Syst. Biol. 67, 216–235 (2018).
 Susko, E. & Roger, A. J. On Reduced Amino Acid Alphabets for Phylogenetic Inference. Molecular Biology and Evolution 24, 2139–2150 (2007).
 Shimodaira, H. An Approximately Unbiased Test of Phylogenetic Tree Selection. Systematic Biology 51, 492–508 (2002).
 ```

### Installation

No installation is needed for WhereDoGGo?. Just download the compressed file, extract it to your desired location, and place the ‘bin’ directory in your PATH.

On Linux, this can be done by adding the following line at the end of the hidden ‘.bashrc’ file:

```
export PATH="/path/to/WhereDoGGo_version/bin/:$PATH"
```

### Citation
The manuscript for WhereDoGo? is currently under preparation. For now, please cite this GitHub repository. Also cite the programs and libraries used by the different modules as above.
In a publication using WhereDoGGo?, for example, to determine the taxonomic placement of a suspected new prokaryotic lineage (as suggested by a run of GTDB-Tk), we would suggest something along these lines for the relevant Methods section (version numbers just examples, let's assume a new class of p__Asgardarchaeota):
```
The taxonomic placement of XYZ was further determined using WhereDoGGo? v20240621 (https://github.com/MEDEAlab/WhereDoGGo) run under Python 3.12.3. A local database of the XYZ genomes was created with doggo_herd (-ext .fna -prj XYZplacement). Open-reading frames (protein sequences) were predicted with Pyrodigal v3.4.1 (Larralde, 2022). To confirm that the GTDB-Tk taxonomy is correct, we created a local genome database using doggo_fetch (-lvl d__Archaea -res o -n 1 -ig) with the parsed metadata from GTDB r220 for Archaea and the included ignore list for WhereDoGGo? v20240621. Genomes were downloaded with ncbi-datasets-cli 16.21.0 (https://github.com/ncbi/datasets), and ORFs predicted with Pyrodigal 3.4.1. Then doggo_sniff was used to create a concatenation of the 120 bacterial marker genes of GTDB. The XYZplacement and d__Archaea_o_1 databases were combined and homologs were searched using HMMER v3.4.1 (Eddy, 2011) with the bitscore cutoffs for this WhereDoGGo? version. Sequences were extracted from the combined database with Pullseq (https://github.com/bcthomas/pullseq), adjacent fragmented sequences fused automatically based on their protein accessions, and the taxa in each marker that still had multiple sequences were removed. Datasets were aligned with MAFFT 7.526 (E-INS-i) (Katoh & Standley, 2013) and the alignments trimmed with BMGE 1.12 (Criscuolo & Gribaldo, 2010) and concatenated. The concatenation was used to run phylogenies with doggo_zoomies (-MFP -C60 -SR4 -D6 -SR4C60 -GHOST -desat MFP -AU MFP) running IQ-TREE 2.3.4 (Minh et al., 2020). The following phylogenies were run: 1) model automatically selected by MODELFINDER (Kalyaanamoorthy et al., 2017) (-mset JTT,WAG,LG -mfreq FU,F,FO), 2) Posterior Mean Site Frequency (PMSF) model (Wang et al., 2018) with the matrix selected by MFP and 10 Free-rate categories (LG+C60+F+R10), 3) recoded concatenation under the Susko-Roger 4-state reduced alphabet (Susko & Roger, 2007) (-mset GTR -mfreq FU,F,FO), 4) recoded concatenation under the Dayhoff 6-state reduced alphabet (Susko & Roger, 2007 and references therein), 5) SR4 recoded alignment with PMSF (C60), 6) a series of phylogenies with progressively desaturated subsets of the original concatenation, under the model automatically selected by Modelfinder as above, 7) testing all possible alternative positions of the XYZ clade on the tree using the implementation of the Approximate Unbiased test (Shimodaira, 2002) in IQ-TREE for the MFP phylogeny. All branch supports were calculated with 1000 ultrafast bootstrap (Hoang et al., 2018) and 1000 aLRT SH-like (Guindon et al., 2010) replicates, and branches with at least 95 for ultrafast bootstraps and 80 for aLRT SH-like were considered strongly supported as per the IQ-TREE manual. Non-standard libraries used by WhereDoGGo? were Biopython 1.83 (Cock et al., 2009), ETE3 3.1.3 (Huerta-Cepas et al., 2016), NumPy 1.26.4 (Harris et al., 2020), and pandas 2.2.2 (McKinney, 2010).
```
