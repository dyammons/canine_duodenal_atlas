## Instructions to obtain processed Seurat objects  

The processed data are available on Zenodo and can be downloaded by visiting the [project repository web page](https://zenodo.org/records/11153534).
Once on the web page scroll down and select download for the file(s) of interest.

Alternatively, use `wget` in a terminal to retrieve the data:
```sh
wget https://zenodo.org/record/11153534/files/AllCells_duod_annotated.rds    # Full dataset
wget https://zenodo.org/record/11153534/files/Myeloid_duod_annotated.rds     # Myeloid cell dataset
wget https://zenodo.org/record/11153534/files/Tcell_duod_annotated.rds       # T cell dataset
wget https://zenodo.org/record/11153534/files/Epithelial_duod_annotated.rds  # Epithelial cell dataset
```

Prefer to use tools in python or R? Check out `zenodo_get` or `inborutils` to download within the respective software. 

<details><summary>Example usage of zenodo_get </summary>
<p>

Below is the code needed to install `zendo_get` using `pip` and the command to download the repositiry specific to this project (this should be completed in an environment with python3 installed).  

Visit the [`zendo_get`](https://github.com/dvolgyes/zenodo_get) page for most up to date instructions.

```sh
#install the python tool using pip
pip3 install zenodo_get

#download the Zenodo repository
zenodo_get 10.5281/zenodo.11153534
```

</p>
</details>

## Instructions to obtain count matrices from NCBI GEO  

### Retrieve the data

To download via command line, navigate to directory in which you wish to download the data and run the following:
```sh
#pull down the data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254005/suppl/GSE254005_RAW.tar

#unpack the tar ball
tar -xvf GSE254005_RAW.tar
```

To download via NCBI webpage, navigate to the [GSE254005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254005) page and download the zip folder containing the count matrices. Once downloaded, you will likely need to extract/unpack the data before using (should be able to do by right clicking on the `GSE254005_RAW.tar` file).

### Rearrange the data

The count matrices were renamed when uploaded to NCBI. If you wish to use the analysis in it's entirety you will need to rename the files accordingly. Feel free to reach out if you have any questions regarding renaming the files.

## Instructions to obtain raw data from SRA
Navigate to directory of interest and run for each file you wish to pull down. Then with SRA toolkit installed run:

NOTE: all raw data files are around 700 GB of data
```sh
prefetch -v --max-size=55000000 SRR27700549 #smallest file for ex
fastq-dump --gzip --split-files SRR27700549
```
From there you will have to modify the file names to be compatible with the Cell Ranger input (if using Cell Ranger). It expects:  
`[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`  

Feel free to reach out (email or create an issue on GitHub) if you have trouble getting the data downloaded.

