## GWRSFinder: Genome Wide Repeat Sequence Finder

<b>Aim:</b> Idetnify multicopy sequences within a genomic dataset (typically all genomes from the same species or genera) that can be used as a target for PCR. 

<b>Description:</b> This Nextflow pipeline uses the software following software to identify multicopy sequences: 
- Tallymer (part of the genometools software, http://genometools.org/pub/) to identify <i>k</i>mers repeated within a sequence, and their respective position within the sequence. 
- Bedtools (https://bedtools.readthedocs.io/en/latest/#) merges overlapping or bookended <i>k</i>mers. 
- BLAST command line tools (https://www.ncbi.nlm.nih.gov/books/NBK279690/) retrieves the copy number for each merged sequence in each genome in the user provided set of genomes.
  
<b>Requirements:<b>
1. Docker.
2. Nextflow.

<b>How to Run:</b>
1. Build the Dockerfile within the repository:
   
   ```
   docker build . -t gwrsfinder
   ```

2. pull ncbi/blast docker image from docker hub

    ```
    docker image pull ncbi/blast
    ```

3. Prepare the params.yml file located within the repository. 

4. Run the nextflow script from the command line:
   
   ```
   nextflow run workflow.nf -params-file params.yml 
   ```


<b>Output Files:</b>
1. PCR_targets.json: Contains the sequence name, total copy number, and copy number per genome in the database for each multicopy sequence found. 
2. PCR_targets.fasta: Contains the sequence of each multicopy sequence found. 
