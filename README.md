# taxMaps
taxMaps is an ultra-efficient, customizable and fully scalable taxonomic classification tool for short-read data designed to deal with large DNA/RNA metagenomics samples. Its performance and comprehensiveness makes it highly suitable for unbiased contamination detection in large-scale sequencing operations, microbiome studies comprising a large number of samples, and for applications where the analysis delivery time is a critical factor, such as pathogen identification from clinical or environmental samples. 

* Version: 0.2
* Author: Andre Corvelo, [New York Genome Center](https://www.nygenome.org)

taxMaps is freely available for academic and non-commercial research purposes ([`LICENSE.txt`](https://github.com/nygenome/taxmaps/blob/master/LICENSE.txt)).  

# Features
* Flexible input - `.bam`, `.fastq` (compressed/uncompressed, interleaved or not) or any custom command that streams `.fastq` to `STDOUT`
* Thorough preprocessing - quality trim, multiple adapter removal, hard trim (5' and 3'), low complexity filter.
* Mapping prioritization - through the use of multiple custom indexes (useful for spike-in and host subtraction).
* Performance - taxMaps takes full advantage of the outstanding features and performance of [GEM](http://www.nature.com/nmeth/journal/v9/n12/full/nmeth.2221.html). When using pre-compressed indexes, this state-of-the-art aligner can map in a very sensitive (up to 20% edit distance) manner several million reads per hour against a single-file nucleotide database encompassing more than 370 Gbases, including NCBI’s nt database, microbial RefSeq and thousands of microbial and viral genomes, on a single 16 CPU machine. 
* Detailed mapping reports - including coverage and edit-distance histograms per taxa, which can be visualized interactively.
* SGE support.  

# Manual
## List of contents
* [Installation](#installation)
   * [Dependencies](#dependencies)
   * [Download](#download)
* [Databases](#databases)
   * [Pre-built indexes](README.md#pre-built-indexes)
   * [Custom indexes](README.md#custom-indexes)
      * [Building the taxonomic table](README.md#building-the-taxonomic-table)
      * [Preprocessing and indexing your custom database](README.md#preprocessing-and-indexing-your-custom-database)
* [Running taxMaps](#running-taxmaps)
   * [Basic usage](#basic-usage)
      * [Output directory and sample prefix](#output-directory-and-sample-prefix)
   * [Input](#input)
      * [FASTQ](README.md#fastq)
      * [BAM](README.md#bam)
      * [Custom command](README.md#custom-command)
      * [Downsampling](README.md#downsampling)
   * [Preprocessing](#preprocessing)
      * [Quality scores encoding](README.md#quality-scores-encoding)
      * [Quality trimming](README.md#quality-trimming)
      * [Adapter removal](README.md#adapter-removal)
      * [Low complexity filtering](README.md#low-complexity-filtering)
      * [Hard-trimming your data](README.md#hard-trimming-your-data)
      * [Minimum read length](README.md#minimum-read-length)
   * [Mapping](#mapping)
      * [Index and length file](README.md#index-and-length-file) 
      * [Maximum edit distance](README.md#maximum-edit-distance)
      * [Number of CPUs](README.md#numbers-of-cpus)
      * [Multiple indexes](README.md#multiple-indexes)
   * [Taxonomic classification](#taxonomic-classification)
   * [Reporting](#reporting)
      * [Coverage histograms](README.md#coverage-histograms)
      * [Minimum reporting evidence](README.md#minimum-reporting-evidence)
   * [Output](#output)
      * [Output directory](README.md#output-directory)
      * [Merged mapping file](README.md#merged-mapping-file)
      * [Summary table](README.md#summary-table)
      * [Report table](README.md#report-table)
      * [Interactive abundance chart](README.md#interactive-abundance-chart)
   * [SGE](#sge)
   * [Dry run](#dry-run)
* [Quick reference](#quick-reference)
   * [./taxMaps](#taxmaps-1)
   * [./taxMaps-taxtbl](#taxmaps-taxtbl)
   * [./taxMaps-index](#taxmaps-index)
   * [./bin](#bin)

## Installation
#### Dependencies
* python 2.7 [www.python.org](https://www.python.org/)
* numpy 1.7.0 or higher [www.numpy.org](http://www.numpy.org/) 
* samtools 0.1.19 or higher [http://sourceforge.net/projects/samtools](http://sourceforge.net/projects/samtools/)
* cutadapt 1.4.1 or higher [code.google.com/p/cutadapt](https://code.google.com/p/cutadapt/)
* PRINSEQ-lite 0.20.3 or higher [prinseq.sourceforge.net](http://prinseq.sourceforge.net/)
* GEM pre-release 3 [sourceforge.net/projects/gemlibrary](http://sourceforge.net/projects/gemlibrary)
* Krona 2.4 [http://sourceforge.net/projects/krona](http://sourceforge.net/projects/krona/)

All of the executables from programs above should be in your `$PATH`  
For Krona, you will need to set `PERL5LIB` with the following command:  

```
$ export PERL5LIB=path_to_krona/lib/
```

After setting all the dependencies, confirm that you have the right version of `python` and `numpy`:

```
$ python --version ; python -c "import numpy; print 'NumPy ' + numpy.__version__"
```

and that you can run the following commands: `samtools`, `cutadapt`, `prinseq-lite.pl`, `gem-indexer`, `gem-mapper` and `ktImportText`.  

#### Download

taxMaps is available through github and can be obtained by using [`git`](http://git-scm.com/downloads) with following command:  

```
$ git clone git://github.com/nygenome/taxmaps.git
```

This will download the taxMaps repository locally into a `taxmaps/` directory.

Alternatively, you can go to the "Releases" tab in https://github.com/nygenome/taxmaps and download the latest version (`tar.gz`). In this case you need to extract taxMaps with the following command:

```
$ tar -zxvf taxmaps.{latest_version}.tar.gz
```

taxMaps doesn't require installation. You can copy `taxmaps/` anywhere you want and, to run `taxMaps` from anywhere, you just need make sure that `{path_to_taxMaps}/` gets added to your `PATH`:

```
$ export PATH=$PATH:{path_to_taxMaps}
```      

Adding line above to `.bash_rc` or `.bash_profile` configuration files will make it permanent.  


## Databases

### Pre-built indexes  
taxMaps uses pre-built indexes that are available via FTP on the [NYGC public FTP](http://tinyurl.com/krn3e7o) site. These indexes have been compressed using a kmer LCA pre-assignation for several different kmer sizes, which greatly improves mapping efficiency against large databases. To ensure maximum mapping fidelity, you  should pick a kmer size just above your read length or hard-trim down your reads to match the kmer size used in compression. For this, you can use the corresponding `taxMaps` option `-L`. You should also use the `taxonomy.tbl` file provided. 

### Custom indexes
We strongly recommend the use of pre-built compressed indexes. While we provide scripts to generate taxMaps-compatible indexes, we do not offer any support on it. 

#### Building the taxonomic table 
taxMaps requires a taxonomic table file for several of its operations. To generate this table you need to download [taxdump.tar.gz](http://tinyurl.com/lj5wkh8) from the NCBI FTP site and extract its contents with the following command:

```
$ tar -zxvf taxdump.tar.gz
```
  
Among other files, you should find `names.dmp` and `nodes.dmp`. The following command:

```
$ taxMaps-taxtbl -n names.dmp -t nodes.dmp > taxonomy.tbl
```
  
will generate a tab-delimited file named `taxonomy.tbl`, which contains one taxonomic entity per line with the following field order: `TaxId`, `Rank`, `Name` and `PathFromRoot`.  
 
#### Preprocessing and indexing your custom database
To generate your indexed database you need to download [gi_taxid_nucl.dmp.gz](http://tinyurl.com/oez34rw) from the NCBI FTP site and extract it using the following command:

```
$ gunzip gi_taxid_nucl.dmp.gz
```
  
into `gi_taxid_nucl.dmp`. This file contains the correspondence between `GI` and `TaxID`.

Ensure that your sequence database is in FASTA format and the FASTA headers are in the following NCBI format: `>gi|1234|...` - the `GI` (1234 in this case) must be in the second field of the `|`-delimited header. Also, your FASTA file should not contain any duplicated `GIs`.

At this point you should have all the required files to run `taxMaps-index` with the following command:

```
$ taxMaps-index -i {your_db}.fa -c gi_taxid_nucl.dmp -t taxonomy.tbl -p {prefix}
```
 
At the end of this process, among others, you should have 4 files:
* A formated `{prefix}.gitax.fa` file where the headers have been changed to `>GI_TaxID`.
* A `{prefix}.discarded.fa` file with sequences for which there was no associated `TaxID` (Due to NCBI's frequent updates, sometimes not all files are in 'sync'). These sequences will not show up in you index.
* A `{prefix}.len` tab-delimited file specifying the length of each sequence in your database containing the following fields: `SeqID` and `Length`.
* A `{prefix}.gem` index file for `gem-mapper`.  

Generating an index can take from a few minutes to several hours, depending on the size of your database. You can speed it up by running `taxMaps-index` on multiple threads by using the option `-T`.  You can also submit your job to a SGE cluster by specifying a partition/queue with the option `-Q`. The number of slots can be specified using the option `-S`.

## Running taxMaps

To run taxMaps, first make sure you have all the required files:
* Taxonomic table (`taxonomy.tbl`)
* Database index (`{db}.gem`)
* Length file (`{db}.len`)
* Short read sequences (`{reads}.fq`, `{reads}.bam` or any other format that can be converted into `.fq` with a single command to `STDOUT`)  

### Basic usage  

A simple taxMaps command should look something like this:

```
$ taxMaps -f reads.fq -d refseq_viral.gem -e 0.16 -c 8 -t taxonomy.tbl -p samp1 -o txM/sample1
```

The command above will use `8` CPUs to map all reads in `reads.fq` to `refseq_viral.gem`, allowing up to an edit distance of 16% (`0.16`). All the results will be written to `$PWD/txM/sample1` with the prefix `samp1`.  

##### Output directory and sample prefix

By default, taxMaps runs on the current working directory. It is **highly recomended** to use the option `-o` to specify an output directory that will contain all the results and log files. Doing it so, will keep your analysis results organized. It is also recommended to provide a sample prefix with `-p`. The following command:

### Input

##### FASTQ  
The most basic type of input for taxMaps is FASTQ (`.fq` or `.fastq`). You can specify a single file with the `-f` option.
In case you have chosen a 'paired-end' [taxonomic classification mode](https://github.com/nygenome/taxmaps#taxonomic-classification) (`-m p` or `-m P`), this file should be interleaved with the sequence identifiers for read 1 and read 2 ending in `/1` and `/2`, respectively. taxMaps can also deal with Casava 1.8 sequence identifiers (see [Illumina sequence identifiers](http://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers)).  

If your read data is separated in two files (one for read 1 and another for read 2), you can specify both files with the `-1` and `-2` options. Again, the sequence identifiers must be as explained above and both files should be in perfect synchrony. 

taxMaps accepts compressed FASTQ files (`.fq.gz` or `.fastq.gz`)

##### BAM  
If you are trying to map previously unmapped reads to a specific host and you have the alignments in `.bam` format, you can use the `-b` option. You need to make sure you have the corresponding `.bai` index file in the same folder. taxMaps will extract the unmapped reads from your `.bam` file. If you have chosen a 'paired-end' classification mode (`-m p|P`), taxMaps will only extract pairs for which both mates are unmapped.

##### Custom command  
Finally, taxMaps can take custom input commands as long as they write fastq to `STDOUT` by using the option `-i`. In this case, you must use the absolute paths to any file in your command. Here's an example:

```
$ taxMaps -i "zcat /data/project_x/sample*.fq.gz" ...
```  
  
##### Downsampling  
If you wish to downsample your data, you can specify the downsample probability by specifying a float number between 0 and 1 through the option `-D`. For instance, if you want taxMaps to classify 10% of your data, you should use: 

```
$ taxMaps -f reads.fq -D 0.1 ...
```  

### Preprocessing

With taxMaps you can process your data prior to mapping. Controlling the quality of the input data usually improves mapping sensitivity and specificity.  

##### Quality scores encoding  
taxMaps expects quality scores to be encoded in Phred+33 format (see [FASTQ encoding](http://en.wikipedia.org/wiki/FASTQ_format#Encoding)). If your quality data is encoded in Phred+64 format, you should add the `--phred64` flag to your command (usually not required for Illumina 1.8+, Sanger, Roche/454, Ion Torrent, PacBio data). 

##### Quality trimming  
Using the option `-q` you can set the threshold for 3' quality-based trimming of your read data. taxMaps will pass this option to `cutadapt`. If you wish to remove low quality bases from both ends of the reads, just specify the quality threshold values for the 5' and 3' ends, separed by `,` (e.g. `-q 20,20`). 

##### Adapter removal  
To remove adapter sequence, you can pass a string of `cutadapt` options by using the option `-a`. For instance, if you want to remove the adapter sequence `AGATCGGAAGAGC` from both ends of your reads in two passes, you should specify `-a "-b AGATCGGAAGAGC -n 2"`. You can also trim multiple adapter sequences - see [cutadapt user manual](https://cutadapt.readthedocs.org/en/latest/guide.html). You **should not** pass any cutadapt options related to input and output as these will likely cause taxMaps to crash (e.g., `-a "-b AGATCGGAAGAGC -o output.fq"`).

##### Low complexity filtering  
taxMaps, by default and through `prinseq`, filters out low-complexity reads defined by entropy values lower than `70` or containing more than `4`% of ambiguous positions (Ns). Both cutoffs can be changed with the options `-C` and `-N`, respectively. 

##### Hard-trimming your data  
To setup an upper-bound for read length, you can use the option `-L`. taxMaps will trim your reads down to the specified length from the 3' end. To trim your read data from the 5' end, you should use the option `-w` and specify the number of bases to be clipped. 

##### Minimum read length  
Finally, with the option `-l` you can control the minimum length a read should have to be eligible for mapping. It defaults to `50`. Very low values will lead to poor mapping specificity.  

### Mapping

##### Index and length file  
In taxMaps, you must specify the database index file with the option `-d`. Just make sure that for every index `{db}.gem` you have the corresponding `{db}.len` in the same directory. 

##### Maximum edit distance  
By default taxMaps allows up to 20% edit distance. This threshold can be decreased in situations where a greater precision is desired. With the option `-e` you can set the maximum edit distance allowed. You must provide a number between 0 and 1 and edit distance will be considered as a fraction of the read length. For instance, if your read length is 150bp, the following command:

```
$ taxMaps -f reads.fq -d ncbi.gem -e 0.2 ...
```

will allow a maximum absolute distance of 30 edits. It is not recomended to use values greater than 20% of the read length as it can severely affects mapping performance and/or specificity.

##### Numbers of CPUs  
Mapping performance is directly dependent on the number of threads used by the mapper (GEM). This is particularly true when mapping against very large databases, such as NCBI's nt. With the `-c` option you can increase the number of threads used for mapping from the default of `1`. The following command:

```
$ taxMaps -f reads.fq -d ncbi.gem -c 16 ... 
```
will use `16` CPUs.  

##### Multiple indexes  
If you want to map against more than one index in a sequential manner, you should provide a `,`-delimited list of index files:

```
$ taxMaps -f reads.fq -d phix.gem,human.gem,ncbi.gem ...
```
The command above will first map all reads against `phix.gem`. Then, unmapped reads will be mapped against `human.gem` and, finally, reads that were not mapped will be mapped against the `ncbi.gem` index. When mapping against more than one index, you can still control the abovementioned parameters in a index-specific manner by providing `,`-delimited lists of values for `-e` and `-c`:

```
$ taxMaps -f reads.fq -d phix.gem,human.gem,ncbi.gem \
          -e 0.1,0.2,0.2 -c 4,16,16 ...
```

### Taxonomic classification

To run taxMaps, you need to specify through the `-t` option, the path to the taxonomic table downloaded from the [NYGC public FTP](http://tinyurl.com/krn3e7o).

For read classification you can choose between three LCA modes:
* 'single-end' (`-m s`) - reads are classified independently. This is the default option.
* 'paired-end' (`-m p`) - reads are classified as pairs. In case one of the mates has been discarded in the filtering step or simply didn't map, the pair will be classified solely based on the other mate. This mode is more specific than the 'single-end' mode.
* 'paired-end strict' (`-m P`) - reads are also classified as pairs but only when both mates are mapped. This is the more specific classification mode.  

For the paired end modes you need to provide either an interleaved FASTQ file, two 'in-sync' FASTQ files, an indexed BAM file where unmapped mates are contigous or a custom command that writes interleaved FASTQ to `STDOUT`. 

### Reporting

##### Coverage histograms  
taxMaps can compute basepair coverage histograms for all the taxa represented in your sample. This option is turned OFF by default but you can turn it ON by using the flag `--cov`. Bare in mind that by setting this option ON will increase the memory requirements of taxMaps.

For sake of computational efficiency, you can exclude specific taxa from the coverage analysis. For that, you need to provide a `,`-delimited list of TaxIds using the option `-x`. This option can be used, for instance, in situations in which you are not interested in coverage over certain sequences, such as spike-ins or hosts. The following command:

```
$ taxMaps -f reads.fq -d phix.gem,human.gem,ncbi.gem --cov -x 9606,374840 ...
```

will not compute any coverage histogram for phiX (`374840`) and human (`9606`) sequences.

##### Minimum reporting evidence  
In taxMaps you can also exclude from your report, taxa with few reads assigned. These are usually non-specific hits or mapping artifacts. By default, taxMaps will not report taxa with less than 10 reads assigned per million reads in your sample (`-z 0.00001`). You can also specify a minimum number of reads regardeless of the proportion. For that you should use integer values:

```
$ taxMaps -f reads.fq -d ncbi.gem -z 15 ...
```

The command  above will report taxa with at least `15` reads classified as such.

### Output

##### Output directory  
By default, taxMaps will run on the current working directory, but it is strongly recommended to specify a different output directory through the option `-o` for every run.  


The  output directory contains the following subdirectories:
* `txM.{prefix}.base/` - contains either links to your read data file(s) or a `.sh` script with custom input command specified using the option `-i`, all the links to the specified indexes and corresponding length files, and the taxonomy table specified by `-t`.
* `txM.{prefix}.map/` - in this directory you will find all the `.map` files and a `.map.lca` file with the mapping results plus the taxonomic classification for every read (see below for format) .
* `txM.{prefix}.out/` - this is where summary results will be written to. It also contains the graphical and tabular reports. 
* `txM.{prefix}.log/` - where you can find the log files for `cutadapt`, `prinseq`, `GEM` and `samtools`(if your input data is in `bam` format).

In the output directory you will also find a `.sh` file with all commands that taxMaps will run. Also, if you used the `-Q` option for SGE, you will find a `.sge` custom specification file and an additional subdirectory will be created upon submission, containing several files related to the SGE run.

##### Merged mapping file  
After mapping, taxMaps will merge all the `.map` files and perform the taxonomic classification of every read, based on the all its best hits. These results can be found in `{prefix}.merged.map.lca` located in `{out_dir}/txM.{prefix}.map/`. The format is an extended [GEM alignment format](http://algorithms.cnag.cat/wiki/FAQ:The_GEM_alignment_format):  

1. `Name` - Read name.  
2. `Sequence` - Read sequence.  
3. `Quality` - Read quality.
4. `MatchSummary` - histogram-like `:`-delimited string of the number of hits per edit distance (e.g., `0:1:0:2` means that there was `1` hit at distance 1 and `2` hits at distance 3).
5. `Alignments` - `,`-separated list of alignments (see [GEM format](http://algorithms.cnag.cat/wiki/FAQ:The_GEM_alignment_format) for complete description).  
6. `TaxId` - NCBI taxon identifier.
7. `TaxRank` - NCBI taxon rank.
8. `TaxName` - NCBI taxon name.  
9. `PairMappingStatus` - two character code mapping status (`F` filtered out, `U` unmapped, `M` mapped - one per mate). **\***
10. `PairedTaxId` - NCBI taxon identifier. **\***
11. `PairedTaxRank` - NCBI taxon rank. **\***
12. `PairedTaxName` - NCBI taxon name. **\***

**\*** The last 4 fields are exclusive for paired-end classification modes. 

##### Summary table  
taxMaps then summarizes all mapping information into `{prefix}.merged.map.lca.summary` located in `{out_dir}/txM.{prefix}.out/`. This file consists of one-line records for every taxonomic entity identified in the sample or with at least one hit. Here are the fields:

1. `TaxId` - NCBI taxon identifier.
2. `TaxRank` - NCBI taxon rank.
3. `TaxName` - NCBI taxon name.
4. `TaxPath` - Path from the root (`1`) of the taxonomic tree, given as a `:`-separated list of `TaxId`s.
5. `ReadsClassified`  - Number of reads classified as `TaxName`.
6.  `ReadsMappedSpec` - Number of reads mapped to `TaxName`in a specific manner.
7.  `ReadsMappedNoSpec` - Number of reads mapped to `TaxName` in a non-specific manner. 
8. `PairsStrict` - Number of read-pairs classified as `TaxName`, where both mates have the same mapping status (`MM`, `UU` or `FF` as `PairMappingStatus`).
9. `PairsNoStrict`Number of read-pairs classified as `TaxName`, for mates with different mapping status (`MF`, `FM`, `MU`, `UM`, `UF` or `FU` as `PairMappingStatus`).
10. `DistanceHistogram` - `:`-delimited distance histogram string. For instance if `TaxName` had two reads mapping at distance `0` and three at distance `2`, the edit histogram string should look like this: `2:0:3`. 
11. `NumberOfAlignments` - Total number of alignments to `TaxName`.
12. `BasesAligned` - Total number of read bases aligned to `TaxName` sequences.
13. `BasesCovered` - Number of `TaxName` bases covered.
14. `TotalBases` - Total number of bases in `TaxName`sequences.
15. `CoverageHistogram` - The coverage histogram string follows a *log* progression, with the first value representing the number of uncovered bases, the second - the number of 1x covered bases, the third - the number of 2-3x covered based, the fourth - the number of 4-7x covered bases and so on... The following string `10:6432:456:0:120:0...` means that there were `10` uncovered bases, `6432` bases at 1x coverage, `456` at 2-3x, `0` between 4-7x and `120` between 8 -15x coverage.

##### Report table  
In this tab-delimited file (`{out_dir}/txM.{prefix}.out/{prefix}.tbl`) you will find the read counts for the classic taxonomic ranks (`kingdom`, `phylum`, `class`, etc...). Read counts for a given taxon are not limited to reads classified as that taxon, but also include all reads mapping to the taxa under it. These are the fields in the report table:

1. `TaxId` - NCBI taxon identifier.
2. `TaxRank` - NCBI taxon rank
3. `TaxLevel` - Number of nodes in `TaxPath`.
4. `TaxPath` - Path from the root (`1`) of the taxonomic tree, given as a `:`-separated list of `TaxId`s.
5. `TaxName` - NCBI taxon name
6. `ReadCounts` - Number of reads classified under `TaxName`
7. `PercentFromTotal` - Corresponding percentage relative to the root of the report.

##### Interactive abundance chart  

Finally, taxMaps generates an independent `.html` document that allows you to explore in an interactive manner your results (`{out_dir/txM.{prefix}.out/{prefix}.krona.html`). For that, you only need a Web browser and internet access. In the interactive HTML5 chart generated using [Krona](http://sourceforge.net/p/krona/home/krona/), taxomic entities are displayed as nested sectors from the top level (center) to the bottom (outward) of the hierarchy. Here is a screenshot:

![image](doc/img/example_krona.png)  

For more information about Krona, see the project [website](http://sourceforge.net/p/krona/home/krona/).

### SGE  

If you work on a SGE environment, you can submit your taxMaps job directly to a SGE cluster by specifying the partition/queue with `-Q`, the number of slots with `-S`. For instance, the command:

```
$ taxMaps -f reads.fq -d ncbi.gem ... -Q research.q -S 20
```

will generate a custom file with SGE specifications and submit your taxMaps job to the `research.q` queue, requesting `20` slots. 

### Dry run

You can do a dry run of taxMaps by using the flag `--dry`. taxMaps will create the output directory, link the input files and generate the master `.sh` script as usual but it will not launch it. In case you used the option `-Q` for a SGE run, taxMaps will also generate the SGE specification file but it will not submitt the job. Dry runs are good for debugging purposes and to ensure that the path to the input files are correct. Moreover, you can hack a taxMaps run by altering the `.sh` script to your convenience and launching independently.

## Quick reference
### ./taxMaps
`taxMaps` is pipeline generating script. Given an input, list of indexes and options it generates an a working directory and .sh script with all commands to be executed. If a SGE partition/queue is specified with the `-Q` option, it submits the job.  

##### Options:  
`--version` show program's version number and exit  
`-h`,`--help` show this help message and exit   

*Input*  
`-i` Input command (absolute paths!). Quoted. Interleaved `fq` on `STDOUT`. [`STR, Mandatory | -f | (-1 & -2) | -b`]  
`-f` Input `fq`, `fastq`, `fq.gz` or `fastq.gz`. Interleaved. [`CSV_FILE_LIST, Mandatory | -i | (-1 & -2) | -b`]  
`-b` Input `.bam` file. Requires `.bam.bai` in the same folder. [`FILE, Mandatory | -i | -f | (-1 & -2)`]  
`-1` Input `fq`, `fastq`, `fq.gz` or `fastq.gz` read 1. In sync with, and of the same sort as `-2`. [`CSV_FILE_LIST, Mandatory | -f | -i | -b`]  
`-2` Input `fq`, `fastq`, `fq.gz` or `fastq.gz` read 2. In sync with, and of the same sort as `-1`. [`CSV_FILE_LIST, Mandatory | -f | -i | -b`]  
`-D` Downsampling probability. [`FLOAT, 1`]  

*Preprocessing*  
`-q` Cutadapt quality cutoff value (default = None). [`INT, `]  
`-a` Adapters string. Quoted.  e.g., "-a AGATCGGAAGAGC -n2 "  (default = None). [`STR, `]  
`-l` Minimum read length for mapping. [`INT, 50`]  
`-L` Maximum read length (hard trimmming). [`INT, `]  
`-w` 5p trim length (hard trimmming). [`INT, `]  
`-C` Entropy cutoff for low complexity filtering. Use 0 for no filtering. [`INT, 70`]  
`-N` Filter reads with more than N% of 'N' characters. Use 100 for no filtering. [`INT, 4`]   
`--phred64` Quality scores in Phred+64 format. [`FLAG, OFF`]

*Mapping*  
`-d` Index file(s). CSV. [`FILE_LIST, Mandatory`]  
`-e` Edit distance(s). CSV. Value(s) between 0 and 1. [`FLOAT_LIST, 0.2`]  
`-c` Number of CPUs. CSV. [`INT_LIST, 1`]  

*Taxonomy*  
`-t` Taxonomic table file. [`FILE, Mandatory`]  
`-m` Taxa determination mode ('s' single-end; 'p' paired-end; 'P' paired-end strict). [`STR, s`]  

*Reporting*  
`--cov` Compute coverage histograms. [`FLAG, OFF`]  
`-x` Excluded taxids. CSV (only if --cov flag ON). [`STR_LIST, `]  
`-z` Reporting cutoff (INT or FLOAT < 1). [`FLOAT, 0.00001`]  

*Miscellaneous*  
`-p` Sample prefix. [`STR, sample`]  
`-o` Output directory. [`DIR, .`]  
`-Q` Cluster queue/partition. [`STR, `]  
`-S` Cluster slots (default = 20). [`INT, 20`]  

`--dry` Dry run. [`FLAG, OFF`]

### ./taxMaps-taxtbl     
`taxMaps-taxtbl` is a script that generates the taxonomic table required by `taxMaps` containg all taxids, corresponding descriptions and paths-from-root. It writes to `stdout`.  

##### Options:  
`--version` show program's version number and exit.  
`-h`,`--help` show this help message and exit.  

`-n` NCBI Taxonomy names.dmp file. [`FILE, Mandatory`]   
`-t` NCBI Taxonomy nodes.dmp file. [`FILE, Mandatory`]  

### ./taxMaps-index
`taxMaps-index` is a script that formats the FASTA headers and generates the GEM indexes required by `taxMaps`.  

##### Options:  
`--version` show program's version number and exit.  
`-h`,`--help` show this help message and exit.  

`-i` Input fasta file. [`FILE, Mandatory`]  
`-c` Correspondence file gi -> tax file. [`FILE, Mandatory`]  
`-t` Taxonomy table file. [`FILE, Mandatory`]  
`-p` Prefix for output files. [`STR, Mandatory`]  

`-T` Number of threads. [`INT, 1`]  
`-Q` Cluster queue/partition. [`STR, `]  
`-S` Cluster slots. [`INT, 1`]  

`--dry` Dry run. [`FLAG, OFF`]  

### ./bin    
`txM_fastalen` Computes the length of `FASTA` sequences.   
`txM_fq2gem` Converts `FASTQ` into GEM mapper format.  
`txM_fqhtrim` Trims from the 3' end of the reads down to the specified length.  
`txM_fqintlv` Interleaves read1 and read2 FASTQ files.  
`txM_fqltrim` Trims the specified number of bases from the 5' end of reads.  
`txM_fqmate` Checks for presence of mates when extracting from `BAM`.   
`txM_fqminlen` Filters reads taht are shorter than the specified read length.  
`txM_fqsample` Downsamples `FASTQ`.  
`txM_gem2fq` Converts GEM mapper format into `FASTQ`.  
`txM_gitax` Formats `FASTA` header for taxMaps.  
`txM_lca` Given a GEM mapper output, computes the LCA for all reads, either in sigle-end or paired-end mode.  
`txM_mapout` Filters mapped reads.  
`txM_mergintlv` Merges and interleaves multiple GEM mapper files assuming the pairing order is consistent across files.   
`txM_report` Generates a tables and plots with the number of reads mapping for all represented taxa.  
`txM_rescore` Rescores alignments from edit to Levenshtein distance.  
`txM_sge` Submits jobs specified in a custom format to a SGE cluster.  
`txM_summary` Generates a summary per taxon with hits, including edit distance and coverage histogram. 


