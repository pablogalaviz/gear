# GEAR, Genomic sEquence AnalyzeR.

GEAR (Genomic sEquence AnalyzeR) is CMRI's high performance multipurpose software. It is a C++/CUDA software designed for whole genome data analysis. GEAR works with second and third generation sequencing data. 

It is developed following best practices including: documentation, unit tests, code reviews, multi-threading  and multi-platform development. 

Currently there are three modules: 

MotifCount: finds and count DNA motifs and regular expressions in fasta, fastq and bam files. For aligned bam files it is possible to define regions of interest in the genome. 

VariantCallAnalysis: creates a summary analysis from vcf files. The analysis can be performed on regions of interest as well. 

TelomereAnalysis: uses a novel DNA sequence network approach to evaluate telomere variants in long sequence reads. 



## Getting Started

git clone repository 
```
> git clone https://cmri_user@bitbucket.cmri.com.au/scm/bioinf/gear.git
> cd gear
```

### Prerequisites

1. Build and install boost libraries >= 1.67.0 (default location /opt/boost) 
    ```
    > mkdir /path/to/boost/source
    > cd /path/to/boost/source
    > wget https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.gz
    > tar -zxvf boost_1_67_0.tar.gz
    > cd boost_1_67_0
    > ./bootstrap.sh --prefix=/opt/boost
    > sudo ./b2 install 
    ```
2. Install HSTlib (default location /usr/)
    ```
   > wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
   > bzip2 -d htslib-1.10.2.tar.bz2
   > tar -xvf htslib-1.10.2.tar
   > cd htslib-1.10.2
   > ./configure --prefix=/where/to/install
   > make
   > (sudo) make install 
    ``` 


## Installation
Edit library paths in CMakeList.txt if necessary. 

Inside GEAR root directory:

```
> mkdir build
> cd build 
> cmake .. 
> make 
> (sudo) make DESTDIR=/some/path/ install
```

## Usage

```
gear [options]
Allowed options::
  -b [ --backup ] arg (=1)       Create a backup of previous output
  -h [ --help ]                  Shows a help message
  -i [ --input ] arg             Input file
  -v [ --validate ] arg (=0)     Validate sequences (slow)
  -p [ --progress ] arg (=0)     Show progress message every X records (0 - 
                                 off)
  -o [ --output ] arg (=output)  Output directory name
  -m [ --motifs ] arg            Motif per region definition in json format
  -s [ --silent ]                Shows only errors
  -d [ --debug ]                 Shows debug messages in log
```
## History

First release 2020.1

## Credits

Author: Pablo Galaviz             
Email:  pgalaviz@cmri.org.au 


**Children’s Medical Research Institute, finding cures for childhood genetic diseases**  

## License

GEAR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

GEAR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GEAR.  If not, see <http://www.gnu.org/licenses/>.
