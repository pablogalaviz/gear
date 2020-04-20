# GEAR, Genomic sEquence AnalyzeR.

Finds motif sequences using regular expressions in fasta, fastq, csv, or bam files. 
For bam files it is possible to search motifs in mapped regions. 

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

    Follow the [instructions](https://github.com/samtools/htslib)

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
