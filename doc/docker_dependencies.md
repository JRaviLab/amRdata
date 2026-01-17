
# For Macs with silicon
## Pull the docker image for CD-HIT
> docker pull weizhongli1987/cdhit:4.8.1 
### run cd-hit docker
> docker run -v `pwd`:/data -w /data weizhongli1987/cdhit:4.8.1 cd-hit ........ 
## Pull the docker image for BV-BRC CLI
> docker pull danylmb/bvbrc:5.3
### run bvbrc docker
> docker run danylmb/bvbrc:5.3 p3-all-genomes --fields
### Panaroo image is not available for arm64 so running the image using emulation
> docker run --platform linux/amd64 staphb/panaroo:latest
## Pull the docker image for Interproscan
> docker pull interpro/interproscan:5.76-107.0
### follow the given instruction to download the background data   
> curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/alt/interproscan-data-5.76-107.0.tar.gz
> 
> curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/alt/interproscan-data-5.76-107.0.tar.gz.md5
>
> md5sum -c interproscan-data-5.76-107.0.tar.gz.md5

#### if md5sum doesn't work in MAC
##### install homebrew (if not available)
> /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
> 
> echo >> /Users/jravilab/.zprofile
> 
> echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> /Users/jravilab/.zprofile
> 
> eval "$(/opt/homebrew/bin/brew shellenv)"  
>  
> brew install coreutils 
> 
> PATH="/opt/homebrew/opt/coreutils/libexec/gnubin:$PATH" 

##### Try md5sum again 
> md5sum -c interproscan-data-5.76-107.0.tar.gz.md5

##### Finally extract the data files
> tar -pxzf interproscan-data-5.76-107.0.tar.gz
>
> docker run --rm -v "$PWD/interproscan-5.74-105.0/data:/opt/interproscan/data" -v "$PWD:/work" interpro/interproscan:5.74-105.0


# For alpine like linux based HPCs pull the dockers to create singularity images!
> apptainer pull docker://danylmb/bvbrc:5.3
> 
> apptainer pull docker://weizhongli1987/cdhit:4.8.1
> 
> apptainer pull docker://interpro/interproscan_5.76-107.0
> 
> curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/alt/interproscan-data-5.76-107.0.tar.gz
> 
> tar -pxzf interproscan-data-5.76-107.0.tar.gz
> 
> apptainer pull docker://staphb/panaroo:latest
