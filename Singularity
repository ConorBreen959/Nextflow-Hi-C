BootStrap: docker
From: ubuntu:latest

%labels
    PROJECT RNA-Seq-1.0

%pre
    apt-get install -y debootstrap

%post
    apt-get update
    apt-get install -y wget
    apt-get install -y gzip
    apt-get install -y bzip2
    apt-get install -y curl
    apt-get install -y unzip

    ## g++
    apt-get install -y build-essential
    
    # install anaconda
    if [ ! -d /usr/local/anaconda ]; then
       wget https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh \
       	    -O ~/anaconda.sh && \
	    bash ~/anaconda.sh -b -p /usr/local/anaconda && \
	    rm ~/anaconda.sh
    fi

    # set anaconda path
    export PATH=$PATH:/usr/local/anaconda/bin
    conda update conda

    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    
    # Let us save some space
    conda clean --packages -y

    # external tools
    echo "Installing external tools ... "
    conda install -y bowtie2
    conda install -y samtools

    # Python and relevant libraries
    echo "Installing Python ... "
    conda install -y -c conda-forge python=3.7
    conda install -y -c anaconda scipy 
    conda install -y -c anaconda numpy
    conda install -y -c bcbio bx-python
    conda install -y -c bioconda pysam 

    # Install R and relevant packages
    echo "Installing R ... "
    conda update readline
    conda install -c r r-base=3.6
    conda install -c r r-ggplot2
    conda install -c r r-rcolorbrewer
    conda install -c r r-gridbase	

    # Install MultiQC
    conda install -c bioconda multiqc 
   
    # Install HiC-pro
    echo "Installing latest HiC-Pro release ..."
    VERSION="devel"
    echo $VERSION".zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && unzip hicpro_latest.zip
    
    cd $(echo HiC-Pro-$VERSION)
    make configure
    make install
 
    # Let us save some space
    conda clean --packages -y
    conda clean --all -y
    rm -rf /usr/local/anaconda/pkgs

%test
    INSTALLED_HICPRO_VERSION=$(find /usr/local/bin -name HiC-Pro | xargs dirname)
    $INSTALLED_HICPRO_VERSION/HiC-Pro -h

%environment
    export PATH=$PATH:/usr/local/anaconda/bin
    INSTALLED_VERSION=$(find /usr/local/bin -name HiC-Pro | xargs dirname)
    export PATH=$PATH:$INSTALLED_VERSION 
    export LANG=C
