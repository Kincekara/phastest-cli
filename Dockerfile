# adapted from https://phastest.ca/download_file/phastest-docker

FROM ubuntu:jammy AS builder

ARG DEBIAN_FRONTEND=noninteractive

# Install base packages
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    locales \
    wget

RUN apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    perl \
    gcc \
    git \
    gnupg \
    make \
    munge \
    curl \
    libpam0g-dev \
    psmisc \
    bash-completion \
    net-tools \
    iputils-ping \
    cpanminus

# Install libraries that BioPerl dependencies depend on
RUN apt-get install --yes \
    libdb-dev \
    zlib1g-dev \
    graphviz

# Install perl modules 
RUN cpanm --notest \
    CPAN::Meta \
    YAML \
    Digest::SHA \
    Module::Build \
    Test::Most \
    Test::Weaken \
    Test::Memory::Cycle \
    Clone

# Install perl modules for network and SSL (and their dependencies)
RUN apt-get install --yes \
    openssl \
    libssl-dev

RUN cpanm --notest \
    IO::Socket::INET6 \
    LWP \
    LWP::Protocol::https

# Install packages for XML processing
RUN apt-get install --yes \
    expat \
    libexpat-dev \
    libxml2-dev \
    libxslt1-dev \
    libgdbm-dev

RUN cpanm --notest \
    XML::Parser \
    XML::Parser::PerlSAX \
    XML::DOM \
    XML::DOM::XPath \
    XML::SAX \
    XML::SAX::Writer \
    XML::Simple \
    XML::Twig \
    XML::Writer \
    XML::LibXML

# Install what counts as BioPerl dependencies
RUN cpanm --notest \
    HTML::TableExtract \
    Algorithm::Munkres \
    Array::Compare \
    Convert::Binary::C \
    Error \
    Graph \
    GraphViz \
    Inline::C \
    PostScript::TextBlock \
    Set::Scalar \
    Sort::Naturally \
    Math::Random \
    Spreadsheet::ParseExcel \
    IO::String \
    JSON \
    Data::Stag

RUN cpanm --notest \
    Bio::DB::GenBank \
    Bio::Perl \
    JSON::XS \
    LWP::Simple

# Install database connectivity packages
RUN apt-get install --yes \
    libdbi-perl \
    libdbd-mysql-perl \
    libdbd-pg-perl \
    libdbd-sqlite3-perl

RUN cpanm --notest \
    DB_File

# Install GD and other graphics dependencies
RUN apt-get install --yes \
    libgd-dev

RUN cpanm --notest \
    GD \
    SVG \
    SVG::Graph

## App ##
FROM ubuntu:jammy AS app

# Set PHASTEST environment
ENV PHASTEST_HOME="/opt/phastest-app"

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="phastest-cli"
LABEL software.version="0.1"
LABEL description="CLI version of PHASTEST"
LABEL website="https://github.com/Kincekara/phastest-cli"
LABEL license="https://github.com/Kincekara/phastest-cli/blob/master/LICENSE"
LABEL maintainer="Kutluhan Incekara"
LABEL maintainer.email="kutluhan.incekara@ct.gov"

# Set up PHASTEST dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    wget \
    ca-certificates \
    perl \
    libdbi-perl \
    libdb-dev \
    zlib1g-dev \
    libssl-dev \
    expat \
    libexpat-dev \
    libxml2-dev \
    libxslt1-dev \
    libgdbm-dev \
    libdbi-perl \
    libdbd-mysql-perl \
    libdbd-pg-perl \
    libdbd-sqlite3-perl \
    infernal \
    trnascan-se \
    ncbi-blast+ \
    hmmer \
    barrnap \
    prodigal \
    diamond-aligner \
    aragorn \
    graphviz &&\
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Copy over the phastest scripts and dependencies
COPY --from=builder /usr/local/ /usr/local/
COPY ./phastest-app/ ${PHASTEST_HOME}/
COPY ./phastest.sh /usr/local/bin/phastest
RUN chmod +x /usr/local/bin/phastest

# Setup PHASTEST databases
RUN makeblastdb -in ${PHASTEST_HOME}/DB/prophage_virus.db -dbtype prot &&\
    # makeblastdb -in ${PHASTEST_HOME}/DB/swissprot.db -dbtype prot &&\
    diamond makedb --in ${PHASTEST_HOME}/DB/swissprot.db -d ${PHASTEST_HOME}/DB/swissprot.dmnd &&\
    rm ${PHASTEST_HOME}/DB/swissprot.db

ENV PHASTEST_CLUSTER_HOME=$PHASTEST_HOME \
    PATH=$PATH:"${PHASTEST_HOME}/sub_programs/FragGeneScan1.20" \
    LC_ALL=C

CMD [ "phastest", "help" ]

WORKDIR /data