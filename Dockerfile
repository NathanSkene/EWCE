# base image: rocker/verse (with a specific version of R)
#   has R, RStudio, tidyverse, devtools, tex, and publishing-related packages
FROM rocker/verse:3.6.0

# required
MAINTAINER Nathan Skene <nskene@imperial.ac.uk>

# copy the repo contents into the docker image at `/ewce`
COPY . /ewce

RUN apt-get update &&\
    apt-get install -y --no-install-recommends\
    libxml2-dev libcurl4-openssl-dev libssl-dev\
    libssh2-1-dev

# install the dependencies of the R package located at `/ewce`
RUN R -e "devtools::install_dev_deps('/ewce', dep = TRUE, quiet=TRUE)"

# install the package
RUN R -e "devtools::install_github('nathanskene/ewce',quiet=TRUE)"
