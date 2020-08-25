FROM gcc:latest

# install necessary packages for compliling
RUN apt-get update && \
        apt-get -y install libhdf5-serial-dev libgsl-dev libcgal-dev
