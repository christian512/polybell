FROM ubuntu:20.04
RUN apt-get update
ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install some packages
RUN apt-get install -y git wget make

# Install Python3
RUN apt-get install -y python3 python3-pip
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install the Python Framework
WORKDIR /home
RUN git clone https://github.com/christian512/polybell.git
WORKDIR /home/polybell
RUN python -m pip install -r requirements.txt
RUN python -m pip install -e .
RUN python -m pip install -i https://pypi.gurobi.com gurobipy

# Install GAP
WORKDIR /home
RUN apt-get install libgmp-dev -y
RUN wget https://github.com/gap-system/gap/releases/download/v4.11.1/gap-4.11.1.tar.gz
RUN tar -xf gap-4.11.1.tar.gz
RUn rm gap-4.11.1.tar.gz
WORKDIR /home/gap-4.11.1
RUN ./configure
RUN make
WORKDIR /home/gap-4.11.1/pkg
RUN ../bin/BuildPackages.sh
RUN ln -s /home/gap-4.11.1/bin/gap.sh /usr/bin/gap
RUN ln -s /home/gap-4.11.1/bin/gap.sh /usr/bin/gap.sh

# # Install RANDA
WORKDIR /home
RUN apt-get install build-essential cmake -y
RUN git clone https://github.com/christian512/randa.git
RUN mkdir /home/randa/build
WORKDIR /home/randa/build
RUN cmake ..
RUN make
RUN make install
WORKDIR /home
CMD /bin/bash

# Install GUROBI

