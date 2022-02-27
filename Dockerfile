FROM ubuntu:20.04
RUN apt-get update
ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install the Python Framework
COPY . /linearbell
CMD cd linearbell
CMD python -m pip install -r requirements.txt
CMD python -m pip install -e
CMD python -m pip install -i https://pypi.gurobi.com gurobipy
CMD cd ..

# Install GAP
RUN apt-get install libgmp-dev -y
CMD wget https://github.com/gap-system/gap/releases/download/v4.11.1/gap-4.11.1.tar.gz
CMD tar -xf gap-4.11.1.tar.gz
CMD cd gap-4.11.1
CMD ./configure
CMD make
CMD cd pkg
CMD ../bin/BuildPackages.sh
CMD cd ../../

# TODO: SETUP PATH


# Install RANDA
RUN apt-get install build-essential cmake -y
CMD git clone https://github.com/christian512/randa.git
CMD cd randa
CMD mkdir build
CMD cd build
CMD cmake ..
CMD make
CMD make install

# Install GUROBI

