# Script to run rivet docker container

# Create new rivet container
#docker run -it --volume /mnt/d/PhD/Work/Majorana/signal_simulation:/work  rivet_2.7.2 bash

docker run -it -v $PWD:$PWD -w $PWD -u `id -u $USER`:`id -g` hepstore/rivet:2.7.2 bash 

# Launch existing rivet container
#docker start c0d3df4e8baa
#docker attach c0d3df4e8baa