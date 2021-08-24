#!/bin/bash
echo Set your workding directory. This is the directory where all your data lives and where all files/code will be saved. \
Please enter full path.
read dir

docker run -d -p 8787:8787 -e PASSWORD=@hackathon2021 -v $dir:/home drbueno/rstudio-single-cell