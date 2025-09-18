# AlphaGenome evaluates the epigenomic impact of non-coding variants.

To get started, first clone the repository and pull the corresponding Docker image:

```git clone https://github.com/DUAN-GAO/alphagenome.git```

```docker pull duangao/alphagenome```

Example

Once the repository and Docker image are prepared, the pipeline can be executed with the following command:

```docker run --rm -it -v $(pwd):/workspace -w /home duangao/alphagenome python main.py --rsid rs5934683```

-v $(pwd):/workspace mounts the current working directory into the container.

-w /home sets the working directory inside the container.

--rsid specifies the variant identifier to be processed.

