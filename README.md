LSC
===

LSC is a long read error correction tool developed by Kin Fai Au and Pegah T Afshar. It offers fast correction with high sensitivity and good accuracy.

### How to run from Docker

#### Running the included example data:

```bash
$ docker run -v $(pwd)/output:/Source/LSC/example/output \
             -v $(pwd)/temp:/Source/LSC/example/temp \
             -w /Source/LSC/example \
             vacation/lsc runLSC.py --long_reads LR.fa \
             --short_reads SR.fa --specific_tempdir temp --output output
```

#### Running the command on your data:

Construct prepare volumes that point to your inputs if everything is staged in your working directory it could be run as:

```bash
$ docker run -v $(pwd):/home vacation/lsc runLSC.py --long_reads LR.fa \
             --short_reads SR.fa --specific_tempdir temp --output output
```

The use of `--specific_tempdir` is pretty import because this folder gets large and it will fill the docker if you don't have it mounted.
