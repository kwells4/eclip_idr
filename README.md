# IDR pipeline

A snakemake pipeline to run idr to test the replicability of two eCLIP samples and to find a conservative peak set.

Follows the framework developed by ENCODE and uses the idr pipeline from the [kundaje lab at Stanford](https://github.com/kundajelab/idr)

This also uses clipper. I recommend running from [this repo](https://github.com/kwells4/clipper) as it fixes an issue with temporary files being generated and not deleted. This is especially necessary if you are running multiple clipper runs in parallel as this repo does.

Additionally, this repo contains a python version of the peak input normalization that is more flexible to experimental design than the original perl version of the script.