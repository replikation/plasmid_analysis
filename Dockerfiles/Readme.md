Dockerfiles
===

* These are the custom `Dockerfiles` for this workflow
* they should be available via `docker pull` and are downloaded automatically via nextflow
* alternatively you can build them locally via:

```
cd psi_docker && docker build -t replikation/psi-cd-hit:latest .
cd abri_docker && docker build -t replikation/abricate_nextflow:latest .
```

# Repositories

Docker images are created from the following repositories:

+ [**cd-hit**](https://github.com/weizhongli/cdhit)
+ [**abricate**](https://github.com/tseemann/abricate)
+ [**prokka**](https://github.com/tseemann/prokka/)

**Cite these tools if you use them**
