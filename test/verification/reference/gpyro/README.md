# Gpyro reference data

Every Gpyro verification comparison is a matched native-input run of V2
revision `bff0126c68f2a35c084ef04490ae327185f281a3`, retained with its inputs,
native outputs, refinements and cost in [`../matched/`](../matched/README.md).
This directory keeps one further Gpyro V2 output, read by the burnout
gradient gates.

## vbo01 — V-BO-01 (MaCFP PMMA gasification, 50 kW/m²)

`vbo01_gpyro_v2_pmma.csv` is the **Gpyro V2** (
<https://github.com/Gedeben/GPYRO-V2>) run made here on the authors' MaCFP
PMMA deck (`PMMA.data` in that repository: 5.9 mm PMMA on 12 mm Kaowool,
50 kW/m², power-law k(T) and c(T) fits), unmodified. Gpyro V2 was built from
that repository without a system gfortran:

```
docker run --rm --user $(id -u):$(id -g) -e HOME=/tmp -e GPYRO_FCOMPL_GNU=gfortran \
  -v $PWD:/gp -w /gp/build/linux gcc:14 \
  bash -c "mkdir -p gpyro && cd gpyro && make -j8 -f ../../Makefile_gpyro gnu_linux"
```

and run from the deck's directory with `OMP_NUM_THREADS=1 gpyro PMMA.data`. The
first header row names the columns; column 2 (`001_MLR`) is the mass loss rate
in **g m⁻² s⁻¹**, column 6 the remaining thickness in m, columns 3–5 the
temperatures (°C) at the surface, 12 mm and 17.9 mm. This file is not a
matched V-BO-01 dataset: native gas transport, burnout and surface
closures differ, and escaped mass loss is not our gas-generation observable.
The retained file does not identify the source commit
or executable hash used for this run, so it must not be attributed to a
specific V2 revision. The source inspected for the comparison audit is
`bff0126c68f2a35c084ef04490ae327185f281a3`; that inspection does not establish
that this revision produced the CSV.
