# pyLeiden

A CLI tool for clustering with the [Leiden](https://www.nature.com/articles/s41598-019-41695-z) algorithm.

## Installation

```bash
# pixi
pixi global -c conda-forge install pyleiden
# pip:
pip install pyleiden
# pipx:
pipx install pyleiden
# conda:
conda install -c conda-forge pyleiden
# mamba:
mamba install -c conda-forge pyleiden
```

## Options

```
usage: pyleiden [-h] [-v] [-o {CPM,modularity}] [-r RESOLUTION_PARAMETER]
                [-b BETA] [-n NUMBER_ITERATIONS] [-s SEED] [-q]
                INPUT OUTPUT

pyLeiden: a CLI tool for clustering with the Leiden algorithm.

positional arguments:
  INPUT                 Tabular file containing the edges of the network.
                        The first two columns should be elelemnts that
                        are connected via an edge in the graph. A third
                        column may be provided with numerical values that
                        represent the weight of the edge.
  OUTPUT                Cluster membership. Each line contains all the
                        members of a given cluster separated by a tab.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -o {CPM,modularity}, --objective-function {CPM,modularity}
                        Whether to use the Constant Potts Model (CPM) or
                        modularity. (default: CPM)
  -r RESOLUTION_PARAMETER, --resolution-parameter RESOLUTION_PARAMETER
                        The resolution parameter to use. Higher
                        resolutions lead to more smaller clusters, while
                        lower resolutions lead to fewer larger clusters.
                        (default: 1.0)
  -b BETA, --beta BETA  Parameter affecting the randomness in the Leiden
                        algorithm. This affects only the refinement step
                        of the algorithm. (default: 0.01)
  -n NUMBER_ITERATIONS, --number-iterations NUMBER_ITERATIONS
                        The number of iterations to iterate the Leiden
                        algorithm. Each iteration may improve the
                        partition further. Using a negative number of
                        iterations will run until a stable iteration is
                        encountered (i.e. the quality was not increased
                        during that iteration). (default: 2)
  -s SEED, --seed SEED  Seed for the random number generator. (default:
                        1953)
  -q, --quiet           Run quietly (will not print output to stout).
                        (default: False)
```

## Example

To showcase the pyLeiden's functionality, we will cluster plasmid sequences retrieved from the [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb/) database into PTUs. These PTUs correspond to distinct groups of plasmids sharing a minimum of 95% average nucleotide identity (ANI) and 85% aligned fraction (AF) among their members.

First we download the plasmid DNA sequences from PLSDB and use [SeqKit](https://github.com/shenwei356/seqkit) to format the FASTA file by leaving only the sequence accessions in the headers:

```bash
curl -L https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta \
    | seqkit seq --remove-gaps --only-id --upper-case --min-len 2500 \
    > plsdb.fna
```

Next, we will estimate the ANI and AF between pairs of plasmids with [skani](https://github.com/bluenote-1577/skani):

```bash
skani triangle -t 16 --diagonal --sparse -i -m 150 -c 30 -s 70 plsdb.fna > skani_output.tsv
```

skani's output will look like this:

```
Ref_file    Query_file   ANI     Align_fraction_ref   Align_fraction_query   Ref_name        Query_name
---------   ----------   -----   ------------------   --------------------   -------------   -------------
plsdb.fna   plsdb.fna    88.13   3.57                 22.72                  NZ_CP047963.1   NZ_LR890753.1
plsdb.fna   plsdb.fna    90.91   1.92                 17.10                  NZ_CP047963.1   NZ_CP135690.1
plsdb.fna   plsdb.fna    97.33   1.03                 18.07                  NZ_CP047963.1   NZ_MN477222.1
plsdb.fna   plsdb.fna    98.82   4.16                 30.04                  NZ_CP047963.1   OW969852.1
plsdb.fna   plsdb.fna    93.74   2.57                 28.76                  NZ_CP047963.1   NZ_CP006919.1
plsdb.fna   plsdb.fna    90.89   0.57                 28.74                  NZ_CP047963.1   NZ_MK973082.1
plsdb.fna   plsdb.fna    94.94   3.04                 30.73                  NZ_CP047963.1   JX843238.1
plsdb.fna   plsdb.fna    98.06   1.72                 31.86                  NZ_CP047963.1   CP067766.1
plsdb.fna   plsdb.fna    89.23   2.02                 19.52                  NZ_CP047963.1   NZ_CP041175.1
```

To use pyLeiden, a tabular file listing pairwise relationships between plasmids (edges) is required as input. This file should have a minimum of two columns indicating the connected plasmids, but a third column indicating the strength of the connection is also accepted. This allows the similarity levels of pairs of plasmids to be considered during clustering.

To generate the input file for pyLeiden from skani's output we can use `awk`:

```bash
awk 'NR>1 && $3>=95 && ($4>=85 || $5>=85) {
    printf("%s\t%s\t%.4f\n", $6, $7, $3 * ($4 > $5 ? $4 : $5) / 10000)
}' skani_output.tsv > edges.tsv
```

This command will:

1. Remove the header from skani's output.
2. Filter out edges that do not meet the minimum ANI or AF criteria (ANI ≥ 95% and AF ≥ 85%).
3. Compute edge weights using the formula: $Weight = \frac{ANI \times \max(AF_{query}, AF_{target})}{10000}$.

The `edges.tsv` file will look like this:

```
NZ_CP047963.1   NZ_CP072419.1   0.8990
NZ_CP047963.1   NZ_CP073001.1   0.9818
NZ_CP047963.1   NZ_CP107371.1   0.9533
NZ_CP047963.1   NZ_CP079649.1   0.8703
NZ_CP047963.1   NZ_CP079814.1   0.9695
NZ_CP047963.1   NZ_CP065442.1   0.9822
NZ_CP047963.1   NZ_CP086290.1   1.0000
NZ_CP047963.1   NZ_MN401418.1   0.9967
NZ_CP047963.1   NZ_CP109990.1   0.9833
NZ_CP079733.1   NZ_CP100083.1   0.9684
```

Finally, we can run pyLeiden to cluster the plasmids into PTUs:

```bash
pyleiden -r 0.9 edges.tsv clusters.txt

    [1/3] Reading input file and building the graph.
    [2/3] Clustering nodes.
    [3/3] Writing output file.
    Total number of nodes: 57,338
    Total number of clusters: 24,579
```

In the output file (`clusters.txt`), plasmids that were clustered together will be in the same line and separated by tab characters.

```
NZ_CP075726.1	NZ_CP074015.1	NZ_CP063491.1	CP134921.1	NZ_CP101870.1	CP099161.1	NZ_CP047474.1	NZ_CP091475.1
NZ_CP114321.1	NZ_CP114165.1	NZ_CP114169.1	NZ_CP114355.1	NZ_CP114362.1	NZ_CP114333.1	NZ_CP114317.1	NZ_CP114329.1
NZ_CP072485.1	NZ_CP114643.1	NZ_CP114696.1	NZ_CP114682.1	NZ_CP114669.1	NZ_CP114744.1	NZ_CP114656.1	NZ_CP114629.1
NZ_CP072487.1	NZ_CP114685.1	NZ_CP114646.1	NZ_CP114699.1	NZ_CP114659.1	NZ_CP114672.1	NZ_CP114633.1	NZ_CP114747.1
NZ_CP072484.1	NZ_CP114668.1	NZ_CP114681.1	NZ_CP114628.1	NZ_CP114743.1	NZ_CP114642.1	NZ_CP114655.1	NZ_CP114695.1
AP028565.1	NZ_AP018455.1	AP028551.1	NZ_CP081885.1	NZ_AP018454.1	LC602700.1	LC635857.1
CP103692.1	NZ_CP135467.1	CP109717.1	NZ_CP103502.1	NZ_CP022826.1	CP020530.1	NZ_CP114357.1
NZ_CP053472.1	NZ_CP033100.1	NZ_MK046687.1	CP052949.1	NZ_CP069220.1	NZ_CP061042.1	NC_001384.1
OW848857.1	OW848811.1	NZ_CP074546.1	NZ_CP117754.1	NC_019083.1	NZ_AP022436.1	NZ_MT039163.1
CP073291.1	NZ_CP102673.1	NZ_CP102676.1	NZ_CP098215.1	NZ_CP067300.1	NZ_CP067304.1	NZ_AP021934.1
```

You can control the granularity of the clustering using `--resolution-parameter`:

```bash
pyleiden --resolution-parameter 0.1 edges.tsv clusters.txt

    [1/3] Reading input file and building the graph.
    [2/3] Clustering nodes.
    [3/3] Writing output file.
    Total number of nodes: 57,338
    Total number of clusters: 18,867
```

When decreasing the `--resolution-parameter` from `0.9` to `0.1`, the number of clusters decreased from 24,579 to 18,867. In practice, lower values of `--resolution-parameter` will result in larger clusters with more loosely connected members.