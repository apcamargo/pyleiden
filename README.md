# pyLeiden

A CLI tool for clustering with the Leiden algorithm.

## Installation

```bash
pip install pyleiden
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
curl -L https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plsdb.fna.bz2 \
    | seqkit seq --only-id --upper-case \
    > plsdb.fna
```

Next, we will estimate the ANI and AF between pairs of plasmids with [skani](https://github.com/bluenote-1577/skani):

```bash
skani dist -t 16 -c 80 -s 70 --qi plsdb.fna --ri plsdb.fna -o skani_output.tsv
```

skani's output will look like this:

```
Ref_file    Query_file   ANI      Align_fraction_ref   Align_fraction_query   Ref_name        Query_name
---------   ----------   ------   ------------------   --------------------   -------------   ----------
plsdb.fna   plsdb.fna    100.00   97.78                97.78                  AB063523.1      AB063523.1
plsdb.fna   plsdb.fna    100.00   99.60                99.60                  AB244976.1      AB244976.1
plsdb.fna   plsdb.fna    99.85    73.12                97.91                  NZ_CP005192.1   AB244976.1
plsdb.fna   plsdb.fna    99.43    17.38                23.08                  NZ_CP005087.1   AB244976.1
plsdb.fna   plsdb.fna    97.62    30.27                15.22                  NC_020563.2     AB244976.1
plsdb.fna   plsdb.fna    97.35    3.90                 17.27                  NZ_CP005190.1   AB244976.1
plsdb.fna   plsdb.fna    96.86    82.16                8.60                   NZ_CP060130.1   AB244976.1
plsdb.fna   plsdb.fna    95.74    5.37                 15.52                  NC_014007.1     AB244976.1
plsdb.fna   plsdb.fna    94.39    48.82                46.12                  NZ_CP009296.1   AB244976.1
```

To use pyLeiden, a tabular file listing pairwise relationships between plasmids (edges) is required as input. This file should have a minimum of two columns indicating the connected plasmids, but a third column indicating the strength of the connection is also accepted. This allows the similarity levels of pairs of plasmids to be considered during clustering.

To generate the input file for pyLeiden from skani's output we can use `awk`:

```bash
awk 'NR>1 && $3>=95 && ($4>=85 || $5>=85) {
    printf("%s\t%s\t%.4f\n", $6, $7, $3 * ($4 > $5 ? $4 : $5) / 10000)
}' skani_output.tsv > edges.tsv
```

This command performs three actions: (1) removes the header that was present in skani's output; (2) filters out edges that do not meet the minimum ANI or AF criteria (ANI ≥ 95% and AF ≥ 85%); (3) computes edge weights using the following formula: $Weight = \frac{ANI \times \max(AF_{query}, AF_{target})}{10000}$.

The `edges.tsv` file will look like this:

```
AB063523.1      AB063523.1   1.0000
AB244976.1      AB244976.1   1.0000
NZ_CP005192.1   AB244976.1   0.9796
NZ_CP060130.1   AB244976.1   0.8365
AB576781.2      AB576781.2   0.9993
JF274991.1      AB576781.2   0.9993
JF274992.1      AB576781.2   0.9999
NZ_CP040565.1   AB576781.2   0.8565
NZ_CP039596.1   AB576781.2   0.8512
NZ_CP043905.1   AB576781.2   0.8544
```

Finally, we can run pyLeiden to cluster the plasmids into PTUs:

```bash
pyleiden edges.tsv clusters.txt

    [1/3] Reading input file and building the graph.
    [2/3] Clustering nodes.
    [3/3] Writing output file.
    Total number of nodes: 34,511
    Total number of clusters: 14,794
```

In the output file (`clusters.txt`), plasmids that were clustered together will be in the same line and separated by tab characters.

```
NZ_CP042577.1	NZ_CP073307.1	NZ_CP056765.1	NZ_CP047764.1	NZ_CP058117.1	NZ_CP068289.1	NZ_CP056154.1
NZ_CP056221.1	NZ_CP056383.1	NZ_CP056363.1	NZ_CP056898.1	NZ_CP056900.1	NZ_CP056353.1	NZ_CP056825.1
NZ_CP073733.1	NZ_CP073724.1	NZ_CP073739.1	NZ_CP073727.1	NZ_CP073730.1	NZ_CP073742.1	NZ_CP073736.1
NZ_LC586266.1	NZ_LC586264.1	NZ_LC586263.1	NZ_LC586268.1	NZ_LC586265.1	NZ_LC586267.2	NZ_LC586262.1
AP018702.1	NZ_CP065814.1	NZ_CP020461.1	NZ_MF383377.1	NZ_LT960792.1	NZ_CP075332.1
NZ_CP075287.1	NZ_CP064176.1	NZ_CP039605.1	NZ_CP070552.1	NZ_MN543574.1	NZ_CP047338.1
AP022387.1	NZ_AP022395.1	NZ_AP022400.1	NZ_AP022514.1	NZ_AP022495.1	NZ_KY014464.1
AP022442.1	NZ_CP067387.1	NZ_AP024180.1	NZ_CP055455.1	LR890443.1	LR890483.1
NZ_CP041424.1	NZ_CP032799.1	NZ_CP041430.1	NZ_CP041434.1	NZ_CP041432.1	CP049610.1
CP049186.1	CP049172.1	CP045525.2	CP049170.1	NZ_CP042338.1	NZ_MW646303.1
```

You can control the granularity of the clustering using `--resolution-parameter`:

```bash
pyleiden --resolution-parameter 0.2 edges.tsv clusters.txt

    [1/3] Reading input file and building the graph.
    [2/3] Clustering nodes.
    [3/3] Writing output file.
    Total number of nodes: 34,511
    Total number of clusters: 13,135
```

When decreasing the `--resolution-parameter` from the default value of `1.0` to `0.2`, the number of clusters decreased from 14,794 to 13,135. In practice, lower values of `--resolution-parameter` will result in larger clusters whose members are more loosely connected.