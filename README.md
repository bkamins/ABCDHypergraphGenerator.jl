# ABCDHypergraphGenerator.jl

Artificial Benchmark for Hypergraphs Community Detection (ABCDH) - A Random Hypergraph Model with Community Structure

Bogumił Kamiński, Paweł Prałat, François Théberge

### Usage

The project is prepared to be run from command line interface.

Installation instructions:
* install Julia and make sure it is in search path
* clone the project repository
* change directory to `/utils` and run `abcdh.jl` file stored there following
  the example instructions given below.

Here is an example session output (assuming you are in `/utils` folder).

```
$ julia --project abcdh.jl --help
julia --project abcdh.jl -n N -d D -c C -x X -q Q -w W -s S [-o O] [-m] [-h] [--stats]

optional arguments:
  -o O        prefix for output file names; the generated file names
              are [prefix]_deg.txt for degree sequence,
              [prefix]_comm.txt for community size sequence,
              [prefix]_assign.txt for assignment of vertices to
              communities, [prefix]_he.txt for hyperedges
  -m          if this flag passed a multi-hypergraph is generated; by
              default a simple hypergraph is generated
  --stats     if this flag passed print generated hypergraph
              statistics
  -h, --help  show this help message and exit

required arguments:
  -n N        number of vertices
  -d D        either a tuple γ,δ,D or file with degree sequence
              (filename cannot contain two commas)
  -c C        either a tuple β,s,S or file with community size
              sequence (filename cannot contain two commas)
  -x X        mixing parameter ξ
  -q Q        either a sequence q₁,q₂,...,qₖ of weights of hyperedges
              of sizes from 1 to k or a file name with such a sequence
              in the same format (if a single number is passed it is
              considered to be a filename)
  -w W        either one of values: ':pure', ':linear', ':equal' or a
              file name with weights for w_cd where in i-th line are
              comma separated weights for hyperedge of size i+2
              ranging from floor(i+2,2)+1 to i; the :pure value
              assumes that only c=d weight is non zero, the :linear
              value assumes weight equal to c, the :equal value
              assumes all weights are 1
  -s S        seed value for generator

$ julia --project abcdh.jl -n 100000 -d 2.5,5,100 -c 1.5,1000,10000 -x 0.5 -q 0.0,0.4,0.3,0.2,0.1 -w :linear -s 1234 -o result --stats
[ Info: failed to sample an admissible community sequence in 1000 draws. Fixing
[ Info: distribution of hyperedge proportions does not add up to 1. Fixing.
[ Info: Moving 1 stumps from community 2 to background graph
[ Info: Moving 1 stumps from community 4 to background graph
[ Info: Moving 1 stumps from community 10 to background graph
[ Info: Moving 1 stumps from community 12 to background graph
[ Info: Moving 1 stumps from community 14 to background graph
[ Info: Moving 1 stumps from community 20 to background graph
[ Info: Moving 1 stumps from community 22 to background graph
[ Info: Moving 1 stumps from community 23 to background graph
[ Info: Moving 1 stumps from community 24 to background graph
[ Info: Moving 1 stumps from community 26 to background graph
[ Info: Moving 1 stumps from community 29 to background graph
[ Info: Moving 1 stumps from community 35 to background graph
[ Info: Moving 1 stumps from community 36 to background graph
[ Info: fixing 291 hyperedges that were multisets
[ Info: fixing 162 hyperedges that were duplicated

[ Info: Degrees
generated degree distribution
Summary Stats:
Length:         100000
Missing Count:  0
Mean:           10.780250
Minimum:        5.000000
1st Quartile:   5.000000
Median:         7.000000
3rd Quartile:   11.000000
Maximum:        100.000000
Type:           Int64


[ Info: Communities
generated community size distribution
Summary Stats:
Length:         38
Missing Count:  0
Mean:           2631.578947
Minimum:        1017.000000
1st Quartile:   1242.000000
Median:         2176.000000
3rd Quartile:   3729.750000
Maximum:        8947.000000
Type:           Int64


[ Info: Hyperedges
5×3 DataFrame
 Row │ x      q        vol
     │ Int64  Float64  Float64?
─────┼─────────────────────────────────
   1 │     1      0.0  missing
   2 │     2      0.4        0.400119
   3 │     3      0.3        0.300013
   4 │     4      0.2        0.199955
   5 │     5      0.1        0.0999142
allunique(hyperedges) = true
all(allunique, hyperedges) = true
(ξ, ξ_emp) = (0.5, 0.4749853975989391)

4×2 DataFrame
 Row │ he_size  ξ_emp
     │ Int64    Float64
─────┼───────────────────
   1 │       2  0.481998
   2 │       3  0.445908
   3 │       4  0.498098
   4 │       5  0.49248
[ Info: w
5×5 Matrix{Float64}:
 1.0  0.0  0.0  0.0       0.0
 0.0  1.0  0.4  0.0       0.0
 0.0  0.0  0.6  0.428571  0.25
 0.0  0.0  0.0  0.571429  0.333333
 0.0  0.0  0.0  0.0       0.416667
[ Info: w_emp
5×5 Matrix{Float64}:
 NaN  0.0  0.0       0.0       0.0
 NaN  1.0  0.444781  0.0       0.0
 NaN  0.0  0.555219  0.417643  0.246044
 NaN  0.0  0.0       0.582357  0.336687
 NaN  0.0  0.0       0.0       0.417269

 julia --project abcdh.jl -n 100 -d 2.5,5,20 -c 1.5,10,30 -x 0.3 -q 0.0,0.4,0.3,0.2,0.1 -w :linear -s 1234 --stats
[ Info: distribution of hyperedge proportions does not add up to 1. Fixing.
[ Info: Moving 1 stumps from community 2 to background graph
[ Info: Moving 1 stumps from community 5 to background graph
[ Info: added degree to the following nodes due to parity issues: [47]
[ Info: fixing 29 hyperedges that were multisets
[ Info: fixing 11 hyperedges that were duplicated
[ Info: skipping saving generated graph

[ Info: Degrees
generated degree distribution
Summary Stats:
Length:         100
Missing Count:  0
Mean:           8.400000
Minimum:        5.000000
1st Quartile:   6.000000
Median:         7.000000
3rd Quartile:   9.250000
Maximum:        20.000000
Type:           Int64
deviations from wanted degrees
1×3 DataFrame
 Row │ degs_wanted  degs_generated  nrow
     │ Int64        Int64           Int64
─────┼────────────────────────────────────
   1 │           7               8      1


[ Info: Communities
generated community size distribution
Summary Stats:
Length:         5
Missing Count:  0
Mean:           20.000000
Minimum:        11.000000
1st Quartile:   16.000000
Median:         17.000000
3rd Quartile:   26.000000
Maximum:        30.000000
Type:           Int64


[ Info: Hyperedges
5×3 DataFrame
 Row │ x      q        vol
     │ Int64  Float64  Float64?
─────┼─────────────────────────────────
   1 │     1      0.0  missing
   2 │     2      0.4        0.42381
   3 │     3      0.3        0.307143
   4 │     4      0.2        0.185714
   5 │     5      0.1        0.0833333
allunique(hyperedges) = true
all(allunique, hyperedges) = true
(ξ, ξ_emp) = (0.3, 0.2870662460567823)

4×2 DataFrame
 Row │ he_size  ξ_emp
     │ Int64    Float64
─────┼───────────────────
   1 │       2  0.353933
   2 │       3  0.139535
   3 │       4  0.282051
   4 │       5  0.357143
[ Info: w
5×5 Matrix{Float64}:
 1.0  0.0  0.0  0.0       0.0
 0.0  1.0  0.4  0.0       0.0
 0.0  0.0  0.6  0.428571  0.25
 0.0  0.0  0.0  0.571429  0.333333
 0.0  0.0  0.0  0.0       0.416667
[ Info: w_emp
5×5 Matrix{Float64}:
 NaN  0.0  0.0       0.0       0.0
 NaN  1.0  0.635135  0.0       0.0
 NaN  0.0  0.364865  0.464286  0.444444
 NaN  0.0  0.0       0.535714  0.222222
 NaN  0.0  0.0       0.0       0.333333

$ julia --project abcdh.jl -n 100 -d 2.5,5,20 -c 1.5,10,30 -x 0.3 -q 0.0,0.4,0.3,0.2,0.1 -w :linear -s 1234 --stats -m
[ Info: distribution of hyperedge proportions does not add up to 1. Fixing.
[ Info: Moving 1 stumps from community 2 to background graph
[ Info: Moving 1 stumps from community 5 to background graph
[ Info: added degree to the following nodes due to parity issues: [47]
[ Info: skipping saving generated graph

[ Info: Degrees
generated degree distribution
Summary Stats:
Length:         100
Missing Count:  0
Mean:           8.400000
Minimum:        5.000000
1st Quartile:   6.000000
Median:         7.000000
3rd Quartile:   9.250000
Maximum:        20.000000
Type:           Int64
deviations from wanted degrees
1×3 DataFrame
 Row │ degs_wanted  degs_generated  nrow
     │ Int64        Int64           Int64
─────┼────────────────────────────────────
   1 │           7               8      1


[ Info: Communities
generated community size distribution
Summary Stats:
Length:         5
Missing Count:  0
Mean:           20.000000
Minimum:        11.000000
1st Quartile:   16.000000
Median:         17.000000
3rd Quartile:   26.000000
Maximum:        30.000000
Type:           Int64


[ Info: Hyperedges
5×3 DataFrame
 Row │ x      q        vol
     │ Int64  Float64  Float64?
─────┼─────────────────────────────────
   1 │     1      0.0  missing
   2 │     2      0.4        0.42381
   3 │     3      0.3        0.307143
   4 │     4      0.2        0.185714
   5 │     5      0.1        0.0833333
allunique(hyperedges) = false
all(allunique, hyperedges) = false
(ξ, ξ_emp) = (0.3, 0.19873817034700317)

4×2 DataFrame
 Row │ he_size  ξ_emp
     │ Int64    Float64
─────┼────────────────────
   1 │       2  0.241573
   2 │       3  0.104651
   3 │       4  0.25641
   4 │       5  0.0714286
[ Info: w
5×5 Matrix{Float64}:
 1.0  0.0  0.0  0.0       0.0
 0.0  1.0  0.4  0.0       0.0
 0.0  0.0  0.6  0.428571  0.25
 0.0  0.0  0.0  0.571429  0.333333
 0.0  0.0  0.0  0.0       0.416667
[ Info: w_emp
5×5 Matrix{Float64}:
 NaN  0.0  0.0       0.0       0.0
 NaN  1.0  0.467532  0.0       0.0
 NaN  0.0  0.532468  0.344828  0.461538
 NaN  0.0  0.0       0.655172  0.230769
 NaN  0.0  0.0       0.0       0.307692
```

The session shows you:
1. How to get help.
2. Generating standard hypergraph.
3. Generating a small "tight" hypergraph (that has significant number of initial collisions).
4. Generating a small "tight" non-simple hypergraph (collisions are accepted).

For steps 3 and 4 results are not saved to disk.

Note that statstics generation with `--stats` option can be time consuming for large hypergraphs.
