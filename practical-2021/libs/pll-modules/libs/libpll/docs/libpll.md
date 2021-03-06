Scope
--------------------------------------------------------------------------------

The goal of this documentation is th explain how to interact with `libpll`,
as long as the reader already knows the broad strokes of phylogenetic inference.
As such, there will be some terms which are present in the text, but explaining
the terms completely is outside of the scope of this documentation.

This documentation is _not_ intended to teach phylogenetic inference, HPC, or
programming in general. Basic knowledge of these topics is assumed. Instead,
this will explain how `libpll` does phylogenetic inference, not how it is done
in general.

High level design
--------------------------------------------------------------------------------

At a high level, `libpll` exists to implement efficient versions of core
functions that take up the majority of the runtime during phylogenetic
inference. To do this, `libpll` has a few important data structures which
contain most of the data required for likelihood computation:

- [`pll_utree_t`](pll_utree_t.md)
- [`pll_partition_t`](pll_partition_t.md)
- [`pllmod_treeinfo_t`](pllmod_treeinfo_t.md)

The `pll_utree_t` data structure contains the information that is relevant to
the tree portion of the model, while `pll_partitition_t` contains the other
model parameters, as well as buffers to store intermediate values called
[CLVs][clvs], and information about the state of computation and the machine.

[clvs]: pll_partition_t.md#clv

For most use cases (especially those involving likelihood calculations) of
`libpll`, both a `pll_utree_t` and a `pll_partition_t` will be required.
Information on how to create and interact with this data structures can be found
in their respective pages.

Likelihood Evaluation
--------------------------------------------------------------------------------

`libpll` evaluates the likelihood of a tree using Felsenstien's Algorithm
[felsenstien]. In summary, the algorithm proceeds like so:

1. Pick a virtual root arbitrarily,
2. Perform a post order traversal from the virtual root. 
3. For each node in the traversal compute the current nodes CLV:
    - If a tip, simply return the assigned CLV
    - Otherwise, compute the CLV using the children's CLVs.

The specific implemenation for this algorithm in `libpll` can be described as
this: allocate a single buffer which will contain all the CLVs. Perform the
steps listed above, but during the traversal, don't actually compute any CLVs.
Instead, assign entries into the CLV buffer for each node when the tree is
created, and record the parent and children's CLVs in records. These records
contain the operations required to calculate all the CLVs. The final step is to
compute the operations. This can be summarized as:

1. Assign Entries into the CLV buffer 
1. Pick a virtual root arbitrarily,
2. Perform a post order traversal from the virtual root. 
3. For each node in the traversal:
    - If its a tip, do nothing
    - Otherwise, record the current CLV index, as well as the children's indices,
      into operations.
4. Using the operations created in the last step, compute the CLVs for the tree.

[felsenstien]: J. Felsenstein, ???Evolutionary trees from DNA sequences: A
maximum likelihood approach,??? Journal of Molecular Evolution, vol. 17, no. 6,
pp. 368???376, Nov. 1981.

Core Tasks
--------------------------------------------------------------------------------

In terms of calculations, `libpll` has optimized routines for 4 main tasks:

- Partials
- Derivatives
- Likelihood
- Probability Matrices

In order to calculate the likelihood of a tree from scratch only Probability
Matrices, Partials, and Likelihood routines are required. The Derivatives
portion of calculation is used for optimization of the branch lengths of the
tree.

A sample tree inference cycle might be:

1. Calculate Probability Matrices,
2. Calculate Partials,
3. Calculate Likelihood,
4. Optimize Numeric Model Parameters/Branch Lengths,
5. Propose Change to Tree,
6. Optimize Numeric Model Parameters/Branch Lengths,
7. Repeat from 5 until done.

### Partials

```
PLL_EXPORT void pll_update_partials(pll_partition_t * partition,
                                    const pll_operation_t * operations,
                                    unsigned int count);
```

This function will compute the CLVs in the nodes specified in the `operations`
array. This array is generated by [`pll_create_operations`][1].

[1]: pll_utree_t.md#Notable-Functions

### Derivatives

`libpll` has the capability to calculate first and second derivatives around a
single branch. To do this, first a sumtable needs to be computed. This can be
done with the function:

```
pll_export int pll_update_sumtable(pll_partition_t * partition,
                                      unsigned int parent_clv_index,
                                      unsigned int child_clv_index,
                                      int parent_scaler_index,
                                      int child_scaler_index,
                                      const unsigned int * params_indices,
                                      double *sumtable);
```

Notable Parameters:

- `parent_clv_index`, `child_clv_index`: The CLVs to the nodes representing the
  two endpoints of the branch in question.
- `param_indices`: A list of the indices for each rate category present in the
  partition.
- `sumtable`: Out parameter for the sumtable. Needs to be allocated with size
  `rates * states_padded`.

Once the sumtable is computed, we can use it to repeatedly call the derivatives
function. This evaluates the derivative of the model _at a point_. The intended
use case is to repeatedly evaluate first and second derivatives of one branch in
order to optimize it using NR-optimization.

```
pll_export int pll_compute_likelihood_derivatives(pll_partition_t * partition,
                                                  int parent_scaler_index,
                                                  int child_scaler_index,
                                                  double branch_length,
                                                  const unsigned int * params_indices,
                                                  const double * sumtable,
                                                  double * d_f,
                                                  double * dd_f);
```

Notable Parameters:

- `branch_length`: branch length to evaluate the derivative at.
- `sumtable`: Sumtable from `pll_update_sumtable`.
- `d_f` : Out parameter which will contain the first derivative.
- `dd_f`: Out parameter which will contain the second derivative.


### Likelihood

There are two ways of computing a likelihood: around a root node, or around a
virtual root node. The difference is the number of CLVs involved. In a rooted
tree, the root CLV is basically finished. All that is required is to multiply
the CLV by the frequency, and compute the product.

In an unrooted tree, there is no single CLV that is nearly done, so we need to
compute it. So, in the case of an unrooted tree, an additional step is required
to make the final CLV.

Accordingly, there are two versions of the function to compute the likelihood of
the tree. The first is for trees with a root, and the second is for those
without.

Common arguments between the two functions:

- `partition`: The partition to calculate the likelihood on.
- `freq_indices`: An array of indices that is `sites` long.
- `persite_lnl`: If not null, then the per-site loglikelihoods are placed in
  this buffer.

```
PLL_EXPORT double pll_compute_root_loglikelihood(pll_partition_t * partition,
                                                 unsigned int clv_index,
                                                 int scaler_index,
                                                 const unsigned int * freqs_indices,
                                                 double * persite_lnl);

PLL_EXPORT double pll_compute_edge_loglikelihood(pll_partition_t * partition,
                                                 unsigned int parent_clv_index,
                                                 int parent_scaler_index,
                                                 unsigned int child_clv_index,
                                                 int child_scaler_index,
                                                 unsigned int matrix_index,
                                                 const unsigned int * freqs_indices,
                                                 double * persite_lnl);
```

The only meaningful difference between the arguments of these two functions is
whether they take 1 or 2 CLV/scalar indicies as arguments. In the case of root
logliklihood calculation, only one CLV is needed, so it takes one. In the case
of edge loglikelihood, we are calculating around and edge, so we need 2 CLVs,
and additionally a matrix index.

### Probability Matrix

To update the probability matrices of the `pll_partition_t`, a list of branch
lengths needs to be obtained. The best way to do this is with
[`pll_utree_create_operations`][2].


```
PLL_EXPORT int pll_update_prob_matrices(pll_partition_t * partition,
                                        const unsigned int * params_index,
                                        const unsigned int * matrix_indices,
                                        const double * branch_lengths,
                                        unsigned int count);
```

- `matrix_indices`: A list of indices to put the matrices into.
- `branch_lengths`: A list of the branch lengths to use to produce the
  probability matrices. Each item in this list needs correspond to the same
  thing in the `matrix_indices` list, as this function will put the probability
  matrices in the location indicated by `matrix_indices`. The upshot of this is
  to just use [`pll_utree_create_operations`][2].
- `count`: The number of matrices to update.

[2]: pll_utree_t.md#Notable-Functions

### Misc

```
PLL_EXPORT unsigned int * pll_compress_site_patterns(char ** sequence,
                                                     const pll_state_t * map,
                                                     int count,
                                                     int * length);
```

Compresses the MSA in place, i.e. the buffer `sequence` is changed to store the
compressed alignment.

- `sequence`: The alignment to compress, should be the one from `pll_msa_t`.
- `map`: Encoding map, ie `pll_map_nt`.
- `count`: number of sequences
- `length`: outparameter of the new length.

Data Structures
-------------------------------------------------------------------------------

Data structures have their own pages:

- [`pll_partition_t`](pll_partition_t.md)
- [`pll_utree_t`](pll_utree_t.md)

Errors
-------------------------------------------------------------------------------

Many functions which don't directly return a value instead return whether or
not the function was successful. In this case, the function was successful, then
`PLL_SUCCESS` is returned. Otherwise, `PLL_FAILURE` is returned. This value is
the same regardless of the error.

The _type_ of error is stored in `pll_errno`, and the message is stored in
`pll_errmsg`. `pll_errno` is an `int`, while `pll_errmsg` is a `char[200]`. Both
of these are present at the global scope.

Here are the possible values for `pll_errno`

- `PLL_ERROR_FILE_OPEN`
- `PLL_ERROR_FILE_SEEK`
- `PLL_ERROR_FILE_EOF`
- `PLL_ERROR_FASTA_ILLEGALCHAR`
- `PLL_ERROR_FASTA_UNPRINTABLECHAR`
- `PLL_ERROR_FASTA_INVALIDHEADER`
- `PLL_ERROR_PHYLIP_SYNTAX`
- `PLL_ERROR_PHYLIP_LONGSEQ`
- `PLL_ERROR_PHYLIP_NONALIGNED`
- `PLL_ERROR_PHYLIP_ILLEGALCHAR`
- `PLL_ERROR_PHYLIP_UNPRINTABLECHAR`
- `PLL_ERROR_NEWICK_SYNTAX`
- `PLL_ERROR_MEM_ALLOC`
- `PLL_ERROR_PARAM_INVALID`
- `PLL_ERROR_TIPDATA_ILLEGALSTATE`
- `PLL_ERROR_TIPDATA_ILLEGALFUNCTION`
- `PLL_ERROR_TREE_CONVERSION`
- `PLL_ERROR_INVAR_INCOMPAT`
- `PLL_ERROR_INVAR_PROPORTION`
- `PLL_ERROR_INVAR_PARAMINDEX`
- `PLL_ERROR_INVAR_NONEFOUND`
- `PLL_ERROR_AB_INVALIDMETHOD`
- `PLL_ERROR_AB_NOSUPPORT`
- `PLL_ERROR_SPR_TERMINALBRANCH`
- `PLL_ERROR_SPR_NOCHANGE`
- `PLL_ERROR_NNI_INVALIDMOVE`
- `PLL_ERROR_NNI_TERMINALBRANCH`
- `PLL_ERROR_STEPWISE_STRUCT`
- `PLL_ERROR_STEPWISE_TIPS`
- `PLL_ERROR_STEPWISE_UNSUPPORTED`
- `PLL_ERROR_EINVAL`
- `PLL_ERROR_MSA_EMPTY`
- `PLL_ERROR_MSA_MAP_INVALID`
- `PLL_ERROR_TREE_INVALID`
