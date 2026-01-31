# Select number of factors via eigenvalue proportion rule

This follows the simple rule used in the nutrition case study: choose K
as the number of eigenvalues explaining \> `cutoff` proportion of total
variance.

## Usage

``` r
select_k_from_sigma(Sigma, cutoff = 0.05)
```

## Arguments

- Sigma:

  A symmetric covariance matrix.

- cutoff:

  Proportion threshold (default 0.05).

## Value

Integer K.
