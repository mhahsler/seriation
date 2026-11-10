# Zoo Data Set

A database containing characteristics of different animals. The database
was created and donated by Richard S. Forsyth and is available from the
UCI Machine Learning Repository (Newman et al, 1998).

## Format

A data frame with 101 observations on the following 17 variables.

- `hair`:

  a numeric vector

- `feathers`:

  a numeric vector

- `eggs`:

  a numeric vector

- `milk`:

  a numeric vector

- `airborne`:

  a numeric vector

- `aquatic`:

  a numeric vector

- `predator`:

  a numeric vector

- `toothed`:

  a numeric vector

- `backbone`:

  a numeric vector

- `breathes`:

  a numeric vector

- `venomous`:

  a numeric vector

- `fins`:

  a numeric vector

- `legs`:

  a numeric vector

- `tail`:

  a numeric vector

- `domestic`:

  a numeric vector

- `catsize`:

  a numeric vector

- `class`:

  a factor with levels `amphibian` `bird` `fish` `insect` `invertebrate`
  `mammal` `reptile`

## Source

David Aha, Patrick Murphy, Christopher Merz, Eamonn Keogh, Cathy Blake,
Seth Hettich, David Newman, Arthur Asuncion, Moshe Lichman, Dheeru Dua,
Casey Graff (2023): UCI Machine Learning Repository,
<https://archive.ics.uci.edu/>, University of California, Irvine.

## See also

Other data:
[`Chameleon`](https://michael.hahsler.net/seriation/reference/Chameleon.md),
[`Irish`](https://michael.hahsler.net/seriation/reference/Irish.md),
[`Munsingen`](https://michael.hahsler.net/seriation/reference/Munsingen.md),
[`Psych24`](https://michael.hahsler.net/seriation/reference/Psych24.md),
[`SupremeCourt`](https://michael.hahsler.net/seriation/reference/SupremeCourt.md),
[`Townships`](https://michael.hahsler.net/seriation/reference/Townships.md),
[`Wood`](https://michael.hahsler.net/seriation/reference/Wood.md),
[`create_lines_data()`](https://michael.hahsler.net/seriation/reference/create_lines_data.md),
[`is.robinson()`](https://michael.hahsler.net/seriation/reference/is.robinson.md)

## Examples

``` r
data("Zoo")
x <- scale(Zoo[, -17])


d <- dist(x)
pimage(d)


order <- seriate(d, method = "tsp")
pimage(d, order)
```
