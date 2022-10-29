# seriate.dist returns expected results

    Code
      str(os[deterMethods])
    Output
      List of 20
       $ BBURCG       : 'ser_permutation_vector' Named int [1:4] 1 2 4 3
        ..- attr(*, "names")= chr [1:4] "a" "b" "d" "c"
        ..- attr(*, "method")= chr "BBURCG"
       $ BBWRCG       : 'ser_permutation_vector' Named int [1:4] 1 2 4 3
        ..- attr(*, "names")= chr [1:4] "a" "b" "d" "c"
        ..- attr(*, "method")= chr "BBWRCG"
       $ GW           :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -4 1 -2 -3 2
        ..$ height     : num [1:3] 1 1 2.24
        ..$ order      : int [1:4] 1 2 4 3
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "complete"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "GW"
       $ GW_average   :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -4 1 -2 -3 2
        ..$ height     : num [1:3] 1 1 1.99
        ..$ order      : int [1:4] 1 2 4 3
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "average"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "GW_average"
       $ GW_complete  :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -4 1 -2 -3 2
        ..$ height     : num [1:3] 1 1 2.24
        ..$ order      : int [1:4] 1 2 4 3
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "complete"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "GW_complete"
       $ GW_single    :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -4 1 -2 -3 2
        ..$ height     : num [1:3] 1 1 1.73
        ..$ order      : int [1:4] 1 2 4 3
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "single"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "GW_single"
       $ GW_ward      :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -4 1 -2 -3 2
        ..$ height     : num [1:3] 1 1 2.65
        ..$ order      : int [1:4] 1 2 4 3
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "ward.D2"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "GW_ward"
       $ HC           :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -3 1 -2 -4 2
        ..$ height     : num [1:3] 1 1 2.24
        ..$ order      : int [1:4] 1 2 3 4
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "complete"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "HC"
       $ HC_average   :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -3 1 -2 -4 2
        ..$ height     : num [1:3] 1 1 1.99
        ..$ order      : int [1:4] 1 2 3 4
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "average"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "HC_average"
       $ HC_complete  :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -3 1 -2 -4 2
        ..$ height     : num [1:3] 1 1 2.24
        ..$ order      : int [1:4] 1 2 3 4
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "complete"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "HC_complete"
       $ HC_single    :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -3 1 -2 -4 2
        ..$ height     : num [1:3] 1 1 1.73
        ..$ order      : int [1:4] 1 2 3 4
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "single"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "HC_single"
       $ HC_ward      :Classes 'ser_permutation_vector', 'hclust'  hidden list of 7
        ..$ merge      : int [1:3, 1:2] -1 -3 1 -2 -4 2
        ..$ height     : num [1:3] 1 1 2.65
        ..$ order      : int [1:4] 1 2 3 4
        ..$ labels     : chr [1:4] "a" "b" "c" "d"
        ..$ method     : chr "ward.D2"
        ..$ call       : language hclust(d = d, method = control$method)
        ..$ dist.method: chr "euclidean"
        ..- attr(*, "method")= chr "HC_ward"
       $ Identity     : 'ser_permutation_vector' Named int [1:4] 1 2 3 4
        ..- attr(*, "names")= chr [1:4] "a" "b" "c" "d"
        ..- attr(*, "method")= chr "Identity"
       $ MDS          : 'ser_permutation_vector' Named int [1:4] 3 4 2 1
        ..- attr(*, "names")= chr [1:4] "c" "d" "b" "a"
        ..- attr(*, "method")= chr "MDS"
       $ MDS_angle    : 'ser_permutation_vector' Named int [1:4] 1 2 4 3
        ..- attr(*, "names")= chr [1:4] "a" "b" "d" "c"
        ..- attr(*, "method")= chr "MDS_angle"
       $ MDS_metric   : 'ser_permutation_vector' Named int [1:4] 3 4 2 1
        ..- attr(*, "names")= chr [1:4] "c" "d" "b" "a"
        ..- attr(*, "method")= chr "MDS_metric"
       $ R2E          : 'ser_permutation_vector' Named int [1:4] 4 3 1 2
        ..- attr(*, "names")= chr [1:4] "d" "c" "a" "b"
        ..- attr(*, "method")= chr "R2E"
       $ Spectral     : 'ser_permutation_vector' Named int [1:4] 3 4 2 1
        ..- attr(*, "names")= chr [1:4] "c" "d" "b" "a"
        ..- attr(*, "method")= chr "Spectral"
       $ Spectral_norm: 'ser_permutation_vector' Named int [1:4] 3 4 2 1
        ..- attr(*, "names")= chr [1:4] "c" "d" "b" "a"
        ..- attr(*, "method")= chr "Spectral_norm"
       $ VAT          : 'ser_permutation_vector' Named int [1:4] 3 4 2 1
        ..- attr(*, "names")= chr [1:4] "c" "d" "b" "a"
        ..- attr(*, "method")= chr "VAT"

# seriate.matrix returns expected results

    Code
      str(os[deterMethods])
    Output
      List of 5
       $ CA       :List of 2
        ..$ : 'ser_permutation_vector' Named int [1:4] 1 2 4 3
        .. ..- attr(*, "names")= chr [1:4] "a" "b" "c" "d"
        .. ..- attr(*, "method")= chr "CA"
        ..$ : 'ser_permutation_vector' Named int [1:5] 2 1 3 4 5
        .. ..- attr(*, "names")= chr [1:5] "A" "B" "C" "D" ...
        .. ..- attr(*, "method")= chr "CA"
        ..- attr(*, "class")= chr [1:2] "ser_permutation" "list"
       $ Identity :List of 2
        ..$ : 'ser_permutation_vector' Named int [1:4] 1 2 3 4
        .. ..- attr(*, "names")= chr [1:4] "a" "b" "c" "d"
        .. ..- attr(*, "method")= chr "Identity"
        ..$ : 'ser_permutation_vector' Named int [1:5] 1 2 3 4 5
        .. ..- attr(*, "names")= chr [1:5] "A" "B" "C" "D" ...
        .. ..- attr(*, "method")= chr "Identity"
        ..- attr(*, "class")= chr [1:2] "ser_permutation" "list"
       $ PCA      :List of 2
        ..$ : 'ser_permutation_vector' Named int [1:4] 1 2 4 3
        .. ..- attr(*, "names")= chr [1:4] "a" "b" "c" "d"
        .. ..- attr(*, "method")= chr "PCA"
        ..$ : 'ser_permutation_vector' Named int [1:5] 2 1 3 4 5
        .. ..- attr(*, "names")= chr [1:5] "A" "B" "C" "D" ...
        .. ..- attr(*, "method")= chr "PCA"
        ..- attr(*, "class")= chr [1:2] "ser_permutation" "list"
       $ PCA_angle:List of 2
        ..$ : 'ser_permutation_vector' Named int [1:4] 1 2 4 3
        .. ..- attr(*, "names")= chr [1:4] "a" "b" "c" "d"
        .. ..- attr(*, "method")= chr "PCA_angle"
        ..$ : 'ser_permutation_vector' Named int [1:5] 4 5 3 1 2
        .. ..- attr(*, "names")= chr [1:5] "A" "B" "C" "D" ...
        .. ..- attr(*, "method")= chr "PCA_angle"
        ..- attr(*, "class")= chr [1:2] "ser_permutation" "list"
       $ Reverse  :List of 2
        ..$ : 'ser_permutation_vector' Named int [1:4] 4 3 2 1
        .. ..- attr(*, "names")= chr [1:4] "a" "b" "c" "d"
        .. ..- attr(*, "method")= chr "Reverse"
        ..$ : 'ser_permutation_vector' Named int [1:5] 5 4 3 2 1
        .. ..- attr(*, "names")= chr [1:5] "A" "B" "C" "D" ...
        .. ..- attr(*, "method")= chr "Reverse"
        ..- attr(*, "class")= chr [1:2] "ser_permutation" "list"

# data.frame seriation works as expected

    Code
      permute(df, o)
    Output
        A B C D E
      a 1 1 0 0 0
      b 1 1 1 0 0
      d 1 0 1 1 1
      c 0 0 1 1 1

