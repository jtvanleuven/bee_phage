library(ggplot2)
library(ggplot2movies)
devtools::unload('UpSetR')
library(ComplexUpset)
head(movies, 3)
genres = colnames(movies)[18:24]
movies[genres] = movies[genres] == 1
t(head(movies[genres], 3))
movies[movies$mpaa == '', 'mpaa'] = NA
movies = na.omit(movies)
upset(movies, genres, name='genre', width_ratio=0.1)

set.seed(0)
upset(
  movies,
  genres,
  annotations = list(
    'Length'=list(
      aes=aes(x=intersection, y=length),
      geom=geom_boxplot()
    ),
    'Rating'=list(
      aes=aes(x=intersection, y=rating),
      geom=list(
        # checkout ggbeeswarm::geom_quasirandom for better results!
        geom_jitter(aes(color=log10(votes))),
        geom_violin(width=1.1, alpha=0.5)
      )
    )
  ),
  min_size=10,
  width_ratio=0.1
)


upset(
  movies,
  genres,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=F,
      aes=aes(fill=mpaa)
    )
  ),
  width_ratio=0.1
)

