

my_ggsave = function(name,
                     width,
                     height,
                     .results.dir = results.dir,
                     .overleaf.dir = overleaf.dir) {
  setwd(.results.dir)
  ggsave( name,
          width = width, 
          height = height)
  
  setwd(.overleaf.dir)
  ggsave( name,
          width = width, 
          height = height)
}

nuni = function(x) length(unique(x))


# internal fn from R package EValue, for use customizing sens_plot
g = Vectorize( function(x) {
  # define transformation in a way that is monotonic over the effective range of B (>1)
  # to avoid ggplot errors in sens_plot
  # helper function for confounded_meta
  if ( is.na(x) ) return(NA)
  if (x < 1) return( x / 1e10 )
  x + sqrt( x^2 - x )
} )