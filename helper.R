

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
