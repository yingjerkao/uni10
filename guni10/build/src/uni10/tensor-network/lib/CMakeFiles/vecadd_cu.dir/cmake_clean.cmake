FILE(REMOVE_RECURSE
  "../../../../lib/libvecadd_cu.pdb"
  "../../../../lib/libvecadd_cu.a"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/vecadd_cu.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
