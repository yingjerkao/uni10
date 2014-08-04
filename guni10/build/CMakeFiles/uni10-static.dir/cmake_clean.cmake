FILE(REMOVE_RECURSE
  "lib/libuni10.pdb"
  "lib/libuni10.a"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/uni10-static.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
