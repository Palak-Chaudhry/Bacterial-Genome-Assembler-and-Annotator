#! usr/bin/env bash

awk '
  # Parse the .fas file and store the virulence factors and their information
  BEGIN { FS = "|" }
  /^>/ {
    id = $1
    gsub(">", "", id)
    sub(/ .*$/, "", id)
    info[id] = $0
  }

  # Parse the TSV file and add the virulence factor information as a new column
  NR == FNR { next }
  FNR == 1 { print $0, "Virulence Info" }
  {
    found = "No Information"
    for (id in info) {
      if (index($0, id) != 0) {
        found = info[id]
        break
      }
    }
    print $0, found
  }
' VFDB_setA_pro.fas $1 > $2
