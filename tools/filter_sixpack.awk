#!/bin/awk -f

!/^>/ {next}
{getline s}
length(s) >=i { print $0 "\n" s}


