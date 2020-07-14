#!/bin/sh
sed -nr '/ORIGIN/,\/\//p' < $1 | awk 'NR > 1 {print $2 $3 $4 $5 $6 $7}'