#!/bin/bash


ls $1 | while read FILE; do
  cat $1/"$FILE" $2/"$FILE" $3/"$FILE" $4/"$FILE" $5/"$FILE" $6/"$FILE" >> $7/"$FILE"
done
echo "done"
