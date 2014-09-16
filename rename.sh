#!/bin/bash

regex=$1

while read fn
do

fn1=$(echo $fn | sed "$regex")

mv $fn $fn1

done