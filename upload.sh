#!/usr/bin/env bash

echo "compile jar"
sbt package
echo "upload to kova server"
scp target/scala-2.10/dnaseqanalyzer_2.10-1.0.jar bswastanto@kova-01.ewi.tudelft.nl:~/LAB_3