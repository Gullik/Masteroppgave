#!/bin/bash

#Run

pdflatex thesis.tex
biber thesis
pdflatex thesis.tex
