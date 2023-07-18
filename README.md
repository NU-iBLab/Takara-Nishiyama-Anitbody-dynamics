# Modeling COVID-19 vaccine booster-elicited antibody response and impact of infection history

Takara Nishiyama, Yuichiro Miyamatsu, Hyeongki Park, Naotoshi Nakamura, Risa Yokokawa Shibata, Shingo Iwami, Yoji Nagasaki

## Description
This code is an analysis code for the paper "Modeling COVID-19 vaccine booster-elicited antibody response and impact of infection history"

## 1 monolix

'sample_observation.csv' are IgG(RBD) data of a randomly generated participant.

'model.txt' is the description for mathematical model of antibody dynamics.

These files are used to estimate parameters of the viral dynamics model in MONOLIX2021R2.

## 2 R_code

'estimatedIndividualParameters.txt' are the estimated parameters of antibody dynamics model.

'model.c' is a file for performing numerical calculations. It is used in 'Analysis_code.R'.

'Analysis_code.R' shows the numerical calculations and statistical analysis performed in the paper.

