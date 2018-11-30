# Worhp-with-ADOL-C
This project presents a way to combine the NLP solver WORHP and the automatic differentiation tool ADOL-C. We have reimplemented the "Getting started" example of WOHRP available [here](https://worhp.de/content/cppexample), but instead of computing the required derivatives by hand, we are using ADOL-C drivers. ADOL-C is available [here](https://www.coin-or.org/download/source/ADOL-C/)

# Motivation
The used case of this combination is rapid prototyping. Modeling a dynamic process is usually an iterative process. The complexity of the model rises gradually. We could either recompute all derivatives analytically every time the model changes or use ADOL-C to do this automatically for us.

# Code Example
ADOL-C makes use of active (adoule) and passive data types (any other type).
