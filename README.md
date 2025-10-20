[![DOI](https://zenodo.org/badge/124278536.svg)](https://zenodo.org/badge/latestdoi/124278536) [![GitHub contributors](https://img.shields.io/github/contributors/michaelwenk/casekit.svg)](https://github.com/michaelwenk/casekit/graphs/contributors/) [![GitHub issues](https://img.shields.io/github/issues/michaelwenk/casekit.svg)](https://github.com/michaelwenk/casekit/issues/) [![GitHub release](https://img.shields.io/github/release/michaelwenk/casekit.svg)](https://github.com/michaelwenk/casekit/releases/)

# The Computer-Assisted-Structure-Elucidation Kit (CASEkit)

## Introduction

This project depends on the Chemistry Development Project (CDK), hosted under https://cdk.github.io/
Please refer to these pages for updated information and the latest version of the CDK. CDK's API documentation is
available though our [Github site](http://cdk.github.io/cdk/).

## Download Source code

This assumes that you have git working on your system and you have initialised your local repository.

Then, downloading casekit is just a matter of

    git clone https://github.com/michaelwenk/casekit

## Install NMRShiftDB Predictor Tool

casekit depends on the prediction tool from NMRShiftDB, the _predictorc.jar_ file. Maven needs it to be locally installed via

    cd casekit && \
    sh install_predictorc_jar.sh

## Compiling

Compiling the library is performed with Apache Maven and requires Java 1.7 or later:

    mvn clean package

will create an all-in-one-jar under ./target which you can use in your Java project.
