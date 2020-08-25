[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# ULTRA-PHAST
One-to-many variant of ULTRA

This framework contains C++ code for ULTRA-PHAST, an algorithmic approach for fast answering of one-to-all and one-to-many queries in public transit networks with a secondary transportation mode (e.g., walking or taxis).
It was developed at [KIT](https://www.kit.edu) in the [group of Prof. Dorothea Wagner](https://i11www.iti.kit.edu/).

## Usage

This framework contains code for the ULTRA preprocessing step (i.e., contraction hierarchies and shortcut computation) as well as two query algorithms (UP-RAPTOR and UP-CSA). Each of the major components of the framework can be used/tested separately using one of the programs in the ``Runnables`` folder:

* ``BuildCoreCH`` takes the public transit network and computes the core-CH required for the ULTRA preprocessing
* ``ComputeShortcuts`` takes the public transit network as well as the core-CH and computes the ULTRA shortcuts
* ``BuildCH`` takes the transfer graph and computes the CH required in the query phase for initial/final transfers
* ``RunRAPTORQueriesToBall`` takes the public transit network, shortcuts, and CH and performs UP-RAPTOR experiments, with the target set randomly chosen from a ball
* ``RunRAPTORQueriesToStops`` takes the public transit network, shortcuts, and CH and performs UP-RAPTOR experiments, with a target set of all stops
* ``RunRAPTORQueriesToVertices`` takes the public transit network, shortcuts, and CH and performs UP-RAPTOR experiments, with a target set of all vertices
* ``RunCSAQueriesToBall`` takes the public transit network, shortcuts, and CH and performs UP-CSA experiments, with the target set randomly chosen from a ball
* ``RunCSAQueriesToStops`` takes the public transit network, shortcuts, and CH and performs UP-CSA experiments, with a target set of all stops
* ``RunCSAQueriesToVertices`` takes the public transit network, shortcuts, and CH and performs UP-CSA experiments, with a target set of all vertices

The ``Makefile`` located in the ``Runnables`` folder contains instructions for building all of the above programs. Simply edit the top part of the Makefile to adjust the compiler and flags available to you and run ``make`` afterwards.

All of the above programs use a custom binary format for loading the public transit network as well as the transfer graph. As an example we provide the public transit network of Switzerland together with a transfer graph extracted from OpenStreetMap in the appropriate binary format at [https://i11www.iti.kit.edu/PublicTransitData/Switzerland/binaryFiles/](https://i11www.iti.kit.edu/PublicTransitData/Switzerland/binaryFiles/). 

## Publications

The algorithms in this framework are based on the following publications:

* *PHAST: Hardware-Accelerated Shortest Path Trees*
  Daniel Delling, Andrew V Goldberg, Andreas Nowatzyk, Renato F Werneck
  In: Journal of Parallel and Distributed Computing, 73(7):940–952, 2013
  [pdf](https://i11www.iti.kit.edu/extra/publications/dgnw-phast-12.pdf)
* *UnLimited TRAnsfers for Multi-Modal Route Planning: An Efficient Solution*
  Moritz Baum, Valentin Buchhold, Jonas Sauer, Dorothea Wagner, Tobias Zündorf
  In: Proceedings of the 27th Annual European Symposium on Algorithms (ESA'19), Leibniz International Proceedings in Informatics, pages 14:1–14:16, 2019
  [pdf](https://drops.dagstuhl.de/opus/volltexte/2019/11135/pdf/LIPIcs-ESA-2019-14.pdf) [arXiv](https://arxiv.org/abs/1906.04832)
* *An Efficient Solution for One-to-Many Multi-Modal Journey Planning*
  Jonas Sauer, Dorothea Wagner, Tobias Zündorf
  To appear in: 20th Workshop on Algorithmic Approaches for Transportation Modelling, Optimization, and Systems (ATMOS 2020), 2020
