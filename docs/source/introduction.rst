.. title:: Introduction

.. _introduction:

============
Introduction
============

**Kai Thin, Tejas Patel - Department of Mechanical Engineering, Michigan State University, East Lansing, MI, 2022**

==============================

vanDANA is a highly efficient FEM Immersed Boundary (IB) based Flow-thermal FSI solver utilizing the `FEniCS <https://fenicsproject.org/>`__ library (version 2019.2.0). The solver is based on the `Distibuted Langrange Multiplier based Fictiious Domain method <https://www.sciencedirect.com/science/article/pii/S0021999105000148>`__ and is extended to deal with `heattransfer <https://www.sciencedirect.com/science/article/pii/S0021999106000167>`__.The interpolation of variables is conducted using the smeared `delta-functions <https://www.sciencedirect.com/science/article/pii/S0021999109004136>`__. Additionally, the flow solver is incompressible and has the option of choosing from various stabilization schemes : SUPG, PSPG and Crosswind; and the structure can be set as either incompressible/compressible.

This manual intends to create a basic understanding of the numerical algorithm and also provides a detailed explanation of the workflow for the reader to set up a user-defined problem. Here we expect the reader to have some basic knowledge of coding partial differential equations in `FEniCS <https://fenicsproject.org>`__.

Note: vanDANA is licensed under the GNU GPL, version 3 or any later version and is Copyright (2022) by the authors.
