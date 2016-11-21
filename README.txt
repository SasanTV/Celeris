This archive contains the source code for Celeris software.
Celeris is developed by Sasan Tavakkol, on top of a previous 
demo project developed by Stephen Thompson.
The Celeris project was advised by Patrick Lynett, PhD.
This readme file is written orginally by Stephen,
and edited later by Sasan.


Contacts for Sasan Tavakkol:
        tavakkol@usc.edu or sasantavakkol@yahoo.com
Contacts for Patrick Lynett:
        plynett@usc.edu
Contacts for Stephen Thompson:
        stephen@solarflare.org.uk


Compiling
=========

A Visual Studio 2008 solution file is provided ("Celeris.sln" in
the "msvc" directory). Compiling should be straightforward, although
you will need to have the latest DirectX SDK correctly installed on
your machine.


Roadmap
=======

The source consists of three projects:

1) Coercri -- This is a simple wrapper library around operating system
functions (graphics, sounds and so on).

2) Guichan -- This is an open source GUI library which was used to
create the buttons and sliders at the right-hand side of the screen.

3) Celeris-- This is the main project. Briefly, the main files are as follows:

 * engine.cpp -- Main "engine" for the simulation, contains all the
   code that drives the GPU. The bulk of the code is found here.

 * compute.hlsl -- Contains shaders for doing the numerical simulation of
   the extended Boussinesq equations on the GPU.

 * graphics.fx -- Contains shaders for creating the graphical
   appearance of the land and water surfaces.

 * gui_manager.cpp -- Creates the GUI at the right-hand side of the
   screen.

 * settings.cpp -- Stores and manages simulation settings.

 * terrain_heightfield.cpp -- Contains formulas for determining the
   shape of the terrain.


Copyright / Legal Information
=============================

The Celeris software source code consists of four open source
components each of which has their own copyright and licence conditions.
The four components are:

1) Celeris itself
2) Coercri
3) Guichan
4) TinyXML

Each source code file identifies clearly (at the top of the file)
which licence applies to that file.