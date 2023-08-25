This archive contains the source code for Celeris software.
Celeris is developed by Sasan Tavakkol, on top of a previous 
demo project developed by Stephen Thompson.
The Celeris project was advised by Patrick Lynett, PhD.
This readme file is written orginally by Stephen,
and edited later by Sasan.


Contacts for Sasan Tavakkol:
        tavakkol [at] usc [.] edu or sasantavakkol [at] yahoo [.] com
Contacts for Patrick Lynett:
        plynett [at] usc [.] edu
Contacts for Stephen Thompson:
        stephen [at] solarflare [.] org [.] uk

Celeris
=========

This repository hosts Celeris Advent, a free wave modeling software. Celeris is also available as Celeris Citius and Celeris Desktop, paid versions with greater capabilities and support:
https://www.celerialabs.com/

Standalone Application
=========

Celeris Advent is a simple executable file that works immediatly for most users. Download the most recent release version (v1.3.1) of Celeris Advent from:
https://www.celeria.org/download

Or download another release (e.g. vX.Y.Z) from the GitHub side-bar. Once downloaded, unzip the folder and navigate to 
"~/Downloads/Celeris.Advent.vX.Y.Z/Celeris Advent (vX.Y.Z)/Celeris Advent/"

Double-click "Celeris.exe" to launch the program. It will prompt you to select a CML file (Celeris Input Script) to set-up the simulation. Examples are contained within folders at:
"~/Downloads/Celeris.Advent.vX.Y.Z/Celeris Advent (vX.Y.Z)/Examples/"

The Celiris user manual provides advice on getting started, it may be found in:
"~/Downloads/Celeris.Advent.vX.Y.Z/Celeris Advent (vX.Y.Z)/Users.Manual.v1.0"

Or at:
https://www.celeria.org/downloads

Requirements for Standalone Application
=========

Operating System (tested): Windows 10 
Operating System (not recently tested): Windows 8, 7, XP
Graphics: DirectX, DirectX 11 Run-Time SDK 
Hardware: GPU or integrated graphics capable of running DirectX 11

For Mac OS and Linux users, there has been some success in running Celeris Advent through emulation tools such as Wine and virtual machines. However, a guide is not currently provided.

Common Issues
=========

If simulations "blow-up" (i.e. the waves go unstable and crash the application), it can be fixed in a few ways
1 - Use a smaller time-step in the CML file, or change the "NOMINAL_CFL" number in the GUI settings tab to be about 0.1.
2 - Use larger "epsilon" number by editing the CML file, but note that the depth is compared to sqrt(sqrt(epsilon)) and the model is not accurate for depth smaller than this value.
3 - Lower the minmod limiter "theta" in the GUI. It can change between 1 (most dissipative) to 2 (least dissipative). Smaller theta results in a more stable simulation, but causes more numerical dissipation.
4 - Smooth the bathymetry with in-built GUI tools. Sudden, sharp changes in the landscape can cause waves to go unstable.

If performance appears very slow, Celeris may be using your systems integrated graphics instead of a dedicated GPU. You can check in task manager, under the "Performance" tab to see if your dedicated GPU is running.

Some features of the DirectX 11 Run-Time SDK are not natively installed on Windows 8/10+. This may crash Celeris when starting a simulation, throwing error "CompileFromFile Failed". You can install the DirectX 11 Run-Time SDK library files for Windows 8/10+ by downloading a Microsoft tool from the following link:
https://www.microsoft.com/en-us/download/details.aspx?displaylang=en&id=35&tduid=(d5c34b75227807aea00cad2fc302b7da)(256380)(2459594)(TnL5HPStwNw-wC6WISUgqQM1hpaZ1Sq5Fg)()

Run "dxwebsetup.exe". When finished, launch Celeris. If the issue is not resolved try the following steps: 
>> Right-click on the dxwebsetup.exe and select "Properties"
>> Open the "Compatibility" tab.
>> Check the "Run this program in compatibility mode" box and select "Windows XP (Service Pack 2)", or whichever Windows XP version is available for you.
>> Click "Apply" and then "OK"
>> Launch dxwebsetup.exe again
>> Double-click Celeris.exe after dxwebsetup.exe completes

If the problem persists, try opening "dxdiag" from the Windows search bar. This application will help diagnose further DirectX issues.

Compiling from Source
=========

For advanced users, a Visual Studio 2008 solution file is provided ("Celeris.sln" in
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
