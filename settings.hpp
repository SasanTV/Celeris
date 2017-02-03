/*
 * FILE:
 *   settings.hpp
 *
 * AUTHORS:
 *	 Sasan Tavakkol <tavakkol@usc.edu> or <sasantavakkol@yahoo.com>
 *   Stephen Thompson <stephen@solarflare.org.uk>
 *
 * LAST UPDATE:
 *	 26-Aug-2016
 * CREATED:
 *   24-Oct-2011
 *
 COPYRIGHT:
 *   Copyright (C) 2016, Sasan Tavakkol.
 *   Copyright (C) 2012, Stephen Thompson.
 *
 *   This file is part of Celeris software. Celeris is an interactive
 *	 Boussinesq-type, coastal wave solver and visualizer. 
 *
 *   Celeris is free software: you can redistribute it
 *   and/or modify it under the terms of the GNU General Public
 *   License as published by the Free Software Foundation, either
 *   version 3 of the License, or (at your option) any later version.
 *
 *   The Shallow Water Demo is distributed in the hope that it will be
 *   useful, but WITHOUT ANY WARRANTY; without even the implied
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *   See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the Shallow Water Demo. If not, see
 *   <http://www.gnu.org/licenses/>.
 *
 */

#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <vector>


enum SettingType {
    S_NULL,
    S_NEW_TAB,
    S_LABEL,
    S_SLIDER,
    S_SLIDER_INT,
    S_SLIDER_MULT_4,
    S_CHECKBOX,
	S_WATER_SURFACE_SHADING,
	S_WATER_SURFACE_SHADING_VARIABLE,
	S_TERRAIN_TEXTURE,
	S_SKYBOX_TEXTURE,
	S_BOUNDARY_DROPDOWN,
	S_BOUNDARY_OPTIONS_DROPDOWN,
	S_LEFT_MOUSE_DROPDOWN,
	S_TEXTFIELD

};


// Sasan: This enum's use case is changed and comments are not correct
enum ResetType {
    R_NONE,
    R_TERRAIN,
    R_MESH,     // remesh, no water
    R_VALLEY,   // remesh, setup water in valley
    R_SEA,      // remesh, fill to sea level
    R_SQUARE    // remesh, fill in a square in the middle
};

enum LeftMouseSettings {
    LM_ADD_WATER,
    LM_REMOVE_WATER,
    LM_STIR_WATER,
    LM_RAISE_TERRAIN,
    LM_LOWER_TERRAIN
};


enum ShadingOptions{
    PHOTOREALISTIC,
    PARULA,
    JET,
    ZEBRA,
    CUSTOM
};

enum BoundaryType{
    B_SOLID,
    B_SPONGE,
    B_SINEWAVE,
	B_IRREGULARWAVE,
	B_UNIFORMTIMESERIES
};

enum BoundarySide{
    B_WEST,
    B_EAST,
    B_SOUTH,
	B_NORTH
};

struct Setting {
    const char * name;
    const char * unit;
    SettingType type;
    ResetType reset_type;
    double min;
    double max;
    double value;
};

struct SineWaveSetting {
	float amplitude;
	float period;
	float theta; //angle in radians to x direction.
	SineWaveSetting (): amplitude(0), period(0), theta(0){}
};

struct FlowParameters {
	float w;
	float hu;
	float hv; 
};


struct UniformTimeSeries{
	std::string fileName;
	std::vector<std::vector<float>> data;
	int iter;
	UniformTimeSeries () : iter(1), fileName("NA"){}
	void reset_iter (){
		iter = 1;
	}
};


struct BoundarySetting {
	std::string type;
	int width; //number of cells in boundary.
	float waterLevel; // this can be sea level, or a control depth at river output. It is measured from datum.
	SineWaveSetting sineWaveSetting;
	UniformTimeSeries uniformTimeSeries;
//	FlowParameters ForcedInput;
	bool hasChanged;
	BoundarySetting (): type("Solid"), width(2), waterLevel(0), hasChanged(false)  {}
};


struct Point {
	int x;
	int y;
	Point ():x(0),y(0){}
	Point (int j, int i):x(j),y(i){}
	Point (const Point &p): x(p.x), y(p.y){}
};

struct Range {
	std::string name;
	Point bottomLeft;
	Point topRight;
	Range (): bottomLeft(0,0), topRight(0,0){}
	Range (Point bl, Point tr): bottomLeft(bl), topRight (tr) {}
	Range (std::string n, Point bl, Point tr): name(n), bottomLeft(bl), topRight (tr) {}
};

// a soliton is a solitary wave
struct Soliton{
	
	float param_h;
	float theta;
	float length;
	float xc;
	float yc;
	Soliton () : param_h(0), theta (0), xc(0), yc(0), length(0){}
	Soliton (float H,float theta_in, float xcenter,float ycenter, float length_in) : param_h(H), theta (theta_in), xc(xcenter), yc(ycenter), length(length_in){}
	void setSoliton(float H, float theta_in, float xcenter, float ycenter, float length_in){
		param_h = H; theta = theta_in; xc = xcenter; yc = ycenter; length = length_in;
	}
};

struct CameraSetting{
	float fov;
	float x, y, z;
	float pitch, yaw;
	CameraSetting() : fov(50), x(0), y(0), z(0), pitch(0), yaw(0){}
};
struct SurfaceShadingSetting{
	ShadingOptions type;
	int shadingVariable;
	bool autoColormap;
	bool autoInundationDepth;
	float colormapMin;
	float colormapMax;
	float drylandDepthOfInundation;
	float minInundation; // for visualization purposes
	float maxInundation; // for visualization purposes


	SurfaceShadingSetting(){
		type = PHOTOREALISTIC;
		shadingVariable = 0; // 0 is eta
		autoColormap = true;
		autoInundationDepth = true;
		colormapMin = -1.0f;	
		colormapMax = 1.0f;	
		drylandDepthOfInundation = 0.1; // 10 cm.
		minInundation = 0.001; // 1 mm
		maxInundation = 1; // 1 m
	}
};

struct TerrainTextureSetting{
	int type;
	bool autoColormap;
	float colormapMin;
	float colormapMax;

	TerrainTextureSetting(){
		type = 0; // 0 is SAND
		autoColormap = true;
		colormapMin = -1.0f;	
		colormapMax = 1.0f;	
	}
};

struct Lighting{
	float ambientLight;
	float sun_altitude; // in degrees
	float sun_azimuth; // in degrees

	Lighting(){
		ambientLight = 1.0f;
		sun_altitude = 45.0f;
		sun_azimuth = 30.0f;
	}
};

struct GraphicsSetting{
	bool autoCam;

	SurfaceShadingSetting surfaceShading;
	TerrainTextureSetting terrainTexture;
	int skyboxType;
	
	Lighting lighting;

	bool gridOn;
	float gridScale;
	float verticalScale;

	float fresnelCoef,refractive_index, att_1, att_2;
	CameraSetting camera;
	GraphicsSetting (){
		autoCam = true;
		skyboxType = 0;
		gridOn = true;
		gridScale = 1.0f;
		verticalScale = 1.0f;
		fresnelCoef = 0.5f; refractive_index = 0.5f; att_1 = 0.5f; att_2 = 0.5f;
	}
};

struct InitSetting {
	
	//name
	std::string project_name;

	//Paths
	std::string exePath;
	std::string initCMLName;
	std::string initWFileName;
	std::string bathymetryFileName;
	std::string westIrrWaveFileName,  eastIrrWaveFileName;
	std::string southIrrWaveFileName, northIrrWaveFileName;

	bool rereadBathy;

	std::string logPath;
	
	//Model
	bool isBoussinesq;
	float epsilon;
	float theta;
	float friction;
	int isManning; // 1 is manning, anything else is quadratic.
	int correctionStepsNum;
	float timestep;

	//Field
	int nx;
	int ny;
	float width; //x direction
	float length; // y direction
	float stillWaterElevation;

	//Boundaries
	BoundarySetting westBoundary;
	BoundarySetting eastBoundary;
	BoundarySetting southBoundary;
	BoundarySetting northBoundary;
	
	//solitary waves
	static const int MAX_NUM_SOLITARY = 100;
	Soliton solitons[MAX_NUM_SOLITARY];
	int countOfSolitons;
	bool is_there_new_solitary_wave;
	
	//Log
	bool doLog;
	int logStep;
	std::string logType;
	static const int MAX_NUM_RANGE = 100;
	Range logRange[MAX_NUM_RANGE];
	int countOfRanges;

	bool saveBathymetry;

	static const int MAX_NUM_GAUGE = 1000;
	std::string gaugesFilename;
	Point logGauges[MAX_NUM_GAUGE];
	int countOfGauges;
	float max_positive_bathy, min_negative_bathy, min_bathy;
	GraphicsSetting graphics;

	InitSetting(){
		project_name = "NA";

		//Paths
		exePath = "NA";
		initCMLName = "NA";
		initWFileName = "NA";
		bathymetryFileName = "NA";
		westIrrWaveFileName = "NA";  eastIrrWaveFileName = "NA";
		southIrrWaveFileName = "NA"; northIrrWaveFileName = "NA";

		rereadBathy = true;

		logPath = "NA";
		
		//Model
		isBoussinesq  = true;
		epsilon = 5e-8f;
		theta = 2;
		friction = 0;
		isManning = 0; // 1 is manning, anything else is quadratic.
		correctionStepsNum = 0;
		timestep = 0.0001;

		//Field
		nx = 100;
		ny = 100;
		width = 100; //x direction
		length = 100; // y direction
		stillWaterElevation = 0;
		

		//Boundaries
		westBoundary = BoundarySetting ();
		eastBoundary = BoundarySetting ();
		southBoundary = BoundarySetting ();
		northBoundary = BoundarySetting ();
		
		//solitary waves
		countOfSolitons = 0;
		is_there_new_solitary_wave = false;
		
		//Log
		doLog = false;
		logStep = 100;
		logType = "NA";
		countOfRanges = 0;

		saveBathymetry = false;

		gaugesFilename = "NA";
		int countOfGauges = 0;
		max_positive_bathy = 0; min_negative_bathy = 0; min_bathy = 0;
	}

	~ InitSetting() {}
};


extern InitSetting & initSetting;
// global array of settings, 'null terminated'
extern Setting g_settings[];
// flag to communicate when the settings have changed
extern ResetType g_reset_type;


// helper functions
float GetSetting(const char * name);
inline float GetSetting(const std::string &s) { return GetSetting(s.c_str()); }
int GetIntSetting(const char *name);

void SetSetting(const char *name, float new_value);
void SetSettingMax(const char *name, float new_value);
void SetSettingMin(const char *name, float new_value);
inline void SetSettingD(const char *name, double val) { SetSetting(name, float(val)); }

#endif
