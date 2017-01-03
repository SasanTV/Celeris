/*
 * FILE:
 *   terrain_heightfield.cpp
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


#include "settings.hpp"
#include "terrain_heightfield.hpp"

#include <fstream>

#include <cmath>
#include <string>

namespace {

    int RoundDown(float x)
    {
        return int(std::floor(x));
    }

    int RoundUp(float x)
    {
        return int(std::ceil(x));
    }

    float Lerp(float x, float down, float up)
    {
        return down + x * (up - down);
    }


    // height calculation.

    // for efficiency we split this API into three separate calls:
    // InitHeight() called once at the beginning
    // InitRow(float y) called once per y-value
    // GetHeight(float x, float &B, float &dB_dx, float &dB_dy) called once per point being sampled.

    float L, W;
    const int NUM_OCTAVES = 8;



    float z0, dz0_dy;
    float C_W, C_D, dCW_dy, dCD_dy;
    float R_C, dRC_dy, R_L, R_R, dRL_dy, dRR_dy;
    float B_L, B_R, dBL_dy, dBR_dy;
    float B_base, dBbase_dy;


}



void decodeLine (std::string* tag, std::string* value, std::string line) {
    unsigned int i = 0;
    bool tagDone = false;
	while (i < line.length() && !tagDone)
    {
        if (line[i] == ' ' || line[i] == '\t' ){ 
            //go
            i++;
        }
        else if (line[i] == '['){
            while (i<line.length() && !tagDone){
                
				*tag += line[i];
                if (line[i] == ']') {
                    tagDone = true;
                    i++;
                    while (i<line.length() && line[i] == ' '){
                        //go
                        i++;
                    }
                    while (i<line.length()) {
                        *value += line[i];
                        i++;
                    }
                }
                i++;
            }
            break;
        }
        else {
            break;
        }
    }

    if (!tagDone) {
        tag[0] = "Invalid";
        value[0] = "Invalid";
    }
}



boost::scoped_array<TerrainEntry> g_terrain_heightfield;
boost::scoped_array<BottomEntry> g_bottom;
boost::scoped_array<float> bathymetry;
int bathymetry_nx;
int bathymetry_ny;
//CAUTION: bathymetry does not have ghost zones
void readBathymetry ()
{

	std::ifstream in ((initSetting.bathymetryFileName).c_str());
	//MISSING catch
	
	if (in){
		std::string line;
		std::string tag = "";
		std::string value = "";

		while (std::getline(in,line))
		{
			if (line[0] == '=')
			{
				break;
			}
			else if (line != ""){
				tag.clear();
				value.clear();
				decodeLine(&tag, &value, line);
				
				if (tag != "Invalid" && value != "Invalid" && value != "") {

					if (tag == "[length]") {
						//go
					}
					else if (tag == "[width]") {
						//go
					}
					else if (tag == "[nx]") {
						bathymetry_nx = atoi(value.c_str());
					}
					else if (tag == "[ny]") {
						bathymetry_ny = atoi(value.c_str());
					}
				}
			}
		}
	}
	else {
 		bathymetry_nx = GetIntSetting("mesh_size_x");
		bathymetry_ny = GetIntSetting("mesh_size_y");
	}
	bathymetry.reset(new float[bathymetry_nx * bathymetry_ny]);

	initSetting.min_bathy = 1e10; // initialized to a very large float!!!
	for (int j = 0; j < bathymetry_ny; ++j) {
		for (int i = 0; i < bathymetry_nx; ++i) {
			float z = 0;
			if (in) {in >> z;}
			static bool min_max_bathy_init = false;
			if (!min_max_bathy_init){
				initSetting.max_positive_bathy = 0;
				initSetting.min_negative_bathy = 0;
				min_max_bathy_init = true;
			}
			bathymetry[j*bathymetry_nx + i] = z;
			initSetting.max_positive_bathy = std::max(initSetting.max_positive_bathy, z);
			initSetting.min_negative_bathy = std::min(initSetting.min_negative_bathy, z);
			initSetting.min_bathy = std::min(initSetting.min_bathy,z);
		}
	}

	float max_celerity = sqrt((initSetting.stillWaterElevation - initSetting.min_bathy)*GetSetting("gravity"));
	const float nx = GetIntSetting("mesh_size_x");
	const float ny = GetIntSetting("mesh_size_y");
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
	const float dx = W / (nx-1);
	const float dy = L / (ny-1);
	SetSetting("nominal_cfl", initSetting.timestep / (std::min(dx,dy)/max_celerity));
}


// must be called after L, W, bathymetry_nx, and bathymetry_ny are set.
float getBathymetry (float x, float y)
{
    // convert x,y to [0,width-1], [0,height-1] scale
    x = (x + W/2) / W * (bathymetry_nx-1);
    y = y / L * (bathymetry_ny-1);

    if (x < 0) x = 0;
    if (x > bathymetry_nx-1) x = float(bathymetry_nx-1);
    if (y < 0) y = 0;
    if (y > bathymetry_ny-1) y = float(bathymetry_ny-1);

    
    const int xdown = RoundDown(x);
    const int xup = RoundUp(x);
    const int ydown = RoundDown(y);
    const int yup = RoundUp(y);

    const float xfrac = x - xdown;
    const float yfrac = y - ydown;

    const int pitch = bathymetry_nx;
    
    const float result_y_down = Lerp(xfrac,
                                     bathymetry[ydown * pitch + xdown],
                                     bathymetry[ydown * pitch + xup]);
    const float result_y_up = Lerp(xfrac,
                                   bathymetry[yup * pitch + xdown],
                                   bathymetry[yup * pitch + xup]);

    return Lerp(yfrac, result_y_down, result_y_up);
}
void UpdateTerrainHeightfield()
{
	if(!initSetting.rereadBathy) return;
    using namespace std;

    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
	L = GetSetting("valley_length");
    W = GetSetting("valley_width");
	const float dx = W / (nx-1);
	const float dy = L / (ny-1);
    const int pitch = nx+4;
    
    g_terrain_heightfield.reset(new TerrainEntry[pitch * (ny+4)]);
    g_bottom.reset(new BottomEntry[pitch * (ny+4)]);
    
    // BX(i,j) = 0.5 * (B(i+1/2, j-1/2) + B(i+1/2, j+1/2))
    // BY(i,j) = 0.5 * (B(i-1/2, j+1/2) + B(i+1/2, j+1/2))
    // BA(i,j) = 0.25 * (B(i-1/2, j-1/2) + B(i-1/2, j+1/2) + B(i+1/2, j-1/2) + B(i+1/2, j+1/2))
	readBathymetry ();

    for (int j = 0; j < ny + 5; ++j) {
        const float y = (float(j)-2.5f) / float(ny-1) * L;    // j - 1/2

		for (int i = 0; i < nx + 5; ++i) {
            const float x = (float(i)-2.5f) / float(nx-1) * W - (W/2);   // i - 1/2

            float z = getBathymetry (x,y);
            //GetHeight(x, h); // change this to getheight (x,y,h)
			//h = 2.5 - 0.01f * y;
			
            if (i < nx+4) {   // (i)F
                if (j < ny+4) {   // (i,j)
                    g_bottom[j*pitch + i].BA = 0.25f * z;
                }
                if (j > 0) {   // (i,j-1)
                    g_bottom[(j-1)*pitch + i].BA += 0.25f * z;
                    g_bottom[(j-1)*pitch + i].BY = 0.5f * z;
                }
            }

            if (i > 0) {  // (i-1)
                if (j < ny+4) {   // (i-1,j)
                    g_bottom[j*pitch + i-1].BA += 0.25f * z;
                    g_bottom[j*pitch + i-1].BX = 0.5f * z;
                }

                if (j > 0) {  // (i-1,j-1)
                    g_bottom[(j-1)*pitch + i-1].BA += 0.25f * z;
                    g_bottom[(j-1)*pitch + i-1].BX += 0.5f * z;
                    g_bottom[(j-1)*pitch + i-1].BY += 0.5f * z;
                }
            }
        }
    }

    // We now use BA instead of B in g_terrain_heightfield.
    // This prevents water "showing through" in steep areas.
    for (int j = 0; j < ny + 4; ++j) {
        for (int i = 0; i < nx + 4; ++i) {
			g_terrain_heightfield[j * pitch + i].B = g_bottom[j * pitch + i].BA + sqrt(sqrt(initSetting.epsilon)); // here I add epsilon, to avoid showing water through land
			g_terrain_heightfield[j * pitch + i].dBdx = (g_bottom[j * pitch + i].BX - g_bottom[j * pitch + i].BA)/(2*dx);
			g_terrain_heightfield[j * pitch + i].dBdy = (g_bottom[j * pitch + i].BY - g_bottom[j * pitch + i].BA)/(2*dy);
        }
    }

	initSetting.rereadBathy = false;
}

float GetTerrainHeight(float x, float y)
{
    const int width = GetIntSetting("mesh_size_x");
    const int height = GetIntSetting("mesh_size_y");    
    
    const float L = GetSetting("valley_length");
    const float W = GetSetting("valley_width");

    // convert x,y to [0,width-1], [0,height-1] scale
    x = (x + W/2) / W * (width-1);
    y = y / L * (height-1);

    if (x < 0) x = 0;
    if (x > width-1) x = float(width-1);
    if (y < 0) y = 0;
    if (y > height-1) y = float(height-1);

    // now add 2 to account for ghost cells
    x += 2;
    y += 2;
    
    const int xdown = RoundDown(x);
    const int xup = RoundUp(x);
    const int ydown = RoundDown(y);
    const int yup = RoundUp(y);

    const float xfrac = x - xdown;
    const float yfrac = y - ydown;

    const int pitch = width + 4;
    
    const float result_y_down = Lerp(xfrac,
                                     g_terrain_heightfield[ydown * pitch + xdown].B,
                                     g_terrain_heightfield[ydown * pitch + xup].B);
    const float result_y_up = Lerp(xfrac,
                                   g_terrain_heightfield[yup * pitch + xdown].B,
                                   g_terrain_heightfield[yup * pitch + xup].B);

    return Lerp(yfrac, result_y_down, result_y_up);
}
