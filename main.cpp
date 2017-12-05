/*
 * FILE:
 *   main.cpp
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

#include "engine.hpp"
#include "gui_manager.hpp"
#include "settings.hpp"
#include "terrain_heightfield.hpp"

#include "coercri/dx11/gfx/dx11_gfx_driver.hpp"
#include "coercri/dx11/gfx/dx11_window.hpp"
#include "coercri/gfx/gfx_context.hpp"
#include "coercri/gfx/window_listener.hpp"
#include "coercri/timer/generic_timer.hpp"

#include "tinyxml_2_6_2/tinyxml/tinyxml.h"
#include "tinyxml_2_6_2/tinyxml/tinystr.h"



   
#include "boost/scoped_ptr.hpp"

#include <cmath>
#include <fstream>

// Following header and library are used to handle relative and absolute paths.
#pragma comment (lib, "Shlwapi.lib")
#include "Shlwapi.h"


const int GUI_WIDTH = 400;
int g_width = 500, g_height = 500;
bool g_quit = false;
bool g_resize = false;


// There is a little trick that you see below! The story behind it is interesting and I, Sasan, am going to share it here.
// InitSetting started as a small structure to encapsulate a few file paths. Then it grew and grew bigger.
// At some point, it was pretty large to keep it on the stack, so I decided to put it on the heap (new InitSetting()).
// But I did not wanted to change all initSetting.foo to initSetting->foo. So, that is the reason behind the below trick :)
InitSetting * p_initSetting = new InitSetting();
InitSetting& initSetting = *p_initSetting;

namespace {
    const float PI = std::atan(1.0f) * 4.0f;
	
	std::string ExePath() {
		char buffer[MAX_PATH];
		GetModuleFileName( NULL, buffer, MAX_PATH );
		std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
		return std::string( buffer ).substr( 0, pos);
	}

	std::string extractPath(std::string path) {
		std::string::size_type pos = path.find_last_of( "\\/" );
		return path.substr( 0, pos);
	}
	
    void UpdateMousePick(int mx, int my, ShallowWaterEngine &engine)
    {
        float wx, wy, wz, d, u, v;
        if (mx >= 0 && mx < g_width && my >= 0 && my < g_height
        && engine.mousePick(mx, my, wx, wy, wz, d, u, v)) {
            SetSetting("x", wx);
            SetSetting("y", wy);
            SetSetting("z", wz);
            SetSetting("depth", d);
            SetSetting("u", u);
            SetSetting("v", v);

            const float spd = std::sqrt(u*u + v*v);
            const float c = std::sqrt(GetSetting("gravity") * d);
            SetSetting("froude", spd / c);
            
        } else {
            SetSetting("x", -1);
            SetSetting("y", -1);
            SetSetting("z", -1);
            SetSetting("depth", -1);
            SetSetting("u", -1);
            SetSetting("v", -1);
            SetSetting("froude", -1);
        }
    }
}

class MyListener : public Coercri::WindowListener {
public:
    MyListener(ShallowWaterEngine &e, GuiManager &g, Coercri::Window &w) : 
        engine(e), gui_manager(g), window(w),
        fwd(false), bwd(false), left(false), right(false), up(false), down(false),
        shift(false),
        cam_x(0), cam_y(0), cam_z(30),
        pitch(0), yaw(0),
        mx(-999), my(-999), left_mouse_down(false), right_mouse_down(false)
    { }

    void onLoseFocus() { fwd = bwd = left = right = up = down = shift = false; left_mouse_down = right_mouse_down = false; }
    
    void onClose() { g_quit = true; }
    void onResize(int w, int h) {
        const int gw = gui_manager.isGuiShown() ? GUI_WIDTH : 0;
        g_width = std::max(0, w - gw);
        g_height = h;
        g_resize = true;
    }

    void onRawKey(bool pressed, Coercri::RawKey rk)
    {
        switch (rk) {
        case Coercri::RK_W: fwd = pressed; break;
        case Coercri::RK_S: bwd = pressed; break;
        case Coercri::RK_A: left = pressed; break;
        case Coercri::RK_D: right = pressed; break;
        case Coercri::RK_PAGE_UP: up = pressed; break;
        case Coercri::RK_PAGE_DOWN: down = pressed; break;
		case Coercri::RK_P: gui_manager.switchPause(); break;
        case Coercri::RK_LEFT_SHIFT: case Coercri::RK_RIGHT_SHIFT: shift = pressed; break;
        }
    }

    void onMouseDown(int x, int y, Coercri::MouseButton mb)
    {
        if (mb == Coercri::MB_RIGHT) {
            mx = x;
            my = y;
            right_mouse_down = true;
            window.captureMouse(true);
			

        } else if (mb == Coercri::MB_LEFT) {
            left_mouse_down = true;
        }
    }

    void onMouseUp(int x, int y, Coercri::MouseButton mb)
    {
        if (mb == Coercri::MB_RIGHT) {
            right_mouse_down = false;
            window.captureMouse(false);
        } else if (mb == Coercri::MB_LEFT) {
            left_mouse_down = false;
        }
    }

    void onMouseMove(int new_x, int new_y) 
    {
        if (right_mouse_down) {
            const int dx = new_x - mx;
            const int dy = new_y - my;

            const float angle_speed = 0.001f;

            yaw += dx * angle_speed;
            pitch -= dy * angle_speed;

            if (pitch > PI/2) pitch = PI/2;
            if (pitch < -PI/2) pitch = -PI/2;
        }

        mx = new_x;
        my = new_y;
    }

    void update(float dt)
    {
        // HACK to communicate camera resets back from the gui.
        gui_manager.getCameraReset(cam_x, cam_y, cam_z, pitch, yaw);
        

        const float W = GetSetting("valley_width");
        const float L = GetSetting("valley_length");
		
        // move the camera
		const float speed_ = std::max(W,L)/10.0f;
        float speed = shift ? 10 * speed_ : speed_;
        speed *= dt;



        if (fwd) {
            cam_x += sin(yaw) * cos(pitch) * speed;
            cam_y += cos(yaw) * cos(pitch) * speed;
            cam_z += sin(pitch) * speed;
        }
        if (bwd) {
            cam_x -= sin(yaw) * speed;
            cam_y -= cos(yaw) * speed;
            cam_z -= sin(pitch) * speed;
        }
        if (left) {
            cam_x -= cos(yaw) * speed;
            cam_y += sin(yaw) * speed;
        }
        if (right) {
            cam_x += cos(yaw) * speed;
            cam_y -= sin(yaw) * speed;
        }
        if (up) {
            cam_z += cos(pitch) * speed;
            cam_x -= sin(pitch) * sin(yaw) * speed;
            cam_y -= sin(pitch) * cos(yaw) * speed;
        }
        if (down) {
            cam_z -= cos(pitch) * speed;
            cam_x += sin(pitch) * sin(yaw) * speed;
            cam_y += sin(pitch) * cos(yaw) * speed;
        }

        // clip camera position
        const float max_ratio = 1.0f;
        float max_clip = max_ratio*std::max(W,L);
        //if (GetIntSetting("clip_camera")) max_clip = 0;
        if (cam_x < -W/2 -max_clip) cam_x = -W/2 - max_clip;
        if (cam_x > W/2 + max_clip) cam_x = W/2 + max_clip;
        if (cam_y < -max_clip) cam_y = -max_clip;
        if (cam_y > L + max_clip) cam_y = L + max_clip;

		const float t_height = (GetTerrainHeight(cam_x, cam_y) + engine.near_plane_dist) * GetSetting("Vertical_Scale");
		cam_z = std::max(t_height, cam_z);
        const float w_height = engine.getWaterHeight(std::min(std::max(-W/2, cam_x), W/2),
                                                   std::min(std::max(0.0f, cam_y), L));
		cam_z = std::min(std::max(W,L) + w_height, cam_z);        
        
        // TODO: only need to call this if the camera actually moved
        // (or the window resized...)
        engine.moveCamera(cam_x, cam_y, cam_z, yaw, pitch, g_width, g_height);

        // do left mouse interaction
        if (left_mouse_down) {
            const float wy = GetSetting("y");
            if (wy >= 0) {
                const float wx = GetSetting("x");
                const float wy = GetSetting("y");
                engine.applyMouseShader(wx, wy, dt);
            }
        }
    }

    int getMX() const { return mx; }
    int getMY() const { return my; }
	void maximizeWindow() {
		window.maximizeWindow();
	}

private:

    ShallowWaterEngine &engine;
    GuiManager &gui_manager;
    Coercri::Window &window;

    bool fwd, bwd, left, right, up, down;
    bool shift;
        
    float cam_x, cam_y, cam_z;
    float pitch, yaw;

    int mx, my;
    bool left_mouse_down, right_mouse_down;
};

void reset_uniforrm_time_series_iters(){
	initSetting.westBoundary.uniformTimeSeries.reset_iter();
	initSetting.eastBoundary.uniformTimeSeries.reset_iter();
	initSetting.southBoundary.uniformTimeSeries.reset_iter();
	initSetting.northBoundary.uniformTimeSeries.reset_iter();
}


int real_main()
{
    // Setup Direct3D

#ifdef _DEBUG
    const bool debug = true;
#else
    const bool debug = false;
#endif
    
    // using texture lookups in the vertex shader, therefore must have feature level >= 10.0.

    std::auto_ptr<Coercri::DX11GfxDriver> gfx_driver(new Coercri::DX11GfxDriver(D3D_DRIVER_TYPE_HARDWARE,
                                                                                D3D11_CREATE_DEVICE_SINGLETHREADED | (debug ? D3D11_CREATE_DEVICE_DEBUG : 0),
                                                                                D3D_FEATURE_LEVEL_10_0));
    boost::shared_ptr<Coercri::Timer> timer(new Coercri::GenericTimer);
#include "coercri/timer/generic_timer.hpp"

    // set up the gui
    std::auto_ptr<MyListener> listener;
    boost::shared_ptr<Coercri::DX11Window> window = 
        boost::static_pointer_cast<Coercri::DX11Window>(
            gfx_driver->createWindow(g_width + GUI_WIDTH, g_height, true, false, initSetting.project_name + " - Celeris Advent (v1.3.3)"));
    GuiManager gui_manager(window, timer, GUI_WIDTH);
	
    // Create the ShallowWaterEngine
    boost::scoped_ptr<ShallowWaterEngine> engine(
        new ShallowWaterEngine(gfx_driver->getDevice(), gfx_driver->getDeviceContext()));

    listener.reset(new MyListener(*engine, gui_manager, *window));
    window->addWindowListener(listener.get());
	
	if (initSetting.graphics.maximized){
		listener->maximizeWindow();
	}
	
    unsigned int last_draw = timer->getMsec();
    bool is_gui_shown = true;
    
    g_resize = true; // make sure it "resizes" first time
	
	unsigned int reference_time = timer->getMsec();
    while (!g_quit) {

        // slight hack to make sure window resizes when gui is shown/hidden
        bool new_is_gui_shown = gui_manager.isGuiShown();
        if (new_is_gui_shown != is_gui_shown) {
            int win_w, win_h;
            window->getSize(win_w, win_h);
            listener->onResize(win_w, win_h);
            is_gui_shown = new_is_gui_shown;
        }
        
        if (g_resize) {
            gui_manager.resize();
            g_resize = false;
        }
		
        // see if updates need to be done
        switch (g_reset_type) {
        
		// R_MESH is currently used to reset simulation 
		case R_MESH:
            engine->remesh(g_reset_type);
			reference_time = timer->getMsec();
			reset_uniforrm_time_series_iters();
            break;
        case R_TERRAIN:
            engine->newTerrainSettings();
            break;
        };

        g_reset_type = R_NONE;        

		float elapsed_time = 0;
        engine->resetTimestep(0.001f,elapsed_time);
        int timestep_count = 0;
        unsigned int timer_at_zero_timesteps = timer->getMsec();
		
        
        while (!g_resize && g_reset_type == R_NONE && !g_quit && gui_manager.isGuiShown() == is_gui_shown) {

            const unsigned int frame_time = 1;  // in msec. acts as fps limiter.

            // see if dt needs to be reset
            const int steps_between_reset = 10;
            if (timestep_count >= steps_between_reset) {
				
                unsigned int timer_now = timer->getMsec();
                float dt = float(timer_now - timer_at_zero_timesteps) / (1000.0f * steps_between_reset);
				elapsed_time = float(timer_now-reference_time)/(1000.0f);
				timestep_count = 0;
                timer_at_zero_timesteps = timer->getMsec();
                engine->resetTimestep(dt,elapsed_time);
				UpdateMousePick(listener->getMX(), listener->getMY(), *engine);
     
            }
	
            
            gfx_driver->pollEvents();
            gui_manager.logic();
		
            if (g_reset_type != R_NONE) break;  // don't continue if the settings are out of date
			
            const unsigned int time_now = timer->getMsec();
            const int time_since_last = int(time_now - last_draw);
            
            bool do_render = true; //int(time_now - last_draw) >= frame_time;
			if (window->needsRepaint()) {
				do_render = true;
			}
			
            
            if (do_render) {

                last_draw = time_now;

                // do update
                listener->update(float(time_since_last) / 1000.0f);
				if (!gui_manager.isPause()){
					// do timestep
					const int timesteps_per_frame = GetIntSetting("timesteps_per_frame");
					for (int i = 0; i < timesteps_per_frame; ++i) {
						engine->timestep();
						++timestep_count;
					}
				}
				engine->afterTimestep();

                // clear the screen
                const float rgba[] = { 0, 0, 0, 1 };
                gfx_driver->getDeviceContext()->ClearRenderTargetView(window->getRenderTargetView(), rgba);

                // draw the landscape & water
                engine->render(window->getRenderTargetView());

                // draw the gui (using coercri 2D routines)
                {
                    std::auto_ptr<Coercri::GfxContext> gc = window->createGfxContext();
					gui_manager.draw(*gc);
                }

                window->cancelInvalidRegion();
            
            } else {
                const int time_to_next_frame = int(last_draw + frame_time - time_now);
                timer->sleepMsec(std::max(1, time_to_next_frame + 1));
            }
        }
    }

    return 0;
}

bool readInputCML();
bool browseInputCML(){

	OPENFILENAME ofn;       // common dialog box structure
	char szFile[260];       // buffer for file name
	// Initialize OPENFILENAME
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	//ofn.hwndOwner = hwnd;
	ofn.lpstrFile = szFile;
	// Set lpstrFile[0] to '\0' so that GetOpenFileName does not 
	// use the contents of szFile to initialize itself.
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = sizeof(szFile);
	ofn.lpstrTitle = "Open Input File";
	ofn.lpstrFilter = "Celeris input file (*.cml)\0*.CML\0All types (*.*)\0*.*\0\0";
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = NULL;
	ofn.nMaxFileTitle = 0;
	ofn.lpstrInitialDir = NULL;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	// Display the Open dialog box. 
	if (GetOpenFileName(&ofn) ==  TRUE){
		initSetting.initCMLName = ofn.lpstrFile;
		return readInputCML();
	} else {
		return 1;
	}
/*
	std::ofstream fileOut;
	fileOut.open ((initSetting.exePath + "setting.init").c_str(),std::ios::out);
	GetOpenFileName(&ofn);
	fileOut << ofn.lpstrFile << "\n";
	GetOpenFileName(&ofn);
	fileOut << ofn.lpstrFile;
	fileOut.close();
*/

}


void readGraphicsInput_helper(TiXmlElement* root);
bool readGraphicsInput();
bool readInputCML()
{
	if(readGraphicsInput())
	{
		return 1;
	}

	TiXmlDocument doc(initSetting.initCMLName);
	bool loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		TiXmlElement* root = doc.FirstChildElement();
		if(root == NULL)
		{
			doc.Clear();
			MessageBox(0, "Invalid CML input file!", "Error!", MB_OK);
			return browseInputCML();
		}
		for(TiXmlElement* elem = root->FirstChildElement(); elem != NULL; elem = elem->NextSiblingElement())
		{
			std::string elemName = elem->Value();
			if(elemName == "fieldDimensions")
			{
				
				float length,width,stillWaterElevation;
				elem->QueryFloatAttribute("length", &length);
				elem->QueryFloatAttribute("width", &width);
				elem->QueryFloatAttribute("stillWaterElevation", &stillWaterElevation);
				initSetting.length = length; //y direction
				initSetting.width  = width; //x direction
				initSetting.stillWaterElevation = stillWaterElevation;
				
			}
			else if(elemName == "name"){
				initSetting.project_name = elem->FirstChild()->ToText()->Value();
			}
			else if(elemName == "model")
			{
				std::string modelType;
				elem->QueryStringAttribute("type", &modelType);
				if (modelType == "BSNQ")
				{
					initSetting.isBoussinesq = true	;
				}
				else if (modelType == "NLSW")
				{
					initSetting.isBoussinesq = false;
				}
				else
				{
					//Catch error or set default
					initSetting.isBoussinesq = true;
				}

				// set shared parameters between models
				float tempFloat;
				elem->FirstChildElement("parameters")->QueryFloatAttribute("epsilon", &tempFloat);
				initSetting.epsilon = tempFloat;

				elem->FirstChildElement("parameters")->QueryFloatAttribute("Theta", &tempFloat);
				initSetting.theta = tempFloat;
		
				
				int tempInt;
				elem->FirstChildElement("parameters")->QueryIntAttribute("correctionStepsNum", &tempInt);
				initSetting.correctionStepsNum = tempInt;

				elem->FirstChildElement("parameters")->QueryFloatAttribute("timestep", &tempFloat);
				initSetting.timestep = tempFloat;

				bool tempBool = false;
				elem->FirstChildElement("parameters")->QueryBoolAttribute("adaptive", &tempBool);
				initSetting.time_scheme = tempBool? predictor_adaptive : predictor;

				std::string frictionType;
				elem->FirstChildElement("friction")->QueryStringAttribute("type", &frictionType);
				initSetting.isManning = (frictionType == "Manning")?1:0;

				elem->FirstChildElement("friction")->QueryFloatAttribute("coef", &tempFloat);
				initSetting.friction = tempFloat;
			}
			else if(elemName == "gridSize")
			{
				int nx,ny;
				elem->QueryIntAttribute("nx", &nx);
				elem->QueryIntAttribute("ny", &ny);
				initSetting.nx = nx;
				initSetting.ny = ny;
			}
			else if(elemName == "bathymetryFilePath")
			{
				std::string tempString = elem->FirstChild()->ToText()->Value();
				const char *pp;
				pp =  tempString.c_str();
				if (PathIsRelative(pp)){
					initSetting.bathymetryFileName = extractPath(initSetting.initCMLName) + "/" + tempString;
				} else {
					initSetting.bathymetryFileName = tempString;
				}
			}
			else if(elemName == "tideSurgeSLR")
			{
				float tideMin = -1.0f;
				float tideMax = +1.0f;
				float tideSet =  0.0f;
				bool autoTide = true;
				elem->QueryBoolAttribute("auto", &autoTide);				
				
				if (!autoTide){
					elem->QueryFloatAttribute("min", &tideMin);
					elem->QueryFloatAttribute("max", &tideMax);
					elem->QueryFloatAttribute("set", &tideSet);
				}
				initSetting.tideSurgeSLR.autoValue = autoTide;
				initSetting.tideSurgeSLR.minValue = tideMax;
				initSetting.tideSurgeSLR.minValue = tideMin;
				initSetting.tideSurgeSLR.setValue = tideSet;
			}
			else if(elemName == "hotStartFilePath")
			{
				std::string tempString = elem->FirstChild()->ToText()->Value();
				const char *pp;
				pp =  tempString.c_str();
				if (PathIsRelative(pp)){
					initSetting.initWFileName = extractPath(initSetting.initCMLName) + "/" + tempString;
				} else {
					initSetting.initWFileName = tempString;
				}
			}
			else if(elemName == "solitaryWave")
			{
				
				static int counterOfSolitons = 0;

				if (counterOfSolitons <= initSetting.MAX_NUM_SOLITARY){
					float tempFloat[5];
					int status;
					elem->QueryFloatAttribute("H", &tempFloat[0]);
					elem->QueryFloatAttribute("theta", &tempFloat[1]);
					elem->QueryFloatAttribute("xc", &tempFloat[2]);
					elem->QueryFloatAttribute("yc", &tempFloat[3]);
					status = elem->QueryFloatAttribute("length", &tempFloat[4]);
					if (status == TIXML_NO_ATTRIBUTE){
						tempFloat[4] = 1e30f; // a very large float! TODO: set this number to a reasonably large float.
					}

					initSetting.solitons[counterOfSolitons] = Soliton (tempFloat[0], PI/180.0f * tempFloat[1], tempFloat[2], tempFloat[3], tempFloat[4]);
					++counterOfSolitons;
				}
				initSetting.countOfSolitons = counterOfSolitons;
			}
			else if(elemName == "westBoundary")
			{
				std::string temp;
				elem->QueryStringAttribute("type", &temp);
				initSetting.westBoundary.type = temp;


				float tempFloat;
				elem->QueryFloatAttribute("seaLevel", &tempFloat);
				initSetting.westBoundary.waterLevel = tempFloat;
				elem->QueryFloatAttribute("widthNum", &tempFloat);
				initSetting.westBoundary.width = tempFloat;


				if (initSetting.westBoundary.type == "SineWave") 
				{
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("amplitude", &tempFloat);
					initSetting.westBoundary.sineWaveSetting.amplitude = tempFloat;
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("period", &tempFloat);
					initSetting.westBoundary.sineWaveSetting.period = tempFloat;
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("theta", &tempFloat);
					initSetting.westBoundary.sineWaveSetting.theta = PI/180.0 * tempFloat;
				} else if (initSetting.westBoundary.type == "IrregularWaves") {
					std::string tempString = elem->FirstChildElement("filePath")->FirstChild()->ToText()->Value();
					const char *pp;
					pp =  tempString.c_str();
					if (PathIsRelative(pp)){
						initSetting.westIrrWaveFileName = extractPath(initSetting.initCMLName) + "/" + tempString;
					} else {
						initSetting.westIrrWaveFileName = tempString;
					}
				} else if (initSetting.westBoundary.type == "UniformTimeSeries") {
					std::string tempString = elem->FirstChildElement("filePath")->FirstChild()->ToText()->Value();
					const char *pp;
					pp =  tempString.c_str();
					if (PathIsRelative(pp)){
						initSetting.westBoundary.uniformTimeSeries.fileName = extractPath(initSetting.initCMLName) + "/" + tempString;
					} else {
						initSetting.westBoundary.uniformTimeSeries.fileName = tempString;
					}
				}

				initSetting.westBoundary.hasChanged = true;

			}
			else if(elemName == "eastBoundary")
			{
				std::string temp;
				elem->QueryStringAttribute("type", &temp);
				initSetting.eastBoundary.type = temp;

				float tempFloat;
				elem->QueryFloatAttribute("seaLevel", &tempFloat);
				initSetting.eastBoundary.waterLevel = tempFloat;
				elem->QueryFloatAttribute("widthNum", &tempFloat);
				initSetting.eastBoundary.width = tempFloat;

				if (initSetting.eastBoundary.type == "SineWave") 
				{
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("amplitude", &tempFloat);
					initSetting.eastBoundary.sineWaveSetting.amplitude = tempFloat;
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("period", &tempFloat);
					initSetting.eastBoundary.sineWaveSetting.period = tempFloat;
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("theta", &tempFloat);
					initSetting.eastBoundary.sineWaveSetting.theta = PI/180.0 * tempFloat;
				} else if (initSetting.eastBoundary.type == "IrregularWaves") {
					std::string tempString = elem->FirstChildElement("filePath")->FirstChild()->ToText()->Value();
					const char *pp;
					pp =  tempString.c_str();
					if (PathIsRelative(pp)){
						initSetting.eastIrrWaveFileName = extractPath(initSetting.initCMLName) + "/" + tempString;
					} else {
						initSetting.eastIrrWaveFileName = tempString;
					}
				} else if (initSetting.eastBoundary.type == "UniformTimeSeries") {
					std::string tempString = elem->FirstChildElement("filePath")->FirstChild()->ToText()->Value();
					const char *pp;
					pp =  tempString.c_str();
					if (PathIsRelative(pp)){
						initSetting.eastBoundary.uniformTimeSeries.fileName = extractPath(initSetting.initCMLName) + "/" + tempString;
					} else {
						initSetting.eastBoundary.uniformTimeSeries.fileName = tempString;
					}
				}

				initSetting.eastBoundary.hasChanged = true;
			}
			else if(elemName == "northBoundary")
			{
				std::string temp;
				elem->QueryStringAttribute("type", &temp);
				initSetting.northBoundary.type = temp;
				
				float tempFloat;
				elem->QueryFloatAttribute("seaLevel", &tempFloat);
				initSetting.northBoundary.waterLevel = tempFloat;
				elem->QueryFloatAttribute("widthNum", &tempFloat);
				initSetting.northBoundary.width = tempFloat;

				if (initSetting.northBoundary.type == "SineWave") 
				{
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("amplitude", &tempFloat);
					initSetting.northBoundary.sineWaveSetting.amplitude = tempFloat;
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("period", &tempFloat);
					initSetting.northBoundary.sineWaveSetting.period = tempFloat;
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("theta", &tempFloat);
					initSetting.northBoundary.sineWaveSetting.theta = PI/180.0 * tempFloat;
				} else if (initSetting.northBoundary.type == "IrregularWaves") {
					std::string tempString = elem->FirstChildElement("filePath")->FirstChild()->ToText()->Value();
					const char *pp;
					pp =  tempString.c_str();
					if (PathIsRelative(pp)){
						initSetting.northIrrWaveFileName = extractPath(initSetting.initCMLName) + "/" + tempString;
					} else {
						initSetting.northIrrWaveFileName = tempString;
					}
				} else if (initSetting.northBoundary.type == "UniformTimeSeries") {
					std::string tempString = elem->FirstChildElement("filePath")->FirstChild()->ToText()->Value();
					const char *pp;
					pp =  tempString.c_str();
					if (PathIsRelative(pp)){
						initSetting.northBoundary.uniformTimeSeries.fileName = extractPath(initSetting.initCMLName) + "/" + tempString;
					} else {
						initSetting.northBoundary.uniformTimeSeries.fileName = tempString;
					}
				}

				initSetting.northBoundary.hasChanged = true;
			}
			else if(elemName == "southBoundary")
			{
				std::string temp;
				elem->QueryStringAttribute("type", &temp);
				initSetting.southBoundary.type = temp;
				float tempFloat;
				elem->QueryFloatAttribute("seaLevel", &tempFloat);
				initSetting.southBoundary.waterLevel = tempFloat;
				elem->QueryFloatAttribute("widthNum", &tempFloat);
				initSetting.southBoundary.width = tempFloat;

				if (initSetting.southBoundary.type == "SineWave") 
				{
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("amplitude", &tempFloat);
					initSetting.southBoundary.sineWaveSetting.amplitude = tempFloat;
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("period", &tempFloat);
					initSetting.southBoundary.sineWaveSetting.period = tempFloat;
					elem->FirstChildElement("sineWave")->QueryFloatAttribute("theta", &tempFloat);
					initSetting.southBoundary.sineWaveSetting.theta = PI/180.0 * tempFloat;
				} else if (initSetting.southBoundary.type == "IrregularWaves") {
					std::string tempString = elem->FirstChildElement("filePath")->FirstChild()->ToText()->Value();
					const char *pp;
					pp =  tempString.c_str();
					if (PathIsRelative(pp)){
						initSetting.southIrrWaveFileName = extractPath(initSetting.initCMLName) + "/" + tempString;
					} else {
						initSetting.southIrrWaveFileName = tempString;
					}
				} else if (initSetting.southBoundary.type == "UniformTimeSeries") {
					std::string tempString = elem->FirstChildElement("filePath")->FirstChild()->ToText()->Value();
					const char *pp;
					pp =  tempString.c_str();
					if (PathIsRelative(pp)){
						initSetting.southBoundary.uniformTimeSeries.fileName = extractPath(initSetting.initCMLName) + "/" + tempString;
					} else {
						initSetting.southBoundary.uniformTimeSeries.fileName = tempString;
					}
				}

				initSetting.southBoundary.hasChanged = true;
			}
			else if(elemName == "logData")
			{
				bool doLog;
				elem->QueryBoolAttribute("doLog", &doLog);
				initSetting.doLog = doLog;

				int logStep;
				elem->QueryIntAttribute("logStep", &logStep);
				initSetting.logStep = logStep;
				int rangeCounter = 0;
				int gaugeCounter = 0;

				for(TiXmlElement* elem2 = elem->FirstChildElement(); elem2 != NULL; elem2 = elem2->NextSiblingElement()){
					std::string elem2Name = elem2->Value() ;
					
					if(elem2Name == "logPath"){
					//initSetting.logPath = elem->FirstChildElement("logPath")->FirstChild()->ToText()->Value();
						TiXmlNode* tempChild = elem2->FirstChild();
						std::string tempString = tempChild? tempChild->ToText()->Value() : "";
						const char *pp;
						pp =  tempString.c_str();
						if (PathIsRelative(pp)){
							initSetting.logPath = extractPath(initSetting.initCMLName) + "/" + tempString;
						} else {
							initSetting.logPath = tempString;
						}
					}
					if (elem2Name == "timeAxis") {
						std::string filename = "";
						elem2->QueryStringAttribute("filename", &filename);
						initSetting.timeAxisFilename = filename;
					}
					else if (elem2Name == "range"){
						int x1,y1,x2,y2;
						std::string filename = "";
						elem2->QueryStringAttribute("filename", &filename);

						elem2->FirstChildElement("bottomLeft")->QueryIntAttribute("x", &x1);
						elem2->FirstChildElement("bottomLeft")->QueryIntAttribute("y", &y1);
						elem2->FirstChildElement("topRight")->QueryIntAttribute("x", &x2);
						elem2->FirstChildElement("topRight")->QueryIntAttribute("y", &y2);
						
						initSetting.logRange[rangeCounter] = Range(filename, Point(x1,y1),Point(x2,y2));
						++rangeCounter;
					}
					else if (elem2Name == "gauges"){
						std::string data;
						std::string filename = "";
						elem2->QueryStringAttribute("filename", &filename);
						initSetting.gaugesFilename = filename;
						data = elem2->FirstChild()->ToText()->Value();
						std::stringstream ss(data);
						char comma;
						int x,y;
						while (ss){
							ss >> x >> comma >> y >> comma;
							initSetting.logGauges[gaugeCounter] = Point(x,y);
							++gaugeCounter;
						}
					}
				}
				initSetting.countOfRanges = rangeCounter;
				initSetting.countOfGauges = gaugeCounter;

			}
			else if(elemName == "graphics")
			{
				readGraphicsInput_helper(elem);
			}
		}
		
		if(initSetting.graphics.autoCam){
			initSetting.graphics.camera.fov = 50.0f;
			initSetting.graphics.camera.x = - initSetting.width / 2;
			initSetting.graphics.camera.y = 0;
			initSetting.graphics.camera.z = sqrt( initSetting.width*initSetting.width + initSetting.length*initSetting.length); //sqrt(initSetting.width * initSetting.length) * 1.2;
			initSetting.graphics.camera.pitch = -1;
			initSetting.graphics.camera.yaw = atan(initSetting.width/initSetting.length);
		}
		return 0;
	}
	else
	{
		//MessageBox(0, doc.ErrorDesc(), "CML input file not found!", MB_OK);
		MessageBox(0, doc.ErrorDesc(), "CML input error!", MB_OK);
		return browseInputCML();
	}

}



void readGraphicsInput_helper(TiXmlElement* root){

	for(TiXmlElement* elem = root->FirstChildElement(); elem != NULL; elem = elem->NextSiblingElement())
	{
		std::string elemName = elem->Value();
		if(elemName == "grid")
		{
			float gridScale;
			bool gridOn;
			elem->QueryBoolAttribute("show", &gridOn);
			elem->QueryFloatAttribute("scale", &gridScale);

			initSetting.graphics.gridOn = gridOn;
			initSetting.graphics.gridScale = gridScale;
		}
		if (elemName == "window")
		{
			bool maximized;
			elem->QueryBoolAttribute("maximized", &maximized);

			initSetting.graphics.maximized = maximized;
		}
		else if(elemName == "surfaceShading")
		{
			int type = 0;
			elem->QueryIntAttribute("type", &type);
			initSetting.graphics.surfaceShading.type = static_cast<ShadingOptions>(type);
			for(TiXmlElement* elem2 = elem->FirstChildElement(); elem2 != NULL; elem2 = elem2->NextSiblingElement()){
				std::string elem2Name = elem2->Value() ;
				if(elem2Name == "colormap"){
					float colormapMin = -1.0f;
					float colormapMax = +1.0f;
					bool autoColormap = true;

					elem2->QueryBoolAttribute("auto", &autoColormap);
					initSetting.graphics.surfaceShading.autoColormap = autoColormap;

					elem2->QueryFloatAttribute("min", &colormapMin);
					initSetting.graphics.surfaceShading.colormapMin = colormapMin;						

					elem2->QueryFloatAttribute("max", &colormapMax);
					initSetting.graphics.surfaceShading.colormapMax = colormapMax;						
				} else if (elem2Name == "shadingVariable"){
					int value = 0;

					elem2->QueryIntAttribute("value", &value);
					initSetting.graphics.surfaceShading.shadingVariable = value;
				}
				else if (elem2Name == "dissipationIntensity") {
					bool show = 0;

					elem2->QueryBoolAttribute("show", &show);
					initSetting.graphics.surfaceShading.showDissipation = show;

					float disspationThreshold = .25f;
					float whiteWaterDecay = 0.1f;

					elem2->QueryFloatAttribute("threshold", &disspationThreshold);
					initSetting.graphics.surfaceShading.dissipationThreshold = disspationThreshold;
					elem2->QueryFloatAttribute("decay", &whiteWaterDecay);
					initSetting.graphics.surfaceShading.whiteWaterDecay = whiteWaterDecay;
				} else if (elem2Name == "drylandDepthOfInundation"){
					
					bool autoInudationDepth = true;
					elem2->QueryBoolAttribute("auto", &autoInudationDepth);
					initSetting.graphics.surfaceShading.autoInundationDepth = autoInudationDepth;

					bool showInundatedArea = false;
					elem2->QueryBoolAttribute("show", &showInundatedArea);
					initSetting.graphics.surfaceShading.showInundatedArea = showInundatedArea;

					float value = 0.01;
					float maxInundation = 1;

					elem2->QueryFloatAttribute("value", &value);
					initSetting.graphics.surfaceShading.drylandDepthOfInundation = value;

					elem2->QueryFloatAttribute("max", &maxInundation);
					initSetting.graphics.surfaceShading.maxInundation = maxInundation;
				}
			}
		}
		else if(elemName == "terrainTexture")
		{
			int type;
			elem->QueryIntAttribute("type", &type);
			initSetting.graphics.terrainTexture.type = type;
			for(TiXmlElement* elem2 = elem->FirstChildElement(); elem2 != NULL; elem2 = elem2->NextSiblingElement()){
				std::string elem2Name = elem2->Value() ;
				if(elem2Name == "colormap"){
					float colormapMin = -1.0f;
					float colormapMax = +1.0f;
					bool autoColormap = true;

					elem2->QueryBoolAttribute("auto", &autoColormap);
					initSetting.graphics.terrainTexture.autoColormap = autoColormap;

					elem2->QueryFloatAttribute("min", &colormapMin);
					initSetting.graphics.terrainTexture.colormapMin = colormapMin;	
					elem2->QueryFloatAttribute("max", &colormapMax);
					initSetting.graphics.terrainTexture.colormapMax = colormapMax;
				}
			}
		}
		else if(elemName == "skybox")
		{
			int type;
			elem->QueryIntAttribute("type", &type);
			initSetting.graphics.skyboxType = type;				
		}
		else if(elemName == "lighting")
		{
			float ambientLight = 1.0f;
			elem->QueryFloatAttribute("ambient", &ambientLight);
			initSetting.graphics.lighting.ambientLight = ambientLight;	

			float sun_altitude = 45.0f;
			elem->QueryFloatAttribute("sunAltitude", &sun_altitude);
			initSetting.graphics.lighting.sun_altitude = sun_altitude;

			float sun_azimuth = 30.0f;
			elem->QueryFloatAttribute("sunAzimuth", &sun_azimuth);
			initSetting.graphics.lighting.sun_azimuth = sun_azimuth;	
		}
		else if(elemName == "vertical")
		{
			float scale;
			elem->QueryFloatAttribute("scale", &scale);
			initSetting.graphics.verticalScale = scale;				
		}
		else if(elemName == "fresnel")
		{
			float fresnelCoef, refractive_index, att_1, att_2;
			elem->QueryFloatAttribute("coef", &fresnelCoef);
			elem->QueryFloatAttribute("refractive_index", &refractive_index);
			elem->QueryFloatAttribute("attenuation_1", &att_1);
			elem->QueryFloatAttribute("attenuation_2", &att_2);
			initSetting.graphics.fresnelCoef = fresnelCoef;
			initSetting.graphics.refractive_index = refractive_index;
			initSetting.graphics.att_1 = att_1;
			initSetting.graphics.att_2 = att_2;
		}
		else if(elemName == "camera")
		{
			bool autoCam;
			elem->QueryBoolAttribute("auto", &autoCam);
			initSetting.graphics.autoCam = autoCam;
			if(!autoCam)
			{ 
				float fov, x, y, z, pitch, yaw;
				elem->QueryFloatAttribute("FOV", &fov);
				elem->QueryFloatAttribute("x", &x);
				elem->QueryFloatAttribute("y", &y);
				elem->QueryFloatAttribute("z", &z);
				elem->QueryFloatAttribute("pitch", &pitch);
				elem->QueryFloatAttribute("yaw", &yaw);

				initSetting.graphics.camera.fov = fov;
				initSetting.graphics.camera.x = x;
				initSetting.graphics.camera.y = y;
				initSetting.graphics.camera.z = z;
				initSetting.graphics.camera.pitch = pitch;
				initSetting.graphics.camera.yaw = yaw;
			}
		}
	}

}
//********
bool readGraphicsInput()
{
	TiXmlDocument doc((initSetting.exePath + "/graphics/graphics.init").c_str());
	bool loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		TiXmlElement* root = doc.FirstChildElement();
		if(root == NULL)
		{	
			doc.Clear();
			MessageBox(0, "Invalid graphics.init file!", "Error!", MB_OK);
			return 1;
		}
		readGraphicsInput_helper(root);

		// Set autocolormap to true, so colormap values can't be set from graphics.init.
		initSetting.graphics.surfaceShading.autoColormap = true;
		initSetting.graphics.terrainTexture.autoColormap = true;
		return 0;
	}else{
		MessageBox(0, doc.ErrorDesc(), "Error!", MB_OK);
	}
}


//********

int initializeSetting()
{

	initSetting.exePath = ExePath();
	std::fstream fileIn ((initSetting.exePath + "/setting.init").c_str());
	std::string line;
	std::getline(fileIn, line);
	initSetting.initCMLName = line; 
	
	if(line == ""){
		return browseInputCML();
	}


	return readInputCML();
}

int main()
{

    try {
        return real_main();
    } catch (std::exception &e) {
        MessageBox(0, e.what(), "Error", MB_ICONEXCLAMATION | MB_OK);
    } catch (...) {
        MessageBox(0, "Unknown exception", "Error", MB_ICONEXCLAMATION | MB_OK);
    }

    return 1;
}

int CALLBACK WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
	int result;
	result = initializeSetting();
	if (!result){
		return main();
	}
	else
	{
		return 1;
	}
}
