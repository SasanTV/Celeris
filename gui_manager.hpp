/*
 * FILE:
 *   gui_manager.hpp
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

#ifndef GUI_MANAGER_HPP
#define GUI_MANAGER_HPP

#include "coercri/gcn/cg_font.hpp"
#include "coercri/gcn/cg_listener.hpp"
#include "coercri/gfx/font.hpp"
#include "coercri/gfx/window.hpp"
#include "coercri/timer/timer.hpp"

#include "boost/shared_ptr.hpp"

class GuiManager : public gcn::ActionListener {
public:

    GuiManager(boost::shared_ptr<Coercri::Window> window,
               boost::shared_ptr<Coercri::Timer> timer_,
               int gui_width_);
    
    ~GuiManager();

    // Caller should call logic() from main loop periodically.
    void logic();

    // Caller should call draw() from main loop, whenever the window is invalid.
    void draw(Coercri::GfxContext &gc);

    // this should be called (by the main loop) whenever the window changes size
    // (technically the GuiManager should create its own listener for this; but the main loop
    // already has a resize listener so it's easier to expose it and do it there)
    void resize();

    bool isGuiShown() const { return gui_shown; }
	bool isPause() const { return pause; }
	void switchPause() { pause = !pause; }

    // communicate camera motion back to main loop (HACK)
    bool getCameraReset(float &x, float &y, float &z, float &pitch, float &yaw);
    
protected:
    void action(const gcn::ActionEvent &e);

private:

    // define our own variant of CGListener that ignores keyboard events.
    // this is because arrow keys changing the GUI sliders can be confusing.
    class MyCGListener : public Coercri::CGListener {
    public:
        MyCGListener(boost::shared_ptr<Coercri::Window> window,
                   boost::shared_ptr<gcn::Gui> gui,
                   boost::shared_ptr<Coercri::Timer> timer) : Coercri::CGListener(window, gui, timer) { }
        //virtual void onCookedKey(Coercri::CookedKey ck, int ch, Coercri::KeyModifier mods) { }
    };
    
    boost::shared_ptr<Coercri::Font> font;
    boost::shared_ptr<Coercri::CGFont> cg_font;
    boost::shared_ptr<Coercri::Timer> timer;

    boost::shared_ptr<gcn::Gui> gui;
    boost::shared_ptr<Coercri::Window> window;
    MyCGListener listener;

    int gui_width;

    boost::shared_ptr<gcn::Container> show_gui_container;
    boost::shared_ptr<gcn::Button> show_gui_button, hide_gui_button, reset_simulation_button, reset_bathymetry_button, save_bathymetry_button, pause_simulation_button;
    boost::shared_ptr<gcn::Button> preset_button_valley, preset_button_valley_hires, preset_button_sea, preset_button_flat;
	boost::shared_ptr<gcn::Button> boundary_set_button, init_condition_set_button, mesh_size_set_button;
    bool gui_shown;
	bool pause;
    
    boost::shared_ptr<gcn::Container> container;
    std::vector<boost::shared_ptr<gcn::Label> > labels;
    std::vector<boost::shared_ptr<gcn::Slider> > sliders;
	std::vector<boost::shared_ptr<gcn::TextField> > textfields;
    std::vector<boost::shared_ptr<gcn::CheckBox> > check_boxes;
    std::vector<boost::shared_ptr<gcn::DropDown> > dropdowns;
    std::vector<boost::shared_ptr<gcn::Label> > values;

    std::vector<boost::shared_ptr<gcn::ListModel> > list_models;

    boost::shared_ptr<gcn::TabbedArea> tabbed_area;
    std::vector<boost::shared_ptr<gcn::Container> > tabs;
    int max_y;

    unsigned int last_time;
    int frame_count;

    bool cam_reset;
    float cam_x, cam_y, cam_z, cam_pitch, cam_yaw;
    
    // private methods
    void createGui();
	void createBoundaryTab(gcn::Container *tab);
	void createInitConditionTab(gcn::Container *tab);
	void createSettingConditionTab(gcn::Container *tab);
	
    void resetSliders();
	void resetColormapSliders();
    
	// private button do
	void boundary_set_button_do();
	void init_condition_set_button_do();
	void mesh_size_set_button_do();
    // prevent copying
    GuiManager(const GuiManager &);
    void operator=(const GuiManager &);
};

#endif
