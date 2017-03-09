/*
 * FILE:
 *   gui_manager.cpp
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
#include "engine.hpp"

#include "gui_manager.hpp"

#include "coercri/gfx/bitmap_font.hpp"
#include "coercri/gfx/load_bmp.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

namespace {
    // Constantly creating and destroying ostringstream objects seems v. bad
    // for performance, so we create a global one and re-use it
    std::ostringstream g_stringstream;
    void SetLabelTxt(gcn::Label &label, double value)
    {
        g_stringstream.str("");
        g_stringstream << value;
        label.setCaption(g_stringstream.str());
        label.adjustSize();
    }
    
    double RealToSlider(const Setting &setting, double value)
    {
        /*
        if (setting.logarithmic) {
            return std::log(value);
        } else {
            return value;
        }
        */
        return value;
    }

    double SliderToReal(const Setting &setting, double value)
    {
        /*
        if (setting.logarithmic) {
            return std::exp(value);
        } else {
            return value;
        }
        */
        return value;
    }

    class MyListModel : public gcn::ListModel {
    public:
        std::string getElementAt(int i) {
            if (i >= 0 && i < int(elts.size())) return elts[i];
            else return "";
        }

        int getNumberOfElements() {
            return int(elts.size());
        }

        void add(const std::string &s) { elts.push_back(s); }

    private:
        std::vector<std::string> elts;
    };



    void SetupFlatPlane()
    {

        SetSettingD("gravity", 9.81);
		SetSettingD("friction", initSetting.friction);
		SetSettingD("Theta", initSetting.theta);
		SetSettingD("nominal_cfl", 0.2);
        SetSettingD("timesteps_per_frame", 10);
        //SetSettingD("time_acceleration", 1);
			
		SetSettingD("Surface Shading",  initSetting.graphics.surfaceShading.type);
		SetSettingD("Shading Variable", initSetting.graphics.surfaceShading.shadingVariable);
		
		SetSettingD("Terrain Texture", initSetting.graphics.terrainTexture.type);
		SetSettingD("Skybox", initSetting.graphics.skyboxType);
		
		SetSettingD("fov", initSetting.graphics.camera.fov);
		SetSettingD("sun_alt", initSetting.graphics.lighting.sun_altitude);
		SetSettingD("sun_az", initSetting.graphics.lighting.sun_azimuth);
		SetSettingD("ambient", initSetting.graphics.lighting.ambientLight); //Default was 0.75
		SetSettingD("Vertical_Scale", initSetting.graphics.verticalScale);
		SetSettingD("fresnel_coeff", initSetting.graphics.fresnelCoef); // Default was 0.983
        //SetSettingD("fresnel_exponent", 5);
        //SetSettingD("specular_intensity", 1);
        //SetSettingD("specular_exponent", 15);
		SetSettingD("refractive_index", initSetting.graphics.refractive_index); // Default was 1.33
		SetSettingD("attenuation_1", initSetting.graphics.att_1); // Default was 0.08
		SetSettingD("attenuation_2", initSetting.graphics.att_2); // Default was 0.08
        SetSettingD("deep_r", 0.05);
        SetSettingD("deep_g", 0.1);
        SetSettingD("deep_b", 0.2);

		SetSettingD("mesh_size_x", initSetting.nx);
		SetSettingD("mesh_size_y", initSetting.ny);
		SetSettingD("valley_length", initSetting.length);
		SetSettingD("valley_width", initSetting.width);
        //SetSettingD("clip_camera", 0);
		SetSettingD("sea_level", initSetting.stillWaterElevation);
		//SetSettingD("use_sea_level", 1);
		SetSettingD("Grid Scale", initSetting.graphics.gridScale);
		SetSettingD("Grid", initSetting.graphics.gridOn);
    }
}

GuiManager::GuiManager(boost::shared_ptr<Coercri::Window> window_,
                       boost::shared_ptr<Coercri::Timer> timer_,
                       int gui_width_)
    : gui(new gcn::Gui),
      window(window_),
      listener(window_, gui, timer_),
      gui_width(gui_width_),
      last_time(0), frame_count(-999),
      timer(timer_),
      gui_shown(true), pause(false),
      cam_reset(false), cam_x(0), cam_y(0), cam_z(0), cam_pitch(0), cam_yaw(0)
{
    // load a font
    std::ifstream str((initSetting.exePath + "/graphics/fonts/knights_sfont.bmp").c_str(), std::ios::in | std::ios::binary);
    boost::shared_ptr<Coercri::PixelArray> parr = Coercri::LoadBMP(str);
    font.reset(new Coercri::BitmapFont(parr));
    cg_font.reset(new Coercri::CGFont(font));
    gcn::Widget::setGlobalFont(cg_font.get());

	
    // add the CGListener
    window->addWindowListener(&listener);
    listener.enableGui();



    // Default settings
    // note: duplicated code (see action())
    SetupFlatPlane();

	// Create the widgets.
    createGui();

    g_reset_type = R_SQUARE;
    cam_reset = true;
	cam_x = initSetting.graphics.camera.x; 
    cam_y = initSetting.graphics.camera.y; 
    cam_z = initSetting.graphics.camera.z; 
	cam_pitch = initSetting.graphics.camera.pitch; 
	cam_yaw = initSetting.graphics.camera.yaw ;
    resetSliders();    
}

GuiManager::~GuiManager()
{
    window->rmWindowListener(&listener);
}

void GuiManager::createGui()
{
    const int value_label_width = 80;

    show_gui_container.reset(new gcn::Container);
    hide_gui_button.reset(new gcn::Button("Hide GUI"));
    show_gui_button.reset(new gcn::Button("Show GUI"));
    hide_gui_button->addActionListener(this);
    show_gui_button->addActionListener(this);
    show_gui_container->add(show_gui_button.get());

    reset_simulation_button.reset(new gcn::Button("Reset Simulation"));
    reset_simulation_button->addActionListener(this);

	reset_bathymetry_button.reset(new gcn::Button("Reset Bathymetry"));
    reset_bathymetry_button->addActionListener(this);

	save_bathymetry_button.reset(new gcn::Button("Save Bathymetry"));
    save_bathymetry_button->addActionListener(this);

	save_inundation_button.reset(new gcn::Button("Save Inundation"));
    save_inundation_button->addActionListener(this);

	pause_simulation_button.reset(new gcn::Button("Pause Simulation"));
    pause_simulation_button->addActionListener(this);
    
	boundary_set_button.reset(new gcn::Button("  Set  "));
	boundary_set_button->addActionListener(this);
	
	init_condition_set_button.reset(new gcn::Button("  Set  "));
	init_condition_set_button->addActionListener(this);

	mesh_size_set_button.reset(new gcn::Button("  Set  "));
	mesh_size_set_button->addActionListener(this);


	container.reset(new gcn::Container);

    tabbed_area.reset(new gcn::TabbedArea);
    container->add(tabbed_area.get());
    tabbed_area->setSize(container->getWidth(), container->getHeight());

    int max_label_width = 0;
    const int pad = 14;
    int y = 12;
	
    boost::shared_ptr<gcn::Container> presets_tab(new gcn::Container);
    presets_tab->setSize(container->getWidth(), container->getHeight());
    tabs.push_back(presets_tab); // TODO: presets_tab does not really exists. but removing this line, causes a crash. figure out why!?
/*	
    //tabbed_area->addTab("Presets", presets_tab.get());
    preset_button_valley.reset(new gcn::Button("Valley"));
	presets_tab->add(preset_button_valley.get(), 5, 10);
    preset_button_valley->addActionListener(this);
    y += pad + preset_button_valley->getHeight();
    
    preset_button_valley_hires.reset(new gcn::Button("Valley (high res)"));
    //presets_tab->add(preset_button_valley_hires.get(), 5, y);
    //preset_button_valley_hires->addActionListener(this);
    //y += pad + preset_button_valley->getHeight();

    preset_button_flat.reset(new gcn::Button("Flat plane"));
    //presets_tab->add(preset_button_flat.get(), 5, y);
    //preset_button_flat->addActionListener(this);
    //y += pad + preset_button_valley->getHeight();
				
				
    preset_button_sea.reset(new gcn::Button("Sea + waves"));
	presets_tab->add(preset_button_sea.get(), 5, 40);
    preset_button_sea->addActionListener(this);

    const int pmaxw = std::max(std::max(std::max(preset_button_valley->getWidth(),
                                                 preset_button_valley_hires->getWidth()),
                                        preset_button_sea->getWidth()),
                               preset_button_flat->getWidth());
    preset_button_valley->setWidth(pmaxw);
    preset_button_valley_hires->setWidth(pmaxw);
    preset_button_sea->setWidth(pmaxw);
    preset_button_flat->setWidth(pmaxw);
    preset_button_valley->setAlignment(gcn::Graphics::LEFT);
    preset_button_valley_hires->setAlignment(gcn::Graphics::LEFT);
    preset_button_sea->setAlignment(gcn::Graphics::LEFT);
    preset_button_flat->setAlignment(gcn::Graphics::LEFT);
*/
    
    for (const Setting *setting = &g_settings[0]; setting->name; ++setting) {


        if (setting->type == S_NEW_TAB && setting->name[0] != 0) {

			boost::shared_ptr<gcn::Container> current_tab(new gcn::Container);
            current_tab->setSize(container->getWidth(), container->getHeight());
            tabs.push_back(current_tab);
            tabbed_area->addTab(setting->name, current_tab.get());
			if(setting->name == "Boundary"){
				createBoundaryTab(current_tab.get());
			} else if(setting->name == "Init.Cond."){
				createInitConditionTab(current_tab.get());
			} else if(setting->name == "settings"){
				createSettingConditionTab(current_tab.get());
			}
			

        }

        boost::shared_ptr<gcn::Label> lab(new gcn::Label(setting->name));
        labels.push_back(lab);

        boost::shared_ptr<gcn::Slider> slider;
		boost::shared_ptr<gcn::TextField> textfield;
        boost::shared_ptr<gcn::CheckBox> checkbox;
        boost::shared_ptr<gcn::DropDown> dropdown;

        switch (setting->type) {
        case S_SLIDER:
        case S_SLIDER_MULT_4:
        case S_SLIDER_INT:
            slider.reset(new gcn::Slider(RealToSlider(*setting, setting->min),
                                         RealToSlider(*setting, setting->max)));
            slider->setValue(RealToSlider(*setting, setting->value));
            slider->addActionListener(this);
            break;
		case S_TEXTFIELD: 
			textfield.reset(new gcn::TextField());
			textfield->addActionListener(this);
			break;


        case S_CHECKBOX:
            checkbox.reset(new gcn::CheckBox("", setting->value == 1));
            checkbox->addActionListener(this);
            break;
		case S_WATER_SURFACE_SHADING:

            {
                boost::shared_ptr<MyListModel> list_model(new MyListModel);
                list_model->add("Photorealistic");
                list_model->add("Parula");
                list_model->add("Jet");
				list_model->add("Zebra");
                list_model->add("Custom");
                list_models.push_back(list_model);
            
                dropdown.reset(new gcn::DropDown(list_model.get()));
                dropdown->addActionListener(this);
            }
            break;

		case S_WATER_SURFACE_SHADING_VARIABLE:

            {
                boost::shared_ptr<MyListModel> list_model(new MyListModel);
                list_model->add("water surface (eta)");
                list_model->add("u velocity");
                list_model->add("v velocity");
				list_model->add("velocity magnitude");
				list_model->add("vorticity");
                list_models.push_back(list_model);
            
                dropdown.reset(new gcn::DropDown(list_model.get()));
                dropdown->addActionListener(this);
            }
            break;
		case S_TERRAIN_TEXTURE:

            {
				//TODO: Make a for loop to load all bmp files in the texture directory.
                boost::shared_ptr<MyListModel> list_model(new MyListModel);
                list_model->add("Sand");
                list_model->add("Grass");
                list_model->add("Concrete");
				list_model->add("Tiles 1");
				list_model->add("Tiles 2");
				list_model->add("Tiles 3");
				list_model->add("Colormap");
				list_model->add("Contours");
				list_model->add("Custom");
                list_models.push_back(list_model);
            
                dropdown.reset(new gcn::DropDown(list_model.get()));
                dropdown->addActionListener(this);
            }
			break;
		case S_SKYBOX_TEXTURE:

            {
				//TODO: Make a for loop to load all bmp files in the texture directory.
                boost::shared_ptr<MyListModel> list_model(new MyListModel);
                list_model->add("Thick Clouds");
                list_model->add("Clouds");
                list_model->add("Ocean");
				list_model->add("Sunset");
				list_model->add("Desert");
				list_model->add("Dark Solid");
				list_model->add("Light Solid");
                list_models.push_back(list_model);
            
                dropdown.reset(new gcn::DropDown(list_model.get()));
                dropdown->addActionListener(this);
            }
            break;
		case S_BOUNDARY_DROPDOWN:

            {
                boost::shared_ptr<MyListModel> list_model(new MyListModel);
                list_model->add("West");
                list_model->add("East");
                list_model->add("South");
				list_model->add("North");
                list_models.push_back(list_model);
            
                dropdown.reset(new gcn::DropDown(list_model.get()));
                dropdown->addActionListener(this);
            }
            break;
		case S_BOUNDARY_OPTIONS_DROPDOWN:

			{
				boost::shared_ptr<MyListModel> list_model(new MyListModel);
				list_model->add("Solid");
				list_model->add("Sponge");
				list_model->add("SineWave");
				list_model->add("IrregularWaves");
				list_model->add("UniformTimeSeries");
				list_models.push_back(list_model);
		    
				dropdown.reset(new gcn::DropDown(list_model.get()));
				dropdown->addActionListener(this);
			}
			break;

        case S_LEFT_MOUSE_DROPDOWN:

            {
                boost::shared_ptr<MyListModel> list_model(new MyListModel);
                list_model->add("Add Water");
                list_model->add("Remove Water");
                list_model->add("Stir Water");
                list_model->add("Raise Terrain");
                list_model->add("Lower Terrain");

                list_models.push_back(list_model);
            
                dropdown.reset(new gcn::DropDown(list_model.get()));
                dropdown->addActionListener(this);
            }
            break;
        }

        sliders.push_back(slider);
        check_boxes.push_back(checkbox);
        dropdowns.push_back(dropdown);
		textfields.push_back(textfield);

        boost::shared_ptr<gcn::Label> val_lbl;
		if (setting->type == S_LABEL || setting->type == S_SLIDER || setting->type == S_SLIDER_MULT_4 || setting->type == S_SLIDER_INT) {
            val_lbl.reset(new gcn::Label);
            SetLabelTxt(*val_lbl, setting->value);
        }
		if (setting->type == S_TEXTFIELD){
			std::ostringstream g_stringstream;
			g_stringstream.str("");
			float temp = GetSetting(setting->name);
			g_stringstream << (temp);
			textfield->setText(g_stringstream.str());
        }

        values.push_back(val_lbl);
        
        const int w = lab->getWidth();
        if (w > max_label_width) max_label_width = w;
    }

    const int slider_width = std::max(10, gui_width - max_label_width - value_label_width);

    y = 0;

    std::vector<boost::shared_ptr<gcn::Container> >::iterator tab_it = tabs.begin();
    max_y = 0;
    
    for (size_t i = 0; i < sliders.size(); ++i) {

        if (g_settings[i].type == S_NEW_TAB) {
            ++tab_it;
            max_y = std::max(max_y, y);
            
            y = 8;
            if (tab_it == tabs.end()) {
                max_y += tabbed_area->getFont()->getHeight() + 30;
                y = max_y + 20;
            }
            
            continue;
        }
        
        if (sliders[i]) {
            sliders[i]->setWidth(slider_width);
            sliders[i]->setHeight(labels[i]->getHeight());
			sliders[i]->setStepLength(0.5); // TODO: Ste step according to maxmin value e.g: (max-min)/20
			
        }

		
        if (textfields[i]) {
            textfields[i]->setWidth(slider_width/2);
			textfields[i]->setHeight(labels[i]->getHeight());
        }

        if (values[i]) {
            values[i]->setWidth(value_label_width);
            values[i]->setHeight(labels[i]->getHeight());
        }

        if (dropdowns[i]) {
            dropdowns[i]->setWidth(slider_width + value_label_width/2);
        }

        gcn::Container & con =
            (tab_it == tabs.end() ? *container : **tab_it);
        
        con.add(labels[i].get(), 0, y);
        if (sliders[i]) con.add(sliders[i].get(), max_label_width, y);
		if (textfields[i]) con.add(textfields[i].get(), max_label_width, y);
        if (check_boxes[i]) con.add(check_boxes[i].get(), max_label_width, y);
        if (dropdowns[i]) con.add(dropdowns[i].get(), max_label_width, y);

        int valx = max_label_width;
        if (g_settings[i].type != S_LABEL) valx += slider_width;
        if (values[i]) con.add(values[i].get(), valx, y);
        
        y += std::max(labels[i]->getHeight(),
                      dropdowns[i] ? dropdowns[i]->getHeight() : 0);
    }

    y += 16;
    container->add(reset_simulation_button.get(), 8, y);
    container->add(hide_gui_button.get(), 8 + reset_simulation_button->getWidth() + 50, y);

	y += 32;
    container->add(reset_bathymetry_button.get(), 8, y);

	y += 32;
    container->add(save_bathymetry_button.get(), 8, y);

	y += 32;
	container->add(pause_simulation_button.get(), 8, y);

	y += 32;
	container->add(save_inundation_button.get(), 8, y);


}

void GuiManager::createBoundaryTab(gcn::Container *tab)
{
	tab->add(boundary_set_button.get(), 300, 110);
}

void GuiManager::createInitConditionTab(gcn::Container *tab)
{
	tab->add(init_condition_set_button.get(), 300, 110);
}

void GuiManager::createSettingConditionTab(gcn::Container *tab)
{
	tab->add(mesh_size_set_button.get(), 300, 16);
}


void GuiManager::resize()
{
    // gui_width is fixed, so we don't have to resize any of the widgets;
    // we just move the container to a new position.

    // (we also shrink the container if the window is so small that
    // gui_width doesn't fit.)
    
    int win_width, win_height;
    window->getSize(win_width, win_height);
    const int width = std::min(win_width, gui_width);
    const int x = win_width - width;
   
    container->setWidth(width);
	
    container->setHeight(win_height);
    container->setX(x);
    show_gui_container->setSize(win_width, win_height);
    show_gui_container->setPosition(0, 0);
    show_gui_container->setOpaque(false);

    tabbed_area->setSize(container->getWidth(), max_y);
    for (std::vector<boost::shared_ptr<gcn::Container> >::iterator it = tabs.begin(); it != tabs.end(); ++it) {
        (*it)->setSize(container->getWidth(), max_y);
    }

    show_gui_button->setPosition(win_width - show_gui_button->getWidth() - 10,
                                 10);
    show_gui_button->adjustSize();

    if (gui_shown) {
        gui->setTop(container.get());
    } else {
        gui->setTop(show_gui_container.get());
    }
}

void GuiManager::logic()
{
    listener.processInput();

    // reset all value labels (every few frames)
    if (frame_count == 0) {
        for (int i = 0; i < int(sliders.size()); ++i) {
            if (values[i]) {
                SetLabelTxt(*values[i], g_settings[i].value);
            }
        }
    }
}

void GuiManager::draw(Coercri::GfxContext &gc)
{
    listener.draw(gc);

    if (frame_count < 0) {
        frame_count = 0;
        last_time = timer->getMsec();

    } else {
        ++frame_count;
        unsigned int time_now = timer->getMsec();
        unsigned int difference = time_now - last_time;
        if (difference > 250) {

            const float fps = (1000.0f * frame_count) / float(difference);
            frame_count = 0;
            last_time = time_now;

            SetSetting("fps", fps);
        }
    }
}

void GuiManager::action(const gcn::ActionEvent &e)
{
    gcn::Slider *slid = static_cast<gcn::Slider*>(e.getSource());
	gcn::TextField *tb = static_cast<gcn::TextField*>(e.getSource());
    gcn::CheckBox *cb = static_cast<gcn::CheckBox*>(e.getSource());
    gcn::DropDown *dd = static_cast<gcn::DropDown*>(e.getSource());

    for (size_t i = 0; i < sliders.size(); ++i) {
        if (sliders[i].get() == slid) {
            g_settings[i].value = SliderToReal(g_settings[i], slid->getValue());
            g_reset_type = std::max(g_reset_type, g_settings[i].reset_type);

            if (g_settings[i].type == S_SLIDER_MULT_4) {
                int val = int(g_settings[i].value + 0.5);
                val = (val + 3) & (~3);
                g_settings[i].value = double(val);
            } else if (g_settings[i].type == S_SLIDER_INT) {
                int val = int(g_settings[i].value + 0.5);
                g_settings[i].value = double(val);
            }
            resetColormapSliders();
            SetLabelTxt(*values[i], g_settings[i].value);
            break;
        }
		 if (textfields[i].get() == tb) {
			g_settings[i].value = std::atof((tb->getText()).c_str());
			g_reset_type = std::max(g_reset_type, g_settings[i].reset_type);
            break;
        }

        if (check_boxes[i].get() == cb) {
            g_settings[i].value = cb->isSelected() ? 1.0 : 0.0;
            g_reset_type = std::max(g_reset_type, g_settings[i].reset_type);
            break;
        }

        if (dropdowns[i].get() == dd) {
            g_settings[i].value = dd->getSelected();
            g_reset_type = std::max(g_reset_type, g_settings[i].reset_type);
			if (g_settings[i].type == S_BOUNDARY_OPTIONS_DROPDOWN) {

				for (size_t j = 0; j < textfields.size(); ++j) { // Seems with the current GUI structure, iteration is the only way to get to a specific object
					gcn::TextField *tempField = textfields[j].get();
					if(tempField != 0){
						if(g_settings[j].name == "wave amplitude" || g_settings[j].name == "wave direction" || g_settings[j].name == "wave period"){
							tempField->setEnabled((GetIntSetting("Boundary Type")==2));
							tempField->setBackgroundColor((GetIntSetting("Boundary Type")==2)? gcn::Color(255,255,255) : gcn::Color(180,180,180));// ((GetIntSetting("Boundary Type")==2));
						}
					}
				}


			}

            break;
        }
    }

    if (e.getSource() == hide_gui_button.get()) {
        gui_shown = false;
        resize();
    } else if (e.getSource() == show_gui_button.get()) {
        gui_shown = true;
        resize();
    } else if (e.getSource() == reset_simulation_button.get()) {
        g_reset_type = R_MESH;
	} else if (e.getSource() == reset_bathymetry_button.get()) {
        g_reset_type = R_MESH;
		initSetting.rereadBathy = true;
	} else if (e.getSource() == save_bathymetry_button.get()) {
		initSetting.saveBathymetry = true;
	} else if (e.getSource() == save_inundation_button.get()) {
		initSetting.saveInundation = true;
	} else if (e.getSource() == pause_simulation_button.get()) {
        pause = !pause ;
		pause_simulation_button.get()->setCaption(pause?"Run Simulation":"Pause Simulation");
	} else if (e.getSource() == init_condition_set_button.get()) {
        init_condition_set_button_do();	
	} else if (e.getSource() == mesh_size_set_button.get()) {
		initSetting.rereadBathy = true;
        mesh_size_set_button_do();	
	} else if (e.getSource() == boundary_set_button.get()) {
		boundary_set_button_do();	
    }
}

void GuiManager::init_condition_set_button_do(){

	float H, theta, xcenter, ycenter;
	const float W = GetSetting("valley_width");
	const float L = GetSetting("valley_length");
	float length_in = W + L; // W + L is just a large enough number.
	
	for (size_t j = 0; j < textfields.size(); ++j) { // Seems with the current GUI structure, iteration is the only way to get to a specific object
		gcn::TextField *tempField = textfields[j].get();

		if(tempField != 0){

			if     (g_settings[j].name == "wave height"){
				H = std::atof(tempField->getText().c_str());
			} else if(g_settings[j].name == "wave direction "){ // Note the white space at the end.
				theta = std::atof(tempField->getText().c_str());
			} else if(g_settings[j].name == "x_0"){
				xcenter = std::atof(tempField->getText().c_str());
			} else if(g_settings[j].name == "y_0"){
				ycenter = std::atof(tempField->getText().c_str());
			}
		}
	}

	initSetting.solitons[initSetting.MAX_NUM_SOLITARY - 1].setSoliton(H, theta, xcenter, ycenter, length_in);
	initSetting.is_there_new_solitary_wave = true;

}



void GuiManager::boundary_set_button_do()
{
		int boundary_side = GetIntSetting("Boundary Side");
		int boundary_type = GetIntSetting("Boundary Type");
		std::string type = "";

		switch(boundary_type){
		case B_SOLID: 
			type = "Solid";
			break;
		case B_SPONGE:
			type = "Sponge";
			g_reset_type = std::max(g_reset_type, R_MESH); // Requires a reset after settin a sponge layer.
			break;
		case B_SINEWAVE:
			type = "SineWave";
			break;
		case B_IRREGULARWAVE:
			type = "IrregularWaves";
			break;
		case B_UNIFORMTIMESERIES:
			type = "UniformTimeSeries";
			break;
		}
		
		int width_num;
		float amplitude, period, theta, seal_level;
		const float PI = std::atan(1.0f) * 4.0f;

		for (size_t j = 0; j < textfields.size(); ++j) { // Seems with the current GUI structure, iteration is the only way to get to a specific object
			gcn::TextField *tempField = textfields[j].get();

			if(tempField != 0){

				if     (g_settings[j].name == "width num"){
					width_num = std::atoi(tempField->getText().c_str());
				} else if(g_settings[j].name == "b. sea level"){
					seal_level = std::atof(tempField->getText().c_str());
				} else if(g_settings[j].name == "wave amplitude"){
					amplitude = std::atof(tempField->getText().c_str());
				} else if(g_settings[j].name == "wave period"){
					period = std::atof(tempField->getText().c_str());
				} else if(g_settings[j].name == "wave direction"){
					theta = std::atof(tempField->getText().c_str());
				}
			}
		}

		switch (boundary_side){
		
		case B_WEST:
			initSetting.westBoundary.type = type;
			initSetting.westBoundary.width = width_num;
			initSetting.westBoundary.waterLevel = seal_level;
			initSetting.westBoundary.sineWaveSetting.amplitude = amplitude;
			initSetting.westBoundary.sineWaveSetting.period = period;
			initSetting.westBoundary.sineWaveSetting.theta = PI/180.0 * theta;
			initSetting.westBoundary.hasChanged = true;
			break;
		case B_EAST:
			initSetting.eastBoundary.type = type;
			initSetting.eastBoundary.width = width_num;
			initSetting.eastBoundary.waterLevel = seal_level;
			initSetting.eastBoundary.sineWaveSetting.amplitude = amplitude;
			initSetting.eastBoundary.sineWaveSetting.period = period;
			initSetting.eastBoundary.sineWaveSetting.theta = PI/180.0 * theta;
			initSetting.eastBoundary.hasChanged = true;
			break;
		case B_SOUTH:
			initSetting.southBoundary.type = type;
			initSetting.southBoundary.width = width_num;
			initSetting.southBoundary.waterLevel = seal_level;
			initSetting.southBoundary.sineWaveSetting.amplitude = amplitude;
			initSetting.southBoundary.sineWaveSetting.period = period;
			initSetting.southBoundary.sineWaveSetting.theta = PI/180.0 * theta;
			initSetting.southBoundary.hasChanged = true;

			break;
		case B_NORTH:
			initSetting.northBoundary.type = type;
			initSetting.northBoundary.width = width_num;
			initSetting.northBoundary.waterLevel = seal_level;
			initSetting.northBoundary.sineWaveSetting.amplitude = amplitude;
			initSetting.northBoundary.sineWaveSetting.period = period;
			initSetting.northBoundary.sineWaveSetting.theta = PI/180.0 * theta;
			initSetting.northBoundary.hasChanged = true;

			break;
		}
}

void GuiManager::mesh_size_set_button_do(){

	const float W = GetSetting("valley_width");
	const float L = GetSetting("valley_length");
	float length_in = W + L; // W + L is just a large enough number.

	
	for (size_t j = 0; j < textfields.size(); ++j) { // Seems with the current GUI structure, iteration is the only way to get to a specific object
		gcn::TextField *tempField = textfields[j].get();

		if(tempField != 0){
			if		 (g_settings[j].name == "mesh_size_x"){
				SetSettingD("mesh_size_x", std::atof(tempField->getText().c_str()));
				g_reset_type = std::max(g_reset_type, g_settings[j].reset_type);
			} else if(g_settings[j].name == "mesh_size_y"){ 
				SetSettingD("mesh_size_y", std::atof(tempField->getText().c_str()));
				g_reset_type = std::max(g_reset_type, g_settings[j].reset_type);
			}
		}
	}
	 
	
}
bool GuiManager::getCameraReset(float &x, float &y, float &z, float &pitch, float &yaw)
{
    if (cam_reset) {
        x = cam_x;
        y = cam_y;
        z = cam_z;
        pitch = cam_pitch;
        yaw = cam_yaw;
        cam_reset = false;
        return true;
    } else {
        return false;
    }
}

void GuiManager::resetSliders()
{
    for (int i = 0; i < int(sliders.size()); ++i) {
        Setting &s = g_settings[i];
        if (s.type == S_SLIDER || s.type == S_SLIDER_MULT_4 || s.type == S_SLIDER_INT) {
            sliders[i]->setValue(RealToSlider(s, s.value));
            SetLabelTxt(*values[i], s.value);
        } else if (s.type == S_CHECKBOX) {
            check_boxes[i]->setSelected(s.value != 0);
        }
    }
}

void GuiManager::resetColormapSliders()
{
    for (int i = 0; i < int(sliders.size()); ++i) {
        Setting &s = g_settings[i];
		if (s.type == S_SLIDER && (s.name == "Colormap Max" || s.name == "Colormap Min" || s.name == "Colormap Min " || s.name == "Colormap Max " || s.name == "Flow Depth")) {
			sliders[i]->setScale(RealToSlider(s, s.min),RealToSlider(s, s.max));	
			sliders[i]->setStepLength((s.max-s.min)/20.0f);
			sliders[i]->setValue(RealToSlider(s, s.value));
            SetLabelTxt(*values[i], s.value);
		} 

    }
}
