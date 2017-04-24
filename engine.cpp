/*
 * FILE:
 *   engine.cpp
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
#include "settings.hpp"
#include "terrain_heightfield.hpp"

#include "coercri/dx11/core/dx_error.hpp"
#include "coercri/gfx/load_bmp.hpp"
#include "coercri/gfx/pixel_array.hpp"

#include "boost/scoped_array.hpp"
#include "coercri/timer/generic_timer.hpp"

#include <d3d11.h>
#include <d3dx11.h>
#include <d3dcompiler.h>
#include <d3dx10math.h>
#include <xnamath.h>

#include <cmath>
#include <boost/math/special_functions/gamma.hpp>

#include <fstream>
#include <string>


#include <iostream>
#include <sstream>
#include <math.h>

#include <time.h>


#ifdef near
#undef near
#endif
#ifdef far
#undef far
#endif

// Turning off USE_KP07 activates an experimental Lax-Wendroff solver.
// Unfortunately this is buggy and unstable currently, so leaving USE_KP07 
// on is recommended.
#define USE_KP07

// Debugging switches
//#define DUMP_TO_FILE
//#define RUN_CHECKS

#if defined(DUMP_TO_FILE) || defined(RUN_CHECKS)
#include <iomanip>
#endif

namespace {

    const float PI = 4.0f * std::atan(1.0f);
	const float BCOEF=1/15.0f;

    float GaussianBrush(float x, float y, float r)
    {
        const float sigma = 1.0f / 3.0f;
        float z = (x*x + y*y) / (r*r);
        if (z > 1) return 0;
        else return 2 * std::exp(-z/(2*sigma*sigma));
    }


	
	void solitaryWave (const float waterDepth, const float param_H, const float theta, const float x_zero, const float y_zero, float x, float y, float *eta, float *hu, float *hv){
		
		const float param_K = sqrt(0.75f * abs(param_H)/pow(waterDepth,3));
		const float param_C = sqrt (GetSetting("gravity") * (param_H+waterDepth));
		
		const float temp_eta = param_H * 1.0f/pow(cosh(param_K * ((x - x_zero) * cos(theta) + (y - y_zero) * sin(theta))),2);
		*eta =  temp_eta; // stillwater elevation is asumed to be zero. this is eta!
		*hu = ( param_C * cos(theta) * temp_eta);
		*hv = ( param_C * sin(theta) * temp_eta);
	}
    
    class MapTexture {
    public:
        MapTexture(ID3D11DeviceContext &cxt, ID3D11Texture2D & tex)
            : context(cxt), texture(tex)
        {
            HRESULT hr = context.Map(&texture, 0, D3D11_MAP_READ, 0, &msr);
            if (FAILED(hr)) {
                throw Coercri::DXError("Map failed", hr);
            }
        }

        ~MapTexture()
        {
            context.Unmap(&texture, 0);
        }

        // allow access to the pointer / pitch values
        D3D11_MAPPED_SUBRESOURCE msr;

    private:
        ID3D11DeviceContext &context;
        ID3D11Texture2D &texture;  // must be a staging texture
    };                   
    
    // VS input for water/terrain rendering
    struct MeshVertex {
        int x, y;
    };

    // Const buffer for water/terrain rendering
    struct MyConstBuffer {


        XMMATRIX tex_to_clip;   // transforms (tex_x, tex_y, world_z, 1) into clip space 
        XMFLOAT3 light_dir;         // in world space 
        float ambient;
        XMFLOAT3 eye_mult;
        float pack6;
        XMFLOAT3 eye_trans;
        float pack7;
        XMFLOAT2 terrain_tex_scale;
        float pack8, pack9;
        XMMATRIX skybox_mtx;
        float world_mult_x, world_mult_y, world_trans_x, world_trans_y;
        float world_to_grass_tex_x, world_to_grass_tex_y;
        float grass_tex_of_origin_x, grass_tex_of_origin_y;
		
        float fresnel_coeff, fresnel_exponent;
        float specular_intensity, specular_exponent;
        float refractive_index;
        float attenuation_1, attenuation_2;
        
        int nx_plus_1, ny_plus_1;
		//float dx, dy;


        float deep_r, deep_g, deep_b;

		float zScale;
		float seaLevel;

		int water_shading;
		int terrain_shading;
		
		float water_colormap_min;
		float water_colormap_max;

		float terrain_colormap_min;
		float terrain_colormap_max;

		int isGridOn; //passing bool as integer.

		int is_dissipation_threshold_on; //passing bool as integer.
		float dissipation_threshold;

		float sqrt_sqrt_epsilon;
		float drylandDepthOfInundation;
    };    

    // Const buffer used for the simulation
    struct SimConstBuffer {
        float two_theta;
        float two_over_nx_plus_four;
        float two_over_ny_plus_four;
        float g;
        float half_g;
		float Bcoef_g;
        float g_over_dx;
        float g_over_dy;
        float one_over_dx;
        float one_over_dy;
        float dt;
        float epsilon;      // usually dx^4
        int nx;
		int ny;
        float friction;  // m s^-2
		int isManning;
		float seaLevel;

		float dissipation_threshold; // For visualization purposes, not simulation.
		float whiteWaterDecayRate; // For visualization purposes, not simulation.
    };
	

    // Const buffer used for boundary conditions
    struct BoundaryConstBuffer {
		float PI;
        float boundary_epsilon;
		int boundary_nx, boundary_ny;  // number of cells, including ghost zones
		float dx, dy;
        int reflect_x, reflect_y;
		int northBoundaryWidth, eastBoundaryWidth, westBoundaryWidth, southBoundaryWidth;
		float northSeaLevel, eastSeaLevel, westSeaLevel, southSeaLevel;
		float eastAmplitude_or_eta, eastPeriod_or_hu, eastTheta_or_hv;
		float westAmplitude_or_eta, westPeriod_or_hu, westTheta_or_hv;
		float northAmplitude_or_eta, northPeriod_or_hu, northTheta_or_hv;
		float southAmplitude_or_eta, southPeriod_or_hu, southTheta_or_hv;
        int solid_wall_flag;
        int inflow_x_min, inflow_x_max;
        float sea_level, inflow_height, inflow_speed;
        float g;
        float total_time;
        float sa1, skx1, sky1, so1;
        float sa2, skx2, sky2, so2;
        float sa3, skx3, sky3, so3;
        float sa4, skx4, sky4, so4;
        float sdecay;
    };

	// Const buffer for setting time integration scheme.
	struct TimeIntegrationConstBuffer{
		int tScheme;
	};

	struct WaveParam{
		float amplitude, period, theta, phase;
	};

	struct IrregularWavesDataConstBuffer{
		WaveParam   wavesWest[MAX_NUM_OF_IRREGULAR_WAVES],
					wavesEast[MAX_NUM_OF_IRREGULAR_WAVES],
					wavesSouth[MAX_NUM_OF_IRREGULAR_WAVES],
					wavesNorth[MAX_NUM_OF_IRREGULAR_WAVES];
		int numberOfWavesWest, numberOfWavesEast, numberOfWavesSouth, numberOfWavesNorth; // This is the number of sinewaves to be superposed.
	}; 

	struct IrregularWavesColumnRowConstBuffer{
		int columnNumber;
		int rowNumber;
	};

	// Const buffer for adding solitary wave.
	struct SolitaryWaveConstBuffer{
		
		float g;
	
		float dx;
		float dy;

		float param_H; 
		float x_zero;
		float y_zero;
		float theta;
	};

    struct LeftMouseConstBuffer {
        float scale_x, scale_y;
        float bias_x, bias_y;
        float two_over_nx_plus_four;
        float two_over_ny_plus_four;
        float disp_A, disp_B;
    };

    float CalcEpsilon()
    {
        const float W = GetSetting("valley_width");
        const float L = GetSetting("valley_length");
        const float nx = GetSetting("mesh_size_x");
        const float ny = GetSetting("mesh_size_y");
        const float dx = W / (nx-1);
        const float dy = L / (ny-1);
        return std::min(initSetting.epsilon, std::pow(std::max(dx, dy), 4));
    }

    float CalcU(float h, float hu)
    {
        float epsilon = CalcEpsilon();
        float h2 = h*h;
        float h4 = h2*h2;
        float divide_by_h = sqrt(2.0f) * h / sqrt(h4 + std::max(h4, epsilon));
        return divide_by_h * hu;
    }
    
    void CompileShader(const char *filename,
                       const char *entry_point,
                       const char *profile,
                       Coercri::ComPtrWrapper<ID3DBlob> &compiled_code)
    {
        DWORD compile_flags = D3DCOMPILE_ENABLE_STRICTNESS | D3DCOMPILE_WARNINGS_ARE_ERRORS;
#ifdef _DEBUG
        compile_flags |= D3DCOMPILE_DEBUG;
#endif

        ID3DBlob *output_blob = 0;
        ID3DBlob *error_blob = 0;
        
        HRESULT hr = D3DX11CompileFromFile(filename,
                                           0,   // defines
                                           0,   // include handler
                                           entry_point,
                                           profile,
                                           compile_flags,
                                           0,  // effect flags (ignored)
                                           0,  // thread pump
                                           &output_blob,
                                           &error_blob,
                                           0);    // used for async compilation only

        compiled_code.reset(output_blob);
        Coercri::ComPtrWrapper<ID3DBlob> error_sentinel(error_blob);
        
        if (FAILED(hr)) {
            std::string msg = "D3DX11CompileFromFile failed";
            if (error_blob) {
                OutputDebugStringA((char*)error_blob->GetBufferPointer());
                msg += "\n\n";
                msg += (char*)error_blob->GetBufferPointer();
            }
            throw Coercri::DXError(msg, hr);
        }
    }

    void CreateVertexShader(ID3D11Device *device,
                            const char *filename,
                            const char *entry_point,
                            Coercri::ComPtrWrapper<ID3DBlob> &bytecode,
                            Coercri::ComPtrWrapper<ID3D11VertexShader> &output)
    {
        const char * vs_profile = "vs_4_0";

        CompileShader(filename, entry_point, vs_profile, bytecode);

        ID3D11VertexShader *vert_shader = 0;
        HRESULT hr = device->CreateVertexShader(bytecode->GetBufferPointer(),
                                                bytecode->GetBufferSize(),
                                                0,
                                                &vert_shader);
        if (FAILED(hr)) {
            throw Coercri::DXError("CreateVertexShader failed", hr);
        }
        output.reset(vert_shader);
    }

    void CreatePixelShader(ID3D11Device *device,
                           const char *filename,
                           const char *entry_point,
                           Coercri::ComPtrWrapper<ID3D11PixelShader> &output)
    {
        const char * ps_profile = "ps_4_0";

        Coercri::ComPtrWrapper<ID3DBlob> bytecode;
        CompileShader(filename, entry_point, ps_profile, bytecode);

        ID3D11PixelShader *pixel_shader = 0;
        HRESULT hr = device->CreatePixelShader(bytecode->GetBufferPointer(),
                                               bytecode->GetBufferSize(),
                                               0,
                                               &pixel_shader);
        if (FAILED(hr)) {
            throw Coercri::DXError("CreatePixelShader failed", hr);
        }
        output.reset(pixel_shader);
    }

    void GetTextureSize(ID3D11Texture2D *tex, int &width, int &height)
    {
        D3D11_TEXTURE2D_DESC td;
        tex->GetDesc(&td);

        width = td.Width;
        height = td.Height;
    }

    void GetRenderTargetSize(ID3D11RenderTargetView *view, int &width, int &height)
    {
        ID3D11Resource *resource;
        view->GetResource(&resource);
        Coercri::ComPtrWrapper<ID3D11Resource> psResource(resource);

        ID3D11Texture2D *texture;
        resource->QueryInterface(__uuidof(ID3D11Texture2D), reinterpret_cast<void**>(&texture));
        Coercri::ComPtrWrapper<ID3D11Texture2D> psTexture(texture);
		
        GetTextureSize(texture, width, height);
    }

    void CreateTextureImpl(ID3D11Device *device,
                           const D3D11_TEXTURE2D_DESC &td,
                           const D3D11_SUBRESOURCE_DATA *srd,
                           Coercri::ComPtrWrapper<ID3D11Texture2D> &out_tex,
                           Coercri::ComPtrWrapper<ID3D11ShaderResourceView> *out_srv,
                           Coercri::ComPtrWrapper<ID3D11RenderTargetView> *out_rtv)
    {
        ID3D11Texture2D *pTexture;
        HRESULT hr = device->CreateTexture2D(&td, srd, &pTexture);
        if (FAILED(hr)) {
            throw Coercri::DXError("Failed to create texture", hr);
        }
        out_tex.reset(pTexture);

        if (out_srv) {
            D3D11_SHADER_RESOURCE_VIEW_DESC sd;
            memset(&sd, 0, sizeof(sd));
            sd.Format = td.Format;
            sd.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
            sd.Texture2D.MipLevels = td.MipLevels;
            
            ID3D11ShaderResourceView *pSRV;
            hr = device->CreateShaderResourceView(pTexture, &sd, &pSRV);
            if (FAILED(hr)) {
                throw Coercri::DXError("Failed to create shader resource view for texture", hr);
            }
            out_srv->reset(pSRV);
        }

        if (out_rtv) {
            D3D11_RENDER_TARGET_VIEW_DESC rd;
            memset(&rd, 0, sizeof(rd));
            rd.Format = td.Format;
            rd.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
            
            ID3D11RenderTargetView *rtv;
            HRESULT hr = device->CreateRenderTargetView(pTexture, &rd, &rtv);
            if (FAILED(hr)) {
                throw Coercri::DXError("Failed to create render target view", hr);
            }
            out_rtv->reset(rtv);
        }
    }
    

    // Creates a texture, with optional initial data,
    // and creates optional shader resource view, and optional render target view
    void CreateTexture(ID3D11Device * device,
                       int width,
                       int height,
                       const D3D11_SUBRESOURCE_DATA *initial_data,
                       DXGI_FORMAT format,
                       bool staging,
                       Coercri::ComPtrWrapper<ID3D11Texture2D> &out_tex,
                       Coercri::ComPtrWrapper<ID3D11ShaderResourceView> *out_srv,
                       Coercri::ComPtrWrapper<ID3D11RenderTargetView> *out_rtv)
    {
        D3D11_TEXTURE2D_DESC td;
        memset(&td, 0, sizeof(td));
        td.Width = width;
        td.Height = height;
        td.MipLevels = 1;
        td.ArraySize = 1;
        td.Format = format;
        td.SampleDesc.Count = 1;
        td.Usage = D3D11_USAGE_DEFAULT;

        if (out_srv) td.BindFlags = D3D11_BIND_SHADER_RESOURCE;
        if (out_rtv) td.BindFlags |= D3D11_BIND_RENDER_TARGET;

        if (staging) {
            td.Usage = D3D11_USAGE_STAGING;
            td.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        }

        CreateTextureImpl(device, td, initial_data, out_tex, out_srv, out_rtv);
    }

    int GetNumMipLevels(int width, int height)
    {
        int num_levels = 1;
        while (width > 1 && height > 1) {
            width /= 2;
            height /= 2;
            ++num_levels;
        }
        return num_levels;
    }

    int TexIndex(int x, int y, int component, int width)
    {
        return (y * width + x)*4 + component;
    }
    
    // creates an IMMUTABLE texture and auto generates the mip maps.
    void CreateTextureWithMips(ID3D11Device *device,
                               int width,
                               int height,
                               const unsigned char * initial_data,  // 32-bit RGBA (8-bit per channel)
                               Coercri::ComPtrWrapper<ID3D11Texture2D> &out_tex,
                               Coercri::ComPtrWrapper<ID3D11ShaderResourceView> &out_view)
    {
        const int num_mip_levels = GetNumMipLevels(width, height);

        D3D11_TEXTURE2D_DESC td;
        memset(&td, 0, sizeof(td));
        td.Width = width;
        td.Height = height;
        td.MipLevels = num_mip_levels;
        td.ArraySize = 1;
        td.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
        td.SampleDesc.Count = 1;
        td.Usage = D3D11_USAGE_IMMUTABLE;
        td.BindFlags = D3D11_BIND_SHADER_RESOURCE;

        boost::scoped_array<D3D11_SUBRESOURCE_DATA> srd(new D3D11_SUBRESOURCE_DATA[num_mip_levels]);
        boost::scoped_array<std::vector<unsigned char> > pixels(new std::vector<unsigned char>[num_mip_levels]);

        srd[0].pSysMem = initial_data;
        srd[0].SysMemPitch = width * 4;

        for (int level = 1; level < num_mip_levels; ++level) {

            const int prev_width = width;
            width = std::max(width/2, 1);
            height = std::max(height/2, 1);

            pixels[level].resize(width*height*4);
            unsigned char * out = &pixels[level][0];
            
            const unsigned char * in;
            if (level == 1) {
                in = initial_data;
            } else {
                in = &pixels[level-1][0];
            }
            
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < height; ++x) {
                    for (int component = 0; component < 4; ++component) {

                        const int in1 = in[TexIndex(2*x,   2*y,   component, prev_width)];
                        const int in2 = in[TexIndex(2*x+1, 2*y,   component, prev_width)];
                        const int in3 = in[TexIndex(2*x,   2*y+1, component, prev_width)];
                        const int in4 = in[TexIndex(2*x+1, 2*y+1, component, prev_width)];

                        const int avg = (in1 + in2 + in3 + in4 + 2) / 4;

                        *out++ = avg;
                    }
                }
            }

            srd[level].pSysMem = &pixels[level][0];
            srd[level].SysMemPitch = width * 4;
        }

        CreateTextureImpl(device, td, &srd[0], out_tex, &out_view, 0);
    }
    
    void CreateSimBuffer(ID3D11Device *device, Coercri::ComPtrWrapper<ID3D11Buffer> &vert_buf, 
                         int i_left, int i_top, int i_right, int i_bottom)
    {
        // assumes viewport of (nx+4) * (ny+4),
        // and renders all pixels from (left,top) (inclusive) to (right,bottom) (exclusive).
         
        // the tex coords sent to the pixel shader are e.g. (0.5f, 0.5f) for the top left pixel in the render target.
        
        const float left = float(i_left);
        const float right = float(i_right);
        const float top = float(i_top);
        const float bottom = float(i_bottom);

        const float vertices[12] = {
            left, top,
            right, top,
            left, bottom,
            right, top,
            right, bottom,
            left, bottom
        };
                
        // create the vertex buffer
        D3D11_BUFFER_DESC bd;
        memset(&bd, 0, sizeof(bd));
        bd.ByteWidth = 12 * sizeof(float);
        bd.Usage = D3D11_USAGE_IMMUTABLE;
        bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;

        D3D11_SUBRESOURCE_DATA sd;
        memset(&sd, 0, sizeof(sd));
        sd.pSysMem = &vertices[0];
        
        ID3D11Buffer *pBuffer;
        HRESULT hr = device->CreateBuffer(&bd, &sd, &pBuffer);
        if (FAILED(hr)) {
            throw Coercri::DXError("Failed to create vertex buffer for the simulation", hr);
        }
        vert_buf.reset(pBuffer);
    }

    int RoundUpTo16(int x)
    {
        return (x + 15) & (~15);
    }

	
	void setTimeIntegrationScheme(ID3D11DeviceContext *context, ShallowWaterEngine::TimeIntegrationScheme *ts, ID3D11Buffer *cbuffer){
		
		
		TimeIntegrationConstBuffer t_cBuffer;
		t_cBuffer.tScheme = (int) (*ts);

		context->UpdateSubresource(cbuffer, 0, 0, &t_cBuffer, 0, 0);
		ID3D11Buffer *t_cst_buf = cbuffer;
	    context->PSSetConstantBuffers(1, 1, &t_cst_buf);
	}
	void st_DumpToFile(ID3D11DeviceContext *context, ID3D11Texture2D *staging, ID3D11Texture2D *tex, int nx, int ny)
	{ // Can get nx and ny using GetTextureSize(...)
		
		context->CopyResource(staging, tex);
		MapTexture m(*context, *staging);

		for (int r = 0; r < initSetting.countOfRanges; ++r){
			std::string fileName = initSetting.logPath + "/" + initSetting.logRange[r].name +".txt";	
			std::ofstream myfile;
			if (!myfile.is_open()){
				myfile.open (fileName.c_str(),std::ios::out | std::ios::app);
			}
			for (int j = initSetting.logRange[r].bottomLeft.y; j <= initSetting.logRange[r].topRight.y; ++j) {
				for (int i = initSetting.logRange[r].bottomLeft.x; i <= initSetting.logRange[r].topRight.x; ++i) {
						const char *q = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
						const float *p = reinterpret_cast<const float*>(q) + i * 4;
						myfile << i << "\t" << j << "\t" << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" << p[3] << "\n";
					
				}
			}
		}
		
		std::string fileName = initSetting.logPath + "/" + initSetting.gaugesFilename +".txt";	
		std::ofstream myfile;
		if (!myfile.is_open()){
			myfile.open (fileName.c_str(),std::ios::out | std::ios::app);
		}
		for (int g = 0; g < initSetting.countOfGauges; ++g){
			int i = initSetting.logGauges[g].x;
			int j = initSetting.logGauges[g].y;
			const char *q = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
			const float *p = reinterpret_cast<const float*>(q) + i * 4;
			myfile << i << "\t" << j << "\t" << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" << p[3] << "\n";
		}
	}



	void st_DumpToFileForDebug(ID3D11DeviceContext *context, ID3D11Texture2D *staging, ID3D11Texture2D *tex, int nx, int ny)
	{ 
		context->CopyResource(staging, tex);
		MapTexture m(*context, *staging);

		std::string fileName = initSetting.exePath + "/" + "debug.txt";	
		std::ofstream myfile;
		if (!myfile.is_open()){
			myfile.open (fileName.c_str(),std::ios::out | std::ios::app);
		}

		for (int j = 0; j < ny + 4; ++j) {
			for (int i = 0; i < nx + 4; ++i) {
					const char *q = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
					const float *p = reinterpret_cast<const float*>(q) + i * 4;
					myfile << i << "\t" << j << "\t" << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" << p[3] << "\n";
			}
		}
		
	}


    void SeaWaveSettings(char c, float &sa, float &skx, float &sky, float &so)
    {
        using std::string;
        sa = GetSetting(string("sa") + c);
        float k = GetSetting(string("sk") + c);
        float kdir = GetSetting(string("sk") + c + "_dir");
        skx = k * std::cos(kdir);
        sky = k * std::sin(kdir);
        so = GetSetting(string("so") + c);
    }
}

ShallowWaterEngine::ShallowWaterEngine(ID3D11Device *device_,
                                       ID3D11DeviceContext *context_)
    : device(device_), context(context_), current_timestep(0), total_time(0)
{



    // create D3D objects
    createShadersAndInputLayout();
    createConstantBuffers();
    loadGraphics();
    //createBlendState();
    loadSkybox(2); // 2 is the default skybox.
    setupMousePicking();
    
    remesh(R_MESH);
    newTerrainSettings();
    //moveCamera(0, 0, .5, 0, 0, 100, 100);
}

void ShallowWaterEngine::remesh(ResetType reset_type)
{
    createMeshBuffers();
    createSimBuffers();
    createTerrainTexture();
    fillTerrainTextureLite();
	//ST_: I added createShadersAndInputLayout to make all the shaders each time we remesh. 
	createShadersAndInputLayout();
	createSimTextures(reset_type);
	
	myTimestep_count = 0;
	total_time = 0;
	
}
int ShallowWaterEngine::st_getCountOfMatrices(int n){   //returns the number of matrices for Cyclic Reduction
    int countOfMatrices = 1;
	while (n > 1) {
		n /= 2;
		++countOfMatrices;
	}
	return countOfMatrices;
}

void ShallowWaterEngine::dumpBathymetryToFile()
{
	const int nx = GetIntSetting("mesh_size_x");
	const int ny = GetIntSetting("mesh_size_y");

	time_t timer = time(0);
	struct tm * now = localtime(&timer);
	std::ostringstream s;
	s << now->tm_year + 1900 << "-" << now->tm_mon + 1 << "-" << now->tm_mday << " " << now->tm_hour << "-" << now->tm_min << "-" << now->tm_sec;
	

	std::string fileName = initSetting.logPath + "/" + s.str() +".cbf";	
	std::ofstream myfile;
	if (!myfile.is_open()){
		myfile.open (fileName.c_str(),std::ios::out | std::ios::app);
	}
	myfile << "[nx] " << nx << "\n";
	myfile << "[ny] " << ny << "\n";
	myfile << "\n" << "\n";
	myfile << "====================================" << "\n";

	for (int j = 2; j < ny + 2; ++j) {
		for (int i = 2; i < nx + 2; ++i) {
			float z = g_bottom[(nx+4) * j + i].BA;
			myfile << z << "\t" ;
		}
		myfile << "\n";
	}
	myfile.close();
}

void ShallowWaterEngine::dumpInundationToFile(ID3D11DeviceContext *context, ID3D11Texture2D *staging, ID3D11Texture2D *tex)
{
 	const int nx = GetIntSetting("mesh_size_x");
	const int ny = GetIntSetting("mesh_size_y");
		
	context->CopyResource(staging, tex);
	MapTexture m(*context, *staging);


	time_t timer = time(0);
	struct tm * now = localtime(&timer);
	std::ostringstream s;
	s << now->tm_year + 1900 << "-" << now->tm_mon + 1 << "-" << now->tm_mday << " " << now->tm_hour << "-" << now->tm_min << "-" << now->tm_sec;
	

	std::string fileName = initSetting.logPath + "/" + s.str() +"-inundation.txt";	
	std::ofstream myfile;
	if (!myfile.is_open()){
		myfile.open (fileName.c_str(),std::ios::out | std::ios::app);
	}

	myfile << "[nx] " << nx << "\n";
	myfile << "[ny] " << ny << "\n";
	myfile << "\n" << "\n";
	myfile << "# This file contains the maximum depth of inundation on each point, in meters.";
	myfile << "\n" << "\n";
	myfile << "====================================" << "\n";

	for (int j = 2; j < ny + 2; ++j) {
		for (int i = 2; i <= nx + 2; ++i) {
			const char *q = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
			const float *p = reinterpret_cast<const float*>(q) + i * 4;
			myfile << p[0] << "\t" ;
		}
		myfile << "\n";
	}
}


void ShallowWaterEngine::newTerrainSettings()
{
    fillTerrainTexture();
}

void ShallowWaterEngine::moveCamera(float x, float y, float z, float yaw_, float pitch_, int vw, int vh)
{
    pitch = pitch_;
    yaw = yaw_;
    camera_x = x;
    camera_y = y;
    camera_z = z;
    
    vp_width = vw;
    vp_height = vh;

    near_plane_dist = std::min(GetSetting("valley_width")/(GetSetting("mesh_size_x")-1)/2.0f, GetSetting("valley_length")/(GetSetting("mesh_size_y")-1)/2.0f);;
    far_plane_dist = 5 * std::max(GetSetting("valley_width"), GetSetting("valley_length"));

    const float fov = GetSetting("fov");  // x fov in degrees
    const float fov_r = fov / 45.0f * std::atan(1.0f);
    xpersp = 1.0f / tan(fov_r * 0.5f);
    ypersp = float(vw)/float(vh) * xpersp;
    //ypersp = 1.0f / tan(fov_r * 0.5f);
    //xpersp = float(vw)/float(vh) * ypersp;
    
    fillConstantBuffers();
}    

void ShallowWaterEngine::timestep()
{


    ///// Allow only a certain number of timesteps for debugging
#if 0

	static bool debug_flag = true;
	
	if(timestep_count == 0){
		debug_flag = true;	
	}
	if (debug_flag) return;
#endif
    /////////////////////////////////////

  
    /*
	boost::shared_ptr<Coercri::Timer> timer(new Coercri::GenericTimer);
	
	float tr=GetSetting("time_ratio");
	float ta=GetSetting("time_acceleration");
	float myAlpha=(tr-ta);
	if (myAlpha<1){
		myAlpha=abs(myAlpha)*myAlpha;
	}
	
	myDelay+=myAlpha*(realworld_timestep-myDelay)/100.0f;
	if(myDelay<0){
		myDelay=0;
	}
	*/
	//timer->sleepMsec(1000.0f*myDelay);
	
    // Get some resource pointers
    ID3D11Buffer * cst_buf = m_psSimConstantBuffer.get();

    ID3D11ShaderResourceView * bottom_tex = m_psBottomTextureView.get();
    ID3D11ShaderResourceView * new_state_or_h_tex = m_psSimTextureView[1 - sim_idx].get();
	ID3D11ShaderResourceView * old_state_tex = m_psSimTextureView[sim_idx].get();
	ID3D11ShaderResourceView * old_gradients_tex = m_psSimTextureView[old_index].get();
	ID3D11ShaderResourceView * old_old_gradients_tex = m_psSimTextureView[old_old_index].get();
	ID3D11ShaderResourceView * predicted_gradients_tex = m_psSimTextureView[predicted_index].get();
	ID3D11ShaderResourceView * scratch_gradients_tex = m_psSimTextureView[scratch_index].get();
	ID3D11ShaderResourceView * F_G_star_old_gradients_tex = m_psSimTextureView[F_G_star_old_index].get();
	ID3D11ShaderResourceView * F_G_star_old_old_gradients_tex = m_psSimTextureView[F_G_star_old_old_index].get();
	ID3D11ShaderResourceView * F_G_star_scratch_gradients_tex = m_psSimTextureView[F_G_star_scratch_index].get();
	ID3D11ShaderResourceView * F_G_star_predicted_gradients_tex = m_psSimTextureView[F_G_star_predicted_index].get();
	ID3D11ShaderResourceView * st_RHS_tex = st_psRHSTextureView.get();
	ID3D11ShaderResourceView * u_tex = m_psSimTextureView[2].get();
    ID3D11ShaderResourceView * v_tex = m_psSimTextureView[3].get();
    ID3D11ShaderResourceView * xflux_tex = m_psSimTextureView[4].get();
    ID3D11ShaderResourceView * yflux_tex = m_psSimTextureView[5].get();
	ID3D11ShaderResourceView * normal_tex = m_psSimTextureView[6].get();     // {nX, nY, nZ, unused}
	ID3D11ShaderResourceView * auxiliary1_tex = m_psAuxiliary1TextureView.get();
	ID3D11ShaderResourceView * auxiliary2_tex = m_psAuxiliary2TextureView.get();
    
	ID3D11RenderTargetView * new_state_or_h_target = m_psSimRenderTargetView[1 - sim_idx].get();
	ID3D11RenderTargetView * old_state_or_h_target = m_psSimRenderTargetView[sim_idx].get();
    ID3D11RenderTargetView * old_gradients_target = m_psSimRenderTargetView[old_index].get();
	ID3D11RenderTargetView * old_old_gradients_target = m_psSimRenderTargetView[old_old_index].get();
    ID3D11RenderTargetView * scratch_gradients_target = m_psSimRenderTargetView[scratch_index].get();
	ID3D11RenderTargetView * predicted_gradients_target = m_psSimRenderTargetView[predicted_index].get();
	ID3D11RenderTargetView * F_G_star_old_gradients_target = m_psSimRenderTargetView[F_G_star_old_index].get();
	ID3D11RenderTargetView * F_G_star_old_old_gradients_target = m_psSimRenderTargetView[F_G_star_old_old_index].get();
    ID3D11RenderTargetView * F_G_star_scratch_gradients_target = m_psSimRenderTargetView[F_G_star_scratch_index].get();
	ID3D11RenderTargetView * F_G_star_predicted_gradients_target = m_psSimRenderTargetView[F_G_star_predicted_index].get();
	ID3D11RenderTargetView * st_RHS_target = st_psRHSRenderTargetView.get();
    ID3D11RenderTargetView * u_target = m_psSimRenderTargetView[2].get();
    ID3D11RenderTargetView * v_target = m_psSimRenderTargetView[3].get();
    ID3D11RenderTargetView * xflux_target = m_psSimRenderTargetView[4].get();
    ID3D11RenderTargetView * yflux_target = m_psSimRenderTargetView[5].get();
    ID3D11RenderTargetView * normal_target = m_psSimRenderTargetView[6].get();
	ID3D11RenderTargetView * auxiliary1_target = m_psAuxiliary1RenderTargetView.get();
	ID3D11RenderTargetView * auxiliary2_target = m_psAuxiliary2RenderTargetView.get();
		              

    // Common Settings
    
    context->ClearState();

    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->IASetInputLayout(m_psSimInputLayout.get());

    context->VSSetShader(m_psSimVertexShader.get(), 0, 0);
    
    D3D11_VIEWPORT vp;
    memset(&vp, 0, sizeof(vp));
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    vp.Width = float(nx + 4);
    vp.Height = float(ny + 4);
    vp.MinDepth = 0;
    vp.MaxDepth = 1;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    context->RSSetViewports(1, &vp);

    context->VSSetConstantBuffers(0, 1, &cst_buf);
    context->PSSetConstantBuffers(0, 1, &cst_buf);
	setTimeIntegrationScheme(context,&timeScheme,m_psTimeIntegrationConstantBuffer.get());

    ID3D11Buffer *vert_buf;
    const UINT stride = 8;
    const UINT offset = 0;

    ID3D11ShaderResourceView *pNULL = 0;

	const int countOfMatricesXdirection = st_getCountOfMatrices(nx);
	const int countOfMatricesYdirection = st_getCountOfMatrices(ny);

    if (bootstrap_needed) {
        // Pass 1
        // read: old_state; bottom
        // write: h, u, v
        
        vert_buf = m_psSimVertexBuffer11.get();
        context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
        
        ID3D11RenderTargetView * p1_tgt[] = {new_state_or_h_target, u_target, v_target, normal_target, auxiliary2_target};
        context->OMSetRenderTargets(5, &p1_tgt[0], 0);
        
        context->PSSetShader(m_psSimPixelShader[0].get(), 0, 0);
        context->PSSetShaderResources(0, 1, &old_state_tex);
        context->PSSetShaderResources(1, 1, &bottom_tex);
		context->PSSetShaderResources(3, 1, &auxiliary1_tex);
	        
        context->Draw(6, 0);    

        bootstrap_needed = false;
    }
	

    // Pass 2
    // read: h, u, v
    // write: xflux, yflux


	// xflux and yflux are used as scratch texture elsewhere. Here we reset them to (0,0,0,0). Note that D1x[0] and D1y[0] are never used in this code, so they are all (0,0,0,0).
	context->CopyResource(m_psSimTexture[4].get(),
						  st_psTriDigTexture_D1x[0].get()); 
	context->CopyResource(m_psSimTexture[5].get(),
						  st_psTriDigTexture_D1y[0].get());
	
    vert_buf = m_psSimVertexBuffer10.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
	
	context->PSSetShaderResources(0, 1, &new_state_or_h_tex);
    context->PSSetShaderResources(1, 1, &u_tex);
    context->PSSetShaderResources(2, 1, &v_tex);
	context->PSSetShaderResources(3, 1, &normal_tex);
	context->PSSetShaderResources(4, 1, &auxiliary2_tex);

	

    ID3D11RenderTargetView * p2_tgt[] = {xflux_target, yflux_target, auxiliary1_target, 0};
    context->OMSetRenderTargets(4, &p2_tgt[0], 0);

    context->PSSetShader(m_psSimPixelShader[1].get(), 0, 0);


    context->Draw(6, 0);


	// Pass 3
    // read: old_state, bottom, xflux, yflux
    // write: new_state

	// reset new_state_or_h to zero
	context->CopyResource(m_psSimTexture[1-sim_idx].get(),
						  st_psTriDigTexture_D1x[0].get()); //st_psTriDigTexture_D1x is all zeros!!!



	vert_buf = m_psSimVertexBuffer00.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);

    context->PSSetShaderResources(0, 1, &pNULL); // unbind new_state_or_h_target so we can set it as output.

    ID3D11RenderTargetView * p3_tgt[] = {new_state_or_h_target, scratch_gradients_target, F_G_star_scratch_gradients_target}; 
    context->OMSetRenderTargets(3, &p3_tgt[0], 0);
	context->PSSetShader(m_psSimPixelShader[2].get(), 0, 0);


	//TODO: remove these, taking care if they are used anywhere else.
	if(firstTimeStep){
		firstTimeStep = false;
		secondTimeStep = true;
	}
	else if (secondTimeStep){
		secondTimeStep = false;
	}

	

	ID3D11ShaderResourceView * ABCx_tex = st_psTriDigTextureView_ABCx[0].get();
	ID3D11ShaderResourceView * ABCy_tex = st_psTriDigTextureView_ABCy[0].get();
    context->PSSetShaderResources(0, 1, &old_state_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);	
    context->PSSetShaderResources(2, 1, &xflux_tex);
    context->PSSetShaderResources(3, 1, &yflux_tex);
	context->PSSetShaderResources(4, 1, &old_gradients_tex);
	context->PSSetShaderResources(5, 1, &old_old_gradients_tex);
	context->PSSetShaderResources(6, 1, &ABCx_tex); 
	context->PSSetShaderResources(7, 1, &ABCy_tex);
	context->PSSetShaderResources(8, 1, &F_G_star_old_gradients_tex);
	context->PSSetShaderResources(9, 1, &F_G_star_old_old_gradients_tex);
	if(timeScheme == corrector)
	{
		context->PSSetShaderResources(10, 1, &predicted_gradients_tex);
		context->PSSetShaderResources(11, 1, &F_G_star_predicted_gradients_tex);
		ID3D11ShaderResourceView * current_state = m_psSimTextureView[CURRENT_STATE].get();
		context->PSSetShaderResources(12, 1, &current_state);
	}

    context->Draw(6, 0);
	
	
	if (timeScheme == predictor){
		context->CopyResource(m_psSimTexture[predicted_index].get(),
							  m_psSimTexture[scratch_index].get());
		context->CopyResource(m_psSimTexture[F_G_star_predicted_index].get(),
							  m_psSimTexture[F_G_star_scratch_index].get());

		context->CopyResource(m_psSimTexture[CURRENT_STATE].get(),
							  m_psSimTexture[sim_idx].get());
	}

	// Run cyclic reduce on D1x (0) and get D1x(1). D1x(0) is read from new_state_or_h_tex
	//Note that here we have already set ABCx_tex = st_psTriDigTextureView_ABCx[0].get();
	
	ID3D11RenderTargetView * D1x_target = st_psTriDigRenderTargetView_D1x[1].get();

    vert_buf = st_psCyclicReductionX_VertexBuffer[1].get();
	context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
	ID3D11RenderTargetView * p4_tgt[] = {D1x_target,0};
	context->OMSetRenderTargets(2, &p4_tgt[0], 0);
	context->PSSetShader(st_psCyclicReduceDx_PixelShader.get(), 0, 0);
	context->PSSetShaderResources(6, 1, &ABCx_tex);
	context->PSSetShaderResources(1, 1, &new_state_or_h_tex);
	
	context->Draw(6, 0);

	// Run cyclic reduce on rest of D1x (i) and get D1x(i+1)

	for(int i=2;i<countOfMatricesXdirection;++i){
		
		//getting some pointers to resources and targets!
		ID3D11ShaderResourceView * ABCx_tex = st_psTriDigTextureView_ABCx[i-1].get();
		ID3D11ShaderResourceView * D1x_tex = st_psTriDigTextureView_D1x[i-1].get();

		ID3D11RenderTargetView * D1x_target = st_psTriDigRenderTargetView_D1x[i].get();

        vert_buf = st_psCyclicReductionX_VertexBuffer[i].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		ID3D11RenderTargetView * p4_tgt[] = {D1x_target,0};
		context->OMSetRenderTargets(2, &p4_tgt[0], 0);
		context->PSSetShader(st_psCyclicReduceDx_PixelShader.get(), 0, 0);
		context->PSSetShaderResources(6, 1, &ABCx_tex);
		context->PSSetShaderResources(1, 1, &D1x_tex);
		
		context->Draw(6, 0);
    }
    
	// Run cyclic reduce on D1y(0) and get D1y(1). D1y(0) is read from new_state_or_h_tex
	// Note that here we have ABCy_tex = st_psTriDigTextureView_ABCy[0].get();

	ID3D11RenderTargetView * D1y_target = st_psTriDigRenderTargetView_D1y[1].get();
	
    vert_buf = st_psCyclicReductionY_VertexBuffer[1].get();
	context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
	ID3D11RenderTargetView * D1y_tgt[] = {D1y_target,0};
	context->OMSetRenderTargets(2, &D1y_tgt[0], 0);
	context->PSSetShader(st_psCyclicReduceDy_PixelShader.get(), 0, 0);
	context->PSSetShaderResources(7, 1, &ABCy_tex);
	context->PSSetShaderResources(1, 1, &new_state_or_h_tex);
	
	context->Draw(6, 0);

	// Run cyclic reduce on rest of D1y (i) and get D1y(i+1)

	for(int i=2;i<countOfMatricesYdirection;++i){
		
		//getting some pointers to resources and targets!
		ID3D11ShaderResourceView * ABCy_tex = st_psTriDigTextureView_ABCy[i-1].get();
		ID3D11ShaderResourceView * D1y_tex = st_psTriDigTextureView_D1y[i-1].get();

		ID3D11RenderTargetView * D1y_target = st_psTriDigRenderTargetView_D1y[i].get();

        vert_buf = st_psCyclicReductionY_VertexBuffer[i].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		ID3D11RenderTargetView * p4_tgt[] = {D1y_target,0};
		context->OMSetRenderTargets(2, &p4_tgt[0], 0);
		context->PSSetShader(st_psCyclicReduceDy_PixelShader.get(), 0, 0);
		context->PSSetShaderResources(7, 1, &ABCy_tex);
		context->PSSetShaderResources(1, 1, &D1y_tex);
		
		context->Draw(6, 0);
    }
    
		
	//Backward Substitution
	//solving for x, the first step. x direction
	//getting some pointers to resources and targets!
	int myIndex=countOfMatricesXdirection-1;
	ABCx_tex = st_psTriDigTextureView_ABCx[myIndex].get();
	ID3D11ShaderResourceView * D1x_tex = st_psTriDigTextureView_D1x[myIndex].get();
	ID3D11ShaderResourceView * Xx_tex = st_psTriDigTextureView_D1x[myIndex+1].get(); //This is the first step, so we dont have X(i+1), but since a=c=0 for this step, it is safe to pass anything. Try passing null.
	ID3D11RenderTargetView * Xx_target = st_psTriDigRenderTargetView_D1xCopy[myIndex].get();

	vert_buf = st_psCyclicReductionX_VertexBuffer[myIndex].get();
	context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
	ID3D11RenderTargetView * p5_tgt[] = {Xx_target,0};
	context->OMSetRenderTargets(2, &p5_tgt[0], 0);
	context->PSSetShader(st_psCyclicSubstituteXx_PixelShader.get(), 0, 0);
	context->PSSetShaderResources(6, 1, &ABCx_tex);
	context->PSSetShaderResources(1, 1, &D1x_tex);
	context->PSSetShaderResources(2, 1, &Xx_tex);
	context->Draw(6, 0);

	// solving for x, the rest of steps.
	for(int i=countOfMatricesXdirection-1-1;i>=0;--i){
	//getting some pointers to resources and targets!
		ABCx_tex = st_psTriDigTextureView_ABCx[i].get();
		if(i!=0){
			D1x_tex = st_psTriDigTextureView_D1x[i].get();
		}
		else{
			D1x_tex = m_psSimTextureView[1 - sim_idx].get();
		}

		Xx_tex = st_psTriDigTextureView_D1xCopy[i+1].get(); // this is i+1, because we pass the Xx from the previous step as the resource. 
		Xx_target = st_psTriDigRenderTargetView_D1xCopy[i].get();

		vert_buf = st_psCyclicReductionX_VertexBuffer[i].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		ID3D11RenderTargetView * p6_tgt[] = {Xx_target,0 };
		context->OMSetRenderTargets(2, &p6_tgt[0], 0);
		context->PSSetShader(st_psCyclicSubstituteXx_PixelShader.get(), 0, 0); // shader is set before, we don't need to do it here. done for clarity.
		context->PSSetShaderResources(6, 1, &ABCx_tex);
		context->PSSetShaderResources(1, 1, &D1x_tex);
		context->PSSetShaderResources(2, 1, &Xx_tex);
		context->Draw(6, 0);
    }
	// Backward Substitution
	// solving for X, the first step. y direction
	//getting some pointers to resources and targets!
	myIndex=countOfMatricesYdirection-1;
	ABCy_tex = st_psTriDigTextureView_ABCy[myIndex].get();
	ID3D11ShaderResourceView * D1y_tex = st_psTriDigTextureView_D1y[myIndex].get();
	ID3D11ShaderResourceView * Xy_tex = st_psTriDigTextureView_D1y[myIndex+1].get(); //This is the first step, so we dont have X(i+1), but since a=c=0 for this step, it is safe to pass anything. Try passing null.
	ID3D11RenderTargetView * Xy_target = st_psTriDigRenderTargetView_D1yCopy[myIndex].get();

	vert_buf = st_psCyclicReductionY_VertexBuffer[myIndex].get();
	context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
	ID3D11RenderTargetView * Xy_tgt[] = {Xy_target,0};
	context->OMSetRenderTargets(2, &Xy_tgt[0], 0);
	context->PSSetShader(st_psCyclicSubstituteXy_PixelShader.get(), 0, 0);
	context->PSSetShaderResources(7, 1, &ABCy_tex);
	context->PSSetShaderResources(1, 1, &D1y_tex);
	context->PSSetShaderResources(3, 1, &Xy_tex);
	context->Draw(6, 0);

	// solving for X, the rest of steps in y direction.
	for(int i=countOfMatricesYdirection-1-1;i>=0;--i){
	//getting some pointers to resources and targets!
		ABCy_tex = st_psTriDigTextureView_ABCy[i].get();

		if(i!=0){
			D1y_tex = st_psTriDigTextureView_D1y[i].get();
		}
		else{
			D1y_tex = m_psSimTextureView[1 - sim_idx].get();
		}

		Xy_tex = st_psTriDigTextureView_D1yCopy[i+1].get(); // this is i+1, because we pass the Xy from the previous step as the resource. 
		Xy_target = st_psTriDigRenderTargetView_D1yCopy[i].get();

		vert_buf = st_psCyclicReductionY_VertexBuffer[i].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		ID3D11RenderTargetView * p6_tgt[] = {Xy_target,0 };
		context->OMSetRenderTargets(2, &p6_tgt[0], 0);
		context->PSSetShader(st_psCyclicSubstituteXy_PixelShader.get(), 0, 0); // shader is set before, we don't need to do it here. done for clarity.
		context->PSSetShaderResources(7, 1, &ABCy_tex);
		context->PSSetShaderResources(1, 1, &D1y_tex);
		context->PSSetShaderResources(3, 1, &Xy_tex);
		context->Draw(6, 0);
    }

//////////////////////////////
	// Copy values found for X (hu,hv) from Xx and Xy back to new_state_or_h_tex.
	Xx_tex = st_psTriDigTextureView_D1xCopy[0].get(); 
	Xy_tex = st_psTriDigTextureView_D1yCopy[0].get(); 
	vert_buf =st_psCyclicReductionX_VertexBuffer[0].get(); // [0] vertexbuffer is the same for x and y directions. so we can use either one.
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
	
	context->CopyResource(st_psRHSTexture.get(),
						  m_psSimTexture[1 - sim_idx].get());  // new_state_or_h. We can not set the target and the resource to new_state_or_h at the same time. Here I make a copy.
	context->OMSetRenderTargets(2, &p3_tgt[0], 0);  //this is the p3_tgt, defined before. p3_tgt[] = { new_state_or_h_target, 0 }
    context->PSSetShader(st_psCopyFromXxandXyPixelShader.get(), 0, 0);
	context->PSSetShaderResources(0, 1, &st_RHS_tex);
	context->PSSetShaderResources(2, 1, &Xx_tex);
	context->PSSetShaderResources(3, 1, &Xy_tex);
	context->Draw(6, 0);

    // Boundary Conditions
	createBoundaryShaders();

    BoundaryConstBuffer cb;
	cb.PI = PI;
    cb.boundary_epsilon = CalcEpsilon();
	cb.boundary_nx = nx + 4;
	cb.boundary_ny = ny + 4;
	cb.dx = GetSetting("valley_width") / (nx-1);
	cb.dy = GetSetting("valley_length") / (ny-1);
	
	cb.northBoundaryWidth = initSetting.northBoundary.width;
	cb.eastBoundaryWidth  = initSetting.eastBoundary.width;
	cb.westBoundaryWidth  = initSetting.westBoundary.width;
	cb.southBoundaryWidth = initSetting.southBoundary.width;

	cb.northSeaLevel = initSetting.northBoundary.waterLevel;
	cb.eastSeaLevel  = initSetting.eastBoundary.waterLevel;
	cb.westSeaLevel  = initSetting.westBoundary.waterLevel;
	cb.southSeaLevel = initSetting.southBoundary.waterLevel;
	

	if (initSetting.eastBoundary.type == "UniformTimeSeries"){
		float temp_eta = 0, temp_hu = 0, temp_hv = 0;
		eastUniformTimeSeries(temp_eta, temp_hu, temp_hv, total_time);

		cb.eastAmplitude_or_eta = temp_eta;
		cb.eastPeriod_or_hu = temp_hu;
		cb.eastTheta_or_hv = temp_hv;
	} else {
		cb.eastAmplitude_or_eta = initSetting.eastBoundary.sineWaveSetting.amplitude;
		cb.eastPeriod_or_hu = initSetting.eastBoundary.sineWaveSetting.period;
		cb.eastTheta_or_hv = initSetting.eastBoundary.sineWaveSetting.theta;
	}

	if (initSetting.westBoundary.type == "UniformTimeSeries"){
		float temp_eta = 0, temp_hu = 0, temp_hv = 0;
		westUniformTimeSeries(temp_eta, temp_hu, temp_hv, total_time);

		cb.westAmplitude_or_eta = temp_eta;
		cb.westPeriod_or_hu = temp_hu;
		cb.westTheta_or_hv = temp_hv;
	} else {
		cb.westAmplitude_or_eta = initSetting.westBoundary.sineWaveSetting.amplitude;
		cb.westPeriod_or_hu = initSetting.westBoundary.sineWaveSetting.period;
		cb.westTheta_or_hv = initSetting.westBoundary.sineWaveSetting.theta;
	}

	if (initSetting.northBoundary.type == "UniformTimeSeries"){
		float temp_eta = 0, temp_hu = 0, temp_hv = 0;
		northUniformTimeSeries(temp_eta, temp_hu, temp_hv, total_time);

		cb.northAmplitude_or_eta = temp_eta;
		cb.northPeriod_or_hu = temp_hu;
		cb.northTheta_or_hv = temp_hv;
	} else {
		cb.northAmplitude_or_eta = initSetting.northBoundary.sineWaveSetting.amplitude;
		cb.northPeriod_or_hu = initSetting.northBoundary.sineWaveSetting.period;
		cb.northTheta_or_hv = initSetting.northBoundary.sineWaveSetting.theta;
	}

	
	if (initSetting.southBoundary.type == "UniformTimeSeries"){
		float temp_eta = 0, temp_hu = 0, temp_hv = 0;
		southUniformTimeSeries(temp_eta, temp_hu, temp_hv, total_time);

		cb.southAmplitude_or_eta = temp_eta;
		cb.southPeriod_or_hu = temp_hu;
		cb.southTheta_or_hv = temp_hv;
	} else {
		cb.southAmplitude_or_eta = initSetting.southBoundary.sineWaveSetting.amplitude;
		cb.southPeriod_or_hu = initSetting.southBoundary.sineWaveSetting.period;
		cb.southTheta_or_hv = initSetting.southBoundary.sineWaveSetting.theta;
	}

    cb.reflect_x = 2*nx+2;
	cb.reflect_y = 2*ny+2;
    //cb.solid_wall_flag = (GetIntSetting("solid_walls") != 0);
	cb.sea_level =  initSetting.isBoussinesq;

    //const float inflow_width = GetSetting("inflow_width");
    //const float inflow_height = GetSetting("inflow_height");
    const float W = GetSetting("valley_width");
    const float dx = GetSetting("valley_width") / (nx-1);

    //cb.inflow_x_min = int((inflow_width + W/2)/dx) + 1;
    //cb.inflow_x_max = int((inflow_width + W/2)/dx) + 3;
    //cb.inflow_height = inflow_height;
    cb.inflow_speed = 0.01f;
    cb.g = GetSetting("gravity");
    cb.total_time = total_time;
    //SeaWaveSettings('1', cb.sa1, cb.skx1, cb.sky1, cb.so1);
    //SeaWaveSettings('2', cb.sa2, cb.skx2, cb.sky2, cb.so2);
    //SeaWaveSettings('3', cb.sa3, cb.skx3, cb.sky3, cb.so3);
    //SeaWaveSettings('4', cb.sa4, cb.skx4, cb.sky4, cb.so4);
    //cb.sdecay = 0.01f / ny * GetSetting("valley_length");


	const int NORTH = 0, EAST = 1, SOUTH = 2, WEST = 3;
	D3D11_BOX src_box;
	src_box.front = 0; //not used
    src_box.back = 1; //not used

    context->UpdateSubresource(m_psBoundaryConstantBuffer.get(), 0, 0, &cb, 0, 0);
    ID3D11Buffer *b_cst_buf = m_psBoundaryConstantBuffer.get();
    context->PSSetConstantBuffers(0, 1, &b_cst_buf);

	for (int i = 0; i < 3; ++i){
		context->PSSetShaderResources(1+i, 1, &pNULL);
	}
	

	
	//east
    // use XFLUX as scratch space, then we'll copy back to the main output 
    context->OMSetRenderTargets(1, &xflux_target, 0);

    // now rebind the input as the output from the previous step (ie the new state).
    context->PSSetShaderResources(0, 1, &new_state_or_h_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);

	if (initSetting.eastBoundary.type == "IrregularWaves"){
		eastIrregularBoundary();
	} else {
		vert_buf = m_psBoundaryVertexBuffer[EAST].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		context->PSSetShader(m_psBoundaryPixelShader[EAST].get(), 0, 0);
		context->Draw(6, 0);
		src_box.top = 0;
		src_box.bottom = ny + 4;
		src_box.left = nx+4 - initSetting.eastBoundary.width;
		src_box.right = nx+4;
		context->CopySubresourceRegion(m_psSimTexture[1-sim_idx].get(), 0, nx+4 - initSetting.eastBoundary.width, 0, 0, m_psSimTexture[4].get(), 0, &src_box);
	}

	//north
	// use XFLUX as scratch space, then we'll copy back to the main output 
	context->OMSetRenderTargets(1, &xflux_target, 0);

	// now rebind the input as the output from the previous step (ie the new state).
	context->PSSetShaderResources(0, 1, &new_state_or_h_tex);
	context->PSSetShaderResources(1, 1, &bottom_tex);

	if (initSetting.northBoundary.type == "IrregularWaves"){
		northIrregularBoundary();
	} else {
		vert_buf = m_psBoundaryVertexBuffer[NORTH].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		context->PSSetShader(m_psBoundaryPixelShader[NORTH].get(), 0, 0);
		context->Draw(6, 0);
		src_box.left = 0;
		src_box.right = nx + 4;
		src_box.top = ny+4 - initSetting.northBoundary.width;
		src_box.bottom = ny+4;
		context->CopySubresourceRegion(m_psSimTexture[1-sim_idx].get(), 0, 0, ny+4 - initSetting.northBoundary.width , 0, m_psSimTexture[4].get(), 0, &src_box);
	}

	//west
	// use XFLUX as scratch space, then we'll copy back to the main output 
	context->OMSetRenderTargets(1, &xflux_target, 0);

	// now rebind the input as the output from the previous step (ie the new state).
	context->PSSetShaderResources(0, 1, &new_state_or_h_tex);
	context->PSSetShaderResources(1, 1, &bottom_tex);

	if (initSetting.westBoundary.type == "IrregularWaves"){
 		westIrregularBoundary();
 	} else { 
		context->PSSetShader(m_psBoundaryPixelShader[WEST].get(), 0, 0); 
		vert_buf = m_psBoundaryVertexBuffer[WEST].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		context->Draw(6, 0);
		src_box.left = 0;
		src_box.right = initSetting.westBoundary.width;
		src_box.top = 0;
		src_box.bottom = ny + 4;
		context->CopySubresourceRegion(m_psSimTexture[1-sim_idx].get(), 0, 0, 0, 0, m_psSimTexture[4].get(), 0, &src_box);
 	}

	
	//south 
	// use XFLUX as scratch space, then we'll copy back to the main output 
	context->OMSetRenderTargets(1, &xflux_target, 0);

	// now rebind the input as the output from the previous step (ie the new state).
	context->PSSetShaderResources(0, 1, &new_state_or_h_tex);
	context->PSSetShaderResources(1, 1, &bottom_tex);

	if (initSetting.southBoundary.type == "IrregularWaves"){
		southIrregularBoundary();
	} else {
		vert_buf = m_psBoundaryVertexBuffer[SOUTH].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		context->PSSetShader(m_psBoundaryPixelShader[SOUTH].get(), 0, 0);
		context->Draw(6, 0);
		src_box.left = 0;
		src_box.right = nx + 4;
		src_box.top = 0;
		src_box.bottom = initSetting.southBoundary.width;
		context->CopySubresourceRegion(m_psSimTexture[1 - sim_idx].get(),
									   0,  // subresource
									   0,  // dest x
									   0,   // dest y
									   0,  // dest z
									   m_psSimTexture[4].get(),  // xflux tex
									   0, // subresource
									   &src_box);
	}	

    // Now do "pass 1" again, this means the H, U, V textures will be ready for the
    // next timestep. 
    // Also, the Normal texture will be created at this point.

    // first need to unbind 'old_state' from the pixel shader (as it is the target for our 'H' texture)
    context->PSSetShaderResources(0, 1, &pNULL);

    context->PSSetConstantBuffers(0, 1, &cst_buf);
    
    vert_buf = m_psSimVertexBuffer11.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
    ID3D11RenderTargetView * p1_tgt[] = { old_state_or_h_target, u_target, v_target, normal_target, auxiliary2_target};
    context->OMSetRenderTargets(5, &p1_tgt[0], 0);
    context->PSSetShader(m_psSimPixelShader[0].get(), 0, 0);
    context->PSSetShaderResources(0, 1, &new_state_or_h_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);
	context->PSSetShaderResources(3, 1, &auxiliary1_tex);
    context->Draw(6, 0);


#if 1
	/*
	static bool bathymetry = true;
	if (initSetting.doLog == true && bathymetry == true){
		bathymetry = false;
		//char temp[100];
		//sprintf(temp, "/state%d.txt",myTimestep_count);
		std::string fileName = initSetting.logPath + "/bathymetry.txt";
		
		std::ofstream myfile;
		if (!myfile.is_open()){
			myfile.open (fileName.c_str(),std::ios::out | std::ios::app);
		}
		st_DumpToFile(myfile, context, m_psFullSizeStagingTexture.get(), m_psBottomTexture.get(),nx,ny);
		//myfile.close();
	}*/
#endif


	if (timestep_count >= 2) {
		if (correction_steps_count == initSetting.correctionStepsNum)
		{
			timeScheme = predictor;
			correction_steps_count = 0;
		} else {
			timeScheme = corrector;
			correction_steps_count++;
		}
		setTimeIntegrationScheme(context,&timeScheme,m_psTimeIntegrationConstantBuffer.get()) ;
	}

	if (timeScheme == predictor || timeScheme == euler) 
	{
		//// Dumping to file ////
		if (initSetting.doLog == true && timestep_count % initSetting.logStep == 0){
			st_DumpToFile(context, m_psFullSizeStagingTexture.get(), m_psSimTexture[1-sim_idx].get(),nx,ny);
			//myfile.close();
		}
		timestep_count++;
		/*
		context->CopyResource(m_psSimTexture[old_index].get(),
							  m_psSimTexture[scratch_index].get());
		context->CopyResource(m_psSimTexture[F_G_star_old_index].get(),
							  m_psSimTexture[F_G_star_scratch_index].get());
		*/
		//swapping indecies for the next time step.
		old_index = old_index + old_old_index + scratch_index;
		old_old_index     = old_index - (old_old_index + scratch_index);
		scratch_index = old_index - (old_old_index + scratch_index);
		old_index         = old_index - (old_old_index + scratch_index);

		F_G_star_old_index = F_G_star_old_index + F_G_star_old_old_index + F_G_star_scratch_index;
		F_G_star_old_old_index     = F_G_star_old_index - (F_G_star_old_old_index + F_G_star_scratch_index);
		F_G_star_scratch_index = F_G_star_old_index - (F_G_star_old_old_index + F_G_star_scratch_index);
		F_G_star_old_index         = F_G_star_old_index - (F_G_star_old_old_index + F_G_star_scratch_index);

		// Keep track of virtual time.	
		total_time += current_timestep;
	}
	// Swap buffers.
	sim_idx = 1 - sim_idx;
}
void ShallowWaterEngine::afterTimestep()
{
	if (initSetting.saveBathymetry){
		dumpBathymetryToFile();
		initSetting.saveBathymetry = false;
	}

	if (initSetting.saveInundation){
		dumpInundationToFile(context, m_psFullSizeStagingTexture.get(), m_psAuxiliary1Texture.get());
		initSetting.saveInundation = false;
	}
	
}
void ShallowWaterEngine::resetTimestep(float realworld_dt, float elapsed_time)
{
    static float old_virtual_time = 1;
	static float old_real_time = 1;
	const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    realworld_timestep = realworld_dt;  


	static float shift_terrain = GetSetting("Tide/Surge/SLR");
	if (shift_terrain != GetSetting("Tide/Surge/SLR")){

		shiftTerrainSlider(- GetSetting("Tide/Surge/SLR") + shift_terrain);

		context->CopyResource(m_psAuxiliary1Texture.get(), st_psTriDigTexture_D1x[0].get());
		context->CopyResource(m_psAuxiliary2Texture.get(), st_psTriDigTexture_D1x[0].get());  

		shift_terrain  =  GetSetting("Tide/Surge/SLR");
	}
    // Run the GetStats pass
    // note: this uses the xflux, yflux textures as scratch space.

    ID3D11Buffer * cst_buf = m_psSimConstantBuffer.get();

    ID3D11ShaderResourceView * bottom_tex = m_psBottomTextureView.get();
    ID3D11ShaderResourceView * old_state_tex = m_psSimTextureView[sim_idx].get();
    ID3D11RenderTargetView * render_targets[] = { m_psSimRenderTargetView[4].get(),     // xflux texture
                                                  m_psSimRenderTargetView[5].get(),     // yflux texture
                                                  m_psGetStatsRenderTargetView.get() }; // dedicated R32_FLOAT texture
    context->ClearState();

    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->IASetInputLayout(m_psSimInputLayout.get());

    context->VSSetShader(m_psSimVertexShader.get(), 0, 0);

    D3D11_VIEWPORT vp;
    memset(&vp, 0, sizeof(vp));
    vp.Width = float(nx + 4);
    vp.Height = float(ny + 4);
    vp.MinDepth = 0;
    vp.MaxDepth = 1;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    context->RSSetViewports(1, &vp);

    context->VSSetConstantBuffers(0, 1, &cst_buf);
    context->PSSetConstantBuffers(0, 1, &cst_buf);

    ID3D11Buffer *vert_buf;
    const UINT stride = 8;
    const UINT offset = 0;

    vert_buf = m_psGetStatsVertexBuffer.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
    
    context->OMSetRenderTargets(3, render_targets, 0);

    context->PSSetShader(m_psGetStatsPixelShader.get(), 0, 0);
    context->PSSetShaderResources(0, 1, &old_state_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);
    
    context->Draw(6, 0);

    // immediately read back the results.
    // (TODO: this will cause a pipeline stall... could run in a background thread to avoid this?)

    D3D11_BOX src_box;
    src_box.left = 0;
    src_box.right = nx/4;
    src_box.top = 0;
    src_box.bottom = ny/4;
    src_box.front = 0;
    src_box.back = 1;

    // copy first target into left part of the staging texture,
    // second target into right part.
    // third target goes into a separate staging texture.
    
    context->CopySubresourceRegion(m_psGetStatsStagingTexture4.get(),
                                   0,
                                   0,
                                   0,
                                   0,
                                   m_psSimTexture[4].get(),
                                   0,
                                   &src_box);

    context->CopySubresourceRegion(m_psGetStatsStagingTexture4.get(),
                                   0,
                                   nx/4,
                                   0,
                                   0,
                                   m_psSimTexture[5].get(),
                                   0,
                                   &src_box);

    context->CopySubresourceRegion(m_psGetStatsStagingTexture1.get(),
                                   0,
                                   0,
                                   0,
                                   0,
                                   m_psGetStatsTexture.get(),
                                   0,
                                   &src_box);
                          
    float mass = 0;
    float x_mtm = 0, y_mtm = 0;
    float ke = 0, pe = 0;
    float max_speed = 0, max_depth = 0;
    float cfl = 0;
    float max_froude = 0.0f;

    {
        MapTexture m(*context, *m_psGetStatsStagingTexture4);

        for (int j = 0; j < ny/4; ++j) {
            const char *row_ptr = reinterpret_cast<const char *>(m.msr.pData) + j * m.msr.RowPitch;
            for (int i = 0; i < nx/4; ++i) {
                const float *col_ptr_1 = reinterpret_cast<const float*>(row_ptr) + i * 4;
                const float *col_ptr_2 = reinterpret_cast<const float*>(row_ptr) + (i+nx/4) * 4;
                
                mass += col_ptr_1[0];   // sum(h)
                pe += col_ptr_1[1];     // sum(B*h + 0.5 * h^2)
                x_mtm += col_ptr_1[2];  // sum(hu)
                y_mtm += col_ptr_1[3];  // sum(hv)
                ke += col_ptr_2[0];     // sum(h*(u2+v2))
                max_speed = std::max(max_speed, col_ptr_2[1]);   // max(u2+v2)
                max_depth = std::max(max_depth, col_ptr_2[2]);   // max(h)
				
				cfl = std::max(cfl, col_ptr_2[3]);  // max((|u|+c)/dx , (|v|+c)/dy)
	
            }
        }
    }

    {
        MapTexture m(*context, *m_psGetStatsStagingTexture1);

        for (int j = 0; j < ny/4; ++j) {
            const char *row_ptr = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
            for (int i = 0; i < nx/4; ++i) {
                const float *p = reinterpret_cast<const float*>(row_ptr) + i;
                max_froude = std::max(max_froude, *p);   // max((u2+v2)/h)
            }
        }
    }
    
    const float DENSITY = 1000;   // kg m^-3
    const float AREA = GetSetting("valley_width") / float(nx-1) 
        * GetSetting("valley_length") / float(ny-1);   // m^2 (area of one cell)

    const float g = GetSetting("gravity");
    mass *= DENSITY * AREA;
    x_mtm *= DENSITY * AREA;
    y_mtm *= DENSITY * AREA;
    ke *= 0.5f * DENSITY * AREA;
    pe *= g * DENSITY * AREA;
    max_speed = std::sqrt(max_speed);
    max_froude /= g;
    max_froude = std::sqrt(max_froude);
        
    // The CFL number is cfl * dt, and this must be less than safety_factor, so dt < safety_factor/cfl
    const float safety_factor = GetSetting("nominal_cfl");
    

	float max_celerity = sqrt((initSetting.stillWaterElevation - initSetting.min_bathy)*g);
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
	const float dx = W / (nx-1);
	const float dy = L / (ny-1);

	//ST_: The following line is changed to avoid variational dt.
	current_timestep = safety_factor * (std::min(dx,dy)/max_celerity); //std::min(dt * GetSetting("time_acceleration"), safety_factor / cfl);
	

    // update the displays
    SetSetting("mass", mass);
    SetSetting("x_momentum", x_mtm);
    SetSetting("y_momentum", y_mtm);
    SetSetting("kinetic_energy", ke);
    SetSetting("potential_energy", pe);
    SetSetting("total_energy", ke + pe);
    SetSetting("max_speed", max_speed);
    SetSetting("max_depth", max_depth);
    SetSetting("max_froude_number", max_froude);
    SetSetting("timestep", current_timestep);
    SetSetting("cfl_number", cfl * current_timestep);
    
	SetSetting("time_ratio", (total_time - old_virtual_time) / (elapsed_time - old_real_time));
	SetSetting("virtual time", total_time);
	SetSetting("elapsed time", elapsed_time);
	old_virtual_time = total_time;
	old_real_time = elapsed_time;
    
    fillConstantBuffers();  // communicate new dt to the simulation.

	if(initSetting.is_there_new_solitary_wave){
		addSolitaryWave();
		initSetting.is_there_new_solitary_wave = false;
	}


	static bool colormap_initialized = false;
	if (!colormap_initialized ){
		float surfaceColormapMin = initSetting.graphics.surfaceShading.autoColormap ? - max_depth : initSetting.graphics.surfaceShading.colormapMin;
		float surfaceColormapMax = initSetting.graphics.surfaceShading.autoColormap ? + max_depth : initSetting.graphics.surfaceShading.colormapMax;

		SetSettingMax("Colormap Max", surfaceColormapMax);
		SetSettingMax("Colormap Min", surfaceColormapMax);
		SetSettingMin("Colormap Max", surfaceColormapMin);
		SetSettingMin("Colormap Min", surfaceColormapMin);
		
		SetSetting("Colormap Max", surfaceColormapMax/2.0f);
		SetSetting("Colormap Min", surfaceColormapMin/2.0f);

		SetSettingMax("Colormap Max ",  initSetting.graphics.terrainTexture.autoColormap ?
										1.5 * initSetting.max_positive_bathy :
										initSetting.graphics.terrainTexture.colormapMax);
		SetSettingMin("Colormap Max ", 0);
		
		SetSettingMax("Colormap Min ", 0);
		SetSettingMin("Colormap Min ",  initSetting.graphics.terrainTexture.autoColormap ?
										1.5 * initSetting.min_negative_bathy :
										initSetting.graphics.terrainTexture.colormapMin);
		
		SetSetting("Colormap Max ", initSetting.max_positive_bathy);
		SetSetting("Colormap Min ", initSetting.min_negative_bathy);

		SetSettingMax("Flow Depth",  initSetting.graphics.surfaceShading.autoInundationDepth ?
										10.0f * sqrt(sqrt(initSetting.epsilon)):
										initSetting.graphics.surfaceShading.maxInundation);
		SetSettingMin("Flow Depth", 0.5f * sqrt(sqrt(initSetting.epsilon)));


		SetSetting("Flow Depth", initSetting.graphics.surfaceShading.autoInundationDepth ? 
										sqrt(sqrt(initSetting.epsilon)) :
										initSetting.graphics.surfaceShading.drylandDepthOfInundation);

		SetSettingMax("Tide/Surge/SLR", initSetting.tideSurgeSLR.maxValue);
		SetSettingMin("Tide/Surge/SLR", initSetting.tideSurgeSLR.minValue);
		SetSetting("Tide/Surge/SLR", initSetting.tideSurgeSLR.setValue);

		if(max_depth != 0) colormap_initialized = true;
	}

}


// BoundaryConstBuffer must be bound to b0 before using this function.
 void ShallowWaterEngine::westIrregularBoundary(){

	ID3D11ShaderResourceView * bottom_tex = m_psBottomTextureView.get();
	context->PSSetShaderResources(1, 1, &bottom_tex);

	// xflux is used as scratch texture elsewhere. Here we reset it to (0,0,0,0). Note that D1x[0] is never used in this code, so it is all (0,0,0,0).
	context->CopyResource(m_psSimTexture[4].get(),
						  st_psTriDigTexture_D1x[0].get()); 
	ID3D11RenderTargetView * xflux_target = m_psSimRenderTargetView[4].get();

	const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
	const int countOfMatricesXdirection = st_getCountOfMatrices(nx);

	ID3D11Buffer *irr_cst_buf = m_psIrregularWavesDataConstantBuffer.get();
	context->PSSetConstantBuffers(1, 1, &irr_cst_buf);

	ID3D11Buffer *vert_buf;
    const UINT stride = 8;
    const UINT offset = 0;


	for(int col = 0; col < initSetting.westBoundary.width; ++col){

		IrregularWavesColumnRowConstBuffer irr_ColRow_cBuffer;
		irr_ColRow_cBuffer.columnNumber = col;
		irr_ColRow_cBuffer.rowNumber = 0;

		context->UpdateSubresource(m_psIrregularWavesColumnRowConstBuffer.get(), 0, 0, &irr_ColRow_cBuffer, 0, 0);

		ID3D11Buffer *irr_col_row_cst_buf = m_psIrregularWavesColumnRowConstBuffer.get();
		context->PSSetConstantBuffers(2, 1, &irr_col_row_cst_buf);

		context->OMSetRenderTargets(1, &xflux_target, 0);
		vert_buf = st_psCyclicReductionX_VertexBuffer[0].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		
		context->PSSetShader(m_psBoundaryPixelShader[3].get(), 0, 0);
		context->Draw(6, 0);

		context->PSSetShader(m_psColumnSumReduceShader.get(), 0, 0);
		for(int i = 1; i < countOfMatricesXdirection; ++i){

			// Target must be set before resource! Because the target in this step, is the resource in the next step
			// and a resource cannot be a target. So we need to set the new target to unbind the previous target, and bind it as resource!
			ID3D11RenderTargetView * next_target = st_psTriDigRenderTargetView_D1x[i].get();
			ID3D11RenderTargetView * _tgt[] = {next_target,0}; //CHECKME
			context->OMSetRenderTargets(2, &_tgt[0], 0);

			ID3D11ShaderResourceView * last_tex = (i == 1 ? m_psSimTextureView[4].get() : st_psTriDigTextureView_D1x[i-1].get());
			context->PSSetShaderResources(0, 1, &last_tex); // not to conflict with bottom texture

			vert_buf = st_psCyclicReductionX_VertexBuffer[i].get();
			context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
			
			context->Draw(6, 0);
		}

		D3D11_BOX src_box;
		src_box.front = 0;
		src_box.back = 1;

		src_box.left = 2; // skipping two ghost cells
		src_box.right = 2 + 1;
		src_box.top = 0;
		src_box.bottom = ny + 4;
	
		context->CopySubresourceRegion(m_psSimTexture[1 - sim_idx].get(), 0, col, 0, 0, st_psTriDigTexture_D1x[countOfMatricesXdirection - 1].get(), 0, &src_box); 
	}
} 

float lerp_coef(float x0, float x, float x1)
{
	return (x - x0)/(x1 - x0);
}

void timeSeriesHelper(float & temp_eta,float & temp_hu,float & temp_hv,float total_time, std::vector<std::vector<float>> & data_v, int & iter){

	const int TIME = 0, ETA = 1, HU = 2, HV = 3;
	temp_eta = 0; temp_hu = 0; temp_hv = 0;

	int num_of_time_steps = data_v.size();

	while (iter < num_of_time_steps && data_v[iter][TIME] < total_time){
		++iter;
	}
	if (iter >= num_of_time_steps){
		return;
	}
	float m = lerp_coef(data_v[iter - 1][TIME], total_time, data_v[iter][TIME]);
	temp_eta = data_v[iter - 1][ETA] + (data_v[iter][ETA] - data_v[iter - 1][ETA]) * m;
	temp_hu  = data_v[iter - 1][HU]  + (data_v[iter][HU]  - data_v[iter - 1][HU] ) * m;
	temp_hv  = data_v[iter - 1][HV]  + (data_v[iter][HV]  - data_v[iter - 1][HV] ) * m;
}

void ShallowWaterEngine::westUniformTimeSeries(float & temp_eta,float & temp_hu,float & temp_hv,float total_time){

	std::vector<std::vector<float>> & data_v = initSetting.westBoundary.uniformTimeSeries.data;
	int & iter = initSetting.westBoundary.uniformTimeSeries.iter;

	timeSeriesHelper(temp_eta, temp_hu, temp_hv, total_time, data_v, iter);
}


void ShallowWaterEngine::eastUniformTimeSeries(float & temp_eta,float & temp_hu,float & temp_hv,float total_time){

	std::vector<std::vector<float>> & data_v = initSetting.eastBoundary.uniformTimeSeries.data;
	int & iter = initSetting.eastBoundary.uniformTimeSeries.iter;

	timeSeriesHelper(temp_eta, temp_hu, temp_hv, total_time, data_v, iter);
}


void ShallowWaterEngine::southUniformTimeSeries(float & temp_eta,float & temp_hu,float & temp_hv,float total_time){

	std::vector<std::vector<float>> & data_v = initSetting.southBoundary.uniformTimeSeries.data;
	int & iter = initSetting.southBoundary.uniformTimeSeries.iter;

	timeSeriesHelper(temp_eta, temp_hu, temp_hv, total_time, data_v, iter);
}


void ShallowWaterEngine::northUniformTimeSeries(float & temp_eta,float & temp_hu,float & temp_hv,float total_time){

	std::vector<std::vector<float>> & data_v = initSetting.northBoundary.uniformTimeSeries.data;
	int & iter = initSetting.northBoundary.uniformTimeSeries.iter;

	timeSeriesHelper(temp_eta, temp_hu, temp_hv, total_time, data_v, iter);
}


 // BoundaryConstBuffer must be bound to b0 before using this function.
void ShallowWaterEngine::eastIrregularBoundary(){

	ID3D11ShaderResourceView * bottom_tex = m_psBottomTextureView.get();
	context->PSSetShaderResources(1, 1, &bottom_tex);

	// xflux is used as scratch texture elsewhere. Here we reset it to (0,0,0,0). Note that D1x[0] is never used in this code, so it is all (0,0,0,0).
	context->CopyResource(m_psSimTexture[4].get(),
						  st_psTriDigTexture_D1x[0].get()); 
	ID3D11RenderTargetView * xflux_target = m_psSimRenderTargetView[4].get();

	const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
	const int countOfMatricesXdirection = st_getCountOfMatrices(nx);

	ID3D11Buffer *irr_cst_buf = m_psIrregularWavesDataConstantBuffer.get();
	context->PSSetConstantBuffers(1, 1, &irr_cst_buf);

	ID3D11Buffer *vert_buf;
    const UINT stride = 8;
    const UINT offset = 0;

	for(int col = nx + 4 - initSetting.eastBoundary.width; col < nx + 4; ++col){

		IrregularWavesColumnRowConstBuffer irr_ColRow_cBuffer;
		irr_ColRow_cBuffer.columnNumber = col; 
		irr_ColRow_cBuffer.rowNumber = 0;

		context->UpdateSubresource(m_psIrregularWavesColumnRowConstBuffer.get(), 0, 0, &irr_ColRow_cBuffer, 0, 0);

		ID3D11Buffer *irr_col_row_cst_buf = m_psIrregularWavesColumnRowConstBuffer.get();
		context->PSSetConstantBuffers(2, 1, &irr_col_row_cst_buf);

		context->OMSetRenderTargets(1, &xflux_target, 0);
		vert_buf = st_psCyclicReductionX_VertexBuffer[0].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		
		context->PSSetShader(m_psBoundaryPixelShader[1].get(), 0, 0); // 1 is east
		context->Draw(6, 0);

		context->PSSetShader(m_psColumnSumReduceShader.get(), 0, 0);
		for(int i = 1; i < countOfMatricesXdirection; ++i){

			// Target must be set before resource! Because the target in this step, is the resource in the next step
			// and a resource cannot be a target. So we need to set the new target to unbind the previous target, and bind it as resource!
			ID3D11RenderTargetView * next_target = st_psTriDigRenderTargetView_D1x[i].get();
			ID3D11RenderTargetView * _tgt[] = {next_target,0}; //CHECKME
			context->OMSetRenderTargets(2, &_tgt[0], 0);

			ID3D11ShaderResourceView * last_tex = (i == 1 ? m_psSimTextureView[4].get() : st_psTriDigTextureView_D1x[i-1].get());
			context->PSSetShaderResources(0, 1, &last_tex); // not to conflict with bottom texture

			vert_buf = st_psCyclicReductionX_VertexBuffer[i].get();
			context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
			
			context->Draw(6, 0);
		}

		D3D11_BOX src_box;
		src_box.front = 0;
		src_box.back = 1;

		src_box.left = 2; // skipping two ghost cells
		src_box.right = 2 + 1;
		src_box.top = 0;
		src_box.bottom = ny + 4;
	
		context->CopySubresourceRegion(m_psSimTexture[1 - sim_idx].get(), 0, col, 0, 0, st_psTriDigTexture_D1x[countOfMatricesXdirection - 1].get(), 0, &src_box); 
	}

	
} 

// BoundaryConstBuffer must be bound to b0 before using this function.
void ShallowWaterEngine::southIrregularBoundary(){

	ID3D11ShaderResourceView * bottom_tex = m_psBottomTextureView.get();
	context->PSSetShaderResources(1, 1, &bottom_tex);

	// xflux is used as scratch texture elsewhere. Here we reset it to (0,0,0,0). Note that D1x[0] is never used in this code, so it is all (0,0,0,0).
	context->CopyResource(m_psSimTexture[4].get(),
						  st_psTriDigTexture_D1y[0].get()); 
	ID3D11RenderTargetView * xflux_target = m_psSimRenderTargetView[4].get();

	const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
	const int countOfMatricesYdirection = st_getCountOfMatrices(ny);

	ID3D11Buffer *irr_cst_buf = m_psIrregularWavesDataConstantBuffer.get();
	context->PSSetConstantBuffers(1, 1, &irr_cst_buf);

	ID3D11Buffer *vert_buf;
    const UINT stride = 8;
    const UINT offset = 0;

	for(int row = 0; row < initSetting.southBoundary.width; ++row){

		IrregularWavesColumnRowConstBuffer irr_ColRow_cBuffer;
		irr_ColRow_cBuffer.columnNumber = 0;
		irr_ColRow_cBuffer.rowNumber = row;

		context->UpdateSubresource(m_psIrregularWavesColumnRowConstBuffer.get(), 0, 0, &irr_ColRow_cBuffer, 0, 0);

		ID3D11Buffer *irr_col_row_cst_buf = m_psIrregularWavesColumnRowConstBuffer.get();
		context->PSSetConstantBuffers(2, 1, &irr_col_row_cst_buf);

		context->OMSetRenderTargets(1, &xflux_target, 0);
		vert_buf = st_psCyclicReductionY_VertexBuffer[0].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		
		context->PSSetShader(m_psBoundaryPixelShader[2].get(), 0, 0);
		context->Draw(6, 0);

		context->PSSetShader(m_psRowSumReduceShader.get(), 0, 0);
		for(int i = 1; i < countOfMatricesYdirection; ++i){

			// Target must be set before resource! Because the target in this step, is the resource in the next step
			// and a resource cannot be a target. So we need to set the new target to unbind the previous target, and bind it as resource!
			ID3D11RenderTargetView * next_target = st_psTriDigRenderTargetView_D1y[i].get();
			ID3D11RenderTargetView * _tgt[] = {next_target,0}; //CHECKME
			context->OMSetRenderTargets(2, &_tgt[0], 0);

			ID3D11ShaderResourceView * last_tex = (i == 1 ? m_psSimTextureView[4].get() : st_psTriDigTextureView_D1y[i-1].get());
			context->PSSetShaderResources(0, 1, &last_tex); // not to conflict with bottom texture	

			vert_buf = st_psCyclicReductionY_VertexBuffer[i].get();
			context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
			
			context->Draw(6, 0);
		}

		D3D11_BOX src_box;
		src_box.front = 0;
		src_box.back = 1;

		src_box.top = 2; // skipping two ghost cells
		src_box.bottom = 2 + 1;
		src_box.left = 0;
		src_box.right = nx + 4;
	
		context->CopySubresourceRegion(m_psSimTexture[1 - sim_idx].get(), 0, 0, row, 0, st_psTriDigTexture_D1y[countOfMatricesYdirection - 1].get(), 0, &src_box); 
	}
} 



// BoundaryConstBuffer must be bound to b0 before using this function.
 void ShallowWaterEngine::northIrregularBoundary(){

	ID3D11ShaderResourceView * bottom_tex = m_psBottomTextureView.get();
	context->PSSetShaderResources(1, 1, &bottom_tex);

	// xflux is used as scratch texture elsewhere. Here we reset it to (0,0,0,0). Note that D1x[0] is never used in this code, so it is all (0,0,0,0).
	context->CopyResource(m_psSimTexture[4].get(),
						  st_psTriDigTexture_D1y[0].get()); 
	ID3D11RenderTargetView * xflux_target = m_psSimRenderTargetView[4].get();

	const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
	const int countOfMatricesYdirection = st_getCountOfMatrices(ny);

	ID3D11Buffer *irr_cst_buf = m_psIrregularWavesDataConstantBuffer.get();
	context->PSSetConstantBuffers(1, 1, &irr_cst_buf);

	ID3D11Buffer *vert_buf;
    const UINT stride = 8;
    const UINT offset = 0;

	for(int row = ny + 4 - initSetting.northBoundary.width; row < ny + 4; ++row){

		IrregularWavesColumnRowConstBuffer irr_ColRow_cBuffer;
		irr_ColRow_cBuffer.columnNumber = 0;
		irr_ColRow_cBuffer.rowNumber = row;

		context->UpdateSubresource(m_psIrregularWavesColumnRowConstBuffer.get(), 0, 0, &irr_ColRow_cBuffer, 0, 0);

		ID3D11Buffer *irr_col_row_cst_buf = m_psIrregularWavesColumnRowConstBuffer.get();
		context->PSSetConstantBuffers(2, 1, &irr_col_row_cst_buf);

		context->OMSetRenderTargets(1, &xflux_target, 0);
		vert_buf = st_psCyclicReductionY_VertexBuffer[0].get();
		context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
		
		context->PSSetShader(m_psBoundaryPixelShader[0].get(), 0, 0);
		context->Draw(6, 0);

		context->PSSetShader(m_psRowSumReduceShader.get(), 0, 0);
		for(int i = 1; i < countOfMatricesYdirection; ++i){

			// Target must be set before resource! Because the target in this step, is the resource in the next step
			// and a resource cannot be a target. So we need to set the new target to unbind the previous target, and bind it as resource!
			ID3D11RenderTargetView * next_target = st_psTriDigRenderTargetView_D1y[i].get();
			ID3D11RenderTargetView * _tgt[] = {next_target,0}; //CHECKME
			context->OMSetRenderTargets(2, &_tgt[0], 0);

			ID3D11ShaderResourceView * last_tex = (i == 1 ? m_psSimTextureView[4].get() : st_psTriDigTextureView_D1y[i-1].get());
			context->PSSetShaderResources(0, 1, &last_tex); // not to conflict with bottom texture

			vert_buf = st_psCyclicReductionY_VertexBuffer[i].get();
			context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
			
			context->Draw(6, 0);
		}

		D3D11_BOX src_box;
		src_box.front = 0;
		src_box.back = 1;

		src_box.top = 2; // skipping two ghost cells
		src_box.bottom = 2 + 1;
		src_box.left = 0;
		src_box.right = nx + 4;
	
		context->CopySubresourceRegion(m_psSimTexture[1 - sim_idx].get(), 0, 0, row, 0, st_psTriDigTexture_D1y[countOfMatricesYdirection - 1].get(), 0, &src_box); 
	}
} 


void ShallowWaterEngine::addSolitaryWave(){

    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
	
	SolitaryWaveConstBuffer sw_cBuffer;
	sw_cBuffer.dx = W / (nx-1);
	sw_cBuffer.dy = L / (ny-1);
	sw_cBuffer.g = GetSetting("gravity");

	Soliton tempSoliton = initSetting.solitons[initSetting.MAX_NUM_SOLITARY - 1];
	sw_cBuffer.param_H = tempSoliton.param_h;
	sw_cBuffer.theta = tempSoliton.theta * PI / 180.0f;
	sw_cBuffer.x_zero = tempSoliton.xc;
	sw_cBuffer.y_zero = tempSoliton.yc;

	context->UpdateSubresource(m_psSolitaryWaveConstantBuffer.get(), 0, 0, &sw_cBuffer, 0, 0);

	ID3D11ShaderResourceView * bottom_tex = m_psBottomTextureView.get();
	ID3D11ShaderResourceView * old_state_tex = m_psSimTextureView[sim_idx].get();
    ID3D11ShaderResourceView * xflux_tex = m_psSimTextureView[4].get();

    ID3D11RenderTargetView * old_state_or_h_target = m_psSimRenderTargetView[sim_idx].get();

	context->CopyResource(m_psSimTexture[4].get(), m_psSimTexture[sim_idx].get());
    
    context->ClearState();

    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->IASetInputLayout(m_psSimInputLayout.get());

    context->VSSetShader(m_psSimVertexShader.get(), 0, 0);
    
    D3D11_VIEWPORT vp;
    memset(&vp, 0, sizeof(vp));
    vp.Width = float(nx + 4);
    vp.Height = float(ny + 4);
    vp.MinDepth = 0;
    vp.MaxDepth = 1;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    context->RSSetViewports(1, &vp);
	

	ID3D11Buffer *sw_cst_buf = m_psSolitaryWaveConstantBuffer.get();
	ID3D11Buffer * cst_buf = m_psSimConstantBuffer.get();
    context->VSSetConstantBuffers(0, 1, &cst_buf);
	context->PSSetConstantBuffers(0, 1, &sw_cst_buf);

    ID3D11Buffer *vert_buf;
    const UINT stride = 8;
    const UINT offset = 0;
    vert_buf = m_psSimVertexBuffer11.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
    
    ID3D11RenderTargetView * rtv = old_state_or_h_target;
    context->OMSetRenderTargets(1, &rtv, 0);
    
    context->PSSetShader(m_psAddSolitaryWaveShader.get(), 0, 0);
    context->PSSetShaderResources(0, 1, &xflux_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);
    
    context->Draw(6, 0);    
}

void ShallowWaterEngine::render(ID3D11RenderTargetView *render_target_view)
{
    // Setup rendering state
    context->ClearState();

    ID3D11ShaderResourceView * heightfield_tex = m_psTerrainTextureView.get();
    ID3D11ShaderResourceView * water_tex = m_psSimTextureView[sim_idx].get();   // {w, hu, hv, unused}
    ID3D11ShaderResourceView * normal_tex = m_psSimTextureView[6].get();     // {nX, nY, nZ, unused}
    ID3D11ShaderResourceView * skybox_tex = m_psSkyboxView.get();
    ID3D11ShaderResourceView * grass_tex = m_psGrassTextureView.get();
	ID3D11ShaderResourceView * auxiliary1_tex = m_psAuxiliary1TextureView.get();
	ID3D11ShaderResourceView * grid_tex = m_psGridTextureView.get();
	ID3D11ShaderResourceView * colormap_tex = m_psColormapTextureView.get();
    ID3D11SamplerState * linear_sampler = m_psLinearSamplerState.get();
    
    
    // input assembler
    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->IASetInputLayout(m_psInputLayout.get());

    ID3D11Buffer *vert_buf = m_psMeshVertexBuffer.get();
    UINT stride = sizeof(MeshVertex);
    const UINT offset = 0;
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
    context->IASetIndexBuffer(m_psMeshIndexBuffer.get(), DXGI_FORMAT_R32_UINT, 0);

    // vertex shader
    context->VSSetShader(m_psWaterVertexShader.get(), 0, 0);

    context->VSSetShaderResources(0, 1, &heightfield_tex);
    context->VSSetShaderResources(1, 1, &water_tex);
    context->VSSetShaderResources(2, 1, &normal_tex);
	context->VSSetShaderResources(4, 1, &auxiliary1_tex);
    
    ID3D11Buffer * cst_buf = m_psConstantBuffer.get();
    context->VSSetConstantBuffers(0, 1, &cst_buf);

    // viewport
    D3D11_VIEWPORT vp;
    memset(&vp, 0, sizeof(vp));
    vp.Width = float(vp_width);
    vp.Height = float(vp_height);
    vp.MinDepth = 0;
    vp.MaxDepth = 1;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    context->RSSetViewports(1, &vp);

    // pixel shader
    context->PSSetShader(m_psWaterPixelShader.get(), 0, 0);
    context->PSSetShaderResources(0, 1, &grass_tex);
    context->PSSetShaderResources(1, 1, &skybox_tex);
	context->PSSetShaderResources(2, 1, &colormap_tex);
	context->PSSetShaderResources(3, 1, &grid_tex);
	

    context->PSSetConstantBuffers(0, 1, &cst_buf);
    context->PSSetSamplers(0, 1, &linear_sampler);
    
    // setup depth buffer if we need to
    int rw, rh;
    int dw = -999, dh;
    GetRenderTargetSize(render_target_view, rw, rh);
    if (m_psDepthStencil.get()) GetTextureSize(m_psDepthStencil.get(), dw, dh);
    if (rw != dw || rh != dh) createDepthStencil(rw, rh);

    // clear the depth buffer
    context->ClearDepthStencilView(m_psDepthStencilView.get(), D3D11_CLEAR_DEPTH, 1.0f, 0);

    // set render target
    context->OMSetRenderTargets(1, &render_target_view, m_psDepthStencilView.get());

    // draw the mesh
    const int mesh_width = GetIntSetting("mesh_size_x");
    const int mesh_height = GetIntSetting("mesh_size_y");
    context->DrawIndexed(6 * (mesh_width - 1) * (mesh_height - 1), 0, 0);

    // now setup for the terrain mesh
    // (it uses the same mesh, but different shaders.)
    context->VSSetShader(m_psTerrainVertexShader.get(), 0, 0);
    context->PSSetShader(m_psTerrainPixelShader.get(), 0, 0);

    // draw the mesh
    context->DrawIndexed(6 * (mesh_width - 1) * (mesh_height - 1), 0, 0);


    // Now draw the skybox
    context->IASetInputLayout(m_psSkyboxInputLayout.get());
    vert_buf = m_psSkyboxVertexBuffer.get();
    stride = 3 * sizeof(float);   // 3 coords per vertex
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
    context->VSSetShader(m_psSkyboxVertexShader.get(), 0, 0);
    context->PSSetShader(m_psSkyboxPixelShader.get(), 0, 0);
    context->Draw(36, 0);
}


void ShallowWaterEngine::createShadersAndInputLayout()
{
    Coercri::ComPtrWrapper<ID3DBlob> bytecode;
	
    
	std::string shallow_water_fx = initSetting.exePath + "/shaders/graphics.fx";
	CreateVertexShader(device, shallow_water_fx.c_str(), "TerrainVertexShader", bytecode, m_psTerrainVertexShader);
    CreatePixelShader(device, shallow_water_fx.c_str(), "TerrainPixelShader", m_psTerrainPixelShader);

    CreateVertexShader(device, shallow_water_fx.c_str(), "WaterVertexShader", bytecode, m_psWaterVertexShader);
    CreatePixelShader(device, shallow_water_fx.c_str(), "WaterPixelShader", m_psWaterPixelShader);

    D3D11_INPUT_ELEMENT_DESC layout[] = {
        { "INT_POSITION", 0, DXGI_FORMAT_R32G32_SINT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
    };

    ID3D11InputLayout *input_layout = 0;
    HRESULT hr = device->CreateInputLayout(&layout[0],
                                           1,
                                           bytecode->GetBufferPointer(),
                                           bytecode->GetBufferSize(),
                                           &input_layout);
    if (FAILED(hr)) {
        throw Coercri::DXError("CreateInputLayout failed", hr);
    }
    m_psInputLayout.reset(input_layout);

	std::string compute_hlsl = initSetting.exePath + "/shaders/compute.hlsl";

	// ColumnSumReduce is used to superpose sinwaves to generate irregular waves.
	CreatePixelShader(device, compute_hlsl.c_str(), "ColumnSumReduce", m_psColumnSumReduceShader);
	CreatePixelShader(device, compute_hlsl.c_str(), "RowSumReduce", m_psRowSumReduceShader);

    CreateVertexShader(device, compute_hlsl.c_str(), "SimVertexShader", bytecode, m_psSimVertexShader);
    CreatePixelShader(device, compute_hlsl.c_str(), "Pass1", m_psSimPixelShader[0]);
    CreatePixelShader(device, compute_hlsl.c_str(), "Pass2", m_psSimPixelShader[1]);
    CreatePixelShader(device, compute_hlsl.c_str(), "Pass3Predictor", m_psSimPixelShader[2]);
		
	CreatePixelShader(device, compute_hlsl.c_str(), "CopyFromXxAndXy", st_psCopyFromXxandXyPixelShader);
	CreatePixelShader(device, compute_hlsl.c_str(), "CyclicReduceDx", st_psCyclicReduceDx_PixelShader);
	CreatePixelShader(device, compute_hlsl.c_str(), "CyclicReduceDy", st_psCyclicReduceDy_PixelShader);
	CreatePixelShader(device, compute_hlsl.c_str(), "CyclicSubstituteXx", st_psCyclicSubstituteXx_PixelShader);
	CreatePixelShader(device, compute_hlsl.c_str(), "CyclicSubstituteXy", st_psCyclicSubstituteXy_PixelShader);
	
	CreatePixelShader(device, compute_hlsl.c_str(), "AddSolitaryWave", m_psAddSolitaryWaveShader);
	
    CreatePixelShader(device, compute_hlsl.c_str(), "GetStats", m_psGetStatsPixelShader);
	
	createBoundaryShadersInit();
    
    D3D11_INPUT_ELEMENT_DESC sim_layout[] = {
        { "TEX_IDX", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
    };

    hr = device->CreateInputLayout(&sim_layout[0],
                                   1,
                                   bytecode->GetBufferPointer(),
                                   bytecode->GetBufferSize(),
                                   &input_layout);
    if (FAILED(hr)) {
        throw Coercri::DXError("CreateInputLayout failed", hr);
    }
    m_psSimInputLayout.reset(input_layout);
}

void ShallowWaterEngine::createBoundaryShadersInit(){


	std::string compute_hlsl = initSetting.exePath + "/shaders/compute.hlsl";

	if(initSetting.northBoundary.hasChanged){
		std::string northBoundary = "NorthBoundary" + initSetting.northBoundary.type;
		CreatePixelShader(device, compute_hlsl.c_str(), northBoundary.c_str(), m_psBoundaryPixelShader[0]);
		initSetting.northBoundary.hasChanged = false;
	}
	if(initSetting.eastBoundary.hasChanged){
		std::string eastBoundary = "EastBoundary" + initSetting.eastBoundary.type;
		CreatePixelShader(device, compute_hlsl.c_str(), eastBoundary.c_str(), m_psBoundaryPixelShader[1]);
		initSetting.eastBoundary.hasChanged = false;
	}
	if(initSetting.southBoundary.hasChanged){
		std::string southBoundary = "SouthBoundary" + initSetting.southBoundary.type;
		CreatePixelShader(device, compute_hlsl.c_str(), southBoundary.c_str(), m_psBoundaryPixelShader[2]);
		initSetting.southBoundary.hasChanged = false;
	}
	if(initSetting.westBoundary.hasChanged){
		std::string westBoundary = "WestBoundary" + initSetting.westBoundary.type;
		CreatePixelShader(device, compute_hlsl.c_str(), westBoundary.c_str(), m_psBoundaryPixelShader[3]);
		initSetting.westBoundary.hasChanged = false;
	}
}

void ShallowWaterEngine::createBoundaryShaders(){

	std::string compute_hlsl = initSetting.exePath + "/shaders/compute.hlsl";

	if(initSetting.northBoundary.hasChanged){
		std::string northBoundary = "NorthBoundary" + initSetting.northBoundary.type;
		CreatePixelShader(device, compute_hlsl.c_str(), northBoundary.c_str(), m_psBoundaryPixelShader[0]);
		
		if (initSetting.northBoundary.type == "IrregularWaves"){
			fillIrregularWavesDataConstantBuffer ("NORTH");  
		}
		//fillUniformTimeSeriesMainMemoryBuffer ();

		initSetting.northBoundary.hasChanged = false;
	}
	if(initSetting.eastBoundary.hasChanged){
		std::string eastBoundary = "EastBoundary" + initSetting.eastBoundary.type;
		CreatePixelShader(device, compute_hlsl.c_str(), eastBoundary.c_str(), m_psBoundaryPixelShader[1]);

		if (initSetting.eastBoundary.type == "IrregularWaves"){
			fillIrregularWavesDataConstantBuffer ("EAST");  
		}
		//fillUniformTimeSeriesMainMemoryBuffer ();

		initSetting.eastBoundary.hasChanged = false;
	}
	if(initSetting.southBoundary.hasChanged){
		std::string southBoundary = "SouthBoundary" + initSetting.southBoundary.type;
		CreatePixelShader(device, compute_hlsl.c_str(), southBoundary.c_str(), m_psBoundaryPixelShader[2]);

		if (initSetting.southBoundary.type == "IrregularWaves"){
			fillIrregularWavesDataConstantBuffer ("SOUTH");  
		}
		//fillUniformTimeSeriesMainMemoryBuffer ();

		initSetting.southBoundary.hasChanged = false;
	}
	if(initSetting.westBoundary.hasChanged){
		std::string westBoundary = "WestBoundary" + initSetting.westBoundary.type;
		CreatePixelShader(device, compute_hlsl.c_str(), westBoundary.c_str(), m_psBoundaryPixelShader[3]);
		
		if (initSetting.westBoundary.type == "IrregularWaves"){
			fillIrregularWavesDataConstantBuffer ("WEST");  
		}
		//fillUniformTimeSeriesMainMemoryBuffer ();		
		
		initSetting.westBoundary.hasChanged = false;
	}
}

// create the heightfield vertex and index buffers.
// (needs to be called again every time the mesh size changes.)
void ShallowWaterEngine::createMeshBuffers()
{
    const int width = GetIntSetting("mesh_size_x");
    const int height = GetIntSetting("mesh_size_y");

    // fill in the array of vertex positions
    // (these are texture indices -- add two to avoid the ghost zones)

    boost::scoped_array<MeshVertex> vertices(new MeshVertex[width*height]);
    
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            MeshVertex &v = vertices[j * width + i];
            v.x = i + 2;
            v.y = j + 2;
        }
    }

    // create data for the index buffer
    boost::scoped_array<int> indices(new int[6 * (width-1) * (height-1)]);
    int *p = &indices[0];
    for (int j = 0; j < height - 1; ++j) {
        for (int i = 0; i < width - 1; ++i) {
            const int tl = (j+1) * width + i;
            const int bl = j * width + i;
            const int tr = (j+1) * width + i+1;
            const int br = j * width + i+1;

            // write triangles in clockwise order.
            *p++ = tl;
            *p++ = br;
            *p++ = bl;

            *p++ = tl;
            *p++ = tr;
            *p++ = br;
        }
    }
            
    // create the vertex buffer
    D3D11_BUFFER_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.ByteWidth = sizeof(MeshVertex) * width * height;
    bd.Usage = D3D11_USAGE_IMMUTABLE;
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;

    D3D11_SUBRESOURCE_DATA sd;
    memset(&sd, 0, sizeof(sd));
    sd.pSysMem = &vertices[0];

    ID3D11Buffer *pBuffer;
    HRESULT hr = device->CreateBuffer(&bd, &sd, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create mesh vertex buffer", hr);
    }
    m_psMeshVertexBuffer.reset(pBuffer);

    // create the index buffer
    bd.ByteWidth = sizeof(int) * 6 * (width-1) * (height-1);
    bd.BindFlags = D3D11_BIND_INDEX_BUFFER;
    sd.pSysMem = &indices[0];
    sd.SysMemPitch = 0;

    hr = device->CreateBuffer(&bd, &sd, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create mesh index buffer", hr);
    }
    m_psMeshIndexBuffer.reset(pBuffer);
}

void ShallowWaterEngine::createSimBuffers()
{
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    
    CreateSimBuffer(device, m_psSimVertexBuffer11, 1, 1, nx + 3, ny + 3);   // single ghost layer around each side
    CreateSimBuffer(device, m_psSimVertexBuffer10, 1, 1, nx + 2, ny + 2);   // west/south ghost layer only
    CreateSimBuffer(device, m_psSimVertexBuffer00, 2, 2, nx + 2, ny + 2);   // interior zones only

    CreateSimBuffer(device, m_psGetStatsVertexBuffer, 0, 0, nx/4, ny/4);

	int northBoundaryWidth = initSetting.northBoundary.width;
	int eastBoundaryWidth  = initSetting.eastBoundary.width;
	int southBoundaryWidth = initSetting.southBoundary.width;
	int westBoundaryWidth  = initSetting.westBoundary.width;

    CreateSimBuffer(device, m_psBoundaryVertexBuffer[0], 0, ny+4 - northBoundaryWidth, nx+4, ny+4);  // north border
    CreateSimBuffer(device, m_psBoundaryVertexBuffer[1], nx+4 - eastBoundaryWidth, 0, nx+4, ny+4);  // east border
    CreateSimBuffer(device, m_psBoundaryVertexBuffer[2], 0, 0, nx+4, southBoundaryWidth);   // south border
    CreateSimBuffer(device, m_psBoundaryVertexBuffer[3], 0, 0, westBoundaryWidth, ny+4);   // west border

	//ST_ creating vertex buffers for Cyclic reduction process
	const int countOfMatricesXdirection = st_getCountOfMatrices(nx);
	int arrayLength = nx;
	for(int i = 0; i < countOfMatricesXdirection; i++){
		CreateSimBuffer(device, st_psCyclicReductionX_VertexBuffer[i], 2, 2, arrayLength + 2, ny+2);   // interior zone in X direction. ghost zones are ignored.
		arrayLength/=2;
	}

	const int countOfMatricesYdirection = st_getCountOfMatrices(ny);
	arrayLength = ny;
	for(int i = 0; i < countOfMatricesYdirection; i++){
		CreateSimBuffer(device, st_psCyclicReductionY_VertexBuffer[i], 2, 2, nx+2, arrayLength + 2);   // interior zone in Y direction. ghost zones are ignored.
		arrayLength/=2;
	}

}

// Creates the terrain texture (leaving it uninitialized)
void ShallowWaterEngine::createTerrainTexture()
{
    // allow space for two ghost zones around each edge (four in total)
    CreateTexture(device,
                  GetIntSetting("mesh_size_x") + 4,
                  GetIntSetting("mesh_size_y") + 4,
                  0,  // initial data
                  DXGI_FORMAT_R32G32B32_FLOAT,
                  false,  // staging
                  m_psTerrainTexture,
                  &m_psTerrainTextureView,
                  0);

    CreateTexture(device,
                  GetIntSetting("mesh_size_x") + 4,
                  GetIntSetting("mesh_size_y") + 4,
                  0,  // initial data
                  DXGI_FORMAT_R32G32B32_FLOAT,
                  false, 
                  m_psBottomTexture,
                  &m_psBottomTextureView,
                  0);

    CreateTexture(device,
				  GetIntSetting("mesh_size_x") + 4,
				  GetIntSetting("mesh_size_y") + 4,
				  0,  // initial data
				  DXGI_FORMAT_R32G32B32A32_FLOAT,
				  false,  // staging
				  m_psAuxiliary1Texture,
				  &m_psAuxiliary1TextureView,
				  &m_psAuxiliary1RenderTargetView);

	CreateTexture(device,
				  GetIntSetting("mesh_size_x") + 4,
				  GetIntSetting("mesh_size_y") + 4,
				  0,  // initial data
				  DXGI_FORMAT_R32G32B32A32_FLOAT,
				  false,  // staging
				  m_psAuxiliary2Texture,
				  &m_psAuxiliary2TextureView,
				  &m_psAuxiliary2RenderTargetView);
}



void ShallowWaterEngine::createTridiagonalCoefTextures()
{

    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
	const int nx_plus_4=nx+4;
	const int ny_plus_4=ny+4;
	const float W = GetSetting("valley_width");
	const float L = GetSetting("valley_length");
	const float delta_x = W / (nx-1);
	const float delta_y = L / (ny-1);
	const float delta_x2 = delta_x * delta_x;
	const float delta_y2 = delta_y * delta_y;
	const float sea_level = initSetting.isBoussinesq ? initSetting.stillWaterElevation : -9999; // is model is not boussinesq set sea level to -9999 to make the coeffs ABC = (0,1,0);

	boost::scoped_array<float> matCoeffs[ST_MAX_CR_MATRIX_NUM_X];
	const int countOfMatricesXdirection = st_getCountOfMatrices(nx);
	int arrayLength=nx;

	for (int i=0; i<ST_MAX_CR_MATRIX_NUM_X;++i){  //ST_:This can be done till countOfMatricesXdirection, to be on the safe side, I intilize all the ST_MAX_CR_MATRIX_NUM_X scoped arrays.
		matCoeffs[i].reset(new float[(arrayLength+4) * (ny_plus_4) * 4]); //ST_: ny could be used instead of ny_plus_4
		arrayLength /= 2;
	}
	
/////////////
	
	//ST_: Making the first (num 0) MatCoeff of ABCx.
	int counter = 0;
	for (int i = 0; i < ny_plus_4; ++i) {
		const int ii = std::max(2, std::min(ny+1, i));

		for (int j = 0; j <= 1; ++j) { // set the Left Margin to (A,B,C)=(0,1,0)
	        //const int ii = std::max(2, std::min(nx+1, i));
			//const float B = g_bottom[(nx+4) * jj + ii].BA;
			matCoeffs[0][counter]=0; // .r=A
			counter++;
			matCoeffs[0][counter]=1; // .g=B
			counter++;
			matCoeffs[0][counter]=0; // .b=C
			counter++;
			matCoeffs[0][counter]=0; // .a=not used.
			counter++;
		}

		for (int j = 2; j <= nx_plus_4-2-1; ++j) {
			const int jj = std::max(2, std::min(nx+1, j));
			float d_here = std::max(0.0f, sea_level - g_bottom[(nx+4) * ii + jj].BA);
			float d_right_right = std::max(0.0f, sea_level - g_bottom[(nx+4) * ii + jj+2].BA);
			float d_right = std::max(0.0f, sea_level - g_bottom[(nx+4) * ii + jj+1].BA);
			float d_left = std::max(0.0f, sea_level - g_bottom[(nx+4) * ii + jj-1].BA);
			float d_left_left = std::max(0.0f, sea_level - g_bottom[(nx+4) * ii + jj-2].BA);

			const float K1x = (BCOEF+1/3.0f) * (d_here * d_here /delta_x2); //delta_x2 is delta_x*delta_x
			const float K2x = (d_here /delta_x2) * (-d_right_right + 8*d_right - 8*d_left + d_left_left); //delta_x2 is delta_x*delta_x
			matCoeffs[0][counter] = -K1x + K2x/72.0f;//.1*(j+1-2); // .r=A
			if (j==2) {
				matCoeffs[0][counter] = 0; // .r=A  set the first meaningful A to zero.
			}
			counter++;
			matCoeffs[0][counter] = 1 + 2*K1x;//j+1-2; // .g=B
			counter++;
			matCoeffs[0][counter] = -(K1x + K2x/72.0f);//.1*(j+1-2)+0.05; // .b=C
			if (j==nx_plus_4-2-1) {
				matCoeffs[0][counter] = 0; // .b=C  set the last meaningful C to zero.
			}
			counter++;
			matCoeffs[0][counter] = 0; // .a=not used.
			counter++;
		}

		for (int j = 0; j <= 1; ++j) { // set the Right Margin to (A,B,C)=(0,1,0)
	        //const int jj = std::max(2, std::min(ny+1, j));
			//const float B = g_bottom[(nx+4) * jj + ii].BA;
			matCoeffs[0][counter]=0; // .r=A
			counter++;
			matCoeffs[0][counter]=1; // .g=B
			counter++;
			matCoeffs[0][counter]=0; // .b=C
			counter++;
			matCoeffs[0][counter]=0; // .a=not used.
			counter++;
		}
	}

	arrayLength = nx;

	//ST_: Buliding the rest of the matrices using CR. 
	for (int k=1; k<countOfMatricesXdirection;++k){
		int counter=0;
		int previousCounter=8;
		int endX=arrayLength+4-2;
		for (int i = 0; i < ny_plus_4; ++i) {
			
			for (int j = 0; j <= 1; ++j) { // set the Left Margin to (A,B,C)=(0,1,0)
				matCoeffs[k][counter]=0; // .r=A
				counter++;
				matCoeffs[k][counter]=1; // .g=B
				counter++;
				matCoeffs[k][counter]=0; // .b=C
				counter++;
				matCoeffs[k][counter]=0; // .a=not used.
				counter++;
			}

			int startX=4; // skip the margin values.
			previousCounter=((startX-1)*4)+((i)*(arrayLength+4))*4;
			
			for (int j = startX; j <= endX ; j=j+2) {

				float k1= matCoeffs[k-1][previousCounter]/matCoeffs[k-1][previousCounter-4+1]; // A_i/B_(i-1)
				float k2= matCoeffs[k-1][previousCounter+1+1]/matCoeffs[k-1][previousCounter+4+1]; // C_i/B_(i+1)
	
				matCoeffs[k][counter]=-matCoeffs[k-1][previousCounter-4]*k1;  // .r=A   % aOut=-aInLeft*k1
				counter++;

				matCoeffs[k][counter]=matCoeffs[k-1][previousCounter+1]-matCoeffs[k-1][previousCounter+2-4]*k1-matCoeffs[k-1][previousCounter+4]*k2;  // .g=B  % bOut=bIn-cInLeft*k1-aInRight*k2
				counter++;

				matCoeffs[k][counter]=-matCoeffs[k-1][previousCounter+2+4]*k2;  // .b=C  % cOut=-cInRight*k2
				counter++;
				matCoeffs[k][counter]=0; // .a=not used.
				counter++;
				previousCounter=previousCounter+2*4; // two jumps to reach the next even/odd equation
			}

			for (int j = 0; j <= 1; ++j) { // set the Right Margin to (A,B,C)=(0,1,0)
				matCoeffs[k][counter]=0; // .r=A
				counter++;
				matCoeffs[k][counter]=1; // .g=B
				counter++;
				matCoeffs[k][counter]=0; // .b=C
				counter++;
				matCoeffs[k][counter]=0; // .a=not used.
				counter++;
			}
		}
		arrayLength/=2;
	}
	
	D3D11_SUBRESOURCE_DATA matCoeffSubData;
    memset(&matCoeffSubData, 0, sizeof(matCoeffSubData));
	arrayLength=nx;
	for (int k=0; k<countOfMatricesXdirection;++k){

		matCoeffSubData.pSysMem = &matCoeffs[k][0];
		matCoeffSubData.SysMemPitch = (arrayLength+4) * 4 * sizeof(float);
		CreateTexture(device,    //The Coeff matrix
					  arrayLength+4,
					  ny_plus_4,
					  &matCoeffSubData,
					  DXGI_FORMAT_R32G32B32A32_FLOAT,
					  false,
					  st_psTriDigTexture_ABCx[k],
					  &st_psTriDigTextureView_ABCx[k],
					  &st_psTriDigRenderTargetView_ABCx[k]);
		arrayLength/=2;
	}

	boost::scoped_array<float> matCoeffs_y[ST_MAX_CR_MATRIX_NUM_Y];
	const int countOfMatricesYdirection = st_getCountOfMatrices(ny);
	arrayLength=ny;

	for (int i=0; i<ST_MAX_CR_MATRIX_NUM_Y;++i){  //ST_:This can be done till countOfMatricesYdirection, to be on the safe side, I intilize all the ST_MAX_CR_MATRIX_NUM_Y scoped arrays.
		matCoeffs_y[i].reset(new float[(arrayLength+4) * (nx_plus_4) * 4]); //ST_: ny can be used instead of ny_plus_4
		arrayLength/=2;
	}
	
	//ST_: Making the first (num 0) MatCoeff of ABCy.
	counter=0;
	for (int i = 0; i < ny_plus_4; ++i) {
		const int ii = std::max(2, std::min(ny+1, i));
		for (int j = 0; j < nx_plus_4; ++j) {
			if ((i==0 || i==1) || (i==ny_plus_4-2 || i==ny_plus_4-1) ){
				matCoeffs_y[0][counter]=0;//.1*(j+1-2); // .r=A
				counter++;
				matCoeffs_y[0][counter]=1;//j+1-2; // .g=B
				counter++;
				matCoeffs_y[0][counter]=0;//.1*(j+1-2)+0.05; // .b=C
				counter++;
				matCoeffs_y[0][counter]=0; // .a=not used.
				counter++;
			}
			else{
				const int jj = std::max(2, std::min(nx+1, j));
				
				float d_here = std::max(0.0f, sea_level - g_bottom[(nx+4) * ii + jj].BA);
				float d_up_up = std::max(0.0f, sea_level - g_bottom[(nx+4) * (ii+2) + jj].BA);
				float d_up = std::max(0.0f, sea_level - g_bottom[(nx+4) * (ii+1) + jj].BA);
				float d_down = std::max(0.0f, sea_level - g_bottom[(nx+4) * (ii-1) + jj].BA);
				float d_down_down = std::max(0.0f, sea_level - g_bottom[(nx+4) * (ii-2) + jj].BA);

	

				const float K1y = (BCOEF + 1/3.0f) * (d_here * d_here /delta_y2); //delta_y2 is delta_y*delta_y
				const float K2y = (d_here / delta_y2) * (-d_up_up + 8*d_up - 8*d_down + d_down_down); //delta_x2 is delta_x*delta_x

				matCoeffs_y[0][counter] = -K1y + K2y/72.0f;//.1*(i+1-2); // .r=A
				if (i==2) {
					matCoeffs_y[0][counter] = 0; // .r=A  set the first meaningful A to zero.
				}
				counter++;
				matCoeffs_y[0][counter] = 1 + 2*K1y;//i+1-2; // .g=B
				counter++;
				matCoeffs_y[0][counter] = -(K1y + K2y/72.0f);//.1*(i+1-2)+0.05; // .b=C
				if (i == ny_plus_4-2-1) {
					matCoeffs_y[0][counter] = 0; // .b=C  set the last meaningful C to zero.
				}
				counter++;
				matCoeffs_y[0][counter]=0; // .a=not used.
				counter++;
			}
		}
	}

	arrayLength=ny/2;

	//ST_: Buliding the rest of the matrices using CR.
	for (int k=1; k<countOfMatricesYdirection;++k){
		int counter=0;
		for (int i = 0; i < arrayLength+4; ++i) {
			for (int j = 0; j < nx_plus_4; ++j) {
				if ((i==0 || i==1) || (i==arrayLength+4-2 || i==arrayLength+4-1) ){
					matCoeffs_y[k][counter]=0;// .r=A
					counter++;
					matCoeffs_y[k][counter]=1; // .g=B
					counter++;
					matCoeffs_y[k][counter]=0; // .b=C
					counter++;
					matCoeffs_y[k][counter]=0; // .a=not used.
					counter++;
				}
				else{
					int counterBefore=((i-1)*2+1-1)*(nx+4)*4+j*4; //aInLeft
					int counterHere=((i-1)*2+1)*(nx+4)*4+j*4; //aIn
					int counterAfter=((i-1)*2+1+1)*(nx+4)*4+j*4; //aInRight
		

					float k1= matCoeffs_y[k-1][counterHere]/matCoeffs_y[k-1][counterBefore+1]; // A_i/B_(i-1)
					float k2= matCoeffs_y[k-1][counterHere+2]/matCoeffs_y[k-1][counterAfter+1]; // C_i/B_(i+1)

					matCoeffs_y[k][counter] = -matCoeffs_y[k-1][counterBefore]*k1;  // .r=A   % aOut=-aInLeft*k1
					counter++;

					matCoeffs_y[k][counter] = matCoeffs_y[k-1][counterHere+1]-matCoeffs_y[k-1][counterBefore+2]*k1-matCoeffs_y[k-1][counterAfter]*k2;  // .g=B  % bOut=bIn-cInLeft*k1-aInRight*k2
					counter++;

					matCoeffs_y[k][counter] = -matCoeffs_y[k-1][counterAfter+2]*k2;  // .b=C  % cOut=-cInRight*k2
					counter++;
					matCoeffs_y[k][counter] = 0; // .a=not used.
					counter++;
				}
			}
		}
		arrayLength/=2;
	}
	
    memset(&matCoeffSubData, 0, sizeof(matCoeffSubData));
	arrayLength=ny;
	for (int k=0; k<countOfMatricesYdirection;++k){

		matCoeffSubData.pSysMem = &matCoeffs_y[k][0];
		matCoeffSubData.SysMemPitch = (nx+4) * 4 * sizeof(float);
		CreateTexture(device,    //The Coeff matrix
					  nx_plus_4,
					  arrayLength+4,
					  &matCoeffSubData,
					  DXGI_FORMAT_R32G32B32A32_FLOAT,
					  false,
					  st_psTriDigTexture_ABCy[k],
					  &st_psTriDigTextureView_ABCy[k],
					  &st_psTriDigRenderTargetView_ABCy[k]);
		arrayLength/=2;
	}
}
// Creates and initializes the simulation textures
// Precondition: terrain heightfield is up to date
void ShallowWaterEngine::createSimTextures(ResetType reset_type)
{
	// *******************************************************************************
	// *******************************************************************************
	// ******************************** IMPORTANT NOTE *******************************
	//
	// Stephen (NLSW developer) used i and j as indices for x and y, respectively.
	// Sasan (Boussinesq developer) used i and j as indices for y and x, respectively.
	// Therefore in each [nested] loop you should double check
	// which index is used in what direction.
	//
	// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ READ THIS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
	const int nx_plus_4=nx+4;
	const int ny_plus_4=ny+4;
	const float W = GetSetting("valley_width");
	const float L = GetSetting("valley_length");
	const float delta_x = W / (nx-1);
	const float delta_y = L / (ny-1);
	const float delta_x2 = delta_x * delta_x;
	const float delta_y2 = delta_y * delta_y;
	const float sea_level = initSetting.isBoussinesq ? initSetting.stillWaterElevation : -9999; // is model is not boussinesq set sea level to -9999 to make the coeffs ABC = (0,1,0);

	//const float dam_pos = GetSetting("dam_position");

	//some random numbers to generate initial condition
    const float xmin = -W/6; 
    const float xmax = W/6;
    const float ymin = L/3;
    const float ymax = 2*L/3;

    float init_w = 0;
	/*
    if (reset_type == R_VALLEY) {
        // find the height of the lowest point along the dam
        float B_min = 99999999.f;
        for (float x = -W/2; x < W/2; x += W/(nx-1)) {
            B_min = std::min(B_min, GetTerrainHeight(x, dam_pos));
        }

        init_w = B_min + 1;
    } else if (reset_type == R_SEA) {
        init_w = initSetting.stillWaterElevation;
    }
	*/
            
    boost::scoped_array<float> ic(new float[(nx+4) * (ny+4) * 4]);
    float *p = &ic[0];

	std::ifstream fileIn ((initSetting.initWFileName).c_str());

	for (int j = 0; j < ny + 4; ++j) {

        for (int i = 0; i < nx + 4; ++i) {
            const int ii = std::max(2, std::min(nx+1, i));
            const int jj = std::max(2, std::min(ny+1, j));
            const float B = g_bottom[(nx+4) * jj + ii].BA;
			const float waterDepth = std::max(0.0f, initSetting.stillWaterElevation - B);
            const float x = (ii-2)*W/(nx-1) - W/2;
            const float y = (jj-2)*L/(ny-1);
            float w = B;
			float hu = 0;
			float hv = 0;

			int xIndex,yIndex;
			float alphaIn;
			
			if (fileIn) //if there is a hotstart file.
			{
				fileIn >> xIndex >> yIndex >> w >> hu >> hv >> alphaIn;
			}
			else
			{
				w = std::max(initSetting.stillWaterElevation,B);
				hu = 0; hv = 0;

			}

			for (int scounter = 0 ; scounter < initSetting.countOfSolitons; scounter++)
			{
				Soliton tempSoliton = initSetting.solitons[scounter];
				if (tempSoliton.param_h != 0 ) {

					
					float vert_dist_from_center = abs(sin(tempSoliton.theta) * (delta_x*i-tempSoliton.xc)) +
												  abs(cos(tempSoliton.theta) * (delta_y*j-tempSoliton.yc));

					float tempState[3];
					if (vert_dist_from_center > tempSoliton.length/2.0f ){
						float temp_param_h = 0;
						/*					
						const float x_for_temp_param_h = tempSoliton.xc + tempSoliton.length/2.0f * sin(tempSoliton.theta)
														 + (vert_dist_from_center - tempSoliton.length/2.0f) * cos(tempSoliton.theta);
						const float y_for_temp_param_h = tempSoliton.yc + tempSoliton.length/2.0f * cos(tempSoliton.theta)
														 - (vert_dist_from_center - tempSoliton.length/2.0f) * sin(tempSoliton.theta);

						solitaryWave (waterDepth, tempSoliton.param_h, tempSoliton.theta, tempSoliton.xc, tempSoliton.yc, x_for_temp_param_h, y_for_temp_param_h, &temp_param_h, &tempState[1], &tempState[2]);	
						*/
						const float smoothing_length = PI / sqrt(0.75f * tempSoliton.param_h/pow(waterDepth,3));
						if (vert_dist_from_center - tempSoliton.length/2.0f < smoothing_length){
							temp_param_h = tempSoliton.param_h * (0.5f * cos(PI/smoothing_length *(vert_dist_from_center - tempSoliton.length/2.0f)) + 0.5f); 
						} else {
							temp_param_h  = 0;
						}
						solitaryWave (waterDepth, temp_param_h, tempSoliton.theta, tempSoliton.xc, tempSoliton.yc, delta_x*i, delta_y*j, &tempState[0], &tempState[1], &tempState[2]);	
					}
					else{
						solitaryWave (waterDepth, tempSoliton.param_h, tempSoliton.theta, tempSoliton.xc, tempSoliton.yc, delta_x*i, delta_y*j, &tempState[0], &tempState[1], &tempState[2]);	
					
					}
					w += tempState[0];
					hu += tempState[1];
					hv += tempState[2];
				}
			}

            // initial condition
            *p++ = w;  // w
            *p++ = hu;  // hu
            *p++ = hv;  // hv
            *p++ = 0;  // unused
		}
    }
	fileIn.close();

    D3D11_SUBRESOURCE_DATA sd;
    memset(&sd, 0, sizeof(sd));
    sd.pSysMem = &ic[0];
    sd.SysMemPitch = (nx+4) * 4 * sizeof(float);

    for (int i = 0; i < NUM_OF_TEXTURES; ++i) {
        CreateTexture(device,
                      nx + 4,
                      ny + 4,
                      i < 2 ? &sd : 0,
                      DXGI_FORMAT_R32G32B32A32_FLOAT,
                      false,
                      m_psSimTexture[i],
                      &m_psSimTextureView[i],
                      &m_psSimRenderTargetView[i]);
    }


//// Initializing right hand side (RHS)	
	boost::scoped_array<float> tempRHS;
	tempRHS.reset(new float[(nx_plus_4) * (ny_plus_4) * 4]);
	
	int counter=0;
	for (int j = 0; j < ny_plus_4; ++j) {
		for (int i = 0; i < nx_plus_4; ++i) {
			tempRHS[counter]=0; // .r=h
			counter++;
			tempRHS[counter]=0; // .g=D1x to find hu
			counter++;
			tempRHS[counter]=0; // .b=D1y to find hv
			counter++;
			tempRHS[counter]=0; // .a=not used
			counter++;
		}
	}

	sd.pSysMem = &tempRHS[0];
    sd.SysMemPitch = (nx_plus_4) * 4 * sizeof(float);
	CreateTexture(device,    //The main Coeff matrix
				  nx_plus_4,
				  ny_plus_4,
				  &sd,
				  DXGI_FORMAT_R32G32B32A32_FLOAT,
				  false,
				  st_psRHSTexture,
				  &st_psRHSTextureView,
				  &st_psRHSRenderTargetView);
	
	createTridiagonalCoefTextures();
	
	const int countOfMatricesXdirection = st_getCountOfMatrices(nx);
	int arrayLength=nx;

	// Making D1x matrices.
	boost::scoped_array<float> D1x[ST_MAX_CR_MATRIX_NUM_X];
	
	arrayLength=nx;
	for (int i=0; i<ST_MAX_CR_MATRIX_NUM_X;++i){  //ST_:This can be done till countOfMatricesXdirection, to be on the safe side, I intilize all the ST_MAX_CR_MATRIX_NUM_X scoped arrays.
		D1x[i].reset(new float[(arrayLength+4) * (ny_plus_4)*4]); //ST_: ny can be used instead of ny_plus_4
		arrayLength/=2;
	}
	arrayLength=nx;
	for (int k=0; k<countOfMatricesXdirection;++k){
		counter=0;
		for (int i = 0; i < ny_plus_4; ++i) {
			
			for (int j = 0; j <=1; ++j) {
				D1x[k][counter]=0; // .r=D1x
				counter++;
				D1x[k][counter]=0; // .r=not used
				counter++;
				D1x[k][counter]=0; // .r=not used
				counter++;
				D1x[k][counter]=0; // .r=not used
				counter++;

			}

			for (int j = 2; j <=arrayLength+4-2-1; ++j) {
				D1x[k][counter]=0;//(j-1)/2.0; // .r=D1x
				counter++;
				D1x[k][counter]=0; // .r=not used
				counter++;
				D1x[k][counter]=0; // .r=not used
				counter++;
				D1x[k][counter]=0; // .r=not used
				counter++;
			}

			for (int j = 0; j <=1; ++j) {
				D1x[k][counter]=0; // .r=D1x
				counter++;
				D1x[k][counter]=0; // .r=not used
				counter++;
				D1x[k][counter]=0; // .r=not used
				counter++;
				D1x[k][counter]=0; // .r=not used
				counter++;
			}
		}
		arrayLength/=2;
	}

	D3D11_SUBRESOURCE_DATA matCoeffSubData;
    memset(&matCoeffSubData, 0, sizeof(matCoeffSubData));
	arrayLength=nx;
	for (int k=0; k<countOfMatricesXdirection;++k){

		matCoeffSubData.pSysMem = &D1x[k][0];
		matCoeffSubData.SysMemPitch = (arrayLength+4) * 4 * sizeof(float); 
		CreateTexture(device,    //The D1x matrix, //Initialization to zero!
				  arrayLength+4,
				  ny_plus_4,
				  &matCoeffSubData,
				  DXGI_FORMAT_R32G32B32A32_FLOAT,
				  false,
				  st_psTriDigTexture_D1x[k],
				  &st_psTriDigTextureView_D1x[k],
				  &st_psTriDigRenderTargetView_D1x[k]);

		CreateTexture(device,    //The D1x matrix copy, used to get the unkown values in Ax=D.
					  arrayLength+4,
					  ny_plus_4,
					  &matCoeffSubData,
					  DXGI_FORMAT_R32G32B32A32_FLOAT,
					  false,
					  st_psTriDigTexture_D1xCopy[k],
					  &st_psTriDigTextureView_D1xCopy[k],
					  &st_psTriDigRenderTargetView_D1xCopy[k]);

	    CreateTexture(device,
              arrayLength+4,
              ny + 4,
              0,
              DXGI_FORMAT_R32G32B32A32_FLOAT,
              true,
              st_psStagingTextureX[k],
              0,
              0);

		arrayLength/=2;
	}

	// Making D1y matrices.
	const int countOfMatricesYdirection = st_getCountOfMatrices(ny);
	boost::scoped_array<float> D1y[ST_MAX_CR_MATRIX_NUM_Y];
	
	arrayLength=ny;
	for (int i=0; i<ST_MAX_CR_MATRIX_NUM_Y;++i){ 
		D1y[i].reset(new float[(arrayLength+4) * (nx_plus_4)*4]); //ST_: ny could be used instead of ny_plus_4
		arrayLength/=2;
	}
	arrayLength=ny;
	for (int k=0; k<countOfMatricesYdirection;++k){
		counter=0;
		for (int i = 0; i < arrayLength+4; ++i) {
			for (int j = 0; j <nx_plus_4; ++j) {
				if ((i==0 || i==1) || (i==arrayLength+4-2 || i==arrayLength+4-1) ){
					D1y[k][counter]=0; // .r=not used
					counter++;
					D1y[k][counter]=0; // .g=not used
					counter++;
					D1y[k][counter]=0; // .b=D1y
					counter++;
					D1y[k][counter]=0; // .a=not used
					counter++;

				}
				else{
					D1y[k][counter]=0; // .r=not used
					counter++;
					D1y[k][counter]=0; // .g=not used
					counter++;
					D1y[k][counter]=0;//(i-1)/2.0f; // .b=D1y
					counter++;
					D1y[k][counter]=0; // .a=not used
					counter++;
				}
			}
		}
		arrayLength/=2;
	}

    memset(&matCoeffSubData, 0, sizeof(matCoeffSubData));
	arrayLength=ny;
	for (int k=0; k<countOfMatricesYdirection;++k){

		matCoeffSubData.pSysMem = &D1y[k][0];
		matCoeffSubData.SysMemPitch = (nx+4) * 4 * sizeof(float); 
		CreateTexture(device,    //The D1y matrix, //Initialization to zero!
				  nx_plus_4,
				  arrayLength+4,
				  &matCoeffSubData,
				  DXGI_FORMAT_R32G32B32A32_FLOAT,
				  false,
				  st_psTriDigTexture_D1y[k],
				  &st_psTriDigTextureView_D1y[k],
				  &st_psTriDigRenderTargetView_D1y[k]);

		CreateTexture(device,    //The D1y matrix copy, used to get the unkown values in Ax=D.
					  nx_plus_4,
					  arrayLength+4,
					  &matCoeffSubData,
					  DXGI_FORMAT_R32G32B32A32_FLOAT,
					  false,
					  st_psTriDigTexture_D1yCopy[k],
					  &st_psTriDigTextureView_D1yCopy[k],
					  &st_psTriDigRenderTargetView_D1yCopy[k]);

	    CreateTexture(device,
              nx + 4,
			  arrayLength+4,
              0,
              DXGI_FORMAT_R32G32B32A32_FLOAT,
              true,
              st_psStagingTextureY[k],
              0,
              0);

		arrayLength/=2;
	}

	// TODO: This has no need to be (nx+4) by (ny+4),
    // we only ever use the top left quarter of it ((nx/4) by (ny/4))...
    // (Perhaps could do two or three separate drawcalls in GetStats pass instead of one multiple-output drawcall.)
    CreateTexture(device,
                  nx + 4,
                  ny + 4,
                  0,
                  DXGI_FORMAT_R32_FLOAT,
                  false,
                  m_psGetStatsTexture,
                  0,
                  &m_psGetStatsRenderTargetView);
    
    sim_idx = 0;
	myDelay = 0;
	old_index = 7;
	old_old_index = 8;
	scratch_index = 9;
	predicted_index = 10;
	F_G_star_old_index = 11;
	F_G_star_old_old_index = 12;
	F_G_star_scratch_index = 13;
	F_G_star_predicted_index = 14;
    bootstrap_needed = true;
	firstTimeStep = true;
	secondTimeStep = false;
	timeScheme = euler; //initialization
	timestep_count = 1;
	correction_steps_count = initSetting.correctionStepsNum;

    
    // Create a staging texture so we can read the data back again when required
    // (e.g. for debugging, or when changing terrain level)
    // TODO: we should be able to get away without using this staging texture; that would save a bit of memory
    CreateTexture(device,
                  nx + 4,
                  ny + 4,
                  0,
                  DXGI_FORMAT_R32G32B32A32_FLOAT,
                  true,
                  m_psFullSizeStagingTexture,
                  0,
                  0);

    // Create another staging texture of size (NX/2) * (NY/4) for GetStats
    CreateTexture(device,
                  nx/2,
                  ny/4,
                  0,
                  DXGI_FORMAT_R32G32B32A32_FLOAT,
                  true,
                  m_psGetStatsStagingTexture4,
                  0,
                  0);

    // and another of size (NX/4) * (NY/4) with only one channel
    CreateTexture(device,
                  nx/4,
                  ny/4,
                  0,
                  DXGI_FORMAT_R32_FLOAT,
                  true,
                  m_psGetStatsStagingTexture1,
                  0,
                  0);
}



// Updates the terrain texture given current settings.
// Also updates the water state texture such that the water depth remains unchanged.
// NOTE: If size has changed then call createTerrainTexture first.
void ShallowWaterEngine::fillTerrainTexture()
{
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");

    // Get the existing water state onto the CPU
    context->CopyResource(m_psFullSizeStagingTexture.get(), m_psSimTexture[sim_idx].get());

    // Now loop through and change w values into h values
    // TODO: it would probably be better to do this on the GPU (would avoid a copy / copy back; and could
    // save memory for the staging texture as well).
    {
        MapTexture m(*context, *m_psFullSizeStagingTexture);
        for (int j = 0; j < ny+4; ++j) {
            char * row_ptr = reinterpret_cast<char*>(m.msr.pData) + j * m.msr.RowPitch;
            const BottomEntry *B_row_ptr = &g_bottom[j * (nx+4)];
            
            for (int i = 0; i < nx+4; ++i) {
                float *col_ptr = reinterpret_cast<float*>(row_ptr) + i * 4;
                const BottomEntry *B_col_ptr = B_row_ptr + i;
                
                (*col_ptr) -= (B_col_ptr->BA);
            }
        }
        
        // Compute the new terrain heightfield
        UpdateTerrainHeightfield();
        
        // Write the new terrain textures
        context->UpdateSubresource(m_psTerrainTexture.get(),
                                   0,   // subresource
                                   0,   // overwrite whole resource
                                   &g_terrain_heightfield[0],
                                   (nx+4) * 12,
                                   0);  // slab pitch
        context->UpdateSubresource(m_psBottomTexture.get(),
                                   0,  // subresource
                                   0,  // overwrite whole resource
                                   &g_bottom[0],
                                   (nx+4) * 12,
                                   0); // slab pitch

		context->CopyResource(m_psAuxiliary1Texture.get(), st_psTriDigTexture_D1x[0].get()); 
		context->CopyResource(m_psAuxiliary2Texture.get(), st_psTriDigTexture_D1x[0].get());

        // Loop through and change h values back into w values
        for (int j = 0; j < ny+4; ++j) {
            char * row_ptr = reinterpret_cast<char*>(m.msr.pData) + j * m.msr.RowPitch;
            const BottomEntry *B_row_ptr = &g_bottom[j * (nx+4)];
            
            for (int i = 0; i < nx+4; ++i) {
                float *col_ptr = reinterpret_cast<float*>(row_ptr) + i * 4;
                const BottomEntry *B_col_ptr = B_row_ptr + i;
                
                (*col_ptr) += (B_col_ptr->BA);
            }
        }
    }
    
    // Copy water texture back to the GPU
    context->CopyResource(m_psSimTexture[sim_idx].get(), m_psFullSizeStagingTexture.get());

    // need to re-bootstrap
    bootstrap_needed = true;
}



// Updates the terrain texture, but does not touch the water texture
void ShallowWaterEngine::fillTerrainTextureLite()
{
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");

    // Compute the new terrain heightfield
    UpdateTerrainHeightfield();
        
    // Write the new terrain textures
    context->UpdateSubresource(m_psTerrainTexture.get(),
                               0,   // subresource
                               0,   // overwrite whole resource
                               &g_terrain_heightfield[0],
                               (nx+4) * 12,
                               0);  // slab pitch
    context->UpdateSubresource(m_psBottomTexture.get(),
                               0,  // subresource
                               0,  // overwrite whole resource
                               &g_bottom[0],
                               (nx+4) * 12,
                               0); // slab pitch

	context->CopyResource(m_psAuxiliary1Texture.get(), st_psTriDigTexture_D1x[0].get());
	context->CopyResource(m_psAuxiliary2Texture.get(), st_psTriDigTexture_D1x[0].get());
		
    
    // need to re-bootstrap
    bootstrap_needed = true;
}



// create the constant buffers, leaving them uninitialized
void ShallowWaterEngine::createConstantBuffers()
{
    D3D11_BUFFER_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.ByteWidth = RoundUpTo16(sizeof(MyConstBuffer));
    bd.Usage = D3D11_USAGE_DEFAULT;
    bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;

    ID3D11Buffer *pBuffer;
    HRESULT hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psConstantBuffer.reset(pBuffer);

    bd.ByteWidth = RoundUpTo16(sizeof(SimConstBuffer));
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psSimConstantBuffer.reset(pBuffer);

    bd.ByteWidth = RoundUpTo16(sizeof(BoundaryConstBuffer));
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psBoundaryConstantBuffer.reset(pBuffer);


	bd.ByteWidth = RoundUpTo16(sizeof(TimeIntegrationConstBuffer));
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psTimeIntegrationConstantBuffer.reset(pBuffer);


	bd.ByteWidth = RoundUpTo16(sizeof(IrregularWavesDataConstBuffer));
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psIrregularWavesDataConstantBuffer.reset(pBuffer);

	bd.ByteWidth = RoundUpTo16(sizeof(IrregularWavesColumnRowConstBuffer));
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psIrregularWavesColumnRowConstBuffer.reset(pBuffer);

	bd.ByteWidth = RoundUpTo16(sizeof(SolitaryWaveConstBuffer));
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psSolitaryWaveConstantBuffer.reset(pBuffer);	
	
	fillIrregularWavesDataConstantBuffer ("ALL");  
	fillUniformTimeSeriesMainMemoryBuffer ();
}


void decodeLine2 (std::string* tag, std::string* value, std::string line) {
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

int getNumberOfWavesFromFileStream (std::ifstream &in){
	int result = -1;
	if (in){
		std::string line = "-";
		std::string tag = "";
		std::string value = "";
		while (std::getline(in,line))
		{
			if (line[0] == '=')
			{
				break;
			} else if (line != "" && line[0] != '#'){
				tag.clear();
				value.clear();
				decodeLine2(&tag, &value, line);
				
				if (tag != "Invalid" && value != "Invalid" && value != "") {
					if (tag == "[NumberOfWaves]") {
						result = std::min(MAX_NUM_OF_IRREGULAR_WAVES,atoi(value.c_str()));
					}
				}
			}
		}
	}
	return result;
}

void skipCommentsInFileStream (std::ifstream &in){
	if (in){
		std::string line = "-";
		while (std::getline(in,line))
		{
			if (line[0] == '=')
			{
				break;
			} 
		}
	}
}


std::string ShallowWaterEngine::generateSpectrumWaves(IrregularWaveSpectrumSetting inSpectrum, std::string boundarySide)
{
	if(inSpectrum.depth_at_source == 0 || inSpectrum.peak_period == 0 || inSpectrum.significant_waveheight == 0){
		return "NA";
	}
	// function file to create random wave input files (directional, TMA) for COULWAVE
	// Hs_o=2.; %Hmo at source depth in m
	// Tp=12; % peak period in s
	// Thetap=0;  %mean wave direction, relative to source line
	float gamma_s = 3.3f;  // frequency spreading factor for TMA spectrum
	float spread_o = 50.0f;  // direction spreading factor
	//float h = 10.0f;  // depth at source region in m
	float del_f = 0.005f; // frequency increment, time series will repeat at 1/del_f seconds
	float del_t = 5;  // direction increment
	
	float f_peak = 1/inSpectrum.peak_period;
	float Hs = inSpectrum.significant_waveheight;
	float g = GetSetting("gravity");
	float h = inSpectrum.depth_at_source;
	Hs = std::min(Hs, h * 0.5f);  // limit wave height to 1/2 of water depth at generation
	float Hs_old = Hs;

   // Shalllow water TMA spectrum
    float beta = 0.0624f/(0.23f+0.033f * gamma_s - 0.185f / (1.9f + gamma_s));
    
    // Calculate Energy Spectrum to Very Large Limits
	float f_start = std::max(del_f, f_peak - 5 * del_f);
    float f_end = f_peak + 5 * del_f;

	std::vector<float> f_array;
	float current_value = f_start;
    while(current_value <= f_end) {
        f_array.push_back(current_value);
        current_value += del_f;
    }

	std::vector<float> E_array;
	std::vector<std::vector<float>> D_matrix;
	for (int i = 0; i < f_array.size(); ++i){
		
		float omega_h = 2.0f * 3.1415f * f_array[i] * sqrt(h / g);
		float phiK = 0;
		if (omega_h > 2.0f){
			phiK = 1.0f;
		} else if (omega_h < 1) {
			phiK = 0.5 * omega_h;
		} else {
			phiK = 1.0f - 0.5f * (2.0f - omega_h) * (2.0f - omega_h);
		}

		float sigma = 0;       
		if (f_array[i] <= f_peak){
            sigma = 0.07f;
		} else {
            sigma = 0.09f;
		}

        float frat = f_array[i] / f_peak;
		E_array.push_back(beta * Hs * Hs / (f_array[i] * pow(frat, 4))
							* exp(-1.25f / pow(frat, 4))
							* pow(gamma_s, (exp(-(frat - 1.0f) * (frat - 1.0f) / (2.0f * sigma * sigma)))));  // energy density
	}
	    
	// Directional Spectrum to Very Large Limits
	std::vector<float> theta_array;
	current_value = -20 + inSpectrum.mean_theta;
	while( current_value <= 20 + inSpectrum.mean_theta) {
        theta_array.push_back(current_value);
        current_value += del_t;
    }

	for (int i = 0; i < f_array.size(); ++i){
		float theta_peak = inSpectrum.mean_theta;  // (degrees) Input some function here to determine the mean wave direction as a function of frequency
		float f_rat = f_array[i] / f_peak; 
		float spread = 0;

		if (f_rat < 1.0f){
			spread = spread_o * pow(f_rat, 5);
		} else {
			spread = spread_o * pow(f_rat, -2.5f);
		}
		
		
		
		float beta_s = pow(2, (2 * spread - 1)) / 3.1415f * exp(2*boost::math::lgamma<float>(spread+1) - boost::math::lgamma<float>(2*spread+1));
		std::vector<float> Temp_array;
		for (int j = 0; j < theta_array.size(); ++j){

			Temp_array.push_back(beta_s * pow(cos(0.5f * ( (theta_array[j] - theta_peak) * 3.1415f / 180.0f)), (2.0f * spread)));  // Directional function
		}
		D_matrix.push_back(Temp_array);
	}   
    
	for (int i = 0; i < f_array.size(); ++i){
		float sum = 0;
		for (int j = 0; j < theta_array.size(); ++j){
			sum = sum + D_matrix[i][j];  // Directional Energy Density function
		}
		for (int j = 0; j < theta_array.size(); ++j){
			D_matrix[i][j] = D_matrix[i][j] / sum;
		}
		
	}
	    
	std::vector<std::vector<float>> E_D_matrix;

	for (int i = 0; i < f_array.size(); ++i){
		std::vector<float> Temp_array;
		for (int j = 0; j < theta_array.size(); ++j){
			Temp_array.push_back(E_array[i] * D_matrix[i][j]);  // Directional Energy Density function
		}
		E_D_matrix.push_back(Temp_array);
	}
	  


	float Hmo = 0;

	for (int i = 0; i < f_array.size(); ++i){
		for (int j = 0; j < theta_array.size(); ++j){
			if (i == 0) {
				del_f = (f_array[2] - f_array[1]);
			} else if (i == f_array.size() - 1){
				del_f = (f_array[f_array.size() - 1] - f_array[f_array.size() - 2]);
			} else {
				del_f = (f_array[i+1] - f_array[i]) / 2.0f + (f_array[i] - f_array[i-1]) / 2.0f;
			}
			Hmo = Hmo + E_D_matrix[i][j] * del_f;
		}
	}

	float Hmo_full_spectrum = sqrt(Hmo) * 4.004f;
	
	// Truncate ends of spectrum at values of 5% of the max
	float trunc = 0.05;
/*
	%DO NOT CHANGE ANYTHING BELOW THIS LINE
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

	int j = floor(theta_array.size() / 2.0f + + 0.5f); // + 0.5 turns the floor function to round.
	std::vector<float> f_num;
	std::vector<std::vector<float>> E_D_num;
	int count = 0;

	if (Hmo == 0){
		f_num.push_back(f_peak);
		float E_num = 0;
		float amp_max = 0;
		count = 1;
	} else {
		float max_E = f_array.size() > 0 ? E_D_matrix[0][j] : 0;
		for (int i = 0; i < f_array.size(); ++i){
			max_E = std::max(max_E, E_D_matrix[i][j]);
		}
		float min_E = trunc * max_E;
		count = 0;
		for (int i = 0; i < f_array.size(); ++i){
			if (E_D_matrix[i][j] > min_E){
				count = count + 1;
				f_num.push_back(f_array[i]);
				E_D_num.push_back (E_D_matrix[i]);
			}
		}
	}



	/////////////

	time_t timer = time(0);
	struct tm * now = localtime(&timer);
	std::ostringstream s;
	s << now->tm_year + 1900 << "-" << now->tm_mon + 1 << "-" << now->tm_mday << " " << now->tm_hour << "-" << now->tm_min << "-" << now->tm_sec;
	

	std::string fileName = initSetting.logPath + "/" + s.str() + "-Irr" + boundarySide + ".txt";	
	std::ofstream myfile;
	if (!myfile.is_open()){
		myfile.open (fileName.c_str(),std::ios::out | std::ios::app);
	}

	myfile << "[Hs] " << inSpectrum.significant_waveheight << "\n";
	myfile << "[Tp] " << inSpectrum.peak_period << "\n";
	myfile << "[Theta_p] " << inSpectrum.mean_theta << "\n";
	myfile << "[h_source] " << inSpectrum.depth_at_source << "\n";
	myfile << "\n" << "\n";
	myfile << "# This file contains the sinewaves generated for the spectrum given above." << "\n\n";
	myfile << "[NumberOfWaves] " << count*theta_array.size() << "\n";	
	myfile << "\n" << "\n";
	myfile << "====================================" << "\n";

	srand (time(0));
	float Hmo_truncated_spectrum = 0;
	float cur_ind = 0;
	for (int i = 0; i < count; ++i){
		if (i == 0){
			del_f = (f_array[1] - f_array[0]);
		} else if (i == count - 1) {
			del_f = (f_array[count - 1] - f_array[count - 2]);
		} else { 
			del_f = (f_array[i + 1] - f_array[i]) / 2.0f + (f_array[i] - f_array[i - 1]) / 2.0f;
		}

		for (int j = 0; j < theta_array.size(); ++j){ 
			float amp = sqrt(2 * E_D_num[i][j] * del_f);
			++cur_ind;
			
			float period = 1 / f_num[i];
			float theta = theta_array[j] * 3.1415 / 180.0f;

			float rand_between_0_and_1 = (rand()%10000) / 10000.0f;
			float phase = rand_between_0_and_1 * 2 * 3.1415;
			myfile << amp << "\t" << period << "\t" << theta << "\t" << phase << "\n";
		} 
	}
	myfile.close();

	return fileName;
}
void ShallowWaterEngine::fillIrregularWavesDataConstantBuffer (std::string boundarySide)
{

	IrregularWavesDataConstBuffer irr_cBuffer;
	bool doThisSide = false; // This could be chosen from an enum.

	doThisSide = boundarySide == "ALL" || boundarySide == "WEST"; 
	if (doThisSide && initSetting.westBoundary.type == "IrregularWaves"){
		std::string irrFilePath = "";
		if(initSetting.westBoundary.IrrWaveSpectrumSetting.useSpectrum){
			irrFilePath = generateSpectrumWaves(initSetting.westBoundary.IrrWaveSpectrumSetting, "West");
		} else {
			irrFilePath = initSetting.westIrrWaveFileName;
		}
		
		std::ifstream in (irrFilePath.c_str());
		irr_cBuffer.numberOfWavesWest = getNumberOfWavesFromFileStream (in);
		int count = 0;
		while (in && count <= irr_cBuffer.numberOfWavesWest) {
			in >> irr_cBuffer.wavesWest[count].amplitude >> irr_cBuffer.wavesWest[count].period >> irr_cBuffer.wavesWest[count].theta >> irr_cBuffer.wavesWest[count].phase;
			++count;
		}
	}

	doThisSide = boundarySide == "ALL" || boundarySide == "EAST"; 
	if (doThisSide && initSetting.eastBoundary.type == "IrregularWaves"){
		std::string irrFilePath = "";
		if(initSetting.eastBoundary.IrrWaveSpectrumSetting.useSpectrum){
			irrFilePath = generateSpectrumWaves(initSetting.eastBoundary.IrrWaveSpectrumSetting, "East");
		} else {
			irrFilePath = initSetting.eastIrrWaveFileName;
		}
		
		std::ifstream in (irrFilePath.c_str());
		irr_cBuffer.numberOfWavesEast = getNumberOfWavesFromFileStream (in);
		int count = 0;
		while (in && count <= irr_cBuffer.numberOfWavesEast) {
			in >> irr_cBuffer.wavesEast[count].amplitude >> irr_cBuffer.wavesEast[count].period >> irr_cBuffer.wavesEast[count].theta >> irr_cBuffer.wavesEast[count].phase;
			++count;
		}
	}

	doThisSide = boundarySide == "ALL" || boundarySide == "SOUTH"; 
	if (doThisSide && initSetting.southBoundary.type == "IrregularWaves"){
		std::string irrFilePath = "";
		if(initSetting.southBoundary.IrrWaveSpectrumSetting.useSpectrum){
			irrFilePath = generateSpectrumWaves(initSetting.southBoundary.IrrWaveSpectrumSetting, "South");
		} else {
			irrFilePath = initSetting.southIrrWaveFileName;
		}
		
		std::ifstream in (irrFilePath.c_str());
		irr_cBuffer.numberOfWavesSouth = getNumberOfWavesFromFileStream (in);
		int count = 0;
		while (in && count <= irr_cBuffer.numberOfWavesSouth) {
			in >> irr_cBuffer.wavesSouth[count].amplitude >> irr_cBuffer.wavesSouth[count].period >> irr_cBuffer.wavesSouth[count].theta >> irr_cBuffer.wavesSouth[count].phase;
			++count;
		}
	}

	doThisSide = boundarySide == "ALL" || boundarySide == "NORTH"; 
	if (doThisSide && initSetting.northBoundary.type == "IrregularWaves"){
		std::string irrFilePath = "";
		if(initSetting.northBoundary.IrrWaveSpectrumSetting.useSpectrum){
			irrFilePath = generateSpectrumWaves(initSetting.northBoundary.IrrWaveSpectrumSetting, "North");
		} else {
			irrFilePath = initSetting.northIrrWaveFileName;
		}
		
		std::ifstream in (irrFilePath.c_str());
		irr_cBuffer.numberOfWavesNorth = getNumberOfWavesFromFileStream (in);
		int count = 0;
		while (in && count <= irr_cBuffer.numberOfWavesNorth) {
			in >> irr_cBuffer.wavesNorth[count].amplitude >> irr_cBuffer.wavesNorth[count].period >> irr_cBuffer.wavesNorth[count].theta >> irr_cBuffer.wavesNorth[count].phase;
			++count;
		}
	}

	context->UpdateSubresource(m_psIrregularWavesDataConstantBuffer.get(), 0, 0, &irr_cBuffer, 0, 0);
}

void ShallowWaterEngine::fillUniformTimeSeriesMainMemoryBuffer ()
{

	if (initSetting.westBoundary.type == "UniformTimeSeries"){
		std::ifstream in ((initSetting.westBoundary.uniformTimeSeries.fileName).c_str());
		skipCommentsInFileStream(in);
		float temp_time, temp_eta, temp_hu, temp_hv;
		while (in) {
			in >> temp_time >> temp_eta >> temp_hu >> temp_hv;
			std::vector<float> tempValues;
			tempValues.push_back(temp_time); tempValues.push_back(temp_eta); tempValues.push_back(temp_hu); tempValues.push_back(temp_hv);
			initSetting.westBoundary.uniformTimeSeries.data.push_back(tempValues); 
		}
	}

	if (initSetting.eastBoundary.type == "UniformTimeSeries"){
		std::ifstream in ((initSetting.eastBoundary.uniformTimeSeries.fileName).c_str());
		skipCommentsInFileStream(in);
		float temp_time, temp_eta, temp_hu, temp_hv;
		while (in) {
			in >> temp_time >> temp_eta >> temp_hu >> temp_hv;
			std::vector<float> tempValues;
			tempValues.push_back(temp_time); tempValues.push_back(temp_eta); tempValues.push_back(temp_hu); tempValues.push_back(temp_hv);
			initSetting.eastBoundary.uniformTimeSeries.data.push_back(tempValues); 
		}
	}

	
	if (initSetting.northBoundary.type == "UniformTimeSeries"){
		std::ifstream in ((initSetting.northBoundary.uniformTimeSeries.fileName).c_str());
		skipCommentsInFileStream(in);
		float temp_time, temp_eta, temp_hu, temp_hv;
		while (in) {
			in >> temp_time >> temp_eta >> temp_hu >> temp_hv;
			std::vector<float> tempValues;
			tempValues.push_back(temp_time); tempValues.push_back(temp_eta); tempValues.push_back(temp_hu); tempValues.push_back(temp_hv);
			initSetting.northBoundary.uniformTimeSeries.data.push_back(tempValues); 
		}
	}
	
	if (initSetting.southBoundary.type == "UniformTimeSeries"){
		std::ifstream in ((initSetting.southBoundary.uniformTimeSeries.fileName).c_str());
		skipCommentsInFileStream(in);
		float temp_time, temp_eta, temp_hu, temp_hv;
		while (in) {
			in >> temp_time >> temp_eta >> temp_hu >> temp_hv;
			std::vector<float> tempValues;
			tempValues.push_back(temp_time); tempValues.push_back(temp_eta); tempValues.push_back(temp_hu); tempValues.push_back(temp_hv);
			initSetting.southBoundary.uniformTimeSeries.data.push_back(tempValues); 
		}
	}
}

inline bool isTerraingTextureColormap(int t){
	return (t == 6 || t == 7) ? 1:0; // 6 and 7 are colormaps now, but always double check!
}
// update the constant buffers given current settings
void ShallowWaterEngine::fillConstantBuffers()
{
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");

    MyConstBuffer cb;

    const float sun_alt = GetSetting("sun_alt") * PI / 180.0f;
    const float sun_az = GetSetting("sun_az") * PI / 180.0f;
    cb.light_dir = XMFLOAT3(cos(sun_alt)*sin(sun_az), cos(sun_alt)*cos(sun_az), sin(sun_alt));
    cb.ambient = GetSetting("ambient");

    cb.eye_mult = XMFLOAT3(-W/(nx-1), -L/(ny-1), -1.0f);
    cb.eye_trans = XMFLOAT3(W/2 + 2*W/(nx-1) + camera_x, 2*L/(ny-1) + camera_y, camera_z);

    cb.world_mult_x = W/(nx-1);
    cb.world_mult_y = L/(ny-1);
    cb.world_trans_x = -W/2 - 2*W/(nx-1);
    cb.world_trans_y = -2*L/(ny-1);

	const float mesh_size = GetSetting("Grid Scale") * 8.0f; // 8 is the number of cells in the gird texture in each directions
    const float terrain_tex_width = W/(nx-1) * mesh_size; // width of the entire BMP in metres
	const float terrain_tex_length = L/(ny-1) * mesh_size; // width of the entire BMP in metres

	cb.isGridOn  = GetIntSetting("Grid");

	cb.is_dissipation_threshold_on  = GetIntSetting("Dissipation Intensity");
	cb.dissipation_threshold  = GetSetting("Dissipation Threshold");
	
	cb.terrain_tex_scale = XMFLOAT2(W / (nx-1) / terrain_tex_width, L / (ny-1) / terrain_tex_length);
    cb.world_to_grass_tex_x = 1.0f / terrain_tex_width;
    cb.world_to_grass_tex_y = 1.0f / terrain_tex_length;
    cb.grass_tex_of_origin_x = W/2 / terrain_tex_width;
    cb.grass_tex_of_origin_y = 0;

/*
	cb.terrain_tex_scale = XMFLOAT2(1.0f / (nx-1), 1.0f / (ny-1));
    cb.world_to_grass_tex_x = 1.0f / W ;
    cb.world_to_grass_tex_y = 1.0f / L;
    cb.grass_tex_of_origin_x = 1.0f/2.0f;
    cb.grass_tex_of_origin_y = 0;
*/
	cb.water_colormap_min = GetSetting("Colormap Min");
	cb.water_colormap_max =  GetSetting("Colormap Max");
		
	cb.terrain_colormap_min = GetSetting("Colormap Min ");
	cb.terrain_colormap_max =  GetSetting("Colormap Max ");


    /*
    cb.fresnel_exponent = GetSetting("fresnel_exponent");
    cb.specular_intensity = GetSetting("specular_intensity");
    cb.specular_exponent = GetSetting("specular_exponent");
    */

    cb.refractive_index = GetSetting("refractive_index");

	// HARDCODED GRAPHICS VALUES
	const float FRESNEL_EXPONENT = 5.0f;
	const float SPECULAR_INTENSITY = 1.0f;
	const float SPECULAR_EXPONENT = 15.0f;

//	const float ATTENUATION_1 = 0.08f;
//	const float ATTENUATION_2 = 0.08f;
	const float DEEP_R = 0.05f;
	const float DEEP_G = 0.1f;
	const float DEEP_B = 0.2f;


	cb.fresnel_exponent = FRESNEL_EXPONENT;
    cb.specular_intensity = SPECULAR_INTENSITY;
    cb.specular_exponent = SPECULAR_EXPONENT;

	cb.attenuation_1 = GetSetting("attenuation_1");
	
    cb.attenuation_2 = GetSetting("attenuation_2");
	/*
    cb.deep_r = GetSetting("deep_r");
    cb.deep_g = GetSetting("deep_g");
    cb.deep_b = GetSetting("deep_b");

	
	cb.attenuation_1 = ATTENUATION_1;
    cb.attenuation_2 = ATTENUATION_2;
	*/
    cb.deep_r = DEEP_R;
    cb.deep_g = DEEP_G;
    cb.deep_b = DEEP_B;

    cb.nx_plus_1 = nx + 1;
    cb.ny_plus_1 = ny + 1;
	cb.sqrt_sqrt_epsilon = sqrt(sqrt(abs(initSetting.epsilon)));
	cb.drylandDepthOfInundation =  GetSetting("Inundated Area")? GetSetting("Flow Depth") : 0;
/*
	cb.dx =  W / (nx - 1);
    cb.dy =  L / (ny - 1);
*/
	cb.zScale = GetSetting("Vertical_Scale");
	cb.seaLevel = initSetting.stillWaterElevation;
		
	static bool firstCall = true;
	static int shading_status = GetIntSetting("Surface Shading");
	if (shading_status != GetIntSetting("Surface Shading") || firstCall)
	{
		shading_status = GetIntSetting("Surface Shading");
		if (shading_status == PHOTOREALISTIC){
			SetSettingD("fresnel_coeff",initSetting.graphics.fresnelCoef);
		} else {
				loadColormap(shading_status);
				SetSettingD("fresnel_coeff",0);
		}
		
	}

	static int terrain_texture = GetIntSetting("Terrain Texture");
	if (terrain_texture != GetIntSetting("Terrain Texture") || firstCall)
	{
		terrain_texture = GetIntSetting("Terrain Texture");
		loadTerrainShading(terrain_texture);
	}

	static int skybox_texture = GetIntSetting("Skybox");
	if (skybox_texture != GetIntSetting("Skybox") || firstCall)
	{
		skybox_texture = GetIntSetting("Skybox");
		loadSkybox(skybox_texture);
		firstCall = false;
	}


	
	cb.water_shading = GetIntSetting("Surface Shading") ? GetIntSetting("Shading Variable") + 1 : 0;
	
	
	cb.terrain_shading = isTerraingTextureColormap(terrain_texture);
	cb.fresnel_coeff = GetSetting("fresnel_coeff");
	
	

    // setup matrix transforming (tex_x, tex_y, world_z, 1) into normalized device coordinates.

    const XMMATRIX tex_to_world( W / (nx - 1), 0,             0,  -W/2 - 2*W/(nx-1),
                                 0,            L / (ny - 1),  0,       - 2*L/(ny-1),
                                 0,            0,             1,  0,
                                 0,            0,             0,  1  );
    
    const XMMATRIX translate_camera_pos( 1, 0, 0, -camera_x,
                                         0, 1, 0, -camera_y,
                                         0, 0, 1, -camera_z,
                                         0, 0, 0, 1 );

    const XMMATRIX rotate_around_z( cos(yaw), -sin(yaw), 0, 0,
                                    sin(yaw), cos(yaw),  0, 0,
                                    0,        0,         1, 0,
                                    0,        0,         0, 1 );

    const XMMATRIX rotate_around_x( 1, 0,           0,          0,
                                    0, cos(pitch),  sin(pitch), 0,
                                    0, -sin(pitch), cos(pitch), 0,
                                    0, 0,           0,          1 );

    const XMMATRIX swap_y_and_z( 1, 0, 0, 0,
                                 0, 0, 1, 0,
                                 0, 1, 0, 0,
                                 0, 0, 0, 1 );

    const float far = far_plane_dist;
    const float near = near_plane_dist;
    const XMMATRIX perspective( xpersp, 0,      0,                  0,
                                0,      ypersp, 0,                  0,
                                0,      0,      far / (far - near), -near * far / (far - near),
                                0,      0,      1,                  0  );

    cb.tex_to_clip = perspective * swap_y_and_z * rotate_around_x
        * rotate_around_z * translate_camera_pos * tex_to_world;

    cb.skybox_mtx = perspective * swap_y_and_z * rotate_around_x * rotate_around_z;

    // Now write it to the const buffer
    context->UpdateSubresource(m_psConstantBuffer.get(),
                               0,     // subresource
                               0,     // overwrite whole resource
                               &cb,
                               0,     // row pitch
                               0);    // slab pitch


    // Now do the sim parameters
    SimConstBuffer sb;

    const float theta = GetSetting("Theta");  // 1=most dissipative, 2=most oscillatory, 1.3 = good default.
    const float g = GetSetting("gravity");
    const float dt = current_timestep;
    
    const float dx = W / (nx-1);
    const float dy = L / (ny-1);
    
    sb.two_theta = 2 * theta;
    sb.two_over_nx_plus_four = 2.0f / (nx+4);
    sb.two_over_ny_plus_four = 2.0f / (ny+4);
    sb.g = g;
    sb.half_g = 0.5f * g;
	sb.Bcoef_g=(BCOEF)*g;
    sb.g_over_dx = g / dx;
    sb.g_over_dy = g / dy;
    sb.one_over_dx = 1.0f / dx;
    sb.one_over_dy = 1.0f / dy;
    sb.dt = dt;
    sb.epsilon = CalcEpsilon();
    sb.nx = nx;
    sb.ny = ny;
	sb.friction = GetSetting("friction");
	sb.isManning= initSetting.isManning;
	sb.seaLevel= initSetting.stillWaterElevation;
	sb.dissipation_threshold  = GetSetting("Dissipation Threshold");
	sb.whiteWaterDecayRate = 1.0f / (1 + 2 * GetSetting("Whitewater Decay"));
    // Now write it to the constant buffer
    context->UpdateSubresource(m_psSimConstantBuffer.get(),
                               0,  // subresource
                               0,  // overwrite whole buffer
                               &sb,
                               0,  // row pitch
                               0); // slab pitch
}

void ShallowWaterEngine::createDepthStencil(int w, int h)
{
    D3D11_TEXTURE2D_DESC td;
    memset(&td, 0, sizeof(td));
    td.Width = w;
    td.Height = h;
    td.MipLevels = 1;
    td.ArraySize = 1;
    td.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
    td.SampleDesc.Count = 1;
    td.Usage = D3D11_USAGE_DEFAULT;
    td.BindFlags = D3D11_BIND_DEPTH_STENCIL;

    ID3D11Texture2D *pDepthStencil = 0;
    HRESULT hr = device->CreateTexture2D(&td, 0, &pDepthStencil);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create depth stencil texture", hr);
    }
    m_psDepthStencil.reset(pDepthStencil);
    
    D3D11_DEPTH_STENCIL_VIEW_DESC vd;
    memset(&vd, 0, sizeof(vd));
    vd.Format = td.Format;
    vd.ViewDimension = D3D11_DSV_DIMENSION_TEXTURE2D;

    ID3D11DepthStencilView *pDepthStencilView = 0;
    hr = device->CreateDepthStencilView(pDepthStencil, &vd, &pDepthStencilView);
    if (FAILED(hr)) {
        throw Coercri::DXError("failed to create depth stencil view", hr);
    }
    m_psDepthStencilView.reset(pDepthStencilView);
}

void ShallowWaterEngine::loadGraphics()
{

    loadTerrainShading(0); // 0 is default shading or sand
	loadTerrainGrid(0); // 0 is default shading.
	D3D11_SAMPLER_DESC sd;
    memset(&sd, 0, sizeof(sd));
    sd.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
    sd.AddressU = sd.AddressV = sd.AddressW = D3D11_TEXTURE_ADDRESS_WRAP;
    sd.ComparisonFunc = D3D11_COMPARISON_NEVER;
    sd.MaxLOD = D3D11_FLOAT32_MAX;

    ID3D11SamplerState * sampler = 0;
    HRESULT hr = device->CreateSamplerState(&sd, &sampler);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create sampler state", hr);
    }
    m_psLinearSamplerState.reset(sampler);
	
	
}

void ShallowWaterEngine::loadTerrainGrid(int shading)
{

	std::ifstream str((initSetting.exePath + "/graphics/textures/grid.bmp").c_str(), std::ios::binary);
    boost::shared_ptr<Coercri::PixelArray> parr = Coercri::LoadBMP(str);
    
    CreateTextureWithMips(device,
                          parr->getWidth(),
                          parr->getHeight(),
                          reinterpret_cast<unsigned char *>(&(*parr)(0,0)),
                          m_psGridTexture,
                          m_psGridTextureView);

}
void ShallowWaterEngine::loadTerrainShading(int shading) 
{
	
	//TODO: Make a for loop to load all bmp files in the texture directory.
	std::string shading_filename = "";

	// With any changes to the hardcoded values, check isTerraingTextureColormap(int)
	switch (shading) {
	case 0:
		shading_filename = "Sand.bmp";
        break;
    case 1:
        shading_filename = "Grass.bmp";
        break;
	case 2:
        shading_filename = "Concrete.bmp";
        break;
	case 3:
        shading_filename = "Tile1.bmp";
        break;
	case 4:
        shading_filename = "Tile2.bmp";
        break;
	case 5:
        shading_filename = "Tile3.bmp";
        break;
	case 6:
        shading_filename = "Colormap.bmp";
        break;
	case 7:
        shading_filename = "Contours.bmp";
        break;
	case 8:
        shading_filename = "Custom.bmp";
        break;
    default:
        shading_filename = "Sand.bmp";
        break;
    }

    std::ifstream str((initSetting.exePath + "/graphics/textures/" + shading_filename).c_str(), std::ios::binary);
    boost::shared_ptr<Coercri::PixelArray> parr = Coercri::LoadBMP(str);
    
    CreateTextureWithMips(device,
                          parr->getWidth(),
                          parr->getHeight(),
                          reinterpret_cast<unsigned char *>(&(*parr)(0,0)),
                          m_psGrassTexture,
                          m_psGrassTextureView);

}



void ShallowWaterEngine::loadColormap(int shading)
{
	std::string shading_filename = "";

	switch (shading) {
	case PARULA:
		shading_filename = "parula.bmp";
        break;
    case JET:
        shading_filename = "jet.bmp";
        break;
    case ZEBRA:
		shading_filename = "zebra.bmp";
        break;
	case CUSTOM:
		shading_filename = "custom.bmp";
        break;
    default:
        shading_filename = "parula.bmp";
        break;
    }
    std::ifstream str((initSetting.exePath + "/graphics/colormaps/" + shading_filename).c_str(), std::ios::binary);
    boost::shared_ptr<Coercri::PixelArray> parr = Coercri::LoadBMP(str);
    
    CreateTextureWithMips(device,
                          parr->getWidth(),
                          parr->getHeight(),
                          reinterpret_cast<unsigned char *>(&(*parr)(0,0)),
                          m_psColormapTexture,
                          m_psColormapTextureView);
}

/*
void ShallowWaterEngine::createBlendState()
{
    D3D11_BLEND_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.RenderTarget[0].BlendEnable = true;
    bd.RenderTarget[0].SrcBlend = D3D11_BLEND_ONE;
    bd.RenderTarget[0].DestBlend = D3D11_BLEND_SRC_ALPHA;
    bd.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
    bd.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ZERO;
    bd.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_ZERO;
    bd.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
    bd.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;

    ID3D11BlendState *blend_state;
    HRESULT hr = device->CreateBlendState(&bd, &blend_state);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create blend state", hr);
    }
    m_psWaterBlendState.reset(blend_state);
}
*/

void ShallowWaterEngine::loadSkybox(int skybox_type)
{

	
	//TODO: Make a for loop to load all bmp files in the skybox directory.
	std::string shading_filename = "";

	switch (skybox_type) {
	case 0:
		shading_filename = "thick_clouds.dds";
        break;
    case 1:
        shading_filename = "clouds.dds";
        break;
	case 2:
        shading_filename = "ocean.dds";
        break;
	case 3:
        shading_filename = "sunset.dds";
        break;
	case 4:
        shading_filename = "desert.dds";
        break;
	case 5:
        shading_filename = "dark.dds";
        break;
	case 6:
        shading_filename = "light.dds";
        break;
    default:
        shading_filename = "clear_sky.dds";
        break;
    }

    // load the cube texture
    ID3D11ShaderResourceView * srv;
    HRESULT hr = D3DX11CreateShaderResourceViewFromFile(device, (initSetting.exePath + "/graphics/skyboxes/" + shading_filename).c_str(), 0, 0, &srv, 0);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to load the skybox", hr);
    }
    m_psSkyboxView.reset(srv);

    // create shaders
    Coercri::ComPtrWrapper<ID3DBlob> bytecode;
	std::string shallow_water_fx = initSetting.exePath + "/shaders/graphics.fx";
    CreateVertexShader(device, shallow_water_fx.c_str(), "SkyboxVertexShader", bytecode, m_psSkyboxVertexShader);
    CreatePixelShader(device, shallow_water_fx.c_str(), "SkyboxPixelShader", m_psSkyboxPixelShader);

    D3D11_INPUT_ELEMENT_DESC layout[] = {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
    };

    // create vertex buffer input layout
    ID3D11InputLayout *input_layout = 0;
    hr = device->CreateInputLayout(&layout[0],
                                   1,
                                   bytecode->GetBufferPointer(),
                                   bytecode->GetBufferSize(),
                                   &input_layout);
    if (FAILED(hr)) {
        throw Coercri::DXError("CreateInputLayout (skybox) failed", hr);
    }
    m_psSkyboxInputLayout.reset(input_layout);

    // create vertex buffer

    const float vertices[6*6*3] = {

        // north (+Y)
        -0.5f, 0.5f, 0.5f,
        0.5f,  0.5f, 0.5f,
        0.5f,  0.5f, -0.5f,
        0.5f,  0.5f, -0.5f,
        -0.5f, 0.5f, -0.5f,
        -0.5f, 0.5f, 0.5f,

        // west (-X)
        -0.5f, -0.5f, 0.5f,
        -0.5f, 0.5f,  0.5f,
        -0.5f, 0.5f,  -0.5f,
        -0.5f, 0.5f,  -0.5f,
        -0.5f, -0.5f, -0.5f,
        -0.5f, -0.5f, 0.5f,

        // south (-Y)
        0.5f, -0.5f, 0.5f,
        -0.5f, -0.5f, 0.5f,
        -0.5f, -0.5f, -0.5f,
        -0.5f, -0.5f, -0.5f,
        0.5f, -0.5f, -0.5f,
        0.5f, -0.5f, 0.5f,

        // east (+X)
        0.5f, 0.5f, 0.5f,
        0.5f, -0.5f, 0.5f,
        0.5f, -0.5f, -0.5f,
        0.5f, -0.5f, -0.5f,
        0.5f, 0.5f, -0.5f,
        0.5f, 0.5f, 0.5f,
        
        // up (+Z)
        -0.5f, -0.5f, 0.5f,
        0.5f, -0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
        -0.5f, 0.5f, 0.5f,
        -0.5f, -0.5f, 0.5f,


		// down (-Z)
        -0.5f, 0.5f, -0.5f,
        0.5f, 0.5f, -0.5f,
        0.5f, -0.5f, -0.5f,
        0.5f, -0.5f, -0.5f,
        -0.5f, -0.5f, -0.5f,
        -0.5f, 0.5f, -0.5f,
    };
   
    D3D11_BUFFER_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.ByteWidth = sizeof(vertices);
    bd.Usage = D3D11_USAGE_IMMUTABLE;
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;

    D3D11_SUBRESOURCE_DATA sd;
    memset(&sd, 0, sizeof(sd));
    sd.pSysMem = &vertices[0];

    ID3D11Buffer *pBuffer;
    hr = device->CreateBuffer(&bd, &sd, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create vertex buffer for skybox", hr);
    }
    m_psSkyboxVertexBuffer.reset(pBuffer);
}


// mouse picking

namespace {
    const int MOUSE_PICK_NUM_TRIANGLES = 12;

    struct LeftMouseVertex {
        float x, y;
    };
}

// converts pixel position into a world point on the water surface

// returns true if there was a hit; in this case world_xyz will contain the hit point.
// otherwise, returns false (and corrupts world_xyz).

bool ShallowWaterEngine::mousePick(int mouse_x, int mouse_y,
                                   float &world_x, float &world_y, float &world_z, float &depth, float &u, float &v)
{
	const float W = initSetting.width;;
	const float L = initSetting.length;
	const int nx = initSetting.nx;
	const int ny = initSetting.ny;

    const XMMATRIX translate_camera_pos( 1, 0, 0, camera_x,
                                         0, 1, 0, camera_y,
                                         0, 0, 1, camera_z,
                                         0, 0, 0, 1 );
    const XMMATRIX rotate_around_z( cos(yaw), sin(yaw), 0, 0,
                                    -sin(yaw), cos(yaw), 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, 1 );
    const XMMATRIX rotate_around_x( 1, 0, 0, 0,
                                    0, cos(pitch), -sin(pitch), 0,
                                    0, sin(pitch), cos(pitch), 0,
                                    0, 0, 0, 1);
    const XMMATRIX swap_y_and_z( 1, 0, 0, 0,
                                 0, 0, 1, 0,
                                 0, 1, 0, 0,
                                 0, 0, 0, 1 );

    // transpose because D3D treats position vectors as row vectors for some reason...
    const XMMATRIX eye_to_world = XMMatrixTranspose(translate_camera_pos * rotate_around_z * 
        rotate_around_x * swap_y_and_z);

    const float ndc_x = 2 * (mouse_x + 0.5f) / vp_width - 1;
    const float ndc_y = -2 * (mouse_y + 0.5f) / vp_height + 1;

    const float stepsize = 1;  // in metres

    bool found = false;
    XMFLOAT4 world_pos;


    // step along the ray in multiples of stepsize.
    // this is a bit inefficient, but it works well enough for what we want to use it for...

    // this assumes resetTimestep has previously been called.
    MapTexture m(*context, *m_psGetStatsStagingTexture4);

    for (float eye_z = 1; eye_z < 1000 && !found; eye_z += stepsize) {

        XMFLOAT4 eye_pos( eye_z * ndc_x / xpersp,
                          eye_z * ndc_y / ypersp,
                          eye_z,
                          1 );
        XMVECTOR eye_pos_vec = XMLoadFloat4(&eye_pos);
        XMVECTOR world_pos_vec = XMVector4Transform(eye_pos_vec, eye_to_world);
        XMStoreFloat4(&world_pos, world_pos_vec);
        
        // clip to the terrain region
        if (world_pos.x >= -W/2 && world_pos.x <= W/2 
        && world_pos.y >= 0 && world_pos.y <= L) {

            // convert world pos to a texture position in the staging texture.
            const float tex_x = (world_pos.x + W/2) / W;
            const float tex_y = (world_pos.y) / L;
            int ix = int(tex_x * (nx/4));
            int iy = int(tex_y * (ny/4));
            if (ix<0) ix=0;
            if (iy<0) iy=0;
            if (ix>nx/4-1) ix=nx/4-1;
            if (iy>ny/4-1) iy=ny/4-1;

            // lookup the "h" value at this point
            const char *row_ptr = reinterpret_cast<const char*>(m.msr.pData) + iy * m.msr.RowPitch;
            const float *col_ptr = reinterpret_cast<const float*>(row_ptr) + ix * 4;
            const float avg_h = (*col_ptr) / 16.0f;

            // lookup terrain height at this position
            const float B = GetTerrainHeight(world_pos.x, world_pos.y);
            const float w = B + avg_h;
            
            if (world_pos.z <= w) {
                world_x = world_pos.x;
                world_y = world_pos.y;
                world_z = w;
                depth = avg_h;

                u = CalcU(avg_h, col_ptr[2]/16);
                v = CalcU(avg_h, col_ptr[3]/16);
                
                found = true;
            }
        }
    }

    return found;
}

void ShallowWaterEngine::applyMouseShader(float world_x, float world_y, float dt)
{
    const int action = GetIntSetting("left_mouse_action");
    if (action == LM_RAISE_TERRAIN || action == LM_LOWER_TERRAIN) {
        raiseLowerTerrain(world_x, world_y, dt);
        return;
    }
    
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    const float brush_radius = GetSetting("left_mouse_radius");
    const float strength = GetSetting("left_mouse_strength") * dt;
    
    // transformation from world coords to texture index.
    const float world_to_idx_x_scale = (nx-1)/W;
    const float world_to_idx_x_bias = (nx+4)/2.0f;
    const float world_to_idx_y_scale = (ny-1)/L;
    const float world_to_idx_y_bias = 2.5f;

    // transformation from "radius 1" coords (in the vertex buf) to texture indices.
    const float x_scale = world_to_idx_x_scale * brush_radius;
    const float x_bias = world_to_idx_x_scale * world_x + world_to_idx_x_bias;
    const float y_scale = world_to_idx_y_scale * brush_radius;
    const float y_bias = world_to_idx_y_scale * world_y + world_to_idx_y_bias;

    // setup the const buffer
    LeftMouseConstBuffer cb;
    cb.scale_x = x_scale;
    cb.scale_y = y_scale;
    cb.bias_x = x_bias;
    cb.bias_y = y_bias;
    cb.two_over_nx_plus_four = 2.0f / (nx+4);
    cb.two_over_ny_plus_four = 2.0f / (ny+4);

    switch (action) {
    case LM_ADD_WATER:
        cb.disp_A = strength;
        cb.disp_B = 0;
        break;

    case LM_REMOVE_WATER:
        cb.disp_A = -strength;
        cb.disp_B = 0;
        break;

    case LM_STIR_WATER:
        cb.disp_A = -strength;
        cb.disp_B = 1.5f * strength;
        break;

    default:
        cb.disp_A = cb.disp_B = 0;
        break;
    }
    
    context->UpdateSubresource(m_psLeftMouseConstantBuffer.get(),
                               0, 0, &cb, 0, 0);

    // prepare to copy state texture to a temporary work area.
    // find the texture indices for the copy region.
    const int centre_ix = int(world_to_idx_x_scale * world_x + world_to_idx_x_bias);
    const int centre_iy = int(world_to_idx_y_scale * world_y + world_to_idx_y_bias);
    const int radius_x = int(world_to_idx_x_scale * brush_radius) + 2;  // safety margin
    const int radius_y = int(world_to_idx_y_scale * brush_radius) + 2;

    const int ix_min = std::max(0, std::min(nx+4, centre_ix - radius_x));   // inclusive
    const int ix_max = std::max(0, std::min(nx+4, centre_ix + radius_x + 1));  // exclusive
    const int iy_min = std::max(0, std::min(ny+4, centre_iy - radius_y));
    const int iy_max = std::max(0, std::min(ny+4, centre_iy + radius_y + 1));
    
    // do the copy
    D3D11_BOX src_box;
    src_box.left = ix_min;
    src_box.right = ix_max;
    src_box.top = iy_min;
    src_box.bottom = iy_max;
    src_box.front = 0;
    src_box.back = 1;
    context->CopySubresourceRegion(m_psSimTexture[4].get(),  // dest texture (xflux -- being used here as scratch space)
                                   0,  // subresource
                                   ix_min,  // dest x
                                   iy_min,  // dest y
                                   0,  // dest z
                                   m_psSimTexture[sim_idx].get(),  // src texture
                                   0,  // subresource
                                   &src_box);

    // now we can set up the pixel shader to input from xflux texture,
    // and write back to current state texture.
    
    context->ClearState();

    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->IASetInputLayout(m_psLeftMouseInputLayout.get());

    ID3D11Buffer *vert_buf = m_psLeftMouseVertexBuffer.get();
    const UINT stride = 8;
    const UINT offset = 0;
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);

    ID3D11Buffer *cst_buf = m_psLeftMouseConstantBuffer.get();
    context->VSSetShader(m_psLeftMouseVertexShader.get(), 0, 0);
    context->VSSetConstantBuffers(0, 1, &cst_buf);

    D3D11_VIEWPORT vp;
    memset(&vp, 0, sizeof(vp));
    vp.Width = float(nx+4);
    vp.Height = float(ny+4);
    vp.MaxDepth = 1;
    context->RSSetViewports(1, &vp);

    ID3D11ShaderResourceView *tex_views[] = { m_psSimTextureView[4].get(), m_psBottomTextureView.get() };
    context->PSSetShader(m_psLeftMousePixelShader.get(), 0, 0);
    context->PSSetConstantBuffers(0, 1, &cst_buf);
    context->PSSetShaderResources(0, 2, &tex_views[0]);

    ID3D11RenderTargetView * rtv = m_psSimRenderTargetView[sim_idx].get();
    context->OMSetRenderTargets(1, &rtv, 0);

    context->Draw(MOUSE_PICK_NUM_TRIANGLES * 3, 0);
}

void ShallowWaterEngine::raiseLowerTerrain(float world_x, float world_y, float dt)
{
    // Do this on CPU as it is the easiest way -- just want to get this code written
    // quickly, don't care about efficiency at this point :)

    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    const float dx = W / (nx-1);
    const float dy = L / (ny-1);    
    const float brush_radius = GetSetting("left_mouse_radius") * 2; // *2 to account for 'gaussian' nature of brush
    float strength = GetSetting("left_mouse_strength");
    const int action = GetIntSetting("left_mouse_action");    
    float displacement;
    if (action == LM_LOWER_TERRAIN) displacement = -strength*dt; else displacement = strength*dt;
    
    const float world_to_idx_x_scale = (nx-1)/W;
    const float world_to_idx_x_bias = (nx+4)/2.0f;
    const float world_to_idx_y_scale = (ny-1)/L;
    const float world_to_idx_y_bias = 2.5f;

    const int centre_ix = int(world_to_idx_x_scale * world_x + world_to_idx_x_bias);
    const int centre_iy = int(world_to_idx_y_scale * world_y + world_to_idx_y_bias);
    const int radius_x = int(world_to_idx_x_scale * brush_radius) + 3;  // safety margin
    const int radius_y = int(world_to_idx_y_scale * brush_radius) + 3;

    const int ix_min = std::max(0, std::min(nx+4, centre_ix - radius_x));   // inclusive
    const int ix_max = std::max(0, std::min(nx+4, centre_ix + radius_x + 1));  // exclusive
    const int iy_min = std::max(0, std::min(ny+4, centre_iy - radius_y));
    const int iy_max = std::max(0, std::min(ny+4, centre_iy + radius_y + 1));    

    D3D11_BOX src_box;
    src_box.left = ix_min;
    src_box.right = ix_max;
    src_box.top = iy_min;
    src_box.bottom = iy_max;
    src_box.front = 0;
    src_box.back = 1;
    context->CopySubresourceRegion(m_psFullSizeStagingTexture.get(),
                                   0,  // subresource
                                   ix_min,  // dest x
                                   iy_min,  // dest y
                                   0,  // dest z
                                   m_psSimTexture[sim_idx].get(),  // src texture -- current state
                                   0,  // subresource
                                   &src_box);

    {
        MapTexture m(*context, *m_psFullSizeStagingTexture);

        for (int iy = iy_min; iy < iy_max; ++iy) {

            // have to be careful about the exact definitions of BX, BY, BA (see terrain_heightfield.cpp)
            // (Water level w must be brought up in line with BA)
            
            const float wy_plus = (iy + 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;
            const float wy_minus = (iy - 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;

            char * row_ptr = static_cast<char*>(m.msr.pData) + m.msr.RowPitch * iy;
            
            for (int ix = ix_min; ix < ix_max; ++ix) {

                const float wx_plus = (ix + 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;
                const float wx_minus = (ix - 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;

                const float plus_plus = GaussianBrush(wx_plus, wy_plus, brush_radius) * displacement;
                const float plus_minus = GaussianBrush(wx_plus, wy_minus, brush_radius) * displacement;
                const float minus_plus = GaussianBrush(wx_minus, wy_plus, brush_radius) * displacement;
                const float minus_minus = GaussianBrush(wx_minus, wy_minus, brush_radius) * displacement;

                float * ptr = reinterpret_cast<float*>(row_ptr) + 4 * ix;

                *ptr += (plus_plus + plus_minus + minus_plus + minus_minus) / 4;
            }
        }
    }

    context->CopySubresourceRegion(m_psSimTexture[sim_idx].get(),  // dest texture
                                   0,  // subresource
                                   ix_min,  // dest x
                                   iy_min,  // dest y
                                   0,  // dest z
                                   m_psFullSizeStagingTexture.get(),  // src texture
                                   0,  // subresource
                                   &src_box);    

    // Now we are going to update the in-memory copy of BottomTexture and
    // upload it to the GPU via UpdateSubresource.

    for (int iy = iy_min; iy < iy_max; ++iy) {

        const float wy_plus = (iy + 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;
        const float wy_minus = (iy - 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;

        BottomEntry * row_ptr = &g_bottom[iy * (nx+4)];
            
        for (int ix = ix_min; ix < ix_max; ++ix) {
            
            const float wx_plus = (ix + 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;
            const float wx_minus = (ix - 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;
            
            const float plus_plus = GaussianBrush(wx_plus, wy_plus, brush_radius) * displacement;
            const float plus_minus = GaussianBrush(wx_plus, wy_minus, brush_radius) * displacement;
            const float minus_plus = GaussianBrush(wx_minus, wy_plus, brush_radius) * displacement;
            const float minus_minus = GaussianBrush(wx_minus, wy_minus, brush_radius) * displacement;

            BottomEntry * ptr = row_ptr + ix;

            ptr->BX += (plus_plus + plus_minus)/2;
            ptr->BY += (plus_plus + minus_plus)/2;
            ptr->BA += (plus_plus + plus_minus + minus_plus + minus_minus)/4;
        }
    }

    context->UpdateSubresource(m_psBottomTexture.get(),  // dest texture
                               0,  // subresource
                               &src_box,  // dest box
                               &g_bottom[iy_min * (nx+4) + ix_min],  // src data
                               (nx+4) * 3 * sizeof(float),   // row pitch
                               0);  // depth pitch (unused)

    // more crappy copy paste code...
    
    for (int iy = iy_min; iy < iy_max; ++iy) {

        const float wy_plus = (iy + 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;
        const float wy_minus = (iy - 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;
        
        TerrainEntry * row_ptr = &g_terrain_heightfield[iy * (nx+4)];
            
        for (int ix = ix_min; ix < ix_max; ++ix) {

            const float wx_plus = (ix + 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;
            const float wx_minus = (ix - 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;

            const float plus_plus = GaussianBrush(wx_plus, wy_plus, brush_radius) * displacement;
            const float plus_minus = GaussianBrush(wx_plus, wy_minus, brush_radius) * displacement;
            const float minus_plus = GaussianBrush(wx_minus, wy_plus, brush_radius) * displacement;
            const float minus_minus = GaussianBrush(wx_minus, wy_minus, brush_radius) * displacement;
            
            TerrainEntry * ptr = row_ptr + ix;

            ptr->B += (plus_plus + plus_minus + minus_plus + minus_minus)/4;
            ptr->dBdx += (plus_plus + plus_minus - minus_plus - minus_minus) / (2*dx);
            ptr->dBdy += (plus_plus + minus_plus - plus_minus - minus_minus) / (2*dy);
        }
    }  

    context->UpdateSubresource(m_psTerrainTexture.get(),  // dest texture
                               0,  // subresource
                               &src_box,  // dest box
                               &g_terrain_heightfield[iy_min * (nx+4) + ix_min],  // src data
                               (nx+4) * 3 * sizeof(float),   // row pitch
                               0);  // depth pitch (unused)    

	createTridiagonalCoefTextures();
}


void ShallowWaterEngine::shiftTerrainSlider(float shift_up)
{
    // Do this on CPU as it is the easiest way -- just want to get this code written
    // quickly, don't care about efficiency at this point :)
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
 
/////////////

    context->CopySubresourceRegion(m_psFullSizeStagingTexture.get(),
                                   0,  // subresource
                                   0,  // dest x
                                   0,  // dest y
                                   0,  // dest z
                                   m_psSimTexture[sim_idx].get(),  // src texture -- current state
                                   0,  // subresource
                                   0);
	MapTexture m(*context, *m_psFullSizeStagingTexture);


    // We are going to update the in-memory copy of BottomTexture and
    // upload it to the GPU via UpdateSubresource.

    for (int iy = 0; iy < ny + 4; ++iy) {

        BottomEntry * row_ptr_b = &g_bottom[iy * (nx+4)];
        TerrainEntry * row_ptr_t = &g_terrain_heightfield[iy * (nx+4)];
		char * row_ptr_w = static_cast<char*>(m.msr.pData) + m.msr.RowPitch * iy;
        for (int ix = 0; ix < nx + 4; ++ix) {
            
			BottomEntry * ptr_b = row_ptr_b + ix;
			float * ptr_w = reinterpret_cast<float*>(row_ptr_w) + 4 * ix;
			TerrainEntry * ptr_t = row_ptr_t + ix;
			
			if ((*ptr_w) - ptr_b->BA < 0.01){
				if (ptr_b->BA < initSetting.stillWaterElevation){
					*ptr_w = initSetting.stillWaterElevation;
				} else{
					*ptr_w += shift_up;
				}
			}
			
			ptr_b->BX += shift_up;
			ptr_b->BY += shift_up;
			ptr_b->BA += shift_up;

			
			ptr_t->B += shift_up;


        }
    }

    context->UpdateSubresource(m_psBottomTexture.get(),  // dest texture
                               0,  // subresource
                               0,  // dest box
                               &g_bottom[0],  // src data
                               (nx+4) * 3 * sizeof(float),   // row pitch
                               0);  // depth pitch (unused) 

    context->UpdateSubresource(m_psTerrainTexture.get(),  // dest texture
                               0,  // subresource
                               0,  // dest box
                               &g_terrain_heightfield[0],  // src data
                               (nx+4) * 3 * sizeof(float),   // row pitch
                               0);  // depth pitch (unused)

    context->CopySubresourceRegion(m_psSimTexture[sim_idx].get(),  // dest texture
                               0,  // subresource
                               0,  // dest x
                               0,  // dest y
                               0,  // dest z
                               m_psFullSizeStagingTexture.get(),  // src texture
                               0,  // subresource
                               0);  

	createTridiagonalCoefTextures();
}



void ShallowWaterEngine::setupMousePicking()
{
    // create shaders
    Coercri::ComPtrWrapper<ID3DBlob> bytecode;
    CreateVertexShader(device, (initSetting.exePath + "/shaders/left_mouse.hlsl").c_str(), "LeftMouseVertexShader", bytecode, m_psLeftMouseVertexShader);
    CreatePixelShader(device, (initSetting.exePath + "/shaders/left_mouse.hlsl").c_str(), "LeftMousePixelShader", m_psLeftMousePixelShader);
        
    // create input layout
    D3D11_INPUT_ELEMENT_DESC layout[] = {
        { "POSITION", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
    };
    
    ID3D11InputLayout *input_layout = 0;
    HRESULT hr = device->CreateInputLayout(&layout[0],
                                           1,
                                           bytecode->GetBufferPointer(),
                                           bytecode->GetBufferSize(),
                                           &input_layout);
    if (FAILED(hr)) {
        throw Coercri::DXError("CreateInputLayout failed", hr);
    }
    m_psLeftMouseInputLayout.reset(input_layout);

    // create a vertex buffer for a circle of radius 1.
    std::vector<LeftMouseVertex> buf(MOUSE_PICK_NUM_TRIANGLES * 3);
    const float d_theta = 2*PI/MOUSE_PICK_NUM_TRIANGLES;
    for (int i = 0; i < MOUSE_PICK_NUM_TRIANGLES; ++i) {
        const float theta = d_theta * i;
        const float theta2 = d_theta * (i+1);

        buf[3*i].x = 0;
        buf[3*i].y = 0;

        buf[3*i+1].x = cos(theta);
        buf[3*i+1].y = sin(theta);

        buf[3*i+2].x = cos(theta2);
        buf[3*i+2].y = sin(theta2);
    }

    D3D11_BUFFER_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.ByteWidth = sizeof(LeftMouseVertex) * MOUSE_PICK_NUM_TRIANGLES * 3;
    bd.Usage = D3D11_USAGE_IMMUTABLE;
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;

    D3D11_SUBRESOURCE_DATA sd;
    memset(&sd, 0, sizeof(sd));
    sd.pSysMem = &buf[0];

    ID3D11Buffer *pBuffer = 0;
    hr = device->CreateBuffer(&bd, &sd, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create left mouse vertex buffer", hr);
    }
    m_psLeftMouseVertexBuffer.reset(pBuffer);

    // create the constant buffer.
    bd.ByteWidth = RoundUpTo16(sizeof(LeftMouseConstBuffer));
    bd.Usage = D3D11_USAGE_DEFAULT;
    bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psLeftMouseConstantBuffer.reset(pBuffer);
}

float ShallowWaterEngine::getWaterHeight(float world_x, float world_y)
{
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    MapTexture m(*context, *m_psGetStatsStagingTexture4);
    const float tex_x = (world_x + W/2) / W;
    const float tex_y = (world_y) / L;
    int ix = int(tex_x * (nx/4));
    int iy = int(tex_y * (ny/4));
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    const char *row_ptr = reinterpret_cast<const char*>(m.msr.pData) + iy * m.msr.RowPitch;
    const float *col_ptr = reinterpret_cast<const float*>(row_ptr) + ix * 4;
    const float avg_h = (*col_ptr) / 16.0f;
    const float B = GetTerrainHeight(world_x, world_y);
    return B + avg_h;
}
