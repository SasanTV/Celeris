/*
 * FILE:
 *   engine.hpp
 *
 * PURPOSE:
 *   "Engine" class for Celeris. Manages all the
 *   Direct3D objects and contains top-level rendering routines.
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

#ifndef ENGINE_HPP
#define ENGINE_HPP

#include "settings.hpp"

#include "coercri/dx11/core/com_ptr_wrapper.hpp"

#include <d3d11.h>
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif



class ShallowWaterEngine {
public:

    ShallowWaterEngine(ID3D11Device *device_, ID3D11DeviceContext *context_);

    void newTerrainSettings();  // call if any terrain settings change.
    void remesh(ResetType rt);   // call if mesh size changes.

	void dumpBathymetryToFile(); // writes bathymetry to file.
	void dumpInundationToFile(ID3D11DeviceContext *context, ID3D11Texture2D *staging, ID3D11Texture2D *tex); // writes inundation to file.

	// call following camera move or window resize  
	void moveCamera(float cam_x, float cam_y, float cam_z, float yaw, float pitch, int vp_width, int vp_height);

    // update the water
    void timestep();

	// Deals with stuff that happens after a timestep, but needs to while software is on pause.
	void afterTimestep();

    // set timestep to the given dt (multiplied by time_acceleration),
    // or to safety_factor * CFL-timestep,
    // whichever is smaller.
    void resetTimestep(float dt, float elapsed_time);
    
    // render a frame to the given render target
    void render(ID3D11RenderTargetView *rtv);


    // mouse picking
    bool mousePick(int mouse_x, int mouse_y,
                   float &world_x, float &world_y, float &world_z, float &depth, float &u, float &v);
    void applyMouseShader(float world_x, float world_y, float dt);

    // get water height (h + B) at the given point
    float getWaterHeight(float world_x, float world_y);

	void createBoundaryShaders();
	void createBoundaryShadersInit();
    
    
private:
    void createShadersAndInputLayout();
	
    void createMeshBuffers();
    void createSimBuffers();
	void addSolitaryWave();

	void westUniformTimeSeries(float & temp_eta,float & temp_hu,float & temp_hv,float total_time);
	void eastUniformTimeSeries(float & temp_eta,float & temp_hu,float & temp_hv,float total_time);
	void southUniformTimeSeries(float & temp_eta,float & temp_hu,float & temp_hv,float total_time);
	void northUniformTimeSeries(float & temp_eta,float & temp_hu,float & temp_hv,float total_time);

	void westIrregularBoundary();
	void eastIrregularBoundary();
	void southIrregularBoundary();
	void northIrregularBoundary();
    
    void createTerrainTexture();
    void fillTerrainTexture();
    void fillTerrainTextureLite();
    void createSimTextures(ResetType reset_type);
	void createTridiagonalCoefTextures();
    void createConstantBuffers();
	std::string ShallowWaterEngine::generateSpectrumWaves(IrregularWaveSpectrumSetting inSpectrum, std::string boundarySide);
	void fillIrregularWavesDataConstantBuffer(std::string boundarySide);
	void fillUniformTimeSeriesMainMemoryBuffer();
    void fillConstantBuffers();
    void createDepthStencil(int w, int h);
    void loadGraphics();
	void loadColormap(int shading);
	void loadTerrainShading(int shading);
	void loadTerrainGrid(int shading);
	
    void loadSkybox(int skybox_type);

    void setupMousePicking();
    void raiseLowerTerrain(float wx, float wy, float dt);
	void shiftTerrainSlider(float shift_up);

	int  st_getCountOfMatrices(int n);

	void setMinMaxSetting();
  
private:
    ID3D11Device *device;
    ID3D11DeviceContext *context;

    // vertex & pixel shaders
    Coercri::ComPtrWrapper<ID3D11VertexShader> m_psTerrainVertexShader, m_psWaterVertexShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psTerrainPixelShader, m_psWaterPixelShader;
    Coercri::ComPtrWrapper<ID3D11VertexShader> m_psSimVertexShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psSimPixelShader[3], m_psGetStatsPixelShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> st_psCopyFromXxandXyPixelShader, st_psCyclicReduceDx_PixelShader, st_psCyclicSubstituteXx_PixelShader;
	Coercri::ComPtrWrapper<ID3D11PixelShader> st_psCyclicReduceDy_PixelShader, st_psCyclicSubstituteXy_PixelShader;
	Coercri::ComPtrWrapper<ID3D11PixelShader> m_psLaxWendroffPixelShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psBoundaryPixelShader[4];
	Coercri::ComPtrWrapper<ID3D11PixelShader> m_psAddSolitaryWaveShader;
	Coercri::ComPtrWrapper<ID3D11PixelShader> m_psColumnSumReduceShader, m_psRowSumReduceShader;

    // input layouts
    Coercri::ComPtrWrapper<ID3D11InputLayout> m_psInputLayout, m_psSimInputLayout;
    
    // vertex and index buffers (for the triangle meshes).
    // can be used for both terrain & water.
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psMeshVertexBuffer;
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psMeshIndexBuffer;

    // vertex buffers for the simulation render-to-texture.
    // (one for each pass.)
    // TODO: Might be better to have one large vertex buffer, with offsets,
    // rather than lots of little ones like this.
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psSimVertexBuffer11, m_psSimVertexBuffer10, m_psSimVertexBuffer00;
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psGetStatsVertexBuffer;
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psBoundaryVertexBuffer[4];
	enum{
		ST_MAX_CR_MATRIX_NUM_X=15, // this is equal to log(MAX_X_GRID)/log(2)
		ST_MAX_CR_MATRIX_NUM_Y=15 // this is equal to log(MAX_Y_GRID)/log(2)
	};
	Coercri::ComPtrWrapper<ID3D11Buffer> st_psCyclicReductionX_VertexBuffer[ST_MAX_CR_MATRIX_NUM_X];
	Coercri::ComPtrWrapper<ID3D11Buffer> st_psCyclicReductionY_VertexBuffer[ST_MAX_CR_MATRIX_NUM_Y];

    // constant buffers
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psConstantBuffer, m_psSimConstantBuffer,
										 m_psBoundaryConstantBuffer, m_psTimeIntegrationConstantBuffer,
										 m_psIrregularWavesDataConstantBuffer, m_psIrregularWavesColumnRowConstBuffer,
										 m_psSolitaryWaveConstantBuffer;
    
    // textures:
    //  -- terrain (contains B, dB/dx, dB/dy)
    //  -- bottom (contains BY, BX, BA; used for simulation)
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psTerrainTexture;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psTerrainTextureView;
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psBottomTexture;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psBottomTextureView;
	Coercri::ComPtrWrapper<ID3D11Texture2D> m_psAuxiliary1Texture, m_psAuxiliary2Texture;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psAuxiliary1TextureView, m_psAuxiliary2TextureView;
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> m_psAuxiliary1RenderTargetView, m_psAuxiliary2RenderTargetView;
	

    // simulation textures:
    // [sim_idx] = state   0 or 1
    // [1-sim_idx] = output state / H   1 or 0
    // [2] = U
    // [3] = V
    // [4] = XFLUX
    // [5] = YFLUX
	// [6] = normals 
	// [7] = old_gradients
	// [8] = old_old_gradients
	// [9] = scratch_gradients
	// [10] = predicted_gradients
	// [11] = F_G_star_old_gradients
	// [12] = F_G_star_old_old_gradients
	// [13] = F_G_star_scratch_gradients
	// [14] = F_G_star_predicted_gradients
	// [15] = current state (used as old state for prediction and correction)
#define NUM_OF_TEXTURES 16
#define MAX_NUM_OF_IRREGULAR_WAVES 1000
#define CURRENT_STATE 15

	Coercri::ComPtrWrapper<ID3D11Texture2D> m_psSimTexture[NUM_OF_TEXTURES], m_psGetStatsTexture;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psSimTextureView[NUM_OF_TEXTURES];
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> m_psSimRenderTargetView[NUM_OF_TEXTURES], m_psGetStatsRenderTargetView;
	
	int sim_idx;
	bool bootstrap_needed;
	bool firstTimeStep;
	bool secondTimeStep;
	TimeIntegrationScheme timeScheme;

	int timestep_count;
	int correction_steps_count;

	float realworld_timestep;
	float myDelay;
	int old_index, old_old_index, scratch_index, predicted_index;

	int F_G_star_old_index, F_G_star_old_old_index, F_G_star_scratch_index, F_G_star_predicted_index;


    // tridiagonal textures:
	// ST_:For each level of reduction in Cyclic Reducation, one m_psTriDigTexture is used. 15 levels is too many! But let's stay on the safe side for now.

    Coercri::ComPtrWrapper<ID3D11Texture2D> st_psTriDigTexture_ABCx[ST_MAX_CR_MATRIX_NUM_X];
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> st_psTriDigTextureView_ABCx[ST_MAX_CR_MATRIX_NUM_X];
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> st_psTriDigRenderTargetView_ABCx[ST_MAX_CR_MATRIX_NUM_X];

	Coercri::ComPtrWrapper<ID3D11Texture2D> st_psTriDigTexture_D1x[ST_MAX_CR_MATRIX_NUM_X]; //D1 is the right hand side of the tridiagonal equation. D1x is in x direction (corresponding to hu)
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> st_psTriDigTextureView_D1x[ST_MAX_CR_MATRIX_NUM_X];
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> st_psTriDigRenderTargetView_D1x[ST_MAX_CR_MATRIX_NUM_X];

	Coercri::ComPtrWrapper<ID3D11Texture2D> st_psTriDigTexture_D1xCopy[ST_MAX_CR_MATRIX_NUM_X]; //D1 is the right hand side of the tridiagonal equation. D1x is in x direction (corresponding to hu)
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> st_psTriDigTextureView_D1xCopy[ST_MAX_CR_MATRIX_NUM_X];
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> st_psTriDigRenderTargetView_D1xCopy[ST_MAX_CR_MATRIX_NUM_X];

	Coercri::ComPtrWrapper<ID3D11Texture2D> st_psTriDigTexture_ABCy[ST_MAX_CR_MATRIX_NUM_Y];
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> st_psTriDigTextureView_ABCy[ST_MAX_CR_MATRIX_NUM_Y];
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> st_psTriDigRenderTargetView_ABCy[ST_MAX_CR_MATRIX_NUM_Y];

	Coercri::ComPtrWrapper<ID3D11Texture2D> st_psTriDigTexture_D1y[ST_MAX_CR_MATRIX_NUM_Y]; //D1 is the right hand side of the tridiagonal equation. D1x is in x direction (corresponding to hu)
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> st_psTriDigTextureView_D1y[ST_MAX_CR_MATRIX_NUM_Y];
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> st_psTriDigRenderTargetView_D1y[ST_MAX_CR_MATRIX_NUM_Y];

	Coercri::ComPtrWrapper<ID3D11Texture2D> st_psTriDigTexture_D1yCopy[ST_MAX_CR_MATRIX_NUM_Y]; //D1 is the right hand side of the tridiagonal equation. D1x is in x direction (corresponding to hu)
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> st_psTriDigTextureView_D1yCopy[ST_MAX_CR_MATRIX_NUM_Y];
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> st_psTriDigRenderTargetView_D1yCopy[ST_MAX_CR_MATRIX_NUM_Y];

	Coercri::ComPtrWrapper<ID3D11Texture2D> st_psRHSTexture;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> st_psRHSTextureView;
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> st_psRHSRenderTargetView;

    // staging textures
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psFullSizeStagingTexture;
	Coercri::ComPtrWrapper<ID3D11Texture2D> m_psGetStatsStagingTexture4, m_psGetStatsStagingTexture1;
	Coercri::ComPtrWrapper<ID3D11Texture2D> st_psStagingTextureX[ST_MAX_CR_MATRIX_NUM_X];
	Coercri::ComPtrWrapper<ID3D11Texture2D> st_psStagingTextureY[ST_MAX_CR_MATRIX_NUM_Y];

    // graphical textures
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psGrassTexture, m_psColormapTexture, m_psGridTexture ;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psGrassTextureView, m_psColormapTextureView, m_psGridTextureView;
    Coercri::ComPtrWrapper<ID3D11SamplerState> m_psLinearSamplerState;

    // skybox
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psSkyboxView;
    Coercri::ComPtrWrapper<ID3D11VertexShader> m_psSkyboxVertexShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psSkyboxPixelShader;
    Coercri::ComPtrWrapper<ID3D11InputLayout> m_psSkyboxInputLayout;
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psSkyboxVertexBuffer;
    
    // depth buffer stuff
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psDepthStencil;
    Coercri::ComPtrWrapper<ID3D11DepthStencilView> m_psDepthStencilView;

    // mouse picking
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psLeftMouseConstantBuffer, m_psLeftMouseVertexBuffer;
    Coercri::ComPtrWrapper<ID3D11VertexShader> m_psLeftMouseVertexShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psLeftMousePixelShader;
    Coercri::ComPtrWrapper<ID3D11InputLayout> m_psLeftMouseInputLayout;
    
    // current timestep
    float current_timestep;
	float old_timestep;
	float old_old_timestep;
    float total_time;
	int myTimestep_count;

    // viewport/camera state
    int vp_width, vp_height;
    float pitch, yaw, camera_x, camera_y, camera_z;
    float xpersp, ypersp;  // perspective coefficients.
	
public:
	float near_plane_dist, far_plane_dist;
};

#endif
