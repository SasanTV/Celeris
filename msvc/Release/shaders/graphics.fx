/* -*-c++-*-
 *
 * FILE:
 *   graphics.fx
 *
 * PURPOSE:
 *   Graphical shaders for shallow water demo
 *   (land, water, and sky)
 *
 * AUTHOR:
 *   Stephen Thompson <stephen@solarflare.org.uk>
 *
 * CREATED:
 *   24-Oct-2011
 *
 * COPYRIGHT:
 *   Copyright (C) 2012, Stephen Thompson.
 *   Sasan Tavakkol <tavakkol@usc.edu> or <sasan.tavakkol@yahoo.com>
 *
 * CREATED:
 *   27-Oct-2016
 *-------------------------------------------------------------------------------------------
 * COPYRIGHT:
 *   Copyright (C) 2012, Stephen Thompson.
 *
 *   This file is part of Stephen Thompson's Shallow Water Demo.
 *
 *   The Shallow Water Demo is free software: you can redistribute it
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
 *
 * COPYRIGHT:
 *   Copyright (C) 2016,  Sasan Tavakkol
 *
 *   This file is part of Celeris Software.
 *
 *   Celeris is free software: you can redistribute it
 *   and/or modify it under the terms of the GNU General Public
 *   License as published by the Free Software Foundation, either
 *   version 3 of the License, or (at your option) any later version.
 *
 *   Celeris is distributed in the hope that it will be
 *   useful, but WITHOUT ANY WARRANTY; without even the implied
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *   See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Celeris. If not, see
 *   <http://www.gnu.org/licenses/>.
 *

 *
 *
 */


// Textures

// heightfield texture: 3 x 32-bit float
// r = B (height)
// g = dB/dx
// b = dB/dy

Texture2D<float3> txHeightfield : register( t0 );

// grass image texture.
Texture2D<float4> txGrass : register( t0 );

// a linear sampler
SamplerState samLinear : register( s0 );

// skybox cube texture
TextureCube txSkybox : register( t1 );


// colormap image texture.
Texture2D<float4> txColormap : register( t2 );

// grid image texture.
Texture2D<float4> txGrid : register( t3 );



// water texture: 4 x 32-bit float
// r = w
// g = hu
// b = hv
// a = (unused)

Texture2D<float4> txWater : register( t1 );

// normal texture

// r = nX
// g = nY
// b = nZ
// a = unused

Texture2D<float4> txNormal : register( t2 );

// r = maxFlowDepth
// g = unused
// b = breaking
Texture2D<float4> txAuxiliary1 : register( t4 );



// Constant Buffers

cbuffer MyConstBuffer : register( b0 )
{
	
    // used by renderer
    row_major float4x4 tex_to_clip;    // transforms (x_tex_idx, y_tex_idx, z_world, 1) to clip space
    float3 light_dir;   // in world space
    float ambient;
    float3 eye_mult, eye_trans;

    float2 terrain_tex_scale;

    // skybox matrix
    row_major float3x3 skybox_mtx;  // transforms world to clip space
    float4 pack;

    // refraction stuff
    float2 world_mult_ie_dx_dy, world_trans;
    float2 world_to_grass_tex;
    float2 grass_tex_of_origin;
	
    // more settings
    float fresnel_coeff, fresnel_exponent;
    float specular_intensity, specular_exponent;
    float refractive_index;
    float attenuation_1, attenuation_2;
    
    int nx_plus_1, ny_plus_1;
//	float dx, dy;
	
    
    float3 deep_col;
	
	float zScale;
	float seaLevel;
	
	//  0 is Photorealisitic,  1 is eta, 2 is u, 3 is v, 4 is speed, 5 is vorticity, 6 is max elevation (eta) .
	int water_shading;
	// 0 is texture, 1 and anything else is colormap.
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


// Structure definitions

struct VS_INPUT {
    int2 pos : INT_POSITION;       // texture indices
};

struct TERRAIN_PS_INPUT {
    float B : BOTTOM_HEIGHT;
	float4 pos : SV_POSITION;      // clip space
    float3 normal : NORMAL;        // world space
    float2 tex_coord : TEXCOORD;
	float4 auxiliary : INUNDATED;
};

struct WATER_PS_INPUT {
    float4 pos : SV_POSITION;      // clip space
    float2 tex_coord : TEXCOORD;
    float3 normal : NORMAL;        // world space
    float3 eye : CAMERA_DIRECTION; // world space
    float water_depth : WATER_DEPTH;   // water depth (h) as float
	float eta : SURFACE_ELEVATION;   // water depth (h) as float
	float u : U_VELOCITY;   // u velocity
	float v : V_VELOCITY;   // v velocity
	float vorticity : VORTICITY;
    float2 world_pos : WORLD_XY_POS;   // world space
    float3 terrain_normal : TERRAIN_NORMAL;
	float B : BOTTOM_HEIGHT;
	float4 auxiliary : AUXILIARY;
};

struct SKYBOX_PS_INPUT {
    float4 pos : SV_POSITION;
    float3 tex : TEXCOORD;
};


// lighting function for terrain.

float3 TerrainColour(float3 tex_colour, float3 terrain_normal)
{
    float light = saturate(dot(light_dir, terrain_normal)) + ambient;
    return light * tex_colour;
}


// Vertex Shader for Terrain

TERRAIN_PS_INPUT TerrainVertexShader( VS_INPUT input )
{
    TERRAIN_PS_INPUT output;
	
	const int3 idx = int3(input.pos.x, input.pos.y, 0);
    const float B = txHeightfield.Load(idx).r;
    output.B = B;
	output.auxiliary = float4(txAuxiliary1.Load(idx).r, 0, txAuxiliary1.Load(idx).b, 0);
	
    // lookup texture values at the input point
    const int3 tpos = int3(input.pos.x, input.pos.y, 0);
    const float3 tex_value = txHeightfield.Load( tpos );

    // compute position in clip space
    const float4 pos_in = float4( input.pos.x, input.pos.y, zScale * tex_value.r, 1 );
    output.pos = mul( tex_to_clip, pos_in );

    const float3 norm = float3(-tex_value.g, -tex_value.b, 1.0f);
    output.normal = normalize(norm);

    // texture coords.
    // TODO: the tex coords might be better input in the vertex buffer rather than calculated here.
    output.tex_coord = terrain_tex_scale * float2(input.pos.x - 2, input.pos.y - 2);
    
    return output;
}


// Pixel Shader for Terrain

float4 TerrainPixelShader( TERRAIN_PS_INPUT input ) : SV_Target
{
    // tex lookup
    float3 tex_colour;
	
	if (terrain_shading == 0){ // 0 is texture, anything else is colormap.
		tex_colour = txGrass.Sample( samLinear, input.tex_coord ).rgb;
	} else {
		const float temp_B = input.B;
		float colormap_u = 0;
		if (temp_B <= 0){
			colormap_u = 0.5 * (input.B - terrain_colormap_min) / (0.0001f - terrain_colormap_min);  
		} else {
			colormap_u = 0.5 + 0.5 * (input.B - 0) / (terrain_colormap_max - (-0.0001));  
		}
		tex_colour = txGrass.Sample(samLinear, float2(clamp(colormap_u,0.01,.99),0.5f)).rgb;   
	}
	
	if(input.auxiliary.r > drylandDepthOfInundation && drylandDepthOfInundation > 0){
		tex_colour = float3(0.9,0,1);
	}

	if(isGridOn){
		tex_colour += (1-txGrid.Sample( samLinear, input.tex_coord ).rgb);
	}

    // Lighting calculation
    return float4(TerrainColour(tex_colour, input.normal), 1);
}



// Vertex Shader for Water

WATER_PS_INPUT WaterVertexShader( VS_INPUT input )
{
    WATER_PS_INPUT output;

	const float h_min = sqrt_sqrt_epsilon; 

    // lookup the terrain level
    const int3 idx = int3(input.pos.x, input.pos.y, 0);
	
    const float3 ground_tex = txHeightfield.Load(idx);
    const float B =  ground_tex.r;
	output.B = B;
	output.auxiliary = float4(txAuxiliary1.Load(idx).r, txAuxiliary1.Load(idx).g, txAuxiliary1.Load(idx).b, txAuxiliary1.Load(idx).a);

    const float3 ground_normal = normalize(float3(-ground_tex.g, -ground_tex.b, 1));
    
    // lookup the water level
    float3 w =  txWater.Load(idx).rgb;

    // simple way to ensure there is no "gap" at the edge of the box
    if (idx.x==2 || idx.y==2 || idx.x==nx_plus_1 || idx.y==ny_plus_1) w.r = B;

    // calculate water depth, this is used in the pixel shader
	float h = (w.r - B);
    output.water_depth = zScale * h; // water depth is used for refraction in pixel shader. therefore we scale it.
	output.eta =  w.r - seaLevel;
	if (h > h_min){
		output.u =  w.g/ h; 
		output.v =  w.b/ h;
	} else {
		//This can be hu/(sqrt_sqrt_h). check, which one is better for visualization. 
		output.u =  0; 
		output.v =  0;
	}

	if (water_shading == 5){ // if colormapping variable is vorticity
		const float3 w_right =  txWater.Load(idx + int3(1, 0, 0)).rgb;
		const float3 w_up =  txWater.Load(idx + int3(0, 1, 0)).rgb;
		const float3 ground_tex_right = txHeightfield.Load(idx + int3(1, 0, 0));
		const float3 ground_tex_up = txHeightfield.Load(idx + int3(0, 1, 0));
		const float B_right =  ground_tex_right.r;
		const float B_up =  ground_tex_up.r;
		const float h_right = (w_right.r - B_right);
		const float h_up = (w_up.r - B_up);
		float v_right = 0;
		float u_up = 0;
		if (h_right > h_min){
			v_right =  w_right.b / h_right; 
		}
		if (h_up > h_min){
			u_up =  w_up.g / h_up; 
		}
		output.vorticity = (v_right - output.v)/ world_mult_ie_dx_dy.x - (u_up - output.u) / world_mult_ie_dx_dy.y;
	} else {
		output.vorticity = 0;
	}


    // compute position in clip space
    const float dmin = sqrt_sqrt_epsilon;
    const float vert_bias = (output.water_depth < dmin ? 2*(output.water_depth - dmin) : 0);
    const float4 pos_in = float4( input.pos.x, input.pos.y, zScale * w.r + vert_bias, 1 );
    output.pos = mul( tex_to_clip, pos_in );
	
	output.tex_coord = terrain_tex_scale * float2(input.pos.x - 2, input.pos.y - 2);

    // now compute the normal
    output.normal = txNormal.Load(idx).rgb;

    // now compute the eye vector (from vertex towards camera)
    output.eye = normalize(eye_mult * pos_in.xyz + eye_trans);
    
    // send the world pos (xy) through, for refraction calculations
    output.world_pos = world_mult_ie_dx_dy * pos_in.xy + world_trans;

    output.terrain_normal = ground_normal;
    
    return output;
}



// simplified gamma correction.

float3 ToLinear(float3 srgb)
{
    return pow(abs(srgb), 2.2f);
}

float3 FromLinear(float3 lin)
{
    return pow(abs(lin), 1.0f / 2.2f);
}


float pseudoRandom( float2 p )
{
  // We need irrationals for pseudo randomness.
  // Most (all?) known transcendental numbers will (generally) work.
  const float2 r = float2(
    23.1406926327792690,  // e^pi (Gelfond's constant)
     2.6651441426902251); // 2^sqrt(2) (Gelfond–Schneider constant)
  return frac( abs(cos( ( 123456789. % 1e-7 + 256. * dot(p,r) ) ) ) );  
}
   

// Pixel Shader for Water

float4 WaterPixelShader( WATER_PS_INPUT input ) : SV_Target
{
    // approximate Fresnel factor
	float3 normal_vec = input.normal;
/*	
	if (is_dissipation_threshold_on){
		int2 tempCoord = int2(input.tex_coord.x, input.tex_coord.y);
		float2 randomSeed = txGrass.Sample(samLinear, tempCoord).rg;
		normal_vec *= (0.5f + (pseudoRandom (float2(randomSeed)) * input.auxiliary.a)* 0.5f);
	//		normal_vec = normalize(normal_vec);
	}
*/
    // (this is the percentage that reflects, as opposed to transmits)
    float fresnel = (1-fresnel_coeff) + fresnel_coeff * pow(max(0, 1 - dot(input.eye, normal_vec)), fresnel_exponent);
	
	// this avoids jittering of water graphics artifact on land.
	if (input.water_depth < min(1e-10, sqrt_sqrt_epsilon)) {
		fresnel = 0;
	}
    // reflected light
    float3 reflection_dir = reflect(-input.eye, normal_vec);
    float3 reflect_col; 
	if (water_shading == 0){
		reflect_col = txSkybox.Sample(samLinear, reflection_dir).rgb;
	} else {
		float colormap_u = 0; //  this u is the u from texture coordinate not velocity (u,v).
		if (water_shading == 1){ //eta
			colormap_u = (input.eta - water_colormap_min) / (water_colormap_max - water_colormap_min);
		} else if (water_shading == 2){ //u
			colormap_u = (input.u - water_colormap_min) / (water_colormap_max - water_colormap_min);
		} else if (water_shading == 3){ //v
			colormap_u = (input.v - water_colormap_min) / (water_colormap_max - water_colormap_min);
		} else if (water_shading == 4){// v mag
			float speed = sqrt(input.v * input.v + input.u * input.u);
			colormap_u = (speed - water_colormap_min) / (water_colormap_max - water_colormap_min);  
		} else if (water_shading == 5){ // vorticity
			colormap_u = (input.vorticity - water_colormap_min) / (water_colormap_max - water_colormap_min);  
		} else if (water_shading == 6){ // max eta
			colormap_u = (input.auxiliary.r + input.B - water_colormap_min) / (water_colormap_max - water_colormap_min);  
		}

		reflect_col = txColormap.Sample(samLinear, float2(clamp(colormap_u,0.01,.99),0.5f)).rgb;   
		if(isGridOn){
			reflect_col += (1 - txGrid.Sample(samLinear, input.tex_coord).rgb);
		}
	}

	// to visualize breaking
	/*if (is_dissipation_threshold_on && input.auxiliary.b > dissipation_threshold){
		reflect_col = float3(1,1,1);
	}*/
	if (is_dissipation_threshold_on){
		float greyScale = input.auxiliary.a ;
		reflect_col += float3(greyScale, greyScale, greyScale);
	}
    //specular
    reflect_col += specular_intensity * pow(max(0,dot(reflection_dir, light_dir)), specular_exponent);

    // refracted light
    float3 refract_col;
    float attenuation_factor;
    float3 refract_dir = refract( -input.eye, normal_vec, 1.0f / refractive_index);

    if (refract_dir.z > 0) {
        // refracted ray goes up into the sky
        // (this can only happen with slanted water surfaces)
        // note: we have no way of knowing the distance through the water that the
        // ray has travelled, so we just hard code the attenuation
        refract_col = txSkybox.Sample(samLinear, refract_dir).rgb;
        attenuation_factor = 0.8f;
        
    } else {
        // refracted ray hits the grass
        float dist = -input.water_depth / refract_dir.z;
        float2 base_pos = input.world_pos.xy + dist * refract_dir.xy;
        float2 refract_tex_coord = world_to_grass_tex * base_pos + grass_tex_of_origin;
        if (terrain_shading == 0){
			refract_col = txGrass.Sample(samLinear, refract_tex_coord).rgb;
		} else {
			const float temp_B = input.B;
			float colormap_u = 0;
			if (temp_B <= 0){
				colormap_u = 0.5 * (input.B - terrain_colormap_min) / (0.0001f - terrain_colormap_min);  
			} else {
				colormap_u = 0.5 + 0.5 * (input.B - 0) / (terrain_colormap_max - (-0.0001));  
			}
			refract_col = txGrass.Sample(samLinear, float2(clamp(colormap_u,0.01,.99),0.5f)).rgb;   
		}
		
		if(isGridOn){
			refract_col += (1 - txGrid.Sample(samLinear, refract_tex_coord).rgb);
		}

		refract_col = TerrainColour(refract_col, input.terrain_normal);  // apply terrain lighting (approximation)
        attenuation_factor = (1-attenuation_1) * exp(-attenuation_2 * dist);
    }
	if (is_dissipation_threshold_on){
		float greyScale = input.auxiliary.a ;
		refract_col += float3(greyScale, greyScale, greyScale);
	}
  
    // combine reflection & refraction
    float3 result = FromLinear(lerp(lerp(ToLinear(deep_col),
                                         ToLinear(refract_col),
                                         attenuation_factor),
                                    ToLinear(reflect_col),
                                    fresnel));
    
	
    return float4(result, 1);
}


// Skybox shaders

SKYBOX_PS_INPUT SkyboxVertexShader( float3 pos : POSITION )
{
    SKYBOX_PS_INPUT output;

    output.pos = mul(skybox_mtx, pos).xyzz;   // set w=z
    output.pos.w *= 1.0001f;  // prevent it colliding with the far clip plane

    output.tex = pos;

    return output;
}

float4 SkyboxPixelShader( SKYBOX_PS_INPUT input ) : SV_TARGET
{
    // tex lookup
    return txSkybox.Sample(samLinear, input.tex);
}
