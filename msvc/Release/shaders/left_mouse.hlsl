/* -*- c++ -*-
 *
 * FILE:
 *   left_mouse.hlsl
 *
 * PURPOSE:
 *   Shaders for the left mouse button effects (add/remove water, etc)
 *
 * AUTHOR:
 *   Stephen Thompson <stephen@solarflare.org.uk>
 *
 * CREATED:
 *   9-Nov-2011
 *
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
 */

cbuffer LeftMouseConstBuffer : register( b0 )
{
    float2 scale;
    float2 bias;

    float two_over_nx_plus_four;
    float two_over_ny_plus_four;

    float disp_A, disp_B;   // control the amount to add or subtract from the water height.
};

// .r = w_bar
// (This is a COPY of the current state, because we can't read from the buffer we are writing to)
Texture2D<float4> txState : register( t0 );

// .b = B_average
Texture2D<float3> txBottom : register( t1 );

struct VS_INPUT {
    float2 pos : POSITION;   // Position in "brush space" (centre (0,0), radius 1)
};

struct VS_OUTPUT {
    float4 pos : SV_POSITION;
    float2 tex_idx : TEX_IDX;  // should be cast to an integer by the pixel shader
    float radius : RADIUS;
};

VS_OUTPUT LeftMouseVertexShader(VS_INPUT input)
{
    VS_OUTPUT output;

    output.tex_idx = scale * input.pos + bias;
    
    output.pos.x = output.tex_idx.x * two_over_nx_plus_four - 1;
    output.pos.y = 1 - output.tex_idx.y * two_over_ny_plus_four;
    output.pos.z = 0.5f;
    output.pos.w = 1.0f;

    output.radius = input.pos.x == 0 && input.pos.y == 0 ? 0.0f : 1.0f;

    return output;
}

float4 LeftMousePixelShader(VS_OUTPUT input) : SV_TARGET
{
    const int3 idx = int3(input.tex_idx.x, input.tex_idx.y, 0);
    float4 current = txState.Load(idx);
    float B = txBottom.Load(idx).b;

    float dh = disp_A + disp_B * input.radius;
    
    
    float old_w = current.r;
    float old_h = old_w - B;
    float new_w = max(B, old_w + dh);
    float new_h = new_w - B;

    // crude velocity calculation. no need to do the full "epsilon" treatment, we will just put a min on
    // old_h, which has the effect of lowering the velocity if h is currently small.
    float limited_old_h = max(0.5f, old_h);
    float old_u = current.g / limited_old_h;
    float old_v = current.b / limited_old_h;
    float new_hu = new_h * old_u;
    float new_hv = new_h * old_v;
    
    return float4(new_w, new_hu, new_hv, 0);
}
