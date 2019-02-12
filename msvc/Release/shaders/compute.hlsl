/* -*- c++ -*-
 * 
 * FILE:
 *   compute.hlsl
 *
 * PURPOSE:
 *   HLSL implementation of the Kurganov-Petrova scheme for solving
 *   the 2D shallow water equations. See:
 *
 *   A. Kurganov and G. Petrova, "A Second-Order Well-Balanced
 *   Positivity Preserving Central-Upwind Scheme for the Saint-Venant
 *   System", Commun. Math. Sci. 5, 133-160, 2007.
 *
 * AUTHOR:
 *   Stephen Thompson <stephen@solarflare.org.uk>
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
 *  Celeris is free software: you can redistribute it
 *   and/or modify it under the terms of the GNU General Public
 *   License as published by the Free Software Foundation, either
 *   version 3 of the License, or (at your option) any later version.
 *
 *  Celeris is distributed in the hope that it will be
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

cbuffer SimConstBuffer : register( b0 )
{
    // THETA = parameter for minmod limiter. between 1 and 2 inclusive.
    // 1 = more dissipative, 2 = more oscillatory. 1.3 is a good default.
    // (Note: we actually store 2*THETA.)
    float TWO_THETA;

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
	float dt_old;
	float dt_old_old;
    float epsilon;      // suggestion: dx^4 or dy^4. but that may be too small for floats?
    int nx; // number of cells in x direction, excluding ghost zones
	int	ny;// number of cells in y direction, excluding ghost zones  
    float friction;              // m s^-1
	int isManning;
	float seaLevel;

	float dissipation_threshold; // For visualization purposes, not simulation.
	float whiteWaterDecayRate; // For visualization purposes, not simulation.
};

cbuffer SolitaryWave : register( b0 )
{
    float sw_g;
	
    float sw_dx;
    float sw_dy;
	
	float sw_param_H; 
	float sw_x_zero;
	float sw_y_zero;
	float sw_theta;
};

cbuffer IrregularWavesConstBuffer : register( b1 )
{
    float4  irregularWavesWest[1000],
			irregularWavesEast[1000],
			irregularWavesSouth[1000],
			irregularWavesNorth[1000]; //  IMPORTANT NOTE: The hardcoded number needs to be the same as MAX_NUM_OF_IRREGULAR_WAVES defined in engine.hpp.
	int numberOfWavesWest, numberOfWavesEast, numberOfWavesSouth, numberOfWavesNorth;
};

cbuffer  IrregularWavesColumnRowConstBuffer: register( b2 ){
		int columnNumber;
		int rowNumber;
}

cbuffer TimeIntegrationBuffer : register( b1 )
{
    int timeScheme; // in the main code enum:int is used. but hlsl does not compile enums, so I use int explicitly here.
};

cbuffer BoundaryConstBuffer : register( b0 )
{
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
    float boundary_g;
	float boundary_dt;
    float total_time;
    float sa1, skx1, sky1, so1;
    float sa2, skx2, sky2, so2;
    float sa3, skx3, sky3, so3;
    float sa4, skx4, sky4, so4;
    float sdecay;
};

// .r = B(j, k+1/2)   "BN" or "BY"
// .g = B(j+1/2, k)   "BE" or "BX"
// .b = B(j, k)       "BA"
// Note: this is the same size as the main textures i.e. includes ghost zones.
Texture2D<float3> txBottom : register( t1 );

// .r = w_bar
// .g = hu_bar
// .b = hv_bar
// .a = (unused)
Texture2D<float4> txState : register( t0 );


// .r = dU1_dt (for w)
// .g = dU2_dt (for hu)
// .b = dU3_dt (for hv)
// .a = (unused)
Texture2D<float4> oldGradients : register( t4 );
Texture2D<float4> oldOldGradients : register( t5 );

// .r = F_star 
// .g = G_star 
// .b = (unused) 
// .a = (unused)
Texture2D<float4> F_G_star_oldGradients : register( t8 );
Texture2D<float4> F_G_star_oldOldGradients : register( t9 );

Texture2D<float4> predictedGradients : register( t10 );
Texture2D<float4> F_G_star_predictedGradients : register( t11 );

Texture2D<float4> current_state : register( t12 );

// .r = h_bar
// .g = D1x (to find hu_bar)
// .b = D1y (to find hv_bar)
// .a = (unused)
// TODO: Remove myRHS
Texture2D<float4> myRHS : register( t0 );

// .r = N
// .g = E
// .b = S
// .a = W
Texture2D<float4> txH : register( t0 );
Texture2D<float4> txU : register( t1 );
Texture2D<float4> txV : register( t2 );


// normal texture
// r = nX
// g = nY
// b = nZ
// a = unused

Texture2D<float4> txNormal : register( t3 );

// r = wasInudated
Texture2D<float4> txAuxiliary1 : register( t3 );

// r = maxFlowDepth, b = breaking
Texture2D<float4> txAuxiliary2 : register( t4 );

// .r = 1 (w-flux)
// .g = 2 (hu-flux)
// .b = 3 (hv-flux)
// .a = (unused)
Texture2D<float4> txXFlux : register( t2 );
Texture2D<float4> txYFlux : register( t3 );


// .r = h_bar
// .g = D1x (to find hu_bar)
// .b = D1y (to find hv_bar)
// .a = (unused)

Texture2D<float4> coefMatx : register( t6 );
Texture2D<float4> coefMaty : register( t7 );
Texture2D<float4> D1x : register( t1 );
Texture2D<float4> Xx : register( t2 );

Texture2D<float4> D1y : register( t1 );
Texture2D<float4> Xy : register( t3 );

// .r = w_bar
// .g = hu_bar
// .b = hv_bar
// .a = (unused)
Texture2D<float4> irregularWavesLastReduction : register( t0 );

struct VS_INPUT {
    float2 tex_idx : TEX_IDX;
};

struct VS_OUTPUT {
    float4 pos : SV_POSITION;
    float2 tex_idx : TEX_IDX; // should be cast to integer.
};

struct PASS_1_OUTPUT {
    float4 h  : SV_TARGET0;    // {hN, hE, hS, hW}
    float4 u  : SV_TARGET1;    // {uN, uE, uS, uW}
    float4 v  : SV_TARGET2;    // {vN, vE, vS, vW}
    float4 n  : SV_TARGET3;    // {nX, nY, nZ, unused}
	float4 aux: SV_TARGET4;    // (maxFlowDepth, unused, breaking, unused)
};

struct PASS_2_OUTPUT {
    float4 xflux : SV_TARGET0;   // {Hx1, Hx2, Hx3, unused}
    float4 yflux : SV_TARGET1;   // {Hy1, Hy2, Hy3, unused}
	float4 auxiliary : SV_TARGET2; //(MaxDepth, unused, unused, unused). The unused variables will be used for max v, u, and |V| later. 
};

struct PASS_3_OUTPUT {
    float4 newState : SV_TARGET0;   // {h, hu, hv, unused}
    float4 dU_by_dt : SV_TARGET1;   // {dU1/dt, dU2/dt, dU3/dt, unused} U1 is w, U2 and U3 are used to get hu and hv.
	float4 F_G_star : SV_TARGET2;
};

struct ADD_SOLITARY_OUTPUT {
    float4 newState : SV_TARGET0;   // {h, hu, hv, unused}
};

VS_OUTPUT SimVertexShader(VS_INPUT input)
{
    VS_OUTPUT output;
    output.pos.x = input.tex_idx.x * two_over_nx_plus_four - 1;
    output.pos.y = 1 - input.tex_idx.y * two_over_ny_plus_four;
    output.pos.z = 0.5f;
    output.pos.w = 1.0f;
    output.tex_idx = input.tex_idx;
    return output;
}

int3 GetTexIdx(VS_OUTPUT input)
{
    return int3(input.tex_idx.x, input.tex_idx.y, 0);
}

float MinMod(float a, float b, float c)
{
    return (a > 0 && b > 0 && c > 0) ? min(min(a,b),c)
        : (a < 0 && b < 0 && c < 0) ? max(max(a,b),c) : 0;
}

void Reconstruct(float west, float here, float east,
                 out float out_west, out float out_east, out float standard_deviation)
{
    // west, here, east = values of U_bar at j-1, j, j+1 (or k-1, k, k+1)
    // out_west, out_east = reconstructed values of U_west and U_east at (j,k)

	//ST_: YOU MAY WANT TO ADD DISSPATION CORRECTION. LOOK INTO YOUR MATLAB CODE!
    
	float z1 = TWO_THETA * (here - west);
	float z2 = (east - west);
	float z3 = TWO_THETA * (east - here);
    float dx_grad_over_two = 0.25f * MinMod( z1, z2, z3);

	float mu = 0.5f * (z1 + z2 + z3) / 3.0f;
	standard_deviation = sqrt(((z1 - mu)*(z1 - mu) + (z2 - mu)*(z2 - mu) + (z3 - mu)*(z3 - mu)) / 3.0f);
    out_east = here + dx_grad_over_two;
    out_west = here - dx_grad_over_two;
}

void CorrectW(float B_west, float B_east, float w_bar,
              inout float w_west, inout float w_east)
{
    if (w_east < B_east) {
        w_east = B_east;
        w_west = max(B_west, 2 * w_bar - B_east);
            
    } else if (w_west < B_west) {
        w_east = max(B_east, 2 * w_bar - B_west);
        w_west = B_west;
    }
}              

void CalcUV(float4 h, float4 hu, float4 hv, out float4 u, out float4 v)
{
    // in:  {hN, hE, hS, hW},  {huN, huE, huS, huW},  {hvN, hvE, hvS, hvW}
    // out: {uN, uE, uS, uW},  {vN, vE, vS, vW}
    const float4 h2 = h * h;
    const float4 h4 = h2 * h2;
    const float4 divide_by_h = sqrt(2.0f) * h / sqrt(h4 + max(h4, epsilon));
	//const float4 divide_by_h =  h / max(h2 ,epsilon);
    u = divide_by_h * hu;
    v = divide_by_h * hv;
}

float CalcUV_Scalar(float h, float hu, float hv, out float u, out float v)
{
    const float h2 = h * h;
    const float h4 = h2 * h2;
    const float divide_by_h = sqrt(2.0f) * h / sqrt(h4 + max(h4, epsilon));
	//const float divide_by_h =  h / max(h2 , epsilon);
    u = divide_by_h * hu;
    v = divide_by_h * hv;
    return divide_by_h;
}

float CalcU_Boundary(float h, float hu)
{
    const float h2 = h * h;
    const float h4 = h2 * h2;

    const float divide_by_h = sqrt(2.0f) * h / sqrt(h4 + max(h4, boundary_epsilon));
	//const float divide_by_h =  h / max(h2 , boundary_epsilon);
    return divide_by_h * hu;
}

float NumericalFlux(float aplus, float aminus, float Fplus, float Fminus, float Udifference)
{
    if (aplus - aminus != 0) {  
        return (aplus * Fminus - aminus * Fplus + aplus * aminus * Udifference) / (aplus - aminus);
    } else {
        return 0;
    }
}


// Pass 1 -- Reconstruct h, u, v at the four edges of each cell.
//
// t0: txState
// t1: txBottom
// Output: txH, txU, txV, txNormal
//
// Runs on bulk + first ghost layer either side
PASS_1_OUTPUT Pass1(VS_OUTPUT input)
{
    // Read in relevant texture values

    const int3 idx = int3(input.tex_idx.x, input.tex_idx.y, 0);

    // in = {w, hu, hv, _} (cell average)
    const float4 in_here = txState.Load(idx);
    const float4 in_south = txState.Load(idx + int3(0, -1, 0));
    const float4 in_north = txState.Load(idx + int3(0, 1, 0));
    const float4 in_west = txState.Load(idx + int3(-1, 0, 0));
    const float4 in_east = txState.Load(idx + int3(1, 0, 0));

    float4 B;     // {BN, BE, BS, BW}
    B.rg = txBottom.Load(idx).rg;
    B.b = txBottom.Load(idx + int3(0, -1, 0)).r;
    B.a = txBottom.Load(idx + int3(-1, 0, 0)).g;
    
    
    // Reconstruct w, hu and hv at the four cell edges (N, E, S, W)
    
    float4 w;    // {wN, wE, wS, wW}
    float4 hu;   // {huN, huE, huS, huW}
    float4 hv;   // {hvN, hvE, hvS, hvW}
    
	float max_sd2;
	float temp_sd2;

    Reconstruct(in_west.r, in_here.r, in_east.r, w.a, w.g, temp_sd2); 
	max_sd2 = temp_sd2 * one_over_dx;
	Reconstruct(in_south.r, in_here.r, in_north.r, w.b, w.r, temp_sd2); 
    max_sd2 = atan(max(max_sd2, temp_sd2 * one_over_dy));
	
    Reconstruct(in_west.g, in_here.g, in_east.g, hu.a, hu.g, temp_sd2); 
    Reconstruct(in_south.g, in_here.g, in_north.g, hu.b, hu.r, temp_sd2); 
    
    Reconstruct(in_west.b, in_here.b, in_east.b, hv.a, hv.g, temp_sd2); 
    Reconstruct(in_south.b, in_here.b, in_north.b, hv.b, hv.r, temp_sd2); 


    // Correct the w values to ensure positivity of h

    CorrectW(B.a, B.g, in_here.r, w.a, w.g);    // wW and wE (from BW, BE, wbar)
    CorrectW(B.b, B.r, in_here.r, w.b, w.r);    // wS and wN (from BS, BN, wbar)


    // Reconstruct h from (corrected) w
    // Calculate u and v from h, hu and hv

    PASS_1_OUTPUT output;

	float4 h = w - B;
    output.h = h;
    CalcUV(output.h, hu, hv, output.u, output.v);
	
	//output.v = 0;

    // Calculate normal 

    float3 normal;
    //normal.x = (w.g - w.a) * one_over_dx;
    //normal.y = (w.r - w.b) * one_over_dy;
    //normal.z = -1;
    normal.x = (in_west.r - in_east.r) * one_over_dx;
    normal.y = (in_south.r - in_north.r) * one_over_dy;
    normal.z = 2;
    normal = normalize(normal);

	float maxInundatedDepth = max((h.r + h.g + h.b + h.a)/4.0f, txAuxiliary1.Load(idx).r);
	
    output.n = float4(normal.x, normal.y, normal.z, 0);
	
	float breaking_white = txAuxiliary1.Load(idx).a;
	if (max_sd2 > dissipation_threshold){
		breaking_white = 1.0f;
	}
    output.aux = float4(maxInundatedDepth, 0, max_sd2, breaking_white);
    return output;
}


// Pass 2 -- Calculate fluxes
//
// t0: txH
// t1: txU
// t2: txV
// t3: normal_tex
// Output: txXFlux and txYFlux
//
// Runs on bulk + first ghost layer to west and south only
PASS_2_OUTPUT Pass2( VS_OUTPUT input )
{
    const int3 idx = GetTexIdx(input);
    
    const float2 h_here = txH.Load(idx).rg;                 // {hN, hE}   evaluated here
    const float hW_east = txH.Load(idx + int3(1,0,0)).a;    // hW evaluated at (j+1, k)
    const float hS_north = txH.Load(idx + int3(0,1,0)).b;   // hS evaluated at (j, k+1)

    const float2 u_here = txU.Load(idx).rg;
    const float uW_east = txU.Load(idx + int3(1,0,0)).a;
    const float uS_north = txU.Load(idx + int3(0,1,0)).b;

    const float2 v_here = txV.Load(idx).rg;
    const float vW_east = txV.Load(idx + int3(1,0,0)).a;
    const float vS_north = txV.Load(idx + int3(0,1,0)).b;
    
    // compute wave speeds
    const float2 cNE = sqrt((g * h_here.rg));    // {cN, cE} evaluated here ?? ST_: I removed a max (0,h) here. Correction step makes sure that h is positive. If it is not, then there is a mistake in the modeling.
    const float cW = sqrt((g * hW_east));        // cW evaluated at (j+1, k)
    const float cS = sqrt((g * hS_north));       // cS evaluated at (j, k+1)

    // compute propagation speeds
    const float aplus  = max(max(u_here.g + cNE.g, uW_east + cW), 0);
    const float aminus = min(min(u_here.g - cNE.g, uW_east - cW), 0);
    const float bplus  = max(max(v_here.r + cNE.r, vS_north + cS), 0);
    const float bminus = min(min(v_here.r - cNE.r, vS_north - cS), 0);

    // compute fluxes
    PASS_2_OUTPUT output = (PASS_2_OUTPUT) 0;
    output.xflux.r = NumericalFlux(aplus,
                                   aminus,
                                   hW_east * uW_east,
                                   h_here.g * u_here.g,
                                   hW_east - h_here.g);
    
    output.xflux.g = NumericalFlux(aplus,
                                   aminus,
                                   hW_east * (uW_east * uW_east + half_g * hW_east),
                                   h_here.g * (u_here.g * u_here.g + half_g * h_here.g),
                                   hW_east * uW_east - h_here.g * u_here.g);

    output.xflux.b = NumericalFlux(aplus,
                                   aminus,
                                   hW_east * uW_east * vW_east,
                                   h_here.g * u_here.g * v_here.g,
                                   hW_east * vW_east - h_here.g * v_here.g);

    output.yflux.r = NumericalFlux(bplus,
                                   bminus,
                                   hS_north * vS_north,
                                   h_here.r * v_here.r,
                                   hS_north - h_here.r);

    output.yflux.g = NumericalFlux(bplus,
                                   bminus,
                                   hS_north * uS_north * vS_north,
                                   h_here.r * u_here.r * v_here.r,
                                   hS_north * uS_north - h_here.r * u_here.r);

    output.yflux.b = NumericalFlux(bplus,
                                   bminus,
                                   hS_north * (vS_north * vS_north + half_g * hS_north),
                                   h_here.r * (v_here.r * v_here.r + half_g * h_here.r),
                                   hS_north * vS_north - h_here.r * v_here.r);
	
	float breaking_white = txAuxiliary2.Load(idx).a;
	
	breaking_white *= pow(abs(whiteWaterDecayRate),dt);
	
	// max_sd2 (i.e. b) varies in every time step, unlike r0
	output.auxiliary = float4(txAuxiliary2.Load(idx).r, 0, txAuxiliary2.Load(idx).b, breaking_white); 
    return output;
}

float FrictionCalc(float hu, float hv, float h)
{	
	float f = friction;
	const float divide_by_h = sqrt(2.0f) * h / sqrt(pow(h,4) + max(pow(h,4), epsilon));
	if(isManning == 1)
	{
		f = g * friction*friction * pow(abs(divide_by_h),0.3333f);
	}
	return    f * sqrt(hu*hu + hv*hv) * divide_by_h * divide_by_h;
}

// Returns the coefficient for Y in our AdamsBashforth third order
// time integration with variable time-step.
// .r = coefficient for Y ^ (n - 0)
// .g = coefficient for Y ^ (n - 1)
// .b = coefficient for Y ^ (n - 2)
float3 AdamsBashforthCoef(float dt_n, float dt_n_minus_1, float dt_n_minus_2) {
	float3 coef;
	coef.r = (dt_n / dt_n_minus_1) * ((2 * dt_n + 6 * dt_n_minus_1 + 3 * dt_n_minus_2) / (dt_n_minus_1 + dt_n_minus_2)) + 6;
	coef.g = -(dt_n / dt_n_minus_1) * ((2 * dt_n + 3 * dt_n_minus_1 + 3 * dt_n_minus_2) / (dt_n_minus_2));
	coef.b = (dt_n / dt_n_minus_2) * ((2 * dt_n + 3 * dt_n_minus_1) / (dt_n_minus_1 + dt_n_minus_2));
	return coef;
}

// Returns the second order backward derivative at dt_n.
float SecondOrderBackward(float y_n, float y_n_minus_1, float y_n_minus_2, float dt_n, float dt_n_minus_1, float dt_n_minus_2) {
	float3 coef;
	coef.r = (2*dt_n_minus_1 + dt_n_minus_2)/(dt_n_minus_1 * (dt_n_minus_1 + dt_n_minus_2));
	coef.g = -(dt_n_minus_1 + dt_n_minus_2)/(dt_n_minus_1 * dt_n_minus_2);
	coef.b = (dt_n_minus_1)/(dt_n_minus_2 * (dt_n_minus_1 + dt_n_minus_2));
	return (coef.r * y_n) + (coef.g * y_n_minus_1) + (coef.b * y_n_minus_2);
}

// Returns the second order central derivative at dt_n_minus_1.
float SecondOrderCentral(float y_n, float y_n_minus_1, float y_n_minus_2, float dt_n, float dt_n_minus_1, float dt_n_minus_2) {
	float3 coef;
	coef.r = (dt_n_minus_2) / (dt_n_minus_1 * (dt_n_minus_1 + dt_n_minus_2));
	coef.g = (dt_n_minus_1 - dt_n_minus_2) / (dt_n_minus_1 * dt_n_minus_2);
	coef.b = -(dt_n_minus_1) / (dt_n_minus_2 * (dt_n_minus_1 + dt_n_minus_2));
	return (coef.r * y_n) + (coef.g * y_n_minus_1) + (coef.b * y_n_minus_2);
}

// Returns the second order forward derivatives at dt_n_minus_2.
float SecondOrderForward(float y_n, float y_n_minus_1, float y_n_minus_2, float dt_n, float dt_n_minus_1, float dt_n_minus_2) {
	float3 coef;
	coef.r = -(dt_n_minus_2) / (dt_n_minus_1 * (dt_n_minus_1 + dt_n_minus_2));
	coef.g = (dt_n_minus_1 + dt_n_minus_2) / (dt_n_minus_1 * dt_n_minus_2);
	coef.b = -(dt_n_minus_1 + 2*dt_n_minus_2) / (dt_n_minus_2 * (dt_n_minus_1 + dt_n_minus_2));
	return (coef.r * y_n) + (coef.g * y_n_minus_1) + (coef.b * y_n_minus_2);
}

// Returns the second order derivatives of hte given variable.
// .r = dy/dt at t_n (backward)
// g. = dy/dt at t_n_minus_1 (central)
// b. = dy/dt at t_n_minus_2 (forward)
float3 SecondOrderDerivatives(float y_n, float y_n_minus_1, float y_n_minus_2, float dt_n, float dt_n_minus_1, float dt_n_minus_2) {
	float3 derivatives;
	derivatives.r = SecondOrderBackward(y_n, y_n_minus_1, y_n_minus_2, dt_n, dt_n_minus_1, dt_n_minus_2);
	derivatives.g = SecondOrderCentral(y_n, y_n_minus_1, y_n_minus_2, dt_n, dt_n_minus_1, dt_n_minus_2);
	derivatives.b = SecondOrderForward(y_n, y_n_minus_1, y_n_minus_2, dt_n, dt_n_minus_1, dt_n_minus_2);
	return derivatives;
}

// Pass 3 -- Do timestep and calculate new w_bar, hu_bar, hv_bar.
//
// t0: txState
// t1: txBottom
// t2: txXFlux
// t3: txYFlux
// Output: New txState
//
// Runs on interior points only
PASS_3_OUTPUT Pass3Predictor( VS_OUTPUT input ) 
{
    const int3 idx = GetTexIdx(input);

	const float3 B_here = txBottom.Load(idx); //north, east, here
	const float BX_west = txBottom.Load(idx + int3(-1,0,0)).g;
    const float BY_south = txBottom.Load(idx + int3(0,-1,0)).r;

	const float3 in_state_here  = txState.Load(idx).rgb;   // w, hu and hv (cell avgs, evaluated here)
	const float3 in_state_right = txState.Load(idx + int3(+1,0,0)).rgb;   // w, hu and hv at right cell (cell avgs, evaluated here)
	const float3 in_state_left  = txState.Load(idx + int3(-1,0,0)).rgb;   // w, hu and hv at left cell (cell avgs, evaluated here)
	const float3 in_state_up    = txState.Load(idx + int3(0,+1,0)).rgb;   // w, hu and hv at up cell (cell avgs, evaluated here)
	const float3 in_state_down  = txState.Load(idx + int3(0,-1,0)).rgb;   // w, hu and hv at down cell (cell avgs, evaluated here)
	const float3 in_state_up_left    = txState.Load(idx + int3(-1,+1,0)).rgb;
	const float3 in_state_up_right   = txState.Load(idx + int3(+1,+1,0)).rgb;
	const float3 in_state_down_left  = txState.Load(idx + int3(-1,-1,0)).rgb;
	const float3 in_state_down_right = txState.Load(idx + int3(+1,-1,0)).rgb;

	const float3 myCoefMatx = coefMatx.Load(idx).rgb;
	const float3 myCoefMaty = coefMaty.Load(idx).rgb;

	const float h = max(0, in_state_here.r - B_here.b);

	//d stencil
	//must pass d as a resource. redundant calculation	
	const float d_here = max(0, seaLevel - B_here.b);
	const float d2_here = d_here * d_here;
	const float d3_here = d2_here * d_here;

	const float d_left  = max(0, seaLevel - txBottom.Load(idx + int3(-1,0,0)).b);
	const float d_right = max(0, seaLevel - txBottom.Load(idx + int3(+1,0,0)).b);
	const float d_down  = max(0, seaLevel - txBottom.Load(idx + int3(0,-1,0)).b);
	const float d_up    = max(0, seaLevel - txBottom.Load(idx + int3(0,+1,0)).b);
	
	const float d_up_left    = max(0, seaLevel - txBottom.Load(idx + int3(-1,+1,0)).b);
	const float d_up_right   = max(0, seaLevel - txBottom.Load(idx + int3(+1,+1,0)).b);
	const float d_down_left  = max(0, seaLevel - txBottom.Load(idx + int3(-1,-1,0)).b);
	const float d_down_right = max(0, seaLevel - txBottom.Load(idx + int3(+1,-1,0)).b);

	const float d_left_left   = max(0, seaLevel - txBottom.Load(idx + int3(-2,0,0)).b);
	const float d_right_right = max(0, seaLevel - txBottom.Load(idx + int3(+2,0,0)).b);
	const float d_down_down   = max(0, seaLevel - txBottom.Load(idx + int3(0,-2,0)).b);
	const float d_up_up       = max(0, seaLevel - txBottom.Load(idx + int3(0,+2,0)).b);

	//eta stencil.
	
	const float eta_here = in_state_here.r - seaLevel;
	
//I can use in_state
	const float eta_left  = txState.Load(idx + int3(-1,0,0)).r - seaLevel;
	const float eta_right = txState.Load(idx + int3(+1,0,0)).r - seaLevel;
	const float eta_down  = txState.Load(idx + int3(0,-1,0)).r - seaLevel;
	const float eta_up    = txState.Load(idx + int3(0,+1,0)).r - seaLevel;

	const float eta_left_left   = txState.Load(idx + int3(-2,0,0)).r - seaLevel;
	const float eta_right_right = txState.Load(idx + int3(+2,0,0)).r - seaLevel;
	const float eta_down_down   = txState.Load(idx + int3(0,-2,0)).r - seaLevel;
	const float eta_up_up       = txState.Load(idx + int3(0,+2,0)).r - seaLevel;

	const float eta_up_left    = txState.Load(idx + int3(-1,+1,0)).r - seaLevel;
	const float eta_up_right   = txState.Load(idx + int3(+1,+1,0)).r - seaLevel;
	const float eta_down_left  = txState.Load(idx + int3(-1,-1,0)).r - seaLevel;
	const float eta_down_right = txState.Load(idx + int3(+1,-1,0)).r - seaLevel;


    const float3 xflux_here = txXFlux.Load(idx).rgb;
    const float3 xflux_west = txXFlux.Load(idx + int3(-1,0,0)).rgb;
    const float3 yflux_here = txYFlux.Load(idx).rgb;
    const float3 yflux_south = txYFlux.Load(idx + int3(0,-1,0)).rgb;
	
    // friction calculations
    float u, v;
    CalcUV_Scalar(h, in_state_here.g, in_state_here.b, u, v);
	//const float3 modified_in_state_here= float3(in_state_here.r,h*u,h*v); //sets hu=h*u and hv=h*v
	const float one_over_d2x = one_over_dx * one_over_dx;
	const float one_over_d3x = one_over_d2x * one_over_dx;
	
	const float one_over_d2y = one_over_dy * one_over_dy;
	const float one_over_d3y = one_over_d2y * one_over_dy;

	const float dd_by_dx = (-d_right_right + 8*d_right - 8*d_left + d_left_left) * one_over_dx/12.0f;
	const float dd_by_dy = (-d_up_up       + 8*d_up    - 8*d_down + d_down_down) * one_over_dy/12.0f;
	const float eta_by_dx_dy = (0.25f*one_over_dx*one_over_dy * (eta_up_right - eta_down_right - eta_up_left + eta_down_left));
	const float eta_by_dx_dx = one_over_d2x * (eta_right - 2*eta_here + eta_left);
	const float eta_by_dy_dy = one_over_d2y * (eta_up - 2*eta_here + eta_down);

	const float Psi1x = Bcoef_g * (d3_here) * ((eta_right_right - 2*eta_right + 2*eta_left - eta_left_left)*(0.5f*one_over_d3x)
						+ (eta_up_right - eta_up_left - 2*eta_right + 2*eta_left + eta_down_right - eta_down_left) * (0.5f*one_over_dx*one_over_d2y));
	const float Psi2x = Bcoef_g * (d2_here) * (dd_by_dx * (2*eta_by_dx_dx + eta_by_dy_dy) + dd_by_dy * eta_by_dx_dy);
	

	const float Psi1y = Bcoef_g * (d3_here) * ((eta_up_up - 2*eta_up + 2*eta_down - eta_down_down)*(0.5f*one_over_d3y)
						+ (eta_up_right + eta_up_left - 2*eta_up + 2*eta_down - eta_down_right - eta_down_left) * (0.5f*one_over_dy*one_over_d2x));
	const float Psi2y = Bcoef_g * (d2_here) * (dd_by_dy * (eta_by_dx_dx + 2*eta_by_dy_dy) + dd_by_dx * eta_by_dx_dy);
	
	
	float friction_ = FrictionCalc(in_state_here.g, in_state_here.b, h);
	const float3 source_term = 
        float3(0,
               -g_over_dx * h * (B_here.g - BX_west)   - in_state_here.g*friction_ + (Psi1x + Psi2x) ,
               -g_over_dy * h * (B_here.r - BY_south)  - in_state_here.b*friction_ + (Psi1y + Psi2y));
        
    const float3 d_by_dt =
        (xflux_west - xflux_here) * one_over_dx
        + (yflux_south - yflux_here) * one_over_dy
        + source_term;
/*
	const float F_star = 0;
	const float G_star = 0;
*/
	const float F_star = (1.0f / 6.0f * d_here *
		(dd_by_dx * ((one_over_dy / 2.0f) * (in_state_up.b - in_state_down.b)) +
		(dd_by_dy * ((one_over_dx / 2.0f) * (in_state_right.b - in_state_left.b)))) +
		0 * (Bcoef_g / g + 1.0f / 3.0f) * d2_here * (one_over_dx*one_over_dy / 4.0f) * (in_state_up_right.b - in_state_down_right.b - in_state_up_left.b + in_state_down_left.b));
    
	const float G_star = (1.0f / 6.0f * d_here *
		(dd_by_dx * ((one_over_dy / 2.0f) * (in_state_up.g - in_state_down.g)) +
		(dd_by_dy * ((one_over_dx / 2.0f) * (in_state_right.g - in_state_left.g)))) +
		0 * (Bcoef_g / g + 1.0f / 3.0f) * d2_here * (one_over_dx*one_over_dy / 4.0f) * (in_state_up_right.g - in_state_down_right.g - in_state_up_left.g + in_state_down_left.g));

    // time stepping (third order)
	
	float wOut = 0, huOut = 0, hvOut = 0;
	
	if (timeScheme == 0)  //if timeScheme is Euler do:
	{
		wOut = in_state_here.r + dt * d_by_dt.r;
		huOut = myCoefMatx.r * in_state_left.g + myCoefMatx.g * in_state_here.g + myCoefMatx.b * in_state_right.g +
			dt * d_by_dt.g; // FG_Star is neglected
		hvOut = myCoefMaty.r * in_state_down.b + myCoefMaty.g * in_state_here.b + myCoefMaty.b * in_state_up.b +
			dt * d_by_dt.b; // FG_Star is neglected
	}
	else if (timeScheme == 1 || timeScheme == 2 || timeScheme == 3) // if time scheme is preditor or corrector
	{
		const float3 oldies = oldGradients.Load(idx).rgb;
		const float3 oldOldies = oldOldGradients.Load(idx).rgb;
		const float3 F_G_star_oldies = F_G_star_oldGradients.Load(idx).rgb;
		const float3 F_G_star_oldOldies = F_G_star_oldOldGradients.Load(idx).rgb;
		 
		if (timeScheme == 1) //if time scheme is predictor
		{
			wOut = in_state_here.r + dt / 12.0f * (23.0f * d_by_dt.r - 16.0f * oldies.r + 5.0f * oldOldies.r);
			huOut = myCoefMatx.r * in_state_left.g + myCoefMatx.g * in_state_here.g + myCoefMatx.b * in_state_right.g +
				dt / 12.0f * (23.0f * d_by_dt.g - 16.0f * oldies.g + 5.0f * oldOldies.g) +
				(2.0f * F_star - 3.0f * F_G_star_oldies.r + F_G_star_oldOldies.r);
			
			hvOut = myCoefMaty.r * in_state_down.b + myCoefMaty.g * in_state_here.b + myCoefMaty.b * in_state_up.b +
				dt / 12.0f * (23.0f * d_by_dt.b - 16.0f * oldies.b + 5.0f * oldOldies.b) +
				(2.0f * G_star - 3.0f * F_G_star_oldies.g + F_G_star_oldOldies.g);
		} 
		else if (timeScheme == 3) //if time scheme is adaptive predictor
		{
			float3 coef = AdamsBashforthCoef(dt, dt_old, dt_old_old);
			float3 derivativesFStar = SecondOrderDerivatives(F_star, F_G_star_oldies.r, F_G_star_oldOldies.r, dt, dt_old, dt_old_old);
			float3 derivativesGStar = SecondOrderDerivatives(G_star, F_G_star_oldies.g, F_G_star_oldOldies.g, dt, dt_old, dt_old_old);

			wOut = in_state_here.r + dt / 6.0f * (coef.r * d_by_dt.r + coef.g * oldies.r + coef.b * oldOldies.r);
			huOut = myCoefMatx.r * in_state_left.g + myCoefMatx.g * in_state_here.g + myCoefMatx.b * in_state_right.g +
				dt / 6.0f * (coef.r * d_by_dt.g + coef.g * oldies.g + coef.b * oldOldies.g) +
				dt / 6.0f * (coef.r * derivativesFStar.r + coef.g * derivativesFStar.g + coef.b * derivativesFStar.b);

			hvOut = myCoefMaty.r * in_state_down.b + myCoefMaty.g * in_state_here.b + myCoefMaty.b * in_state_up.b +
				dt / 6.0f * (coef.r * d_by_dt.b + coef.g * oldies.b + coef.b * oldOldies.b) +
				dt / 6.0f * (coef.r * derivativesGStar.r + coef.g * derivativesGStar.g + coef.b * derivativesGStar.b);
		}
		else if (timeScheme == 2) // if time scheme is corrector
		{
			const float3 current_state_here  = current_state.Load(idx).rgb;   // w, hu and hv (cell avgs, evaluated here)
			const float3 current_state_right = current_state.Load(idx + int3(+1,0,0)).rgb;   // w, hu and hv at right cell (cell avgs, evaluated here)
			const float3 current_state_left  = current_state.Load(idx + int3(-1,0,0)).rgb;   // w, hu and hv at left cell (cell avgs, evaluated here)
			const float3 current_state_up    = current_state.Load(idx + int3(0,+1,0)).rgb;   // w, hu and hv at up cell (cell avgs, evaluated here)
			const float3 current_state_down  = current_state.Load(idx + int3(0,-1,0)).rgb;   // w, hu and hv at down cell (cell avgs, evaluated here)
			
			const float3 predicted = predictedGradients.Load(idx).rgb;
			const float3 F_G_star_predicted = F_G_star_predictedGradients.Load(idx).rgb;
			wOut = current_state_here.r + dt/24.0f * (9*d_by_dt.r + 19*predicted.r  - 5*oldies.r  + 1*oldOldies.r);
			huOut = myCoefMatx.r * current_state_left.g + myCoefMatx.g * current_state_here.g + myCoefMatx.b * current_state_right.g +
				dt/24.0f * (9*d_by_dt.g + 19*predicted.g - 5*oldies.g + 1*oldOldies.g) + 
				F_star - F_G_star_predicted.r;
			hvOut = myCoefMaty.r * current_state_down.b + myCoefMaty.g * current_state_here.b + myCoefMaty.b * current_state_up.b +
				dt/24.0f * (9*d_by_dt.b + 19*predicted.b - 5*oldies.b + 1*oldOldies.b) + 
				G_star - F_G_star_predicted.g;
		}
	}

	PASS_3_OUTPUT output = (PASS_3_OUTPUT) 0;
	float4 result1;
	result1 = float4(wOut, huOut, hvOut, 0);
	output.newState = result1;
	float4 result2 = float4(d_by_dt.r, d_by_dt.g, d_by_dt.b, 0);
	output.dU_by_dt = result2;
	float4 result3 = float4(F_star, G_star, 0, 0);
	//float4 result3 = float4(in_state_up_right.g,in_state_down_right.g,in_state_up_left.g,-1+eta_down_left);
	output.F_G_star = result3;
	
    return output;
}

float4 CopyFromXxAndXy( VS_OUTPUT input ) : SV_Target
{
    const int3 idx = GetTexIdx(input);
    float4 result;
    const float huOut= Xx.Load(idx).g;   
	const float hvOut= Xy.Load(idx).b;
	result=float4(txState.Load(idx).r,huOut, hvOut, 0); 

    return result;
}

// CyclicReduceDx -- cyclic reduction in x direction for rhs 

// t4 or t5: coefMat: ABC matrix.
// t1: D1x: RHS in x direction 
// Output: RHS for the reduced matrix in x direction
// Runs on interior points only


//coefMat: .r=A, .g=B, .b=C, .a=not used
//D1x: .r=D
float4 CyclicReduceDx( VS_OUTPUT input ) : SV_Target  
{
    float4 result;
	const int3 idx2 = GetTexIdx(input);

	const int3 idx1 = int3((idx2.x-2)*2+2+1,idx2.y,0); //Note: -2 to skip ghost zone, +2 to add it back
	
	const float dIn = D1x.Load(idx1).g;
	const float dInLeft = D1x.Load(idx1+int3(-1,0,0)).g;
	const float dInRight = D1x.Load(idx1+int3(1,0,0)).g;
	const float aIn = coefMatx.Load(idx1).r;
	const float cIn = coefMatx.Load(idx1).b;
	const float bInLeft = coefMatx.Load(idx1+int3(-1,0,0)).g;
	const float bInRight = coefMatx.Load(idx1+int3(+1,0,0)).g;
	const float k1 = aIn/bInLeft;
	const float k2 = cIn/bInRight;
	
    const float dOut = dIn-dInLeft*k1-dInRight*k2;   // RHS coming from a cyclic reduction step.
	result =float4(0,dOut,0,0);
    return result;

}

float4 CyclicReduceDy( VS_OUTPUT input ) : SV_Target  
{
    float4 result;
	const int3 idx2 = GetTexIdx(input); //idx2 is the index for output matrix

	const int3 idx1 = int3(idx2.x,(idx2.y-2)*2+2+1,0); //idx1 is the index for input matrix. Note: -2 to skip ghost zone, +2 to add it back
	
	const float dIn = D1y.Load(idx1).b;
	const float dInLeft = D1y.Load(idx1+int3(0,-1,0)).b;
	const float dInRight = D1y.Load(idx1+int3(0,+1,0)).b;
	const float aIn = coefMaty.Load(idx1).r;
	const float cIn = coefMaty.Load(idx1).b;
	const float bInLeft = coefMaty.Load(idx1+int3(0,-1,0)).g;
	const float bInRight = coefMaty.Load(idx1+int3(0,+1,0)).g;
	const float k1 = aIn/bInLeft;
	const float k2 = cIn/bInRight;
	
    const float dOut = dIn-dInLeft*k1-dInRight*k2;   // RHS coming from a cyclic reduction step.
	result =float4(0,0,dOut,0);
    return result;

}

// CyclicSubstitute -- cyclic substitute to find the unknown (hu) in x direction

// t4 or t5: coefMat: ABC matrix.
// t1: D1x: RHS in x direction 
// t2: Xx: X (hu) in x direction 
// Output: X in the next level of substitute
// Runs on interior points only
//coefMat: .r=A, .g=B, .b=C, .a=not used
//D1x: .r=D
//Xx: .r=x

float4 CyclicSubstituteXx( VS_OUTPUT input ) : SV_Target  
{
    float4 result;
	const int3 idx1 = GetTexIdx(input);
	uint x = idx1.x;
	
	if (x%2==1){
		const int3 idx2 = int3(x/2+1,idx1.y,0); 
		result= float4(0,Xx.Load(idx2).g,0,0);;
		//result= D1x.Load(idx1).r; // Activate this line to ignore Cyclic reduction (Developer only)
	}
	else{
		const int3 idx2_1 = int3(x/2,idx1.y,0); 
		const int3 idx2_2 = int3(x/2+1,idx1.y,0); 
		const float dIn = D1x.Load(idx1).g;
		const float aIn = coefMatx.Load(idx1).r;
		const float bIn = coefMatx.Load(idx1).g;	
		const float cIn = coefMatx.Load(idx1).b;
		const float x1 = Xx.Load(idx2_1).g;
		const float x2 = Xx.Load(idx2_2).g;
		
		const float xOut = (dIn-aIn*x1-cIn*x2)/bIn;   // RHS coming from a cyclic reduction step.

		result = float4(0,xOut,0,0);
	}
	return result;
}


float4 CyclicSubstituteXy( VS_OUTPUT input ) : SV_Target  
{
    float4 result;
	const int3 idx1 = GetTexIdx(input);
	uint y=idx1.y;
	
	if (y%2==1){
		const int3 idx2 = int3(idx1.x,y/2+1,0); 
		result= float4(0,0,Xy.Load(idx2).b,0);;
		//result= D1y.Load(idx1).r; // Activate this line to ignore Cyclic reduction (Developer only)
	}
	else{
		const int3 idx2_1 = int3(idx1.x,y/2,0); 
		const int3 idx2_2 = int3(idx1.x,y/2+1,0); 
		const float dIn = D1y.Load(idx1).b;
		const float aIn = coefMaty.Load(idx1).r;
		const float bIn = coefMaty.Load(idx1).g;	
		const float cIn = coefMaty.Load(idx1).b;
		const float x1 = Xy.Load(idx2_1).b;
		const float x2 = Xy.Load(idx2_2).b;
		
		const float xOut = (dIn-aIn*x1-cIn*x2)/bIn;   // RHS coming from a cyclic reduction step.

		result = float4(0,0,xOut,0);
	}
	return result;
}

float3 solitaryWave (const float waterDepth, const float param_H, const float theta,
					const float x_zero, const float y_zero,
					float x, float y){

		const float param_K = sqrt(0.75f * abs(param_H)/pow(waterDepth,3));
		const float param_C = sqrt (sw_g * (param_H+waterDepth));
		
		const float temp_eta = param_H * 1.0f/pow(cosh(param_K * ((x - x_zero) * cos(theta) + (y - y_zero) * sin(theta))),2);
		return float3(temp_eta, (param_C * cos(theta) * temp_eta), (param_C * sin(theta) * temp_eta));
}

ADD_SOLITARY_OUTPUT AddSolitaryWave(VS_OUTPUT input)
{

    const int3 idx = int3(input.tex_idx.x, input.tex_idx.y, 0);

    // in = {w, hu, hv, _} (cell average)
    const float4 in_here = txState.Load(idx);
    float B_here = txBottom.Load(idx).b;
    
	float3 solit_wave = solitaryWave (max(in_here.r - B_here,0), sw_param_H, sw_theta, sw_x_zero, sw_y_zero, idx.x*sw_dx, idx.y*sw_dy);
	
    ADD_SOLITARY_OUTPUT output;
    output.newState = float4(in_here.r + solit_wave.r, in_here.g + solit_wave.g , in_here.b + solit_wave.b, 0);
	//output.newState = float4(in_here.r + 5.0f, in_here.g + 5 , in_here.b + 5, 0);
	//output.newState = float4(5.0f, 5 ,  5, 0);
    return output;
}


float CalcSeaLevel( float2 tex_idx)
{
    return sea_level; 
}

// Boundary condition shaders
float4 WestBoundarySolid( VS_OUTPUT input ) : SV_TARGET
{
    // mirror on real input
	const  float3 in_state_real = txState.Load(int3(4 - input.tex_idx.x, input.tex_idx.y, 0)).rgb;  
	return float4(in_state_real.r , -in_state_real.g, in_state_real.b, 0);
	//return float4(0, 0, 0, 0);
}
float4 SouthBoundarySolid( VS_OUTPUT input ) : SV_TARGET
{
    // mirror on real input
	const  float3 in_state_real = txState.Load(int3(input.tex_idx.x, 4 - input.tex_idx.y, 0)).rgb;  
	return float4(in_state_real.r, in_state_real.g, -in_state_real.b, 0);
}
float4 EastBoundarySolid( VS_OUTPUT input ) : SV_TARGET
{
    // mirror on real input
	const  float3 in_state_real = txState.Load(int3(reflect_x - input.tex_idx.x,input.tex_idx.y,0)).rgb;  
	return float4(in_state_real.r, -in_state_real.g, in_state_real.b, 0);
}
float4 NorthBoundarySolid( VS_OUTPUT input ) : SV_TARGET
{
    // mirror on real input
	const  float3 in_state_real = txState.Load(int3(input.tex_idx.x,reflect_y - input.tex_idx.y,0)).rgb;  
	return float4(in_state_real.r, in_state_real.g, -in_state_real.b, 0);
}

float SpongeFunction(float L, float x) {
	
	return 0.5 + 0.5 * cos(PI * max(L - x, 0) / L);
}

#define MIN_DEPTH 0.0001
float4 WestBoundarySponge( VS_OUTPUT input ) : SV_TARGET
{
	const int3 idx = GetTexIdx(input);
	
	const float d_here = max(MIN_DEPTH, westSeaLevel - txBottom.Load(idx).b);
	float gamma = SpongeFunction(float(westBoundaryWidth) * dx, float(idx.x) * dx)
				/ SpongeFunction(float(westBoundaryWidth) * dx, float(idx.x) * dx + sqrt(boundary_g * d_here) * boundary_dt);
	const  float4 new_state = txState.Load(idx);
	return float4 (gamma * (new_state.r-westSeaLevel) + westSeaLevel, gamma * new_state.g, gamma * new_state.b, 0);
}

float4 EastBoundarySponge( VS_OUTPUT input ) : SV_TARGET
{
	const int3 idx = GetTexIdx(input);

	const float d_here = max(MIN_DEPTH, eastSeaLevel - txBottom.Load(idx).b);
	float gamma = SpongeFunction(float(eastBoundaryWidth) * dx, float(boundary_nx - idx.x) * dx)
				/ SpongeFunction(float(eastBoundaryWidth) * dx, float(boundary_nx - idx.x) * dx + sqrt(boundary_g * d_here) * boundary_dt);
	const  float4 new_state = txState.Load(idx);
	return float4 (gamma * (new_state.r - eastSeaLevel) + eastSeaLevel, gamma * new_state.g, gamma * new_state.b, 0);
}

float4 SouthBoundarySponge( VS_OUTPUT input ) : SV_TARGET
{
	const int3 idx = GetTexIdx(input);

	const float d_here = max(MIN_DEPTH, southSeaLevel - txBottom.Load(idx).b);
	float gamma = SpongeFunction(float(southBoundaryWidth) * dy, float(idx.y) * dy)
		/ SpongeFunction(float(southBoundaryWidth) * dy, float(idx.y) * dy + sqrt(boundary_g * d_here) * boundary_dt);

	const  float4 new_state = txState.Load(idx);
	return float4 (gamma * (new_state.r-southSeaLevel) + southSeaLevel, gamma * new_state.g, gamma * new_state.b, 0);
}

float4 NorthBoundarySponge( VS_OUTPUT input ) : SV_TARGET
{
	const int3 idx = GetTexIdx(input);

	const float d_here = max(MIN_DEPTH, northSeaLevel - txBottom.Load(idx).b);
	float gamma = SpongeFunction(float(northBoundaryWidth) * dy, float(boundary_ny - idx.y) * dy)
				/ SpongeFunction(float(northBoundaryWidth) * dy, float(boundary_ny - idx.y) * dy + sqrt(boundary_g * d_here) * boundary_dt);
	const  float4 new_state = txState.Load(idx);
	return float4 (gamma * (new_state.r - northSeaLevel) + northSeaLevel, gamma * new_state.g, gamma * new_state.b, 0);
}

//  This function is no longer used. I keep it just for reference.
float calc_wavelength(float T, float d){ // The recursive dispersion relation is troublesome on the GPU, therefore we use an approx formula instead of this function.
	const float eps = 0.01f;
	float delta_l = 2.0f * eps;
	const float half_g = 9.81f/2.0f;
	
	float L1 = T * sqrt(2.0f * half_g * d);
	if(d > eps){ // checks if d != 0
		while(delta_l > eps){
		//for (int j = 0; j<50 ; j++){
			float L2 = half_g * T * T / PI * tanh(2.0f*PI*d/L1);
			delta_l = abs(L2 - L1);
			L1 = L2;
		}
	}
	return L1;
}

float calc_wavenumber_approx(float omega, float d){

	const float k = omega*omega / (boundary_g * sqrt(tanh(omega*omega * d / boundary_g ))); // Eckart's Approx.
	return k;
}

float3 sineWave (float x, float y, float t, float d, float amplitude, float period, float theta, float phase) {

	const float omega = 2*PI/period;

	const float k = calc_wavenumber_approx(omega, d);
	
	const float c = omega/k;
	const float kx = cos (theta) * x * k;
	const float ky = sin (theta) * y * k;

	const float eta = amplitude * sin(omega * t - kx - ky + phase); 
	const float hu = c * cos(theta) * eta;
	const float hv = c * sin(theta) * eta;

	return float3(eta,hu,hv);
}


float4 EastBoundarySineWave( VS_OUTPUT input ) : SV_TARGET
{
	const int3 idx = GetTexIdx(input);
	const float B_here = txBottom.Load(idx).b; // r = north, g = east, b = here
	const float d_here = max(0, eastSeaLevel - B_here);
	const float x = (idx.x - boundary_nx)* dx; // CHECKME. boundary_nx or nx?
	const float y = idx.y * dy;
	float3 result;
	if (d_here > MIN_DEPTH) // if not zero
	{
		result = sineWave(x,y,total_time,d_here,eastAmplitude_or_eta, eastPeriod_or_hu, eastTheta_or_hv, 0);
	} else {
		result = float3(0.0f, 0.0f, 0.0f);
	}
		
	return float4(result.r + eastSeaLevel, result.g, result.b, 0);
	
}

float4 WestBoundarySineWave( VS_OUTPUT input ) : SV_TARGET
{
	const int3 idx = GetTexIdx(input);
	const float B_here = txBottom.Load(idx).b; // r = north, g = east, b = here
	const float d_here = max(0, westSeaLevel - B_here);
	const float x = idx.x * dx;
	const float y = idx.y * dy;
	float3 result;
	if (d_here > MIN_DEPTH) // if not zero
	{
		result = sineWave(x,y,total_time,d_here,westAmplitude_or_eta,westPeriod_or_hu,westTheta_or_hv, 0);
	} else {
		result = float3(0.0f, 0.0f, 0.0f);
	}
	
	return float4(result.r + westSeaLevel, result.g, result.b, 0);
}

float4 NorthBoundarySineWave( VS_OUTPUT input ) : SV_TARGET
{
	const int3 idx = GetTexIdx(input);
	const float B_here = txBottom.Load(idx).b; // r = north, g = east, b = here
	const float d_here = max(0, northSeaLevel - B_here);
	const float x = idx.x * dx;
	const float y = (idx.y - boundary_ny) * dy; // CHECKME. boundary_ny or ny?
	float3 result;
	if (d_here > MIN_DEPTH) // if not zero
	{
		result = sineWave(x, y, total_time, d_here, northAmplitude_or_eta, northPeriod_or_hu, northTheta_or_hv, 0);
	} else {
		result = float3(0.0f, 0.0f, 0.0f);
	}
	return float4(result.r + northSeaLevel, result.g, result.b, 0);	
}

float4 SouthBoundarySineWave( VS_OUTPUT input ) : SV_TARGET
{
	const int3 idx = GetTexIdx(input);
	const float B_here = txBottom.Load(idx).b; // r = north, g = east, b = here
	const float d_here = max(0, southSeaLevel - B_here);
	const float x = idx.x * dx;
	const float y = idx.y * dy;
	float3 result;
	if (d_here > MIN_DEPTH) // if not zero
	{
		result = sineWave(x, y, total_time, d_here, southAmplitude_or_eta, southPeriod_or_hu, southTheta_or_hv, 0);
	} else {
		result = float3(0.0f, 0.0f, 0.0f);
	}	
	return float4(result.r + southSeaLevel, result.g, result.b, 0);	
}


float4 WestBoundaryUniformTimeSeries( VS_OUTPUT input ) : SV_TARGET
{
	return float4(westAmplitude_or_eta + westSeaLevel, westPeriod_or_hu, westTheta_or_hv, 0);
}

float4 EastBoundaryUniformTimeSeries( VS_OUTPUT input ) : SV_TARGET
{
	return float4(eastAmplitude_or_eta + eastSeaLevel, eastPeriod_or_hu, eastTheta_or_hv, 0);
}

float4 NorthBoundaryUniformTimeSeries( VS_OUTPUT input ) : SV_TARGET
{
	return float4(northAmplitude_or_eta + northSeaLevel, northPeriod_or_hu, northTheta_or_hv, 0);
}

float4 SouthBoundaryUniformTimeSeries( VS_OUTPUT input ) : SV_TARGET
{
	return float4(southAmplitude_or_eta + southSeaLevel, southPeriod_or_hu, southTheta_or_hv, 0);
}

float4 ColumnSumReduce( VS_OUTPUT input ) : SV_TARGET
{	
	const int3 idx_out = GetTexIdx(input);
	const int3 idx_in_left  = int3((idx_out.x - 2) * 2 + 2    ,idx_out.y,0); //Note: -2 to skip ghost zone, +2 to add it back
	const int3 idx_in_right = int3((idx_out.x - 2) * 2 + 2 + 1,idx_out.y,0); //Note: -2 to skip ghost zone, +2 to add it back
	
	return irregularWavesLastReduction.Load(idx_in_left) + irregularWavesLastReduction.Load(idx_in_right);
}

float4 RowSumReduce( VS_OUTPUT input ) : SV_TARGET
{
	const int3 idx_out = GetTexIdx(input);
	const int3 idx_in_down = int3(idx_out.x, (idx_out.y - 2) * 2 + 2    ,0); //Note: -2 to skip ghost zone, +2 to add it back
	const int3 idx_in_up   = int3(idx_out.x, (idx_out.y - 2) * 2 + 2 + 1,0); //Note: -2 to skip ghost zone, +2 to add it back
	
	return irregularWavesLastReduction.Load(idx_in_down) + irregularWavesLastReduction.Load(idx_in_up);
}

float4 WestBoundaryIrregularWaves( VS_OUTPUT input ) : SV_TARGET
{
	
	const float x = columnNumber * dx;

	const int rowNum = GetTexIdx(input).y;
	const float y = rowNum * dy;
	
	const float B_here = txBottom.Load(int3(columnNumber, rowNum, 0)).b; // r = north, g = east, b = here
	const float d_here = max(0, westSeaLevel - B_here);

	const int waveNum = GetTexIdx(input).x - 2; //  - 2 to skip gost cells. 
	
	float3 result;
	if (d_here > MIN_DEPTH && waveNum < numberOfWavesWest) // if d not zero
	{	
		float4 waveData = irregularWavesWest[waveNum];
		result = sineWave(x, y, total_time, d_here, waveData.r, waveData.g, waveData.b, waveData.a);
	} else {
		result = float3(0.0f, 0.0f, 0.0f);
	}

	return float4(result.r + westSeaLevel, result.g, result.b, 0);
}


float4 EastBoundaryIrregularWaves( VS_OUTPUT input ) : SV_TARGET
{
	const float x = (columnNumber - boundary_nx) * dx;

	const int rowNum = GetTexIdx(input).y;
	const float y = rowNum * dy;
    
	
	const float B_here = txBottom.Load(int3(columnNumber, rowNum, 0)).b; // r = north, g = east, b = here
	const float d_here = max(0, eastSeaLevel - B_here);

	const int waveNum = GetTexIdx(input).x - 2; //  - 2 to skip gost cells. 
	
	float3 result;
	if (d_here > MIN_DEPTH && waveNum < numberOfWavesEast) // if d not zero
	{	
		float4 waveData = irregularWavesEast[waveNum];
		result = sineWave(x, y, total_time, d_here, waveData.r, waveData.g, waveData.b, waveData.a);
	} else {
		result = float3(0.0f, 0.0f, 0.0f);
	}
	return float4(result.r + eastSeaLevel, result.g, result.b, 0);
}

float4 SouthBoundaryIrregularWaves( VS_OUTPUT input ) : SV_TARGET
{
	const int colNum = GetTexIdx(input).x;
	const float x = colNum * dx;

	const float y = rowNumber * dy;
	
	const float B_here = txBottom.Load(int3(colNum, rowNumber, 0)).b; // r = north, g = east, b = here
	const float d_here = max(0, southSeaLevel - B_here);

	const int waveNum = GetTexIdx(input).y - 2; //  - 2 to skip gost cells. 
	
	float3 result;
	if (d_here > MIN_DEPTH && waveNum < numberOfWavesSouth) // if d not zero
	{	
		float4 waveData = irregularWavesSouth[waveNum];
		result = sineWave(x, y, total_time, d_here, waveData.r, waveData.g, waveData.b, waveData.a);
	} else {
		result = float3(0.0f, 0.0f, 0.0f);
	}
	return float4(result.r + southSeaLevel, result.g, result.b, 0);
}

float4 NorthBoundaryIrregularWaves( VS_OUTPUT input ) : SV_TARGET
{
	const int colNum = GetTexIdx(input).x;
	const float x = colNum * dx;

	const float y = (rowNumber - boundary_ny)* dy;

	
	const float B_here = txBottom.Load(int3(colNum, rowNumber, 0)).b; // r = north, g = east, b = here
	const float d_here = max(0, northSeaLevel - B_here);

	const int waveNum = GetTexIdx(input).y - 2; //  - 2 to skip gost cells. 
	
	float3 result;
	if (d_here > MIN_DEPTH && waveNum < numberOfWavesNorth) // if d not zero
	{	
		float4 waveData = irregularWavesNorth[waveNum];
		result = sineWave(x, y, total_time, d_here, waveData.r, waveData.g, waveData.b, waveData.a);
	} else {
		result = float3(0.0f, 0.0f, 0.0f);
	}
	return float4(result.r + northSeaLevel, result.g, result.b, 0);
}

// GetStats shader
// This runs on a block of 4*4 pixels and returns aggregate stats
// back to the CPU
// The CPU then does the final aggregation over all the blocks.

// input: txState (t0), txBottom (t1)
// output: two RGBA float targets, containing the aggregate information.

struct GET_STATS_OUTPUT {
    float4 target0 : SV_TARGET0;    // sum(h),  sum(Bh + 0.5*h^2),  sum(h*u),  sum(h*v)
    float4 target1 : SV_TARGET1;    // sum(h*(u^2 + v^2)),  max(u^2+v^2),  max(h),  max((|u|+c)/dx, (|v|+c)/dy)
    float target2 : SV_TARGET2;    // max((u^2+v^2)/h)
};

GET_STATS_OUTPUT GetStats( VS_OUTPUT input )
{
    const int3 idx = int3(input.tex_idx.x, input.tex_idx.y, 0) * 4;

    float sum_h = 0;
    float sum_Bhh2 = 0;
    float sum_hu = 0;
    float sum_hv = 0;
    float sum_hu2v2 = 0;
    float max_u2v2 = 0;
    float max_h = 0;
    float max_cfl = 0;
    float max_f2 = 0;
    
    for (int j = 2; j < 6; ++j) {   // add 2 to avoid ghost zones
        for (int i = 2; i < 6; ++i) {
            const int3 idx2 = idx + int3(i, j, 0);
            
            const float4 state = txState.Load(idx2);
            const float w = state.r;
            const float hu = state.g;
            const float hv = state.b;

            const float3 B_here = txBottom.Load(idx2);
            const float B = B_here.b;

            const float h = max(0, w - B);
            const float c = sqrt(g * h);
            
            float u, v;
            float divide_by_h = CalcUV_Scalar(h, hu, hv, u, v);

            sum_h += h;
            sum_Bhh2 += h*(B + 0.5f * h);
            sum_hu += hu;
            sum_hv += hv;
            sum_hu2v2 += (hu * u + hv * v);

            float u2v2 = u*u + v*v;
            max_u2v2 = max(max_u2v2, u2v2);
            max_h = max(max_h, h);
            max_cfl = max(max_cfl, (abs(u) + c) * one_over_dx);
            max_cfl = max(max_cfl, (abs(v) + c) * one_over_dy);
            max_f2 = max(max_f2, u2v2 * divide_by_h);
        }
    }

    GET_STATS_OUTPUT output;
    output.target0 = float4(sum_h, sum_Bhh2, sum_hu, sum_hv);
    output.target1 = float4(sum_hu2v2, max_u2v2, max_h, max_cfl);
    output.target2 = max_f2;

    return output;
}
