///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2016, Intel Corporation
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of
// the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * MegaMol
 * Copyright (c) 2021, MegaMol Dev Team
 * All rights reserved.
 */

layout(local_size_x = 8, local_size_y = 8) in;

// edge-ignorant blur in x and y directions, 9 pixels touched (for the lowest quality level 0)
//vec2 PSNonSmartBlur( in vec4 inPos : SV_POSITION, in vec2 inUV : TEXCOORD0 ) : SV_Target
void main()
{
    vec3 inPos = gl_GlobalInvocationID;
    vec2 inUV = inPos.xy * g_ASSAOConsts.HalfViewportPixelSize;

    vec2 halfPixel = g_ASSAOConsts.HalfViewportPixelSize * 0.5f;

    vec2 centre = textureLod(g_BlurInput, inUV, 0.0 ).xy;

    vec4 vals;
    vals.x = textureLod(g_BlurInput, inUV + vec2( -halfPixel.x * 3, +halfPixel.y ), 0.0 ).x;
    vals.y = textureLod(g_BlurInput, inUV + vec2( +halfPixel.x, +halfPixel.y * 3 ), 0.0 ).x;
    vals.z = textureLod(g_BlurInput, inUV + vec2( -halfPixel.x, -halfPixel.y * 3 ), 0.0 ).x;
    vals.w = textureLod(g_BlurInput, inUV + vec2( +halfPixel.x * 3, -halfPixel.y ), 0.0 ).x;

    vec2 returnVal = vec2(dot( vals, 0.2.xxxx ) + centre.x * 0.2, centre.y);
    imageStore(g_PingPongHalfResultB, ivec2(inPos.xy), vec4(returnVal, 0.f, 0.f));
}
