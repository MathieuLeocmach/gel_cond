// POVRay file exported by vtkPOVExporter
//
// +W1348 +H857

#macro LAB2RGB(Color)
    #local LAB = color Color;
    #local L = (LAB.red);
    #local A = (LAB.green);
    #local B = (LAB.blue);
    #local Y = ( L + 16 ) / 116.0;
    #local X = A / 500.0 + Y;
    #local Z = Y - B / 200.0;
    
    #local Y = ( pow(Y,3) > 0.008856 ? pow(Y,3) : ( Y - 16.0 / 116.0 ) / 7.787);
    #local X = ( pow(X,3) > 0.008856 ? pow(X,3) : ( X - 16.0 / 116.0 ) / 7.787 );
    #local Z = ( pow(Z,3) > 0.008856 ? pow(Z,3) : ( Z - 16.0 / 116.0 ) / 7.787 );
    
    #local ref_X =  95.047;
    #local ref_Y = 100.000;
    #local ref_Z = 108.883;
    #local X = ref_X * X;   //ref_X =  95.047  Observer= 2 deg Illuminant= D65
    #local Y = ref_Y * Y;   //ref_Y = 100.000
    #local Z = ref_Z * Z;   //ref_Z = 108.883

    #local X = X / 100.0;        //X = From 0 to ref_X
    #local Y = Y / 100.0;        //Y = From 0 to ref_Y
    #local Z = Z / 100.0;        //Z = From 0 to ref_Z

    #local R = X *  3.2406 + Y * -1.5372 + Z * -0.4986;
    #local G = X * -0.9689 + Y *  1.8758 + Z *  0.0415;
    #local B = X *  0.0557 + Y * -0.2040 + Z *  1.0570;

    #local R = ( R > 0.0031308 ? 1.055 * ( pow(R, ( 1 / 2.4 )) ) - 0.055 : 12.92 * R);
    #local G = ( G > 0.0031308 ? 1.055 * ( pow(G ,( 1 / 2.4 )) ) - 0.055 : 12.92 * G);
    #local B = ( B > 0.0031308 ? 1.055 * ( pow(B, ( 1 / 2.4 )) ) - 0.055 : 12.92 * B);

    //clip colors. ideally we would do something different for colors
    //out of gamut, but not really sure what to do atm.
    <min(max(R,0),1), min(max(G,0),1), min(max(B,0),1), (LAB.filter),(LAB.transmit)>
#end

#macro MSH2RGB(Color)
    #local MSH = color Color;
    #local M = (MSH.red);
    #local S = (MSH.green);
    #local H = (MSH.blue);
    #local L = M*cos(S);
    #local A = M*sin(S)*cos(H);
    #local B = M*sin(S)*sin(H);
    LAB2RGB(<L, A, B, (MSH.filter),(MSH.transmit)>)
    //LAB2RGB(<37.7, 61.9, 33.8>)
    //<L/100,A/100,B/100>
#end

#macro COLORSCALE (X)
    //#local X = mod(X, 1);
    //#local X = (X < 0 ? X+1 : X);
    #local M = (X<0.5 ? 80 + 16*X : 80 +16*(X-0.5) );
    #local S = 1.08 * 2*abs(0.5-X);
    #local H = (X<0.5 ? -1.1 - 2*X*0.561 : 1.061 - 2*(X-0.5)*0.561 );
    MSH2RGB(<M, S, H>)
#end

global_settings {
	ambient_light color rgb <1.0, 1.0, 1.0>
	assumed_gamma 2
}

background { color rgb <1, 1, 1>}

camera {
	perspective
	location <1222.115571, -816.162493, 775.509510>
	sky <-0.141295, 0.367414, 0.919262>
	right <-1, 0, 0>
	angle 30
	look_at <255.519025, 255.487210, 198.617790>
}

light_source {
	<1222.115571, -816.162493, 775.509510>
	color <0.999800, 0.999800, 0.999800>*0.5
	parallel
	point_at <255.519025, 255.487210, 198.617790>
}

light_source {
	<0.111619, 0.766044, 0.633022>
	color <1.000000, 0.972320, 0.902220>*1.5
	parallel
	point_at <0, 0, 0>
}

light_source {
	<-0.044943, -0.965926, 0.254887>
	color <0.908240, 0.933140, 1.000000>*0.5
	parallel
	point_at <0, 0, 0>
}

light_source {
	<0.939693, 0.000000, -0.342020>
	color <0.999800, 0.999800, 0.999800>*0.428372
	parallel
	point_at <0, 0, 0>
}

light_source {
	<-0.939693, 0.000000, -0.342020>
	color <0.999800, 0.999800, 0.999800>*0.428372
	parallel
	point_at <0, 0, 0>
}
union
{
cylinder
{
  <0, 0, 0>, <0, 0, 384>, 1}

cylinder
{
  <0, 512, 0>, <0, 512, 384>, 1
}
cylinder
{
  <512, 0, 0>, <512, 0, 384>, 1
}

cylinder
{
  <512, 512, 0>, <512, 512, 384>, 1
}

cylinder
{
  <0, 0, 0>, <0, 512, 0>, 1
}

cylinder
{
  <512, 0, 0>, <512, 512, 0>, 1
}

cylinder
{
  <0, 0, 384>, <0, 512, 384>, 1
}

cylinder
{
  <512, 0, 384>, <512, 512, 384>, 1
}
cylinder
{
  <0, 0, 0>, <512, 0, 0>, 1
}

cylinder
{
  <0, 512, 0>, <512, 512, 0>, 1
}

cylinder
{
  <0, 0, 384>, <512, 0, 384>, 1
}

cylinder
{
  <0, 512, 384>, <512, 512, 384>, 1
}
texture
  {
    pigment
    {
      color Black
    }
  }
}
