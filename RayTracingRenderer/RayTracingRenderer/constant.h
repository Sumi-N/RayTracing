////////////////////////////////////////////
//This is a place to set a constant values//
////////////////////////////////////////////

/*****These are values for soft shadow*****/
#define RAYPERSHADOW 1

/*****These are values for refrection and refraction*****/
#define ENABLEREFLECTION
#define ENABLEREFRACTION
#define REFLECTIONBOUNCE 2
#define REFRACTIONBOUNCE 3

/*****These are values for glossiness of the object's surface*****/
//#define ENABLEREFLECTIONGLOSSINESS
//#define ENABLEREFRACTIONGLOSSINESS
#define RAYGLOSSINESS 16

/*****These are values for global illumination*****/
#define ENABLEGI
#define GIBOUNCE 2
#define RAYPERGI 16

/*****This is a value for gamma correction*****/
#define ENABLEGAMMA

/*****This is a value for multiple importance sampling*****/
#define ENABLESAMPLING

/*****These are values for bluer effect*****/
//#define BLUREFFECT
#define RAYPERPIXELFORBLUREFFECT 64

/*****These are values for antialiasing*****/
#define ENABLEAA
#define SAMPLEVARIENCE 0.001f

/*****These are values for path tracing*****/
#define ENABLEPT
#define RAYPERPT 256

//These are the constant values which is highly unlikely to be changed
#define MAXSAMPLECOUNT 3
#define SHADOWBIAS 0.0005f
#define RAYPERSAMPLING 4
#define HALF 0.5f

//Below are script macro. Do not touch
#ifndef ENABLEAA
#define NOANTIALIASING
#endif