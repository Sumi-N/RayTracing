////////////////////////////////////////////
//This is a place to set a constant values//
////////////////////////////////////////////

/*****These are values for refrection and refraction*****/
//#define ENABLEREFLECTION
#define ENABLEREFRACTION
//#define REFLECTIONBOUNCE 2
#define REFRACTIONBOUNCE 3

/*****This is a value for gamma correction*****/
#define ENABLEGAMMA

/*****These are values for bluer effect*****/
//#define BLUREFFECT
#define RAYPERPIXELFORBLUREFFECT 1024

/*****These are values for path tracing*****/
#define ENABLEPT
#define RAYPERPT 256
#define GIBOUNCE 2

//These are the constant values which is highly unlikely to be changed
#define SHADOWBIAS 0.0005f
#define HALF 0.5f