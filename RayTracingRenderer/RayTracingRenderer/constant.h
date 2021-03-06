////////////////////////////////////////////
//This is a place to set a constant values//
////////////////////////////////////////////

/*****This is a value for gamma correction*****/
#define ENABLEGAMMA

/*****These are values for path tracing*****/
#define ENABLEPT
#define RAYPERPT 1024
#define GIBOUNCE 20

//These are the constant values which is highly unlikely to be changed
#define SHADOWBIAS 0.0005f
#define HALF 0.5f

/*****These are values for bluer effect*****/
//#define BLUREFFECT
#define RAYPERPIXELFORBLUREFFECT 16