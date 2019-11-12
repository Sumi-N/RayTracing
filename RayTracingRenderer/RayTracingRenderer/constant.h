////////////////////////////////////////////
//This is a place to set a constant values//
////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//These are the constant values which be changed often for an optimization purpose
//////////////////////////////////////////////////////////////////////////////////
#define BOUNCINGTIME 2
#define RAYPERPIXELFORSHADOW 1
#define RAYPERPIXELFORGLOSSINESS 1
#define MONTECARLOGI 64

/*****These are values for bluer effect*****/
//#define BLUREFFECT
#define RAYPERPIXELFORBLUREFFECT 64

/*****These are values for antialiasing*****/
//#define ANTIALIASING
#define SAMPLEVARIENCE 0.01f

//These are the constant values which is highly unlikely to be changed
#define MAXSAMPLECOUNT 3
#define SHADOWBIAS 0.0005f
#define RAYPERSAMPLING 4
#define HALF 0.5f





//Below are script macro. Do not touch
#ifndef ANTIALIASING
#define NOANTIALIASING
#endif