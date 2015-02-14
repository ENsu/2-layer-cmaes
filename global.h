#ifndef GLOBAL_H
#define GLOBAL_H

int NFE = 0;

double best[25] = {
          -450.0 
        , -450.0 
	, -450.0 
	, -450.0
	, -310.0  //5
	,  390.0
	, -180.0
	, -140
	, -330
	, -330    //10
	,  90
	, -460
	, -130
	, -300
	,  120	  //15	    
	,  120
	,  120
	,  10
	,  10
	,  10	  //20
	,  360
	,  360
	,  360
	,  260
	,  260
};

double domainupbound[25] = {
		 100
		,100
		,100
		,100
		,100 // 5
		,100
		,600
		,32
		,5
		,5   // 10
		,0.5
		,3.14159265359
		,5
		,100
		,5   //15
		,5
		,5
		,5
		,5
		,5   //20  
		,5
		,5
		,5
		,5
		,5
};

double domainlowbound[25] = {
    		-100
		,-100
		,-100
		,-100
		,-100 // 5
		,-100
		,0
		,-32
		,-5
		,-5  // 10
		,-0.5
		,-3.14159265359
		,-5
		,-100
		,-5  //15
		,-5
		,-5
		,-5
		,-5
		,-5  //20
		,-5
		,-5
		,-5
		,-5
		, 2
};


double solupbound[25] = {
		 100
		,100
		,100
		,100
		,100
		,100
		,99999999
		,32
		,5
		,5   // 10
		,0.5
                ,3.14159265359
                ,5
                ,100
                ,5   //15
                ,5
                ,5
                ,5
                ,5
                ,5   //20
                ,5
                ,5
                ,5
                ,5
                ,99999999
		
};

double sollowbound[25] = {
    		-100
		,-100
		,-100
		,-100
		,-100
		,-100
		,-99999999
		,-32
		,-5
		,-5
		,-0.5
                ,-3.14159265359
                ,-5
                ,-100
                ,-5  //15
                ,-5
                ,-5
                ,-5
                ,-5
                ,-5  //20
                ,-5
                ,-5
                ,-5
                ,-5
                ,-99999999


};

double accuracy[25] = {
    	         1e-6
    		,1e-6
    		,1e-6
    		,1e-6
    		,1e-6
    		,1e-2 //6
		,1e-2
		,1e-2
		,1e-2
		,1e-2
		,1e-2
		,1e-2
		,1e-2
		,1e-2
		,1e-2
		,1e-2 //16
		,1e-1
		,1e-1
		,1e-1
		,1e-1
		,1e-1
		,1e-1
		,1e-1
		,1e-1
		,1e-1
};
#endif
