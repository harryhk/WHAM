/*
use two versions to calculate F
histogram 

metadata format:
pullx  biasedlocation springk 


point

metadata format:
biasedfile 

use one generic method to unbias probability 

trajFiles. 

GPU 
use memory copy to copy data to device. 

*/


#include"wham.h"
#include<fstream>
#include<cstdio>
using namespace std;


int main(int argc, char ** argv){
	if ( argc != 15 ){
		cerr<<"usage! ./prog.exe(0) metadata(1) \n"
			"hist_min(2) hist_max(3) hist_bins(4)\n"
			"tol(5) max_iter(6) temperature(7) output(8) \n" 
			"input_traj(9) h_p(10)(histogram / point) \n"
			"pmf_min(11) pmf_max(12) pmf_bins(13) [no]gpu(14)"<<endl;
		return -1;
	}
    
	ifstream input(argv[1]);

	// things changed a little here. 
	// metadata file is pointed to files of bias. 

	double hist_min=atof(argv[2]);
	double hist_max=atof(argv[3]);
	int hist_bins=atoi(argv[4]);
	double tol=atof(argv[5]);
	int max_iteration= atoi( argv[6]);
	double temperature= atof(argv[7]);
	string calcF_method(argv[10]);

	double pmf_min=atof( argv[11] );
	double pmf_max=atof( argv[12] );
	int pmf_bins = atoi( argv[13] );

    string gpu_flag(argv[14]);

	double kt_target= temperature * k_b;
	
	hist_group  hist_g(input, tol, temperature, hist_bins, hist_min, hist_max, calcF_method, pmf_bins, pmf_min, pmf_max);
	
	if( hist_g.calcF_Method == method_hist ){
		hist_g.build_histogram();	
	}
	
	if( hist_g.calcF_Method == method_point ) {
		hist_g.read_bias();		

	}


	ifstream input_pos(argv[9]);

	hist_g.read_pos(input_pos);

	//GPU :  copy bias to device memory. 
    
    if(gpu_flag == "gpu"){
        hist_g.dev_alloc();
        hist_g.dev_cpy_data();
    }

	
	//hist_g.build_histogram( hist_min, hist_max, hist_bins);
	
	int iteration=0; 
	
	while( ! hist_g.is_converged() || iteration == 0 ) {
		hist_g.save_free();
		
		if( hist_g.calcF_Method == method_hist ){
			hist_g.wham_iteration_hist();
		}
		if( hist_g.calcF_Method == method_point) {
			if( gpu_flag == "nogpu"){
                hist_g.wham_iteration_point();
            }
            if( gpu_flag == "gpu" ){
                hist_g.dev_wham_iteration();
            }
        }


		iteration++;

		if ( hist_g.calcF_Method == method_hist && iteration % 100 == 0      ){
			printf("# Iteration  %d:   %f\n", iteration, hist_g.average_diff() );
		}

		if ( hist_g.calcF_Method == method_point && iteration      ){
			printf("# Iteration  %d:   %f\n", iteration, hist_g.average_diff() );
		}

		if(iteration == max_iteration){
			cout<< "Too many iterations !! "<<endl;
			break; 
		}

		//if( iteration % 1000 == 0 ) {
	//		hist_g.writeF();
		//}

	}

	//hist_g.calc_free( kt_target ); 

	// write F to file 
	
	hist_g.calc_free2( kt_target);

	FILE *  output;
	output = fopen( argv[8] , "w" );

	hist_g.writefile( output); 

    if(gpu_flag == "gpu"){
        hist_g.dev_free();
    }


	return 0;
}

