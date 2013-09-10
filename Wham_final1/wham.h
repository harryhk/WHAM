#include<iostream>
#include<vector>
#include<string>
#include"gsl/gsl_histogram.h"


using namespace std;

const double k_b=0.008314;
const size_t DATA_MAXN= 100000;
const size_t WIN_MAXN=40; 
const string method_hist="histogram";
const string method_point="point";

class histogram{
	public:
		histogram(): num_points(0){}

		double get_histval(int index);
	
	
		vector<double> data;
		//int first;
		//int last;
		int num_points;
};

class trajectory{
    public:

        trajectory():num_points(0){}


        vector<double> data;
        int num_points;
};


class traj_bias{
	public:

		traj_bias():num_points(0){
			data = new double * [DATA_MAXN];

			for( size_t i=0; i!= DATA_MAXN; i++ ) {
				data[i] = new double [WIN_MAXN];
			}
		}

		~traj_bias(){
			for(size_t i=0; i!= DATA_MAXN ; i++ ) {
				delete [] data[i];
			}

			delete [] data; 
		}
		
			
		
		double **  data; 
		int num_points;

	private:
		traj_bias( const traj_bias &   );

};


class hist_group{
	public:
		
		hist_group(ifstream & input, double tolerance, double temp, int histn, double histmin, double histmax, string calcF_M, int pmfn,  double pmfmin, double pmfmax );
		/* input file is very similar to grossfield	
        *    		
	    *	filepath  refz  springk
	    *	
		*/
		
		

		~hist_group(){
			if(calcF_Method == method_point ){
				delete [] traj_bias_n;
			}
			
			gsl_histogram_free( gsl_hist);
		}

		
		void build_histogram();
		void read_bias(); // read bias stored in the file_ns; 
		void read_pos_F( ifstream & input ); // read calculated F and the trajectory projected on reaction coor. // depricated, no longer this is version. 
		void read_pos( ifstream & input ); // read position without F. 


		bool is_converged(); // return true when it is converged. 
		double calc_bias(int win_b, int win_traj, int traj_p  ); // calculate the bias potential. 
		double calc_bias(int win, int bin  );
		void save_free();
		double average_diff();
		void wham_iteration_hist();
		void wham_iteration_point(); // version 2 of wham iteration, using previsouly stored bias info. 
		void norm_prob();
		void norm_gsl_hist();
		void calc_prob();
		void calc_prob2(); // calc prob with bias different than the reaction coor. 
		int  calc_free( double kt_target);
		int  calc_free2( double kt_target);
		void writefile( FILE * output);
		void writeF( FILE * output) ; 
		void writeF( ); 
		
		int num_windows;
		double tol;
		
		
		int num_bins;
		double hist_max;
		double hist_min;  // don't confuse with pmf set of parameters. 
		
		
		int pmf_num_bins; 
		double pmf_max;
		double pmf_min;
		
		vector<string> file_ns;
		vector<string> file_pos;  // file names for position files. // format trajectory is simple. single column and each file for each window.  
		vector<double> bias_locations;
		vector<double> spring_constants;
		vector<double> F;
		vector<double> F_old;
		vector<double> kT;
		
		// calcF histogram method
		vector<double> prob; 
		vector<histogram> hists;
		vector<double> hist_range; // histogram method to caclulate F. 
		
		// calcF point method 
		traj_bias *  traj_bias_n;  // point method to calculate F . 
		
		
		// output 
		vector<trajectory> trajs;  // trajectory of interested property.
		vector<trajectory> win_trajs;
		vector<double> pmf;
		vector<double> pmf_range;
		gsl_histogram * gsl_hist ; // point method to unbias data.  
		
		
		string calcF_Method; // specify using histogram for all point method to calculate F. 
			
};

