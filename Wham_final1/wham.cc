#include"wham.h"
#include<string>
#include<vector>
#include<sstream>
#include<fstream>
#include<cmath>

//double histogram::get_histval(int index){
//	if( index < first || index > last  ){
//		return 0.0;
//	}
//	else{
//		return data[index-first];
//	}
//}

hist_group::hist_group(ifstream & input, double tolerance, double temp, int histn, double histmin, double histmax, string calcF_M, int pmfn,  double pmfmin, double pmfmax ):tol(tolerance), num_bins(histn), hist_max(histmax), hist_min(histmin), pmf_num_bins(pmfn), pmf_max(pmfmax), pmf_min(pmfmin),  calcF_Method(calcF_M) {
	
	if(calcF_Method == method_hist || calcF_Method == method_point ){
		cout<< "calculation F method: "<< calcF_Method << endl ;
	}
	else{
		cerr<<"calculation F method not supported, terminate program !" <<endl;
		exit(-1);
	}

	
	string line;
	while( getline( input, line  )) {
		
		istringstream str_s(line);
		
		string file;
		double bias;
		double springk;
		
		
		str_s>> file;
		file_ns.push_back(file);
		
		if(calcF_Method == method_hist ){
			str_s>>bias;
			str_s>>springk;
			bias_locations.push_back(bias);
			spring_constants.push_back(springk);
		}
	}

	num_windows=file_ns.size();

	if ( num_windows > WIN_MAXN ) {
		cerr<< "num_windows > WIN_MAXN "<<endl;
		exit(1);
	}
	

	double kt_t= k_b * temp;  // initialize kT
	
	F.resize(num_windows, 0.0);
	F_old.resize(num_windows, 0.0);
	kT.resize(num_windows, kt_t);
	file_pos.resize(num_windows);
	trajs.resize(num_windows);
	win_trajs.resize(num_windows);
	
	if( calcF_Method == method_hist ) {
		
		prob.resize(num_bins, 0.0);
		hists.resize(num_windows);
		hist_range.resize(num_bins, 0);
		
		//pmf.resize(num_bins, 0.0);
	}

	
	//traj_bias_n.resize( num_windows ); 
	
	if( calcF_Method == method_point ) {

		traj_bias_n = new traj_bias [num_windows ];
		
	}


    gsl_hist = gsl_histogram_alloc(pmf_num_bins);
    gsl_histogram_set_ranges_uniform( gsl_hist, pmf_min, pmf_max );
    // set the range ;

    /*for(int i=0; i!= num_bins; i++ ){
        hist_range[i]= 0.5*( gsl_hist->range[i] + gsl_hist->range[i+1] );
    }
	*/
	
}	

void hist_group::build_histogram(){
	
	//int data_c=3; // data column was 3. 

	for( int i=0; i!=num_windows; i++){
		
		// open data file 
		ifstream input( file_ns[i].c_str() ) ; 
		
		int data_n=0;
		string line; 
		
		// set up histogram 
		gsl_histogram * hist= gsl_histogram_alloc(num_bins);
		gsl_histogram_set_ranges_uniform(hist, hist_min, hist_max);
		
		while( getline( input, line) ) {
			istringstream line_str(line);
			double dump; 
			double dz; 

			line_str>> dump; 
			line_str>> dump;  // 3rd column is the useful data
			line_str>> dz ;
			
			//cout<<dz <<endl;

			gsl_histogram_accumulate( hist, dz, 1.0);
			win_trajs[i].data.push_back(dz);
			//trajs[i].data.push_back(dz);

			data_n++;
		}

		win_trajs[i].num_points=data_n;

		// copy histogram  and store trajectories. 
		//cout<<"window "<< i <<endl;
		
		
		hists[i].num_points=data_n;
		//trajs[i].num_points=data_n;
		
		for(size_t j=0; j!= hist->n; j++){
			hists[i].data.push_back( hist->bin[j] );
		}

		if(i==0){
			for(size_t j=0; j!=hist->n; j++){
				hist_range[j]= 0.5* ( hist->range[j] + hist->range[j+1]  ) ;
			}
		}

		// free histogram

		gsl_histogram_free(hist);
	}
}

void hist_group::read_bias(){
	
	size_t traj_pos=0; 

	for(int i=0; i!= num_windows; i++ ) {
		 ifstream input(file_ns[i].c_str());
		 string line; 
		 traj_pos=0;
		 
		 while( getline(input, line) ){
			
			istringstream line_str(line);
			for(int j=0; j!=num_windows; j++){
				line_str >> traj_bias_n[i].data[traj_pos][j];
			}

			traj_pos ++ ;
		 }

		traj_bias_n[i].num_points = traj_pos; 

		if( traj_pos > DATA_MAXN ){
			cerr<< " traj_pos  > DATA_MAXN " <<endl;
			exit(1);
		}
	}
}
		
void hist_group::read_pos_F( ifstream & input ){
    string line;
    for( int i=0; i!= num_windows; i++ ) {
        getline( input, line ) ;
        istringstream str_s(line);
        string file;
        double free_energy;

        str_s >> file;
        str_s >> free_energy;

        file_pos[i] = file;
        F[i] = free_energy; 
    }

    // read trajectory data. 
    for(int i=0; i!=num_windows; i++ ) {
        ifstream data_input( file_pos[i].c_str() );
        
        size_t data_n =0; 
        while( getline( data_input, line) ) {
            istringstream line_str(line);
            double dz; 

            line_str >> dz ; 

            trajs[i].data.push_back(dz);

            data_n ++ ; 
        }

        trajs[i].num_points = data_n; 
    

        if( trajs[i].num_points != traj_bias_n[i].num_points || trajs[i].num_points > DATA_MAXN  ){
            cerr<< "traj and traj bias number of points no match or num_points > DATA_MAXN !! "<<endl;
		}
        

    }
}


void hist_group::read_pos( ifstream & input ){
    string line;
    for( int i=0; i!= num_windows; i++ ) {
        getline( input, line ) ;
        istringstream str_s(line);
        string file;

        str_s >> file;

        file_pos[i] = file;
    }

    // read trajectory data. 
    for(int i=0; i!=num_windows; i++ ) {
        ifstream data_input( file_pos[i].c_str() );
        
        size_t data_n =0; 
        while( getline( data_input, line) ) {
            istringstream line_str(line);
            double dz; 

            line_str >> dz ; 

            trajs[i].data.push_back(dz);

            data_n ++ ; 
        }

        trajs[i].num_points = data_n; 
    

        if( calcF_Method == method_point &&  trajs[i].num_points != traj_bias_n[i].num_points  ){
            cerr<< "traj and traj bias or hists number of points no match !! "<<endl;
        }
        
		if( calcF_Method == method_hist && trajs[i].num_points != hists[i].num_points ){
            cerr<< "traj and hists number of points no match !! "<<endl;
        }

    }
}

bool hist_group::is_converged(){
	double error=0.0;

	for(int i=0; i!=num_windows; i++){
		error=fabs( F[i] - F_old[i] ) ; 
		if ( error > tol ) return false ;
	}
	return true;
}

double hist_group::average_diff(){
	double error=0.0; 
	for(int i=0; i!=num_windows; i++){
		error += fabs( F[i] - F_old[i]  );
	}
	error /= num_windows; 
	return error;
}


double hist_group::calc_bias( int win_b, int win_traj, int traj_p){
	double loc= bias_locations[win_b];
	double spring = spring_constants[win_b]; 
	double coor= win_trajs[win_traj].data[traj_p];
	double dx= coor - loc;

	return 0.5 * dx * dx * spring; 
}

double hist_group::calc_bias( int win, int bin){
    double loc= bias_locations[win];
    double spring = spring_constants[win];
    double coor= hist_range[bin];
    double dx= coor - loc;

    return 0.5 * dx * dx * spring;
}


void hist_group::save_free(){
	for(int i=0; i!=num_windows; i++){
		F_old[i]=F[i];
		F[i]=0.0;
		
	}
	
	if( calcF_Method == method_hist ){
		for(int i=0; i!= num_bins; i++){
			prob[i]=0.0; 
		}
	}

}


/*void hist_group::wham_iteration(){
	double num, denom, bias, bf;
	
	num=0.0; 
	denom=0.0;

	for(int k=0; k!= num_windows; k++){
		for(int i=0; i!=num_windows; i++){
			for(int l=0; l!=trajs[i].num_points; l++){
				
			
				bias=calc_bias(k, i, l );
				num=exp( - bias /kT[i] ) ;
				
				denom=0.0; 
				for ( int j=0; j!= num_windows ; j++){
					bias=calc_bias( j, i, l ); 
					bf=exp( ( F_old[j] - bias ) / kT[j] ); 

					denom += bf * trajs[j].num_points; 
				}


				F[k] += num / denom; 
			}
		}
	}

	for( int i=0; i!= num_windows; i++ ) {
		F[i] = -kT[i] * log( F[i]);
	}

	double F0= F[0];

	for( int i=0; i!= num_windows; i++ ) {
		F[i] -= F0 ; 
	}
}
*/


void hist_group::wham_iteration_hist(){
    double num, denom, bias, bf ;

    for ( int i=0; i!= num_bins; i++   ){
        num=0.0;
        denom=0.0;

        for( int j=0; j!=num_windows; j++){
            num+=hists[j].data[i];

            bias = calc_bias( j , i );
            bf=exp( ( F_old[j] -bias) /kT[j]  ) ;

            denom += hists[j].num_points * bf ;
        }
        
        prob[i] = num / denom ;  

        // calculate the F 
    
        for ( int j=0; j!= num_windows; j++ ){
            bias= calc_bias( j, i  );
            bf= exp( -bias / kT[j] ) * prob[i];
            F[j] += bf;
        }
    
    }
    
    for(int j=0; j!= num_windows; j++){
        F[j] = -kT[j] * log ( F[j] );
    }
 
	//double F0=F[0];
	double F0=F[num_windows-1];

    for ( int j=0; j !=num_windows  ; j++ ){
        F[j] -= F0;
    }

    //for ( int j=num_windows-1 ; j !=num_windows  ; j++){
    //  F[j] -= F[0];
    //}
}


void hist_group::wham_iteration_point(){
	double num, denom, bias, bf; 

	int i , j , k , l ; 

	num=0.0; 
	denom=0.0;

	double ** bias_il_p=NULL;
	double * bias_ilk_p=NULL; 
	double * bias_ilj_p=NULL;
	traj_bias * traj_bias_pi = NULL; 
	traj_bias * traj_bias_pj = NULL; 

	vector<double>::iterator F_iter;
	vector<double>::iterator kT_iter; 

	for( k=0; k!= num_windows ; k++ ) {
		
		traj_bias_pi = traj_bias_n; 
		
		for( i=0; i!=num_windows; i++ ) {

			bias_il_p = traj_bias_pi ->data;   //traj_bias_n[i].data;
			
			for(  l=0; l != traj_bias_pi->num_points; l++ ) {
			
				bias_ilk_p =  (* bias_il_p )+k   ; // traj_bias_n[i].data[l]
				
				num = * bias_ilk_p;     // traj_bias_n[i].data[l][k]; 

				denom =0.0; 

				bias_ilj_p =  * bias_il_p   ; 
				traj_bias_pj = traj_bias_n;
				
				for( j=0,  F_iter=F_old.begin(), kT_iter=kT.begin()  ; j!=num_windows ; j++, F_iter ++, kT_iter ++ ) {
					

					bias =  *bias_ilj_p;  // traj_bias_n[i].data[l][j];
					
					bf = exp ( * F_iter / * kT_iter ) * bias; 
					denom += bf * ( traj_bias_pj->num_points );

					bias_ilj_p ++ ;
					traj_bias_pj ++; 
				}

				F[k] += num / denom ; 
				
				bias_il_p ++; 
			}

			traj_bias_pi ++; 
		}
	}

	for( i=0, F_iter = F.begin(), kT_iter=kT.begin(); i!= num_windows ; i++, F_iter ++, kT_iter++ ) {
		*F_iter = - *kT_iter * log( *F_iter  );
	}

	double F0=F[0];

	for( i=0, F_iter = F.begin() ; i!=num_windows ; i++ , F_iter++ ) {
		*F_iter -= F0 ;
	}
}

				


void hist_group::norm_prob(){
	double sum=0.0;
	for(int i=0; i!= num_bins; i++){
		sum+= prob[i];
	}

	for(int i=0; i!=num_bins; i++){
		prob[i] /= sum; 
	}
}

void hist_group::norm_gsl_hist(){
	double sum=gsl_histogram_sum(gsl_hist);
	gsl_histogram_scale( gsl_hist, 1.0/sum);
}


int hist_group::calc_free(double kt_target){
	int bin_min=0;
	double offset=0.0;
	double min = 1e50; 

	// normalize prob 
	calc_prob();
	norm_prob();

	for(int i=0; i!=num_bins; i++){
		pmf[i]= -kt_target * log( prob[i] ); 

		if ( pmf[i] < min ){
			min= pmf[i];
			bin_min= i; 
		}
	}
	

	offset=min; 

	for(int i=0; i!= num_bins; i++){
		pmf[i] -= offset; 
	}

	return bin_min; 
}

int hist_group::calc_free2(double kt_target){
	int bin_min=0;
	double offset=0.0;

	// normalize prob 
	calc_prob2();
	norm_gsl_hist();


	// find the none zero range of probility distribution 
	
    int index_first=0;
    int index_last=0;
    
    for(int i=0; i<=pmf_num_bins; i++){
        if( gsl_hist->bin[i] != 0 ){
            index_first = i;
            break;
        }
    }

    for( int i=pmf_num_bins-1; i>=0; i--){
        if( gsl_hist->bin[i] != 0 ){
            index_last =i;
            break;
        }
    }

    int index_n=index_last - index_first +1;

    pmf.resize( index_n, 0.0 );
    pmf_range.resize( index_n, 0.0 );




	// use the last bin as 0 profile .
	
	double min=1e50;

	for(int i=0; i!=index_n; i++){
		pmf[i]= -kt_target * log( gsl_hist->bin[i + index_first] ); 

        if ( pmf[i] < min ){
            min= pmf[i];
            bin_min= i;
        }
	}
	

	offset = min;  //pmf[num_bins-1];

	for(int i=0; i!= index_n ; i++){
		pmf[i] -= offset; 
	}

	for(int i=0; i!=index_n; i++){
		pmf_range[i]=0.5* ( gsl_hist->range[ i+ index_first ] +  gsl_hist->range[ i+ index_first +1 ] ) ;
	}


	return bin_min; 
}

void hist_group::writefile( FILE * output ) {
	for(int i=0; i!= pmf.size() ; i++){
		fprintf( output, "%8.3f    %8.3f  \n", pmf_range[i], pmf[i]  );
	}

	// out F
	fprintf(output, "# F \n  ");
	for( int i=0; i!= num_windows ; i++  ){
		fprintf(output, "#  %d   %8.3f \n ", i, F[i]);
	}

}

void hist_group::writeF( FILE * output ) {
	fprintf(output, "#   F \n");

	for(int i=0; i!=num_windows; i++ ) {
		fprintf( output, "%d    %10.5f \n", i, F[i] ); 
	}
}


void hist_group::writeF( ) {
	printf("#   F \n");

	for(int i=0; i!=num_windows; i++ ) {
		printf("%d    %10.5f \n", i, F[i] ); 
	}
}

void hist_group::calc_prob(){
	
	double num, denom , bias, bf; 

	for ( int i=0; i!= num_bins; i++   ){
        num=0.0;
        denom=0.0;

        for( int j=0; j!=num_windows; j++){
            num+=hists[j].data[i];

            bias = calc_bias( j , i ); 
            bf=exp( ( F_old[j] -bias) /kT[j]  ) ;

            denom += hists[j].num_points * bf ;
        }
        
        prob[i] = num / denom ;	
	}
}

void hist_group::calc_prob2(){

    double num, denom, bias, bf, dz, weight;

    num=1.0;

    for( int i= 0; i!= num_windows ; i++ ){

        denom = 0.0;

        for( int l=0; l != trajs[i].num_points ; l ++ ) {

            dz = trajs[i].data[l];

            for(int j=0; j!=num_windows; j++) {
                
				if( calcF_Method == method_hist ) {
					bias = calc_bias( j, i, l ); 
					bias = exp( -bias / kT[j]);
				}
				
				if( calcF_Method == method_point ) {
					bias = traj_bias_n[i].data[l][j] ;
				}
                
				
				bf = exp ( F[j] / kT[j] ) * bias;
                denom += bf * trajs[j].num_points;
            }

            weight= num/ denom;

            gsl_histogram_accumulate( gsl_hist, dz, weight);
        }
    }

    /* copy the gsl data to prob

	for(int i=0; i!= pmf_num_bins; i++){
		prob[i]= gsl_hist->bin[i];
	}
	*/
}

