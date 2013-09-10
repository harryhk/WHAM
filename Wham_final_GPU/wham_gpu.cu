#include"wham.h"
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<new>

void GPU_error( cudaError_t err, const char * file, int line  ){
	if( err != cudaSuccess  ){
		printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line);
		exit(-1);
	}
}

/*void Mem_error( (void *) ptr ){
	if( ptr == NULL){
		printf("Not enough memory \n");
		exit(-1);
	}
}
*/

void hist_group::dev_alloc(){
	// calc total number of bias
	
	n_bias_tot=0;

	for(int i=0; i!=num_windows; i++){
		n_bias_tot += traj_bias_n[i].num_points;
	}

	n_bias_tot *= num_windows;

	GPU_error( cudaMalloc( (void **) & dev_bias, n_bias_tot* sizeof(float) ), __FILE__, __LINE__ );

	// allocate dev_F to store Free energy 
	GPU_error( cudaMalloc( (void **) & dev_F, num_windows * sizeof( float) ), __FILE__, __LINE__ );
	GPU_error( cudaMalloc( (void **) & dev_F_old, num_windows * sizeof( float) ), __FILE__, __LINE__ );
	
	GPU_error( cudaMalloc( (void **) & dev_ni, num_windows * sizeof(int)  ) , __FILE__, __LINE__);

}

void hist_group::dev_free(){
	cudaFree(dev_bias);
	cudaFree(dev_F);
	cudaFree(dev_F_old);
	cudaFree(dev_ni);
}


void hist_group::dev_cpy_data(){
	
	if(dev_bias == NULL){
		printf("device memory not allocated in function dev_cpy_data !\n");
		exit(-1);
	}
	
	//float * dev_cpy_ptr = dev_bias; 

    float * dev_data_temp= ( float *) malloc( n_bias_tot * sizeof(float) );
    float * dev_data_temp_ptr = dev_data_temp;

	// copy the three dimensional bias[i][l][k] into a linear float array; then copy to GPU 
	for(int i=0; i!=num_windows; i++){
		
		for(int l=0; l!=traj_bias_n[i].num_points; l++){
           for(int k=0; k!=num_windows; k++){ 
                *dev_data_temp_ptr = (float ) traj_bias_n[i].data[l][k];
                dev_data_temp_ptr++;
            }
			
            
                //float * traj_il=  traj_bias_n[i].data[l];
                //GPU_error( cudaMemcpy( dev_cpy_ptr, traj_il, num_windows*sizeof(float), cudaMemcpyHostToDevice), __FILE__, __LINE__  );
                //dev_cpy_ptr += num_windows; 
		}
	}

    GPU_error( cudaMemcpy( dev_bias, dev_data_temp , n_bias_tot*sizeof(float), cudaMemcpyHostToDevice), __FILE__, __LINE__  );

    free(dev_data_temp);

	// copy number of points of each trajectory into GPU

	//int  dev_ni_temp[num_windows];
    int * dev_ni_temp = (int *) malloc( num_windows * sizeof(int)  );
	//Mem_error((void *) dev_ni_temp);

	for(int i=0; i!=num_windows; i++){
		dev_ni_temp[i]=traj_bias_n[i].num_points; 
	}

	GPU_error( cudaMemcpy(dev_ni, dev_ni_temp, num_windows * sizeof(int) , cudaMemcpyHostToDevice ) , __FILE__, __LINE__  );

    free(dev_ni_temp);
	
}


__global__ void kernel_wham( float * dev_bias, float * dev_F, float * dev_F_old, int * dev_ni, size_t n_bias_tot , int num_windows ){
	
	size_t dev_idx=blockIdx.x + blockIdx.y * gridDim.x; 
	
	if(dev_idx < n_bias_tot){
		// calculate index i, l , k in the bias ; stricly follow the paper definition of i,l,k


		int i=0;
		size_t pre_sum=0;
		for(i=0; i!=num_windows; i++){
			pre_sum += dev_ni[i] * num_windows;
			if(dev_idx < pre_sum){
				break;
			}
		}

		/*if(i == num_windows){
			//printf("wrong \n");
			exit(-1);
		} */

		// calc l 
		pre_sum -= dev_ni[i];
		size_t indx_win_i= dev_idx -  pre_sum ;
		size_t l=  indx_win_i / num_windows;
		size_t k= indx_win_i - num_windows * l ; 
		
		// numerator
		float num = dev_bias[dev_idx];

		// denumerator
		float denom=0.0;
		size_t denom_idx= pre_sum + l * num_windows; 
		for(int j=0 ; j!=num_windows; j++){
			denom += dev_bias[ denom_idx + j  ] * dev_F_old[j] * dev_ni[i]; 
		}

		dev_F[k] += num/ denom; 
	}
}

void hist_group::dev_wham_iteration(){
	
	// calculate total number of threads needed. 
	
	int block_dim=50000;
	
	if( block_dim * block_dim < n_bias_tot ){
		cerr<<"number of blocks in GPU not enough "<<endl;
		exit(-1);
	}

	dim3 grid(block_dim, block_dim);

	// copy the old F 
	float * F_temp_s =(float *) malloc ( num_windows * sizeof(float));
	float * F_init_s =(float *) malloc ( num_windows * sizeof(float));

	//Mem_error(F_temp_s);
	//Mem_error(F_init_s);

	for(int i=0; i!=num_windows; i++){
		F_temp_s[i]=(float ) F_old[i];
		F_temp_s[i]= exp( F_temp_s[i] / kT[i]  );  // exp( \beta * Fi )
		F_init_s[i]=0.0;
	}

	GPU_error( cudaMemcpy(dev_F_old, F_temp_s, num_windows * sizeof(float), cudaMemcpyHostToDevice), __FILE__, __LINE__  );
	GPU_error( cudaMemcpy(dev_F, F_init_s, num_windows * sizeof(float), cudaMemcpyHostToDevice) , __FILE__, __LINE__ );



	kernel_wham<<< grid ,1  >>>(dev_bias, dev_F, dev_F_old, dev_ni, n_bias_tot, num_windows);

	GPU_error( cudaMemcpy(F_temp_s, dev_F, num_windows * sizeof( float) , cudaMemcpyDeviceToHost), __FILE__, __LINE__ ) ; 

	for(int i=0; i!=num_windows; i++){
		F[i]=(double) F_temp_s[i];
		F[i]= - kT[i] * log(F[i]);
	}
	double F0=F[0];
	for(int i=0; i!=num_windows; i++){
		F[i] -= F0;
	}

}


