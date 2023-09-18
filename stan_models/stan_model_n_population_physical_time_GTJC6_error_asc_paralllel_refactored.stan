functions{
   vector felsenstein_block(vector p_matrices, vector partials,
               array[] real genotype_tip_partials_r, array[] int total_sample_size_topology_i) {

    int total_sample_size = total_sample_size_topology_i[1];
    int topology[total_sample_size-1,3];
    int col_per_shards =  total_sample_size_topology_i[3*total_sample_size-1];
    real one_over_16 = 1.0/16.0;
    real log_lik_shard= 0.0;
    matrix[16,col_per_shards] partials_internal_nodes[total_sample_size];
    matrix[16,col_per_shards] left;
    matrix[16,col_per_shards] right;
    row_vector[16] all_ones = to_row_vector(rep_vector(1, 16));
    int site_weights_current_shard[col_per_shards] =  total_sample_size_topology_i[(3*total_sample_size):(3*total_sample_size+col_per_shards-1)];

    topology[1:(total_sample_size-1), 1] = total_sample_size_topology_i[2:(total_sample_size)];
    topology[1:(total_sample_size-1), 2] = total_sample_size_topology_i[ (total_sample_size+1):(2*total_sample_size-1)];
    topology[1:(total_sample_size-1), 3] = total_sample_size_topology_i[(2*total_sample_size):(3*total_sample_size-2)];
   

    for( n in 1:(total_sample_size-1) ) {//rows
          
          int from_left = (topology[n,2]-1)*16*col_per_shards+1;
          int from_right = (topology[n,3]-1)*16*col_per_shards+1;

          int left_matrix_idx = 16*16*(2*n-2)+1;
          int right_matrix_idx = 16*16*(2*n-1)+1;
          
          if (topology[n,2] > total_sample_size){
            left  =  partials_internal_nodes[topology[n,2]-total_sample_size];

          }
          else{

             left = to_matrix(genotype_tip_partials_r[from_left:(from_left+16*col_per_shards-1)],16,col_per_shards);

          }
      
          if (topology[n,3] > total_sample_size){

               right = partials_internal_nodes[topology[n,3]-total_sample_size];
          }
          else{
             right = to_matrix(genotype_tip_partials_r[from_right:(from_right+16*col_per_shards-1)],16,col_per_shards);

          }
          
          partials_internal_nodes[topology[n,1]-total_sample_size] = (to_matrix(p_matrices[left_matrix_idx:(left_matrix_idx+255)],16,16)*left) .* (to_matrix(p_matrices[right_matrix_idx:(right_matrix_idx+255)],16,16)*right);
        }
       
          partials_internal_nodes[total_sample_size] =  one_over_16 *partials_internal_nodes[topology[total_sample_size-1,1]-total_sample_size];
   
   log_lik_shard=  sum(to_row_vector(site_weights_current_shard).* log(all_ones* partials_internal_nodes[total_sample_size])); 
   
   return [log_lik_shard]';
   //return (to_vector(log(all_ones* partials_internal_nodes[total_sample_size])));
 }
   vector tip_clv(int observed_genotype, real seq_amp_error, real allele_dropout_error){

     vector[16] result= rep_vector(0.0, 16);
     ///*AA*    //*CC*    //*GG*/  /*TT*/    /*AC*/    /*AG*/     /*AT*/    /*CG*/     /*CT*/   /*GT*/       /*CA*/        /*GA*/    /*TA*/     /*GC*/    /*TC*/     /*TG*/
     //"A"= 1,   "C"= 2,  "G"= 3,  "T"= 4, "M"= AC=5, "R"=AG= 6,"W"=AT= 7,"S"=CG= 8,"Y"=CT= 9,"K"= GT=10, "M"= CA=5, "R"= GA=6 "W"=TA=7, "S"=GC= 8, "Y"=TC= 9, "K"=TG= 10
     //1           2         3         4         5         6         7         8          9        10         11            12        13         14          15       16
     real obs_homo_same_homo=1-seq_amp_error+0.5*seq_amp_error*allele_dropout_error;
     real obs_homo_diff_homo=(1.0/6.0)*seq_amp_error*allele_dropout_error;
     real obs_hetero_overlapped_homo=(1.0/6.0)*(1-allele_dropout_error)*seq_amp_error;
     real obs_hetero_not_overlapped_homo=0.0;
     
     
     real obs_homo_overlapped_hetero=0.5*allele_dropout_error+(1.0/6.0)*seq_amp_error-(1.0/3.0)* seq_amp_error*allele_dropout_error;
     real obs_homo_not_overlapped_hetero=(1.0/6.0)*seq_amp_error*allele_dropout_error;
     real obs_hetero_overlapped_hetero=(1.0/6.0)*(1-allele_dropout_error)*seq_amp_error;
     real obs_hetero_same_hetero=(1-allele_dropout_error)*(1-seq_amp_error);
     
     real obs_hetero_not_overlapped_hetero =0.0;

     if(observed_genotype<=4 )//homozigous
     {
           result[1:4]= rep_vector(obs_homo_diff_homo, 4);
           result[observed_genotype]=obs_homo_same_homo;
          //fill all with obs_homo_overlapped_hetero
           result[5:16] = rep_vector(obs_homo_overlapped_hetero, 12);
           if (observed_genotype==1){
             result[8] = obs_homo_not_overlapped_hetero;
             result[9] = obs_homo_not_overlapped_hetero;
             result[10] = obs_homo_not_overlapped_hetero;
             result[14] = obs_homo_not_overlapped_hetero;
             result[15] = obs_homo_not_overlapped_hetero;
             result[16] = obs_homo_not_overlapped_hetero;
           }
          else if(observed_genotype==2){
              result[6] = obs_homo_not_overlapped_hetero;
             result[7] = obs_homo_not_overlapped_hetero;
             result[10] = obs_homo_not_overlapped_hetero;
             result[12] = obs_homo_not_overlapped_hetero;
             result[13] = obs_homo_not_overlapped_hetero;
             result[16] = obs_homo_not_overlapped_hetero;
          }
          else if(observed_genotype==3){
              result[5] = obs_homo_not_overlapped_hetero;
             result[7] = obs_homo_not_overlapped_hetero;
             result[9] = obs_homo_not_overlapped_hetero;
             result[11] = obs_homo_not_overlapped_hetero;
             result[13] = obs_homo_not_overlapped_hetero;
             result[15] = obs_homo_not_overlapped_hetero;
            
          }
          else{//observed_genotype==4
              result[5] = obs_homo_not_overlapped_hetero;
             result[6] = obs_homo_not_overlapped_hetero;
             result[8] = obs_homo_not_overlapped_hetero;
             result[11] = obs_homo_not_overlapped_hetero;
             result[12] = obs_homo_not_overlapped_hetero;
             result[14] = obs_homo_not_overlapped_hetero;
            
          }
     }
     else{//heterozigous
       
       result[1:4]= rep_vector(0.0, 4);
       result[observed_genotype]=obs_hetero_same_hetero;
       
       if (observed_genotype==5)
       {
            result[1] = obs_hetero_overlapped_homo;
            result[2] = obs_hetero_overlapped_homo;
            
            result[6] = obs_hetero_overlapped_hetero;
            result[7] = obs_hetero_overlapped_hetero;
            result[14] = obs_hetero_overlapped_hetero;
            result[15] = obs_hetero_overlapped_hetero;
       }
      else if (observed_genotype==6){
            result[1] = obs_hetero_overlapped_homo;
            result[3] = obs_hetero_overlapped_homo;
            
            result[5] = obs_hetero_overlapped_hetero;
            result[7] = obs_hetero_overlapped_hetero;
            result[8] = obs_hetero_overlapped_hetero;
            result[16] = obs_hetero_overlapped_hetero;
      }
      else if (observed_genotype==7){
            result[1] = obs_hetero_overlapped_homo;
            result[4] = obs_hetero_overlapped_homo;
            
            result[5] = obs_hetero_overlapped_hetero;
            result[6] = obs_hetero_overlapped_hetero;
            result[9] = obs_hetero_overlapped_hetero;
            result[10] = obs_hetero_overlapped_hetero;
      }
      else if (observed_genotype==8){
            result[2] = obs_hetero_overlapped_homo;
            result[3] = obs_hetero_overlapped_homo;
            
            result[6] = obs_hetero_overlapped_hetero;
            result[9] = obs_hetero_overlapped_hetero;
            result[11] = obs_hetero_overlapped_hetero;
            result[16] = obs_hetero_overlapped_hetero;
      }
      else if (observed_genotype==9){
            result[2] = obs_hetero_overlapped_homo;
            result[4] = obs_hetero_overlapped_homo;
            
            result[7] = obs_hetero_overlapped_hetero;
            result[8] = obs_hetero_overlapped_hetero;
            result[10] = obs_hetero_overlapped_hetero;
            result[11] = obs_hetero_overlapped_hetero;
      }
      else if (observed_genotype==10){
            result[3] = obs_hetero_overlapped_homo;
            result[4] = obs_hetero_overlapped_homo;
            
            result[7] = obs_hetero_overlapped_hetero;
            result[9] = obs_hetero_overlapped_hetero;
            result[12] = obs_hetero_overlapped_hetero;
            result[14] = obs_hetero_overlapped_hetero;
      }
      else if (observed_genotype==11){
            result[1] = obs_hetero_overlapped_homo;
            result[2] = obs_hetero_overlapped_homo;
            
            result[8] = obs_hetero_overlapped_hetero;
            result[9] = obs_hetero_overlapped_hetero;
            result[12] = obs_hetero_overlapped_hetero;
            result[13] = obs_hetero_overlapped_hetero;
      }
      else if (observed_genotype==12){
            result[1] = obs_hetero_overlapped_homo;
            result[3] = obs_hetero_overlapped_homo;
            
            result[10] = obs_hetero_overlapped_hetero;
            result[11] = obs_hetero_overlapped_hetero;
            result[13] = obs_hetero_overlapped_hetero;
            result[14] = obs_hetero_overlapped_hetero;
      }
      else if (observed_genotype==13){
            result[1] = obs_hetero_overlapped_homo;
            result[4] = obs_hetero_overlapped_homo;
            
            result[11] = obs_hetero_overlapped_hetero;
            result[12] = obs_hetero_overlapped_hetero;
            result[15] = obs_hetero_overlapped_hetero;
            result[16] = obs_hetero_overlapped_hetero;
      }
      else if (observed_genotype==14){
            result[2] = obs_hetero_overlapped_homo;
            result[3] = obs_hetero_overlapped_homo;
            
            result[5] = obs_hetero_overlapped_hetero;
            result[10] = obs_hetero_overlapped_hetero;
            result[12] = obs_hetero_overlapped_hetero;
            result[15] = obs_hetero_overlapped_hetero;
      }
      else if (observed_genotype==15){
            result[2] = obs_hetero_overlapped_homo;
            result[4] = obs_hetero_overlapped_homo;
            
            result[5] = obs_hetero_overlapped_hetero;
            result[13] = obs_hetero_overlapped_hetero;
            result[14] = obs_hetero_overlapped_hetero;
            result[16] = obs_hetero_overlapped_hetero;
      }
      else{//observed_genotype==16
            result[3] = obs_hetero_overlapped_homo;
            result[4] = obs_hetero_overlapped_homo;
            
            result[6] = obs_hetero_overlapped_hetero;
            result[8] = obs_hetero_overlapped_hetero;
            result[13] = obs_hetero_overlapped_hetero;
            result[15] = obs_hetero_overlapped_hetero;
      }
     }
     return(result);
 }
  vector rep_vector_times(vector x, int[] times) {
    int N = rows(x);
    int Ntimes = sum(times);
    int sizetimes = size(times);

    vector[Ntimes] y;
    int pos = 1;
    for (n in 1:sizetimes) {
   
      y[pos:(pos+times[n]-1)]= rep_vector(x[n],times[n]);
      pos= pos+times[n];
   
    }
    return y;
  }
vector rep_array_times(real[] x, int[] times) {
    int N = size(x);
    int Ntimes = sum(times);
    int sizetimes = size(times);

    vector[Ntimes] y;
    int pos = 1;
    for (n in 1:sizetimes) {
  
      y[pos:(pos+times[n]-1)]= rep_vector(x[n],times[n]);
      pos= pos+times[n];
   
    }
    return y;
  }
matrix[] calculate_gtr_p_matrices(vector freqs, vector rates, vector blens){
    int bcount = rows(blens);
    matrix[16,16] pmats[bcount]; // probability matrices
    
    matrix[16,16] Q; // rate matrix
    matrix[16,16] P2 = diag_matrix(sqrt(freqs));        // diagonal sqrt frequencies
    matrix[16,16] P2inv = diag_matrix(1.0 ./ sqrt(freqs)); // diagonal inverse sqrt frequencies
    matrix[16,16] A; // another symmetric matrix
    vector[16] eigenvalues;
    matrix[16,16] eigenvectors;
    matrix[16,16] m1;
    matrix[16,16] m2;
    real s ;
    int index;
  

      matrix[16,16] R=[[0.0,  0.0,  0.0,  0.0, rates[1], rates[2], rates[3],      0.0,      0.0,     0.0,   rates[1], rates[2], rates[3],      0.0,     0.0,       0.0], /* AA */
                       [0.0,  0.0,  0.0,  0.0, rates[1],      0.0,      0.0,      1.0, rates[4],     0.0,  rates[1],      0.0,      0.0,      1.0, rates[4],      0.0],  /* CC */
                       [0.0,  0.0,  0.0,  0.0,      0.0, rates[2],      0.0,      1.0,      0.0, rates[5],      0.0, rates[2],      0.0,      1.0,      0.0, rates[5]],  /* GG */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0, rates[3],      0.0, rates[4], rates[5],      0.0,      0.0, rates[3],      0.0, rates[4], rates[5]],  /* TT */
                       [0.0,  0.0,  0.0,  0.0,      0.0,     1.0, rates[4], rates[2], rates[3],      0.0,      0.0,      0.0,      0.0, rates[2], rates[3],      0.0],   /* AC */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0, rates[5], rates[1],      0.0, rates[3],      0.0,      0.0,      0.0,      0.0,      0.0, rates[3]],  /* AG */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0, rates[1], rates[2],      0.0,      0.0,      0.0,      0.0,      0.0,      0.0],  /* AT */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0, rates[5], rates[4], rates[2],      0.0,      0.0,      0.0,      0.0, rates[4]],  /* CG */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      1.0, rates[3],      0.0,      0.0,      0.0,      0.0,      0.0],  /* CT */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[3],      0.0, rates[4],      0.0,      0.0],  /* GT */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      1.0, rates[4],      0.0,      0.0,      0.0],  /* CA */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[5], rates[1],      0.0,      0.0],  /* GA */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[1], rates[2]],  /* TA */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[5],      0.0],  /* GC */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      1.0],  /* TC */
                       [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0]]  /* TG */
                       ;
             
    R= R+R';
    s = 0.0;
    index = 1;
    Q = R * diag_matrix(freqs);
    for (i in 1:16) {
      Q[i,i] = 0.0;
      Q[i,i] = -sum(Q[i,1:16]);
      s -= Q[i,i] * freqs[i];
    }
    Q /= s;
    A = P2 * Q * P2inv;
    eigenvalues = eigenvalues_sym(A);
    eigenvectors = eigenvectors_sym(A);
    m1 = P2inv * eigenvectors;
    m2 = eigenvectors' * P2;
      for( b in 1:bcount ){
        pmats[index] = m1 * diag_matrix(exp(eigenvalues*blens[b])) * m2;
        index += 1;
      }
    return pmats;
  }
  int categorical_test_rng(vector theta){
    int selected_index;
    selected_index=categorical_rng(to_vector(theta));
    return(selected_index);
  }
  real categorical_test_lpmf(int y, vector theta){
    real result;
    result= categorical_lpmf( y |  theta);
    return(result);
  }
  real oneOnX_log(real x){
		return -log(x);
	}
  int count_valid_elements_vector(vector values){
    int count =1;
    int result;
    if (rows(values)==0){
       result=0;
    }
    else{
       while(count>=1 && count<= rows(values) && !is_nan(values[count])){
        count = count +1;
        }
     result = count -1;
     }
    return(result);
  }
  int find_integer(int x, int[] values){
     int i=1;
     int result= -1;
     while (i <= num_elements(values) && (abs(values[i] -x) >0.00001))
     {
          i=i+1;
     }
     if ((num_elements(values)==0) ||  (i>num_elements(values)))
        result=-1;
      else
         result=i;
      return(result);
  }
  int find(real x, vector values){
     int i=1;
     int result= -1;
     while (i <= rows(values) && (fabs(values[i] -x) >0.00001))
     {
          i=i+1;
     }
     if ((rows(values)==0) ||  (i>rows(values)))
        result=-1;
      else
         result=i;
      return(result);
  }
  int find_sorted(real x, vector sorted_values){
     int i=1;
     int result= -1;
     while (i <= rows(sorted_values)  && (fabs(sorted_values[i] -x) >0.00001) && sorted_values[i]< x)
     {
          i=i+1;
     }
     if ((rows(sorted_values)==0) ||  (i>rows(sorted_values)) || (i>= 1 && i<=rows(sorted_values) && sorted_values[i]> x) )
        result=-1;
      else
         result=i;
      return(result);
  }
  int[]  cumulative_sum_integer(int[] values, int N){
       int result[N];
       int current_sum=0;
       for(i in 1:N) {
            current_sum+=values[i];
            result[i]=current_sum;
       }
       return(result);
  }
  vector  difference_sorted_lists(vector sorted_list, vector sorted_sublist){
       vector[0] result;
       int idx_list=1;
       int idx_sublist=1;
       while(idx_list <= rows(sorted_list)) {
            if (sorted_list[idx_list]< sorted_sublist[idx_sublist]){
              result =append_row(result,sorted_list[idx_list] );
              idx_list= idx_list+1;
            }
            else if (sorted_list[idx_list]> sorted_sublist[idx_sublist]){
              result =append_row(result,sorted_list[idx_list] );
              idx_sublist= idx_sublist+1;
            }
            else{
              idx_list= idx_list+1;
              idx_sublist= idx_sublist+1;
            }
       }
       return(result);
  }
   int get_index_coal_time_above_MRCA_in_oldest_population_time_less_than(real time_in_oldest_population_time, vector sorted_coal_times_in_oldest_population_time){
    int idx;
    idx=1;
    while(idx <=rows(sorted_coal_times_in_oldest_population_time) && sorted_coal_times_in_oldest_population_time[idx]<time_in_oldest_population_time){
            idx= idx+1;
       }
    return(idx);
  }
  vector get_coal_times_in_oldest_population_time_less_than(real time_in_oldest_population_time, vector sorted_coal_times_in_oldest_population_time){

    int idx = get_index_coal_time_above_MRCA_in_oldest_population_time_less_than( time_in_oldest_population_time,  sorted_coal_times_in_oldest_population_time);
    vector[idx-1] result= segment(sorted_coal_times_in_oldest_population_time,1, idx-1);
    return(result);
  }

  vector get_inter_coal_times_in_oldest_population_time_less_than(real time_in_oldest_population_time, vector sorted_coal_times_in_oldest_population_time){
    vector[0] result;
    int idx;
    real inter_coal_time;
    idx=1;
    while(idx <=rows(sorted_coal_times_in_oldest_population_time) && sorted_coal_times_in_oldest_population_time[idx]<time_in_oldest_population_time){
            if (idx==1){
              inter_coal_time= sorted_coal_times_in_oldest_population_time[idx];
            }
            else{
              inter_coal_time=sorted_coal_times_in_oldest_population_time[idx]-sorted_coal_times_in_oldest_population_time[idx-1];
            }
            result =append_row(result,sorted_coal_times_in_oldest_population_time[idx] );
            idx= idx+1;
       }
    return(result);
  }

  vector get_candidate_time_MRCAs( real time_origin_in_oldest_population_time, vector sorted_coal_times_in_oldest_population_time,
                   int[,] topology, int[] number_of_tips_below, int idx_first_coal_time_above_torigin){

    real result;
    int number_candidate_MRCAs=0;
    int idx_left;
    int idx_right;
    int idx_left_below;
    int idx_right_below;
    int idx_MRCA;
    int current_pos=1;
    //int idx_first_coal_time_above_torigin=get_index_coal_time_above_MRCA_in_oldest_population_time_less_than( time_origin_in_oldest_population_time,  sorted_coal_times_in_oldest_population_time);
    vector [idx_first_coal_time_above_torigin-1] indexes_candidate_MRCAs;
    vector [idx_first_coal_time_above_torigin-1] coal_time_candidate_MRCAs;
   // print("idx_first_coal_time_above_torigin=", idx_first_coal_time_above_torigin);


    for(i in idx_first_coal_time_above_torigin:rows(sorted_coal_times_in_oldest_population_time)){
                           idx_left= topology[i, 2];
                           idx_right= topology[i, 3];
                           idx_left_below = find_integer(idx_left , topology[1:(idx_first_coal_time_above_torigin-1),1]);
                           idx_right_below = find_integer(idx_right , topology[1:(idx_first_coal_time_above_torigin-1),1]);

                           if (idx_left_below!=-1 ){
                              number_candidate_MRCAs=number_candidate_MRCAs+1;
                              indexes_candidate_MRCAs[current_pos]=idx_left;
                              coal_time_candidate_MRCAs[current_pos] =  sorted_coal_times_in_oldest_population_time[idx_left_below];
                              current_pos=current_pos+1;
                            }
                            else{

                               // if (number_of_tips_below[idx_left]>=1){
                               //    number_candidate_MRCAs=number_candidate_MRCAs+1;
                               //    indexes_candidate_MRCAs[current_pos]=idx_left;
                               //    coal_time_candidate_MRCAs[current_pos] =  0.0;
                               //    current_pos=current_pos+1;
                               //
                               // }
                            }
                          if (idx_right_below!=-1 ){
                            number_candidate_MRCAs=number_candidate_MRCAs+1;
                            indexes_candidate_MRCAs[current_pos]=idx_right ;
                            coal_time_candidate_MRCAs[current_pos] = sorted_coal_times_in_oldest_population_time[idx_right_below];
                            current_pos=current_pos+1;
                           }
                           else{

                               // if (number_of_tips_below[idx_right]>=1){
                               //    number_candidate_MRCAs=number_candidate_MRCAs+1;
                               //    indexes_candidate_MRCAs[current_pos]=idx_right;
                               //    coal_time_candidate_MRCAs[current_pos] =  0.0;
                               //    current_pos=current_pos+1;
                               //
                               // }
                            }
                       }

    //idx_MRCA = categorical_rng(rep_vector(1.0/number_candidate_MRCAs, number_candidate_MRCAs));
    // print("number_candidate_MRCAs 2=", number_candidate_MRCAs);
    // print("indexes_candidate_MRCAs 2=", indexes_candidate_MRCAs);
    // print("coal_time_candidate_MRCAs 2=", coal_time_candidate_MRCAs);

    return(coal_time_candidate_MRCAs);
  }
  vector get_list_coal_times_in_oldest_population_time_below_time_origin(int idx_pop, int order, int total_sample_size,
                       int number_internal_nodes_under_MRCA, real coal_time_MRCA_in_oldest_population,
                       vector inmigrants_MRCA_times_in_model_time_oldest_population,
                       real time_origin_in_oldest_population_time, vector sorted_coal_times_in_oldest_population_time,
                       int[,] topology, int[] number_of_coalescent_times_below){
    vector[rows(sorted_coal_times_in_oldest_population_time)] subtree_coal_times_below;
    int idx;
    int idx_subtree = find(coal_time_MRCA_in_oldest_population, sorted_coal_times_in_oldest_population_time);
    int nodes_to_visit[idx_subtree];
    int idx_MRCA;
    int left_child;
    int right_child;
    real coal_time_above;
    int pos_subtree_coal=1;
    int pos_subtree=1;

    //idx_subtree = find(coal_time_MRCA_in_oldest_population, sorted_coal_times_in_oldest_population_time);
    nodes_to_visit[pos_subtree] = topology[idx_subtree, 1];
    pos_subtree=pos_subtree+1;
    //print("outside idx_subtree =", idx_subtree);
    //print("outside nodes_to_visit=", nodes_to_visit);

    while(idx_subtree>=1 &&  pos_subtree_coal <= number_internal_nodes_under_MRCA){
      if (order > 1)
      {
          if (rows(inmigrants_MRCA_times_in_model_time_oldest_population) >0   &&
               find(sorted_coal_times_in_oldest_population_time[idx_subtree], inmigrants_MRCA_times_in_model_time_oldest_population)== -1
               && find_integer(topology[idx_subtree, 1],nodes_to_visit)!=-1)
               {
               subtree_coal_times_below[pos_subtree_coal]= sorted_coal_times_in_oldest_population_time[idx_subtree];
               if ( topology[idx_subtree, 2] > total_sample_size){
                  nodes_to_visit[pos_subtree] = topology[idx_subtree, 2];
                  pos_subtree=pos_subtree+1;
               }
               if ( topology[idx_subtree, 3] > total_sample_size){
                   nodes_to_visit[pos_subtree] = topology[idx_subtree, 3];
                   pos_subtree=pos_subtree+1;
               }
               pos_subtree_coal= pos_subtree_coal+1;
               }
      }
      else
      {
         if ( find_integer(topology[idx_subtree, 1],nodes_to_visit)!=-1)
               {
               subtree_coal_times_below[pos_subtree_coal]= sorted_coal_times_in_oldest_population_time[idx_subtree];
                if ( topology[idx_subtree, 2] > total_sample_size){
                  nodes_to_visit[pos_subtree] = topology[idx_subtree, 2];
                  pos_subtree=pos_subtree+1;
               }
               if ( topology[idx_subtree, 3] > total_sample_size){
                   nodes_to_visit[pos_subtree] = topology[idx_subtree, 3];
                   pos_subtree=pos_subtree+1;
               }
               pos_subtree_coal= pos_subtree_coal+1;
               }
      }

      //print("nodes_to_visit=", nodes_to_visit);
      idx_subtree= idx_subtree-1;
      }

    //for(i in (pos_subtree):rows(sorted_coal_times_in_oldest_population_time))// padding
    //   subtree_coal_times_below[i]=positive_infinity();

    // print("subtree_coal_times_below=",segment(subtree_coal_times_below, 1, pos_subtree_coal-1));

    //subtree_coal_times_below = sort_asc(subtree_coal_times_below);

    return(sort_asc(segment(subtree_coal_times_below, 1, pos_subtree_coal-1)));
  }
    int[] get_tips_positions_below_time(int order, int number_tips_under_MRCA, real coal_time_MRCA_in_oldest_population,
                     vector inmigrants_MRCA_times_in_model_time_oldest_population,
                     real time_origin_in_oldest_population_time, vector sorted_coal_times_in_oldest_population_time,
                     int[,] topology, int total_sample_size){

    int result[total_sample_size];
    int idx_left;
    int idx_right;
    int idx_subtree= find(coal_time_MRCA_in_oldest_population, sorted_coal_times_in_oldest_population_time);
    int tips_positions[idx_subtree+1];
    int nodes_to_visit[idx_subtree+1];
    int idx_MRCA;
    int left_child;
    int right_child;
    real coal_time_above;
    int pos_subtree=1;
    int indexes_nodes_below_torigin[idx_subtree-1] = topology[1:(idx_subtree-1),1];
    int pos=1;

    nodes_to_visit[pos_subtree] = topology[idx_subtree, 1];
    pos_subtree=pos_subtree+1;

    while(idx_subtree>=1 && pos <= number_tips_under_MRCA){
      if (order > 1)
      {
          if (rows(inmigrants_MRCA_times_in_model_time_oldest_population) >0   &&
               find(sorted_coal_times_in_oldest_population_time[idx_subtree], inmigrants_MRCA_times_in_model_time_oldest_population)== -1
               && find_integer(topology[idx_subtree, 1],nodes_to_visit)!=-1)
            {
               if ( topology[idx_subtree, 2] > total_sample_size){
                  nodes_to_visit[pos_subtree] = topology[idx_subtree, 2];
                  pos_subtree=pos_subtree+1;
               }
               else{
                  tips_positions[pos]=topology[idx_subtree, 2];
                  pos=pos+1;
               }
               if ( topology[idx_subtree, 3] > total_sample_size){
                   nodes_to_visit[pos_subtree] = topology[idx_subtree, 3];
                   pos_subtree=pos_subtree+1;
               }
               else{
                  tips_positions[pos]=topology[idx_subtree, 3];
                  pos=pos+1;
               }
            }
      }
      else
      {
         if ( find_integer(topology[idx_subtree, 1],nodes_to_visit)!=-1)
            {
              if ( topology[idx_subtree, 2] > total_sample_size){
                  nodes_to_visit[pos_subtree] = topology[idx_subtree, 2];
                  pos_subtree=pos_subtree+1;
               }
               else{
                  tips_positions[pos]=topology[idx_subtree, 2];
                  pos=pos+1;
               }
               if ( topology[idx_subtree, 3] > total_sample_size){
                   nodes_to_visit[pos_subtree] = topology[idx_subtree, 3];
                   pos_subtree=pos_subtree+1;
               }
               else{
                  tips_positions[pos]=topology[idx_subtree, 3];
                  pos=pos+1;
               }
            }
      }

      //print("nodes_to_visit=", nodes_to_visit);
      idx_subtree= idx_subtree-1;
      //idx_subtree= max(idx_left, idx_right);
    }
    result= append_array(tips_positions[1:(pos-1)], rep_array(-1, total_sample_size-pos+1));
    return(result);
  }
  vector get_torigins_of_inmigrant_populations(int pos,int order, int N, vector sorted_previous_torigins_in_model_time_oldest_population,
                  vector torigins_in_model_time_oldest_population,
                  int[] indexes_father_populations){

      vector[order-1] result;
      int k=1;
      for(i in 1:N){
         if (i!=pos){//an inmigrant is another population
              if (indexes_father_populations[i]==pos){
                   if (k>(order-1)){
                     reject("index out bounds; found k=", k);
                   }
                  result[k]=torigins_in_model_time_oldest_population[i];
             }
             else{
                if (k>(order-1)){
                     reject("index out bounds; found k=", k);
                   }
               result[k]=-1.0;
             }
           k=k+1;
         }
      }
     return(result);
  }
  int population_of_coalescent_event(real x, vector grouped_coalescent_times, vector cumulative_number_coal_times){
     int i=1;
     int result= -1;
     int idx_population =1;
     while ((i <= rows(grouped_coalescent_times)) && (grouped_coalescent_times[i] !=x)){
          if (i==cumulative_number_coal_times[idx_population]){
            idx_population = idx_population +1;
          }
          i=i+1;
     }
     if ((rows(grouped_coalescent_times)==0) || (i>rows(grouped_coalescent_times)))
        result=-1;
      else
         result=idx_population;
      return(result);
  }
  real model_time_to_standard_time(real t,
         real torigin,
         real delta,
         real K)
         {
          real model_time=0.0;
          real a = exp(delta * t) - 1.0;
          real b = 1.0 - exp(-1.0 * delta * torigin);
          real c = 1.0 - exp(-1.0 * delta * (torigin - t));
          real d;
          real numerator;
          real  denominator;
          if (K==0){
            model_time = a * b / (delta * c);
            return(2.0*model_time/ delta);
          }
         else{
           c = exp(-1.0 * delta * torigin);
           d = 1.0 -c;
           a = ((K / delta) * d) - c;
           b = 1.0 - (K/ delta)*d;
           numerator =(a*exp(delta*t)+b);
           denominator = 1-exp(-delta*(torigin-t));
           numerator =(a*exp(delta*t)+b);
           denominator = 1-exp(-delta*(torigin-t));

           model_time = (2.0 / K )*log(numerator/denominator);
           return(model_time);
            }
         }
    real standard_time_to_model_time(real V,
         real torigin,
         real delta,
         real K)
         {
        real a;
        real b;
        real c;
        real d;
        real e;
        real x;
        real firstTerm;
        real secondTerm = 1.0 / delta;
        real thirdTerm;
        real numerator;
        real denominator;
        real StandardTimeG;
        if (K==0){
           a =  exp(-1.0 * delta * torigin);
           b = (1 - a) * (1 - a) * (1.0 / a);
           c = 1 - a;
           d = (1 - a) * (1.0 / a);
           e = V + d;
           thirdTerm = log(1 - b / e);
           thirdTerm = log(1 - ((1 - a) * (1 - a) * (1.0 / a)) / (V * delta + (1 - a) * (1.0 / a)));
           thirdTerm = log(1 + delta * V - a) - log(1 + (delta * V - 1) * a);
            StandardTimeG = secondTerm * thirdTerm;

        if ( (1 + delta * V - a) <= 0 ||   (1 + (delta * V - 1)*a ) <= 0 ) // do approximation if required
        {
            StandardTimeG = 0.0;
            firstTerm = 0.0;
            secondTerm = 0.0;
            thirdTerm = 0.0;
            a = 0.0;
            b = 0.0;
            c = 0.0;
            d = 0.0;
            e = 0.0;
            a = 1 / delta;
            b = log(1 + delta * V);
            firstTerm = a * b;
            d = (V * V * delta * exp(-1.0 * delta * torigin)) / (1 + V);
            secondTerm =  d;
            StandardTimeG = firstTerm - secondTerm;
           }
       }
       else {
          x = exp(K*V /2.0);
          c = exp(-1.0 * delta * torigin);
          d= 1.0 -c;
          a = ((K/delta) * d ) - c;
          b = 1.0 - (K/ delta)*d;
          numerator = x -b;
          denominator = (x*c) +a;
          thirdTerm= log(numerator/denominator);
          StandardTimeG = secondTerm * thirdTerm;
          }
        return StandardTimeG;
         }
      vector  from_model_time_to_Kingman_coalescent_time(vector coal_times_model_time, real torigin,
         real delta,
         real K){

        vector [rows(coal_times_model_time)] result;
        for(i in 1:rows(coal_times_model_time)){
           if ( !is_nan(coal_times_model_time[i]))
               result[i]= model_time_to_standard_time(coal_times_model_time[i], torigin, delta, K );
           else
              result[i]=-1;
        }
        return(result);
      }
  real log_lik_no_coalescent_between_times(real from,  real to, int number_alive_ind, real torigin, real delta, real  K){
    real  result =  0.0;
    real cumulative_to;
    real cumulative_from;
    if (number_alive_ind>1){
      cumulative_to = model_time_to_standard_time(to, torigin, delta, K);
      cumulative_from = model_time_to_standard_time(from, torigin, delta, K);
      result= log(choose(number_alive_ind,2)) -1.0 *choose(number_alive_ind,2)*(cumulative_to-cumulative_from);
    }
    return(result);
  }
  real conditionalDensityTOrigin_rng(real delta, int sample_size){
        real number_ancestors_population_when_sample_minus_one;
        real U;
        real torigin;
        real torigin_in_model_time_oldest_population;
        real term;
        number_ancestors_population_when_sample_minus_one = 2*sample_size;
        U = gamma_rng(number_ancestors_population_when_sample_minus_one+1, 1);
        term = 1.0- delta /(U+delta);
        if (term <= 0){
          reject("term must not be negative; found term=", term);
        }
        else{
            torigin= -(1.0 / delta)*log(term);
        }
        return(torigin);
  }
  real conditionalDensityTOrigin_pdf( real y, real delta,  int sample_size) {
      real result;
      real term1;
      real term2;
      real term3;
      real partial;
      term1 = exp(-1.0*delta*y);
      term2 = delta * term1;
      term3 = 1.0-term1;
      partial = delta * term2 /(term3 * term3);
      partial = partial * exp( -1 * term2/term3);

      result  = partial;
      return(result);
   }
  real conditionalDensityTOrigin_lpdf( real y, real delta,  int sample_size) {
      real result;
      result = log(conditionalDensityTOrigin_pdf(y, delta,   sample_size));
      return(result);
   }
  real conditionalDensityTOrigin_cdf(real y, real delta,  int sample_size) {
     real term1;
     real term2;
     real term3;
     real result;
     term1 = exp(-1.0 * delta * y);
     term2 = delta * term1;
     term3 = 1.0-term1;
     result = exp( -1.0 * term2/term3);

     return(result);
  }
    real conditionalDensityTOrigin_lcdf(real y, real delta,  int sample_size) {
     real term1;
     real term2;
     real term3;
     real result;
     term1 = exp(-1.0 * delta * y);
     term2 = delta * term1;
     term3 = 1.0-term1;
     result =  -1.0* term2/term3;
    return(log(result));
 }
  real conditionalDensityTOrigin_lccdf(real y, real delta,  int sample_size) {
    real term1;
    real term2;
    real term3;
    real result;
    result = conditionalDensityTOrigin_cdf(y, delta, sample_size );
    return(log(1-result));
 }
  real log_h(real t,
         real torigin,
         real delta,
         real K){

    real a = 1.0 - exp(-1.0 * delta * (torigin - t));
    real first_term = 2.0 * log(a);
    real second_term = -1.0 * delta * t;
    real third_term = exp(delta * t);
    real above_term = first_term + second_term;
    real b = 1.0 - exp(-1.0 * delta * torigin);
    real  below_term = 2.0 * log(b);
    real extra_term = log(1+ (K /delta)*b*(third_term -1.0) / a );
    real logH;
    if (K==0)
        logH = above_term - below_term;
    else
        logH =  above_term - below_term + extra_term;
    return(logH);
  }
  real log_prob_father_population(int order, int index_oldest_population, real father_torigin, real father_delta,
                                 real x_father, real torigin, real delta, real x, vector deltas,
                                 vector torigins_in_model_time,
                                vector torigins_in_model_time_oldest_population_par, real K, int N, vector pop_sizes_proportion){
                                 // , simplex pop_sizes_proportion
         real log_numerator=0.0;
         real log_denominator=0.0;
         real logh;
         real temp;
         real result;
         //print("log_prob_father_population");
         //vector[N] ordered_torigins_in_model_time_oldest_population=sort_asc(torigins_in_model_time_oldest_population_par);
        // order_father_pop = find((x_father/pop_sizes_proportion[index_oldest_population]) * father_torigin, ordered_torigins_in_model_time_oldest_population );
         if (order != N ){

                for(i in (order+1):N){
                  int pos = i;
                  log_denominator = log_denominator + log( pop_sizes_proportion[pos]);
                  temp= torigin * x / pop_sizes_proportion[pos];
                  logh=log_h(temp, torigins_in_model_time[pos], deltas[pos], K);
                  log_denominator = log_denominator  + temp;
                  
                  //print("log_denominator", log_denominator);
                }
         }

         log_numerator = log_numerator + log( x_father);
         temp= torigin * x / x_father;
         logh=log_h(temp, father_torigin, father_delta, K);
         log_numerator = log_numerator  + temp;
         
         //print("log_numerator", log_numerator);
         result= log_numerator - log_denominator;
         //print("result", result);
        return(result);
  }
   real log_likelihood_subpopulation(int order, int sample_size, real K, real delta, real torigin, int[] positions,
                                     vector sorted_coal_time_with_previous_torigins_in_model_time, vector child_torigins_in_model_time){

    real log_lik = 0.0;
    real current_time;
    int  alive_cells;
    real term_only_after_first_coal_event;
    real zero =0.0;
    int currentCoalescentEventInThisEpoch=-1;
    vector[rows(sorted_coal_time_with_previous_torigins_in_model_time)+1] padded_coalescent_times =
                                                  append_row(zero, sorted_coal_time_with_previous_torigins_in_model_time);
    real last_event_time_before_migration=0.0;

    alive_cells = sample_size;//alive_cells goes from sample_size downto 2
    for(j in 1:rows(padded_coalescent_times)){
      current_time = padded_coalescent_times[j];

      //print("current_time=",current_time);

      if (fabs(current_time - torigin)<0.00000001)
        return log_lik;

      if (current_time < 0.0) continue;//those are not real inmigrants

      if (rows(child_torigins_in_model_time) > 0 && find(current_time, child_torigins_in_model_time)!= -1){//= -1
        
    
        if (alive_cells >1){
            log_lik+= log_lik_no_coalescent_between_times( last_event_time_before_migration,  current_time,  alive_cells,  torigin,  delta,  K);
            alive_cells= alive_cells+1;
            currentCoalescentEventInThisEpoch=1;
            last_event_time_before_migration= current_time;
            }
      }
      else{
     
        if (alive_cells >1) {
            log_lik = log_lik + log(alive_cells*(alive_cells-1)/2.0);
            log_lik = log_lik + log(2.0)-log_h(current_time,torigin, delta, K);
            if (fabs(current_time)<0.00001){
              term_only_after_first_coal_event=0.0;
            }
            else{
                   term_only_after_first_coal_event= model_time_to_standard_time(last_event_time_before_migration,torigin, delta, K);
              }
             log_lik = log_lik - (alive_cells*(alive_cells-1)/2.0)*(model_time_to_standard_time(current_time,torigin, delta, K)-
                  term_only_after_first_coal_event);
              alive_cells= alive_cells-1;
              currentCoalescentEventInThisEpoch= currentCoalescentEventInThisEpoch+1;
              last_event_time_before_migration= current_time;
            }
         }
      }
      //print("log_lik=",log_lik);
      return(log_lik);
    }
    vector structured_coalescent_rng(vector deltas,vector torigins_in_model_time, vector torigins_in_model_time_oldest_population,
                                   vector pop_sizes_proportion, int[] sample_sizes, int[] number_internal_nodes,
                                   int[] indexes_father_populations, int[] pos_parent_child_pop_MRCA_in_parent_pop, int[] cum_internal_nodes, 
                                   int[] number_child_pops,
                                   int[] child_pop_indexes,
                                   real K ){
    int total_sample_size=sum(sample_sizes);
    int N = rows(deltas);
    int current_pos=1;
    int current_pos_migrations=1;
    int number_inmigrants=1;
    int index_next_immigrants=1;
    int global_index_next_immigrants=1;
    //counting torigin as immigrant
    int  indexes_inmigrants = 1;
    vector[total_sample_size-1] coal_event_times_model_time_oldest_population;

    coal_event_times_model_time_oldest_population=rep_vector(torigins_in_model_time_oldest_population[1]+1,total_sample_size-1);

    print("inside structured_coalescent_rng=",sample_sizes);
    for(i in 1:num_elements(sample_sizes)){
        int sample_size= sample_sizes[i];
        int current_sample_size= sample_size;
        int index_father_population=0;
        current_pos=1;
        if (i==1){
            vector[sample_size-1] number_ancestors_sample;
            vector[sample_size-1] rate_exp;
            vector[sample_size-1] times;
            //vector[sample_size-1] reverse_times;
            vector[sample_size-1] cum_sum_times;

            for(j in 2:sample_size){
               number_ancestors_sample[j-1]=j;
               rate_exp[j-1]=0.5*j*(j-1);
               times[j-1] = exponential_rng(rate_exp[j-1]);
            }
            //print("times=",times);
            for(j in 1:(sample_size-1)){
              cum_sum_times[j] = sum(times[j:(sample_size-1)]);
              coal_event_times_model_time_oldest_population[ sample_size-j]= standard_time_to_model_time(cum_sum_times[j],
                                                                                                 torigins_in_model_time[i],
                                                                                                 deltas[i],
                                                                                                 K);

            }
            //print("i=1 coal_event_times_model_time=", coal_event_times_model_time_oldest_population);
            coal_event_times_model_time_oldest_population[1:(sample_size-1)] = (pop_sizes_proportion[i]/pop_sizes_proportion[N]) *
                                                     (coal_event_times_model_time_oldest_population[1:(sample_size-1)]);
            coal_event_times_model_time_oldest_population[sample_size]= torigins_in_model_time_oldest_population[i];
            current_pos = sample_size+1;
            //print("coal_event_times_model_time_oldest_population=", coal_event_times_model_time_oldest_population);
            index_father_population=2;
        }
    }
    return(coal_event_times_model_time_oldest_population);
  }
  real structured_coalescent_lpdf(vector coal_event_times_model_time, vector deltas, vector torigins_in_model_time,
                              vector torigins_in_physical_time,  vector pop_sizes_proportion, int[] sample_sizes, int[] number_internal_nodes,
                               int[] indexes_father_populations, int[] pos_parent_child_pop_MRCA_in_parent_pop, int[] cum_internal_nodes, 
                               int[] number_child_pops,
                               int[] child_pop_indexes,
                               real K){

    real log_lik=0.0;
    int cum_sample_size=0;
    int cum_number_coal=0;
    int first_position_coal_event=1;
    int N= rows(deltas);
    int pos_next_torigin = 0;
    int cum_sum_number_child_pops[N] = cumulative_sum_integer(number_child_pops, N);
    int k=1;
    int pos_idx_child_pop=1;

    while(cum_sum_number_child_pops[k]==0){
      k=k+1;

    }

    cum_sample_size=0;

    for(i in 1:num_elements(sample_sizes)){
         int order=i;
         int index_father_population = indexes_father_populations[i];
         int number_events= number_internal_nodes[i];
         int positions[sample_sizes[i]] = rep_array(1,sample_sizes[i] );
         vector[number_events+number_child_pops[i]] current_coal_event_times_with_immigrants;
         vector[number_events+number_child_pops[i]] sorted_current_coal_event_times_with_immigrants;
         vector[number_child_pops[i]] child_pops_inmigration ;
         vector[number_child_pops[i]] sorted_child_pops_inmigration ;
         int start_idx=1;
         int start_idx_number_pop=1;
         pos_idx_child_pop=1;
         
         //print("k=",k);
         if (i>=2){
          start_idx = cum_internal_nodes[i-1]+1;
          start_idx_number_pop= k;
         }
         
         current_coal_event_times_with_immigrants[1:number_events] = (coal_event_times_model_time[ start_idx:(cum_internal_nodes[i])]) ;

         if (number_child_pops[i]>0){

   
            for(j in 1:num_elements(indexes_father_populations)){

               if(indexes_father_populations[j]==i){

    
                  
                   child_pops_inmigration[pos_idx_child_pop]= (1.0/pop_sizes_proportion[i])*torigins_in_physical_time[j];
                   
                    //print(" child_pops_inmigration[pos_idx_child_pop] con physical time",child_pops_inmigration[pos_idx_child_pop]);
                   child_pops_inmigration[pos_idx_child_pop]= (pop_sizes_proportion[j]/pop_sizes_proportion[i])*torigins_in_model_time[j];

                     
                     if(child_pops_inmigration[pos_idx_child_pop]> current_coal_event_times_with_immigrants[pos_parent_child_pop_MRCA_in_parent_pop[j]]){

                    //print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!wrong position of time of origin of child pop");
                      }
                    //}
                 // }
                  

                   pos_idx_child_pop=pos_idx_child_pop+1;
               }

            }
           
    
         }
           current_coal_event_times_with_immigrants[(number_events+1):(number_events+number_child_pops[i])]=child_pops_inmigration;
           
           sorted_current_coal_event_times_with_immigrants= sort_asc(current_coal_event_times_with_immigrants);
           sorted_child_pops_inmigration  = sort_asc(child_pops_inmigration);
         

           log_lik= log_lik+ log_likelihood_subpopulation(i, sample_sizes[i], K, deltas[i], torigins_in_model_time[i], positions,
                                      sorted_current_coal_event_times_with_immigrants,
                                      sorted_child_pops_inmigration
                                      //rep_vector(1,1)
                                      );

        

             if (i<(N-1) && index_father_population>0 && N>2){
              log_lik+= log_prob_father_population(order, N, torigins_in_model_time[indexes_father_populations[i]], deltas[indexes_father_populations[i]],
                                           pop_sizes_proportion[indexes_father_populations[i]],
                                           torigins_in_model_time[i], deltas[i], pop_sizes_proportion[i], deltas,  torigins_in_model_time, torigins_in_physical_time,K,N,  pop_sizes_proportion);
            }

    }
    return(log_lik);
  }
/**
   * Return the simplex corresponding to the specified free vector.
   * A simplex is a vector containing values greater than or equal
   * to 0 that sum to 1.  A vector with (K-1) unconstrained values
   * will produce a simplex of size K.
   *
   * The transform is based on a centered stick-breaking process.
   */
  vector simplex_constrain(vector y) {

    vector[rows(y)+1] x;
    int Km1;
    real stick_len;


    Km1 = rows(y);
    stick_len = 1.0;
    for (k in 1:Km1) {
      real z_k;
      z_k = inv_logit(y[k] - log(Km1 - k + 1));
      x[k] = stick_len * z_k;
      stick_len = stick_len - x[k];
    }
    x[Km1+1] = stick_len;

    return x;
  }

  /**
   * Return the log absolute Jacobian determinant of the simplex
   * transform defined in simplex_constrain().
   */
  real simplex_constrain_lj(vector y) {

    real lj;
    int Km1;
    real stick_len;
    
    lj = 0.0;
    Km1 = rows(y);
    stick_len = 1.0;
    for (k in 1:Km1) {
      real adj_y_k;
      adj_y_k = y[k] - log(Km1 - k + 1);
      lj = lj + log(stick_len) - log1p_exp(-adj_y_k) - log1p_exp(adj_y_k);
      stick_len = stick_len * (1.0 - inv_logit(adj_y_k));
    }
    
    return lj;
  } 
   vector simplex_constrain_new(vector y) {
    real lj = 0;
    vector[rows(y)+1] x;
    int Km1 = rows(y);
    real stick_len = 1.0;



    for (k in 1:Km1) {
      real eq_share = -log(Km1 - k+1);
      real adj_y_k = y[k] + eq_share;      
      real z_k = inv_logit(adj_y_k);
      x[k] = stick_len * z_k;
      lj = lj + log(stick_len) - log1p_exp(-adj_y_k) - log1p_exp(adj_y_k);      
      stick_len = stick_len - x[k];
    }
    x[Km1+1] = stick_len;

    return append_row(lj,x);
  }
}
data{
  real<lower=0.0, upper=2.0> K;//parameter for the model M_K
  int<lower=2> N;//number of populations
  int<lower=2> total_sample_size;
  int<lower=0> L;// alignment length
  int<lower=0> ids_MRCA[N];
  int<lower=0> sample_sizes[N];
  int<lower=1> number_internal_nodes[N];
  int indexes_father_populations[N];
  //for example [0,1,1,2] means the population 2 and 3 have father population 1(oldest population from present time)
  //and population 4 has parent population 2
  int number_child_pops[N];
  int child_pop_indexes[sum(number_child_pops)];
  int  pos_parent_child_pop_MRCA_in_parent_pop[N];
  
  int<lower=0,upper=16> genotype_tipdata[total_sample_size,L];
  int<lower=1> site_weights[L];
 

  int<lower=0,upper=2*total_sample_size> topology[total_sample_size-1,3];
 
  int<lower=0,upper=total_sample_size> tip_association[total_sample_size];


  matrix[2*total_sample_size-2, total_sample_size-1 ] coal_times_to_branch_lengths;
  real<lower=0.0> seq_amp_error;
  real<lower=0.0> allele_dropout_error;
  int<lower=1> number_invariant_sites;
  int<lower=1> n_cores;

}
transformed data{
    int number_branches = 2*total_sample_size-2;
    int number_of_tips_below[2*total_sample_size-1];
    int remainder = (L) % n_cores;
    int additional_n_cols = (remainder != 0)? n_cores - (L) % n_cores: 0;
    int<lower=0> col_per_shards = (L+additional_n_cols) / n_cores ;
    
    int number_of_coalescent_times_below[2*total_sample_size-1] ;
    vector[16] genotype_tip_partials[2*total_sample_size,L];

    int shard_begin = 1;
    real sum_w = sum(to_row_vector(site_weights));
    int site_weights_extended[L+additional_n_cols]=append_array(site_weights, rep_array(1,additional_n_cols));

    array[n_cores, 16*col_per_shards*total_sample_size] real genotype_tip_partials_r;
    array[n_cores, 16*col_per_shards*total_sample_size] int total_sample_size_topology_i;
    array[n_cores] vector[0] empty_array;

    vector[16] e1= rep_vector(0.0,16);
    
    //vector[16*(L+1+additional_n_cols)*2*total_sample_size] genotype_partials;//the last  columns are needed to correct for the ASC

    matrix[16,1] genotype_partials_asc[total_sample_size];
    //matrix[16,16] genotype_partials_asc[total_sample_size];


    matrix[total_sample_size-1, 2*total_sample_size-1] indexes_nodes_below= rep_matrix(0, total_sample_size-1,2*total_sample_size-1 );
    int map_internal_node_topology_row[2*total_sample_size-1]= rep_array(0,2*total_sample_size-1);
   
    vector [3*total_sample_size-4] w;//Return non-zero values in matrix 
    int v[3*total_sample_size-4];
    int u[number_branches+1];
    int max_sample_size= max(sample_sizes);
    int cum_sample_sizes[N];
    int cum_internal_nodes[N];
    int cum_number_child_pops[N];

    e1[1]= 1.0;

    
    cum_sample_sizes= cumulative_sum_integer(sample_sizes, N);
    cum_internal_nodes= cumulative_sum_integer(number_internal_nodes, N);
    cum_number_child_pops= cumulative_sum_integer(number_child_pops, N);
 
    for(i in 1:(total_sample_size)){
       number_of_tips_below[i]=1;
       number_of_coalescent_times_below[i]=0;
    }
    for(i in 1:(total_sample_size-1)){
        int parent_node_idx = topology[i,1];
        int left_node_idx = topology[i,2];
        int right_node_idx = topology[i,3];
        int current_row = i;
        number_of_tips_below[parent_node_idx] = number_of_tips_below[left_node_idx]+ number_of_tips_below[right_node_idx];
        number_of_coalescent_times_below[parent_node_idx] = number_of_coalescent_times_below[left_node_idx]+ number_of_coalescent_times_below[right_node_idx]+1;
        
       
    }
   
    w=csr_extract_w(coal_times_to_branch_lengths);
    v=csr_extract_v(coal_times_to_branch_lengths);
    u=csr_extract_u(coal_times_to_branch_lengths);

   for (k in 1:n_cores) {
        
         
         int shard_end = shard_begin + (total_sample_size*col_per_shards*16) -1;
         int begin_tip = 1;

         total_sample_size_topology_i[k] = rep_array(1, 16*col_per_shards*total_sample_size);
         total_sample_size_topology_i[k, 1] = total_sample_size;
         total_sample_size_topology_i[k, 2:(total_sample_size)] = topology[1:(total_sample_size-1), 1];
         total_sample_size_topology_i[k, (total_sample_size+1):(2*total_sample_size-1)] = topology[1:(total_sample_size-1), 2];
         total_sample_size_topology_i[k, (2*total_sample_size):(3*total_sample_size-2)] = topology[1:(total_sample_size-1), 3];
         total_sample_size_topology_i[k, 3*total_sample_size-1] = col_per_shards;

         total_sample_size_topology_i[k, (3*total_sample_size):(3*total_sample_size+col_per_shards-1)]= site_weights_extended[((k-1)*col_per_shards+1):(k*col_per_shards)];
         for( n in 1:total_sample_size ) {//rows
          
          int end_tip = begin_tip + 16*col_per_shards -1;
          
          int tip_index = tip_association[n];

          int begin_site = begin_tip;

          for( i in ((k-1)*col_per_shards+1):(k*col_per_shards) ) {//columns
               
              int end_site = begin_site  +16 -1 ;

              if(i<=L){

                  if (genotype_tipdata[tip_index,i]!= 0){
                  //by columns
                  // genotype_tip_partials[n][1:16, i ] = tip_clv(genotype_tipdata[tip_index,i], seq_amp_error, allele_dropout_error);
                   genotype_tip_partials_r[k, begin_site:end_site ] = to_array_1d(tip_clv(genotype_tipdata[tip_index,i], seq_amp_error, allele_dropout_error));
                  
                  }
                 else{//N genotype
                   //genotype_tip_partials[n][1:16, i ]  = rep_vector(1.0, 16);
                   genotype_tip_partials_r[k, begin_site:end_site  ]  = rep_array(1.0, 16);
                 }   

              }
              else{//i>=(L+1)


                 genotype_tip_partials_r[k, begin_site:end_site ] = rep_array(1.0, 16);

              }   
              begin_site += 16;

          }

          begin_tip += 16*col_per_shards;
      }

      shard_begin += (total_sample_size*col_per_shards*16);

    }
    

    for( n in 1:total_sample_size ) {//rows
          genotype_partials_asc[n]= to_matrix(e1);
          //genotype_partials_asc[n]= identity_matrix(16);
      }

      for (i in 1:(total_sample_size-1)){

            map_internal_node_topology_row[topology[i,1]]=i;
            if (topology[i,2] > total_sample_size){
              indexes_nodes_below[i,]= indexes_nodes_below[map_internal_node_topology_row[topology[i,2]],];
            }
           if (topology[i,3] > total_sample_size){
                 indexes_nodes_below[i,]= indexes_nodes_below[i,] + indexes_nodes_below[map_internal_node_topology_row[topology[i,3]],];

           }
           indexes_nodes_below[i,topology[i,1]]=1;
      }

}
parameters{
  real<lower=0.0> hyper_parameters_exponential[N];
  //real<lower=0.0> hyper_parameters_exponential_theta;
  vector<lower=0.001>[N] gammas;
  //the last gamma corresponds to the oldest population
  
  real<lower=0.0> theta;
  real<lower=0.0>  torigin_oldest_population;
  
  simplex[N] pop_sizes_proportion;
  //simplex[total_sample_size] population_simplexes[N];
  
  vector<lower=0.0>[sum(number_internal_nodes)] population_numbers;

  vector<lower=0, upper=1>[N-1] positions_torigins_in_intervals; 

}
transformed parameters {
  //declarations

  vector[sum(number_internal_nodes)] coal_event_times_physical_time;
  vector[sum(number_internal_nodes)] coal_event_times_model_time;
  vector<lower=0.0>[N] torigins_in_model_time;
  vector<lower=0.0>[N] torigins_in_physical_time;
  vector<lower=0>[sum(number_internal_nodes)+N] pi;
  #vector<lower=0>[sum(sample_sizes)+N] pi;

  torigins_in_model_time[N] = torigin_oldest_population;
  torigins_in_physical_time[N] = torigins_in_model_time[N]* pop_sizes_proportion[N];
 

   for(i in 1:N){
       int start_idx=1;
       int end_idx=1;
       int start_coal_idx=1;
       int end_coal_idx=1;
       int no_internals ;
       
         if (i>=2){
          start_idx =cum_sample_sizes[i-1]+1+cum_number_child_pops[i-1];
          start_coal_idx= cum_internal_nodes[i-1]+1;
        }
        end_idx = cum_sample_sizes[i]+cum_number_child_pops[i];

        end_coal_idx= cum_internal_nodes[i];

       no_internals = number_internal_nodes[i];
     
       //pi[start_idx:end_idx] = append_row(head(population_simplexes[i], no_internals), sum(tail(population_simplexes[i], total_sample_size - no_internals  )));
      
       pi[start_idx:end_idx] = simplex_constrain( population_numbers[start_coal_idx:end_coal_idx]);
      
   }

   for (k in 1:(N-1)) {
      int j= N-k;
      int idx_parent_pop = indexes_father_populations[j];
      vector [number_internal_nodes[idx_parent_pop]+1] pi_father_pop;
      vector[pos_parent_child_pop_MRCA_in_parent_pop[j]] head_cum_sum;
      real lower_coal_time=0;
      real upper_coal_time;
      real unif;
      int end_idx=1;
      int start_coal_idx=1;
      int end_coal_idx=1;
      int no_internals = number_internal_nodes[idx_parent_pop-1];
      //real torigins_in_model_time_parent_population;
      int start_idx=1;
     
        if (idx_parent_pop>=2){
          start_idx = cum_sample_sizes[idx_parent_pop-1]+1+cum_number_child_pops[idx_parent_pop-1];
          start_coal_idx= cum_internal_nodes[idx_parent_pop-1]+1;
        }
       end_idx = cum_sample_sizes[idx_parent_pop]+cum_number_child_pops[idx_parent_pop];

        end_coal_idx= cum_internal_nodes[idx_parent_pop];

    
     

       pi_father_pop = pi[start_idx:end_idx];

      head_cum_sum= head(cumulative_sum(pi_father_pop), pos_parent_child_pop_MRCA_in_parent_pop[j]);
 
      
      if (pos_parent_child_pop_MRCA_in_parent_pop[j]==1){
         lower_coal_time = 0;

      }else{
        lower_coal_time = 0;// head_cum_sum[pos_parent_child_pop_MRCA_in_parent_pop[j]-1]; 
      }
   


      upper_coal_time = head_cum_sum[pos_parent_child_pop_MRCA_in_parent_pop[j]]; 
      unif= positions_torigins_in_intervals[j];

      torigins_in_physical_time[j] =  torigins_in_physical_time[idx_parent_pop] *unif*upper_coal_time ;

    
      torigins_in_model_time[j] = torigins_in_physical_time[j] / pop_sizes_proportion[j];
    
   }

   for (j in 1:N) {
        int start_idx=1;
        int end_idx=1;
        int start_coal_idx=1;
        int end_coal_idx=1;
        
        if (j>=2){
          start_idx =cum_sample_sizes[j-1]+1+cum_number_child_pops[j-1];
          start_coal_idx= cum_internal_nodes[j-1]+1;
        }
        end_idx = cum_sample_sizes[j]+cum_number_child_pops[j];

        end_coal_idx= cum_internal_nodes[j];

         
      
       coal_event_times_model_time[start_coal_idx:end_coal_idx] = head(cumulative_sum(pi[start_idx:end_idx]), number_internal_nodes[j]) *  torigins_in_model_time[j];
    
       coal_event_times_physical_time[start_coal_idx:end_coal_idx] = head(cumulative_sum(pi[start_idx:end_idx]), number_internal_nodes[j]) *  torigins_in_physical_time[j]  ;
     
    }

}
model{
  real b;
  real a;
  real c;
  real exp_2_3;
  real exp_4_3;
  real a_divide_by_3;
  real asc_correction;
  real log_lik_constant_site;
  real lik_constant_site;
  real one_over_16 = 1.0/16.0;

  //matrix[16,16] left;
  //matrix[16,16] right;
  matrix[16,1] left;
  matrix[16,1] right;
  matrix[16,1] genotype_internal_partials_asc[total_sample_size];
  //matrix[16,16] genotype_internal_partials_asc[total_sample_size];
  vector[n_cores] log_lik_per_shard;

  vector[16] partials[total_sample_size,L];

  vector[16*16*number_branches] p_matrices_vector; 
  matrix[16,16] p_matrices[number_branches]; 
  vector[number_branches] branch_lengths;
  vector[total_sample_size-1] J_diag;
  vector[number_branches] J_diag_theta;
  row_vector[16] all_ones = to_row_vector(rep_vector(1, 16));
  
  int pos = 1;
    
  hyper_parameters_exponential~gamma(rep_vector(0.001,N), rep_vector(0.001,N));
  gammas ~ exponential(hyper_parameters_exponential);


  //theta ~ exponential(1);


  //pop_sizes_proportion~ dirichlet(to_vector(sample_sizes));

  
  torigin_oldest_population~conditionalDensityTOrigin(gammas[N], sample_sizes[N]);

   // for (j in 1:N) {
   //   population_simplexes[j]~ dirichlet(append_row(rep_vector(1, number_internal_nodes[j] ), 
   //                   rep_vector(1.0/(total_sample_size - number_internal_nodes[j]), total_sample_size  - number_internal_nodes[j])));
     
   // }

  J_diag = rep_vector_times(torigins_in_model_time, number_internal_nodes);
  //this is a Jacobian from the simplexes to the  coalescent event times 
 

  target += sum(log(J_diag));



  for (j in 1:N){ 
      int start_coal_idx=1;
      int end_coal_idx=1;
      
      if (j>=2){
          start_coal_idx= cum_internal_nodes[j-1]+1;
        }

        end_coal_idx= cum_internal_nodes[j];
    
     
     target += simplex_constrain_lj( population_numbers[start_coal_idx:end_coal_idx]);

  }


 for (j in 1:(N-1)) {
    int index_father_pop = indexes_father_populations[j];
    int number_internal_nodes_parent_pop = number_internal_nodes[index_father_pop];
    vector[pos_parent_child_pop_MRCA_in_parent_pop[j]] head_cum_sum;
    int start_idx = cum_sample_sizes[index_father_pop-1]+1+cum_number_child_pops[index_father_pop-1];
    int end_idx = cum_sample_sizes[index_father_pop]+cum_number_child_pops[index_father_pop];
    int start_coal_idx=1;
    int end_coal_idx=1;
    vector[number_internal_nodes_parent_pop+1] father_pop_simplex;

   
      if (index_father_pop>=2){
          start_idx = cum_sample_sizes[index_father_pop-1]+1+cum_number_child_pops[index_father_pop-1];
          start_coal_idx= cum_internal_nodes[index_father_pop-1]+1;
        }
       end_idx = cum_sample_sizes[index_father_pop]+cum_number_child_pops[index_father_pop];

        end_coal_idx= cum_internal_nodes[index_father_pop];

   // vector[number_internal_nodes_parent_pop] father_pop_simplex = pi[index_father_pop][1:number_internal_nodes_parent_pop];
  
   
    father_pop_simplex = pi[start_idx:end_idx];
    head_cum_sum= head(cumulative_sum(father_pop_simplex), pos_parent_child_pop_MRCA_in_parent_pop[j]);

    //this if ro the change from torigins_in_model_time[N] to torigins_in_model_time_parent_population[N]

    target += log(torigins_in_model_time[index_father_pop])+log( head_cum_sum[pos_parent_child_pop_MRCA_in_parent_pop[j]]);
  
    //jacobian of the transformation from u to times origin child populations in the same time unit of the father T
    target +=  log(pop_sizes_proportion[index_father_pop])-log(pop_sizes_proportion[j]);
    //jacobian of the transformation to convert to the unit of time origin in the child population

    torigins_in_model_time[j]~conditionalDensityTOrigin(gammas[j], sample_sizes[j]);
}



coal_event_times_model_time~structured_coalescent(gammas, torigins_in_model_time,
                    torigins_in_physical_time,
                    pop_sizes_proportion, sample_sizes, number_internal_nodes, indexes_father_populations, pos_parent_child_pop_MRCA_in_parent_pop, cum_internal_nodes, 
                    number_child_pops,
                    child_pop_indexes,
                    K);

branch_lengths =   csr_matrix_times_vector(number_branches, total_sample_size-1, w,  v, u, coal_event_times_physical_time);

//p_matrices= calculate_gtr_p_matrices(frequencies_genotypes, rates, branch_lengths*theta);
for( k in 1:number_branches){
      int begin_matrix_idx = 16*16*(k-1)+1;
      int end_matrix_idx = 16*16*k;

      exp_4_3 = exp(-4.0*branch_lengths[k]*theta/3.0);
      exp_2_3 = exp(-2.0*branch_lengths[k]*theta/3.0);
      b = (1.0/16.0)*exp_4_3*(exp_2_3 - 1)*(exp_2_3 -1);
      
      p_matrices_vector[begin_matrix_idx:end_matrix_idx ] = rep_vector(b, 256);

      c = (1.0/16.0)*(1-3*exp_4_3+2*exp_2_3);

      
     for(i in 1:16){
          
          int begin_col_idx = begin_matrix_idx+16*(i-1);

          if(i==1){
            //rows 5, 6, 7, 11, 12, 13
             p_matrices_vector[begin_col_idx+4]= c;
             p_matrices_vector[begin_col_idx+5]= c;
             p_matrices_vector[begin_col_idx+6]= c;
             p_matrices_vector[begin_col_idx+10]= c;
             p_matrices_vector[begin_col_idx+11]= c;
             p_matrices_vector[begin_col_idx+12]= c;

          }
          else if(i==2){
           //rows 5, 8, 9, 11, 14, 15
             p_matrices_vector[begin_col_idx+4]= c;
             p_matrices_vector[begin_col_idx+7]= c;
             p_matrices_vector[begin_col_idx+8]= c;
             p_matrices_vector[begin_col_idx+10]= c;
             p_matrices_vector[begin_col_idx+13]= c;
             p_matrices_vector[begin_col_idx+14]= c;
          }
          else if(i==3){
            //rows 6, 8, 10, 12, 14, 16
             p_matrices_vector[begin_col_idx+5]= c;
             p_matrices_vector[begin_col_idx+7]= c;
             p_matrices_vector[begin_col_idx+9]= c;
             p_matrices_vector[begin_col_idx+11]= c;
             p_matrices_vector[begin_col_idx+13]= c;
             p_matrices_vector[begin_col_idx+15]= c;
          }
          else if(i==4){
            //rows 7, 9, 10, 13, 15, 16
             p_matrices_vector[begin_col_idx+6]= c;
             p_matrices_vector[begin_col_idx+8]= c;
             p_matrices_vector[begin_col_idx+9]= c;
             p_matrices_vector[begin_col_idx+12]= c;
             p_matrices_vector[begin_col_idx+14]= c;
             p_matrices_vector[begin_col_idx+15]= c;
    
            
          }
          else if(i==5){
              //rows 1, 2, 6, 7, 14, 15
             p_matrices_vector[begin_col_idx+0]= c;
             p_matrices_vector[begin_col_idx+1]= c;
             p_matrices_vector[begin_col_idx+5]= c;
             p_matrices_vector[begin_col_idx+6]= c;
             p_matrices_vector[begin_col_idx+13]= c;
             p_matrices_vector[begin_col_idx+14]= c;
      
            
          }
          else if(i==6){
               //rows 1, 3, 5, 7, 8, 16
             p_matrices_vector[begin_col_idx+0]= c;
             p_matrices_vector[begin_col_idx+2]= c;
             p_matrices_vector[begin_col_idx+4]= c;
             p_matrices_vector[begin_col_idx+6]= c;
             p_matrices_vector[begin_col_idx+7]= c;
             p_matrices_vector[begin_col_idx+15]= c;
  
            
          }
          else if(i==7){
             //rows 1, 4, 5, 6, 9, 10
             p_matrices_vector[begin_col_idx+0]= c;
             p_matrices_vector[begin_col_idx+3]= c;
             p_matrices_vector[begin_col_idx+4]= c;
             p_matrices_vector[begin_col_idx+5]= c;
             p_matrices_vector[begin_col_idx+8]= c;
             p_matrices_vector[begin_col_idx+9]= c;

            
          }
          else if(i==8){
            //rows 2, 3, 6, 9, 11, 16
             p_matrices_vector[begin_col_idx+1]= c;
             p_matrices_vector[begin_col_idx+2]= c;
             p_matrices_vector[begin_col_idx+5]= c;
             p_matrices_vector[begin_col_idx+8]= c;
             p_matrices_vector[begin_col_idx+10]= c;
             p_matrices_vector[begin_col_idx+15]= c;
          
            
          }
          else if(i==9){
            //rows 2, 4, 7, 8, 10, 11
             p_matrices_vector[begin_col_idx+1]= c;
             p_matrices_vector[begin_col_idx+3]= c;
             p_matrices_vector[begin_col_idx+6]= c;
             p_matrices_vector[begin_col_idx+7]= c;
             p_matrices_vector[begin_col_idx+9]= c;
             p_matrices_vector[begin_col_idx+10]= c;
    
            
          }
          else if(i==10){
            //rows 3, 4, 7, 9, 12, 14
             p_matrices_vector[begin_col_idx+2]= c;
             p_matrices_vector[begin_col_idx+3]= c;
             p_matrices_vector[begin_col_idx+6]= c;
             p_matrices_vector[begin_col_idx+8]= c;
             p_matrices_vector[begin_col_idx+11]= c;
             p_matrices_vector[begin_col_idx+13]= c;
            
          }
          else if(i==11){
            //rows 1, 2, 8, 9, 12, 13
             p_matrices_vector[begin_col_idx+0]= c;
             p_matrices_vector[begin_col_idx+1]= c;
             p_matrices_vector[begin_col_idx+7]= c;
             p_matrices_vector[begin_col_idx+8]= c;
             p_matrices_vector[begin_col_idx+11]= c;
             p_matrices_vector[begin_col_idx+12]= c;
          
            
          }
          else if(i==12){
            //rows 1, 3, 10, 11, 13, 14
             p_matrices_vector[begin_col_idx+0]= c;
             p_matrices_vector[begin_col_idx+2]= c;
             p_matrices_vector[begin_col_idx+9]= c;
             p_matrices_vector[begin_col_idx+10]= c;
             p_matrices_vector[begin_col_idx+12]= c;
             p_matrices_vector[begin_col_idx+13]= c;
      
            
          }
          else if(i==13){
            //rows 1, 4, 11, 12, 15, 16
             p_matrices_vector[begin_col_idx+0]= c;
             p_matrices_vector[begin_col_idx+3]= c;
             p_matrices_vector[begin_col_idx+10]= c;
             p_matrices_vector[begin_col_idx+11]= c;
             p_matrices_vector[begin_col_idx+14]= c;
             p_matrices_vector[begin_col_idx+15]= c;
            
          }
          else if(i==14){
            //rows 2, 3, 5, 10, 12, 15
             p_matrices_vector[begin_col_idx+1]= c;
             p_matrices_vector[begin_col_idx+2]= c;
             p_matrices_vector[begin_col_idx+4]= c;
             p_matrices_vector[begin_col_idx+9]= c;
             p_matrices_vector[begin_col_idx+11]= c;
             p_matrices_vector[begin_col_idx+14]= c;
            
          }
          else if(i==15){
            //rows 2, 4, 5, 13, 14, 16
             p_matrices_vector[begin_col_idx+1]= c;
             p_matrices_vector[begin_col_idx+3]= c;
             p_matrices_vector[begin_col_idx+4]= c;
             p_matrices_vector[begin_col_idx+12]= c;
             p_matrices_vector[begin_col_idx+13]= c;
             p_matrices_vector[begin_col_idx+15]= c;
    
            
          }
          else{//i==16
            //rows 3, 4, 6, 8, 13, 15
             p_matrices_vector[begin_col_idx+2]= c;
             p_matrices_vector[begin_col_idx+3]= c;
             p_matrices_vector[begin_col_idx+5]= c;
             p_matrices_vector[begin_col_idx+7]= c;
             p_matrices_vector[begin_col_idx+12]= c;
             p_matrices_vector[begin_col_idx+14]= c;
          
          }
          
        p_matrices_vector[begin_col_idx+i-1] = (1.0/16.0)*exp_4_3*(exp_2_3 +3)*(exp_2_3 +3);
        p_matrices[k] = to_matrix(p_matrices_vector[begin_matrix_idx:end_matrix_idx],16,16);
      }
    }


  log_lik_per_shard = map_rect(felsenstein_block, p_matrices_vector, empty_array, genotype_tip_partials_r, total_sample_size_topology_i);
    
    
  for( n in 1:(total_sample_size-1) ) {//rows

           left  = (topology[n,2] > total_sample_size)? genotype_internal_partials_asc[topology[n,2]-total_sample_size] :  genotype_partials_asc[topology[n,2]];
           right = (topology[n,3] > total_sample_size)? genotype_internal_partials_asc[topology[n,3]-total_sample_size] :  genotype_partials_asc[topology[n,3]];
          
           genotype_internal_partials_asc[topology[n,1]-total_sample_size] = (p_matrices[2*n-1]*left) .* (p_matrices[2*n]*right);
    }
       
  genotype_internal_partials_asc[total_sample_size] =  one_over_16 *genotype_internal_partials_asc[topology[total_sample_size-1,1]-total_sample_size];
  lik_constant_site = sum(all_ones*genotype_internal_partials_asc[total_sample_size]);
  log_lik_constant_site = log(lik_constant_site);
  //asc_correction = -1.0* sum_w * log(1.0-lik_constant_site);
  asc_correction = number_invariant_sites* log_lik_constant_site;
  target += sum(log_lik_per_shard)+asc_correction;

}
generated quantities {
  vector<lower=0.001>[N] unscaled_gammas;
 

  for (i in 1:N) {
    unscaled_gammas[i] = gammas[i]/pop_sizes_proportion[i] ;
  }

}
