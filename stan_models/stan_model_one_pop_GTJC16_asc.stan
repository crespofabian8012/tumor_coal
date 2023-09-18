functions{
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
  matrix[] calculate_gtr_p_matrices_gt16(vector freqs, vector rates, vector blens){
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
		// symmetric rate matrix
		                  /*AA*//*CC*//*GG*/ /*TT*/   /*AC*/    /*AG*/     /*AT*/  /*CG*/     /*CT*/   /*GT*/    /*CA*/    /*GA*/    /*TA*/     /*GC*/    /*TC*/     /*TG*/
		matrix[16,16] R = [[0.0,  0.0,  0.0,  0.0, rates[1], rates[2], rates[3],     0.0,       0.0,     0.0,  rates[1], rates[2], rates[3],     0.0,       0.0,      0.0],  /* AA */
						           [0.0,  0.0,  0.0,  0.0, rates[1],      0.0,      0.0, rates[4], rates[5],     0.0,  rates[1],      0.0,      0.0, rates[4], rates[5],      0.0],  /* CC */
						           [0.0,  0.0,  0.0,  0.0,      0.0, rates[2],      0.0, rates[4],      0.0, rates[6],      0.0, rates[2],      0.0, rates[4],      0.0, rates[6]],  /* GG */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0, rates[3],      0.0, rates[5], rates[6],      0.0,      0.0, rates[3],      0.0, rates[5], rates[6]],  /* TT */
						           [0.0,  0.0,  0.0,  0.0,      0.0, rates[4], rates[5], rates[2], rates[3],      0.0,      0.0,      0.0,      0.0, rates[2], rates[3],      0.0],  /* AC */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0, rates[6], rates[1],      0.0, rates[3],      0.0,      0.0,      0.0,      0.0,      0.0, rates[3]],  /* AG */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0, rates[1], rates[2],      0.0,      0.0,      0.0,      0.0,      0.0,      0.0],  /* AT */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0, rates[6], rates[5], rates[2],      0.0,      0.0,      0.0,      0.0, rates[5]],  /* CG */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[4], rates[3],      0.0,      0.0,      0.0,      0.0,      0.0],  /* CT */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[3],      0.0, rates[5],      0.0,      0.0],  /* GT */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[4], rates[5],      0.0,      0.0,      0.0],  /* CA */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[6], rates[1],      0.0,      0.0],  /* GA */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[1], rates[2]],  /* TA */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[6],      0.0],  /* GC */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, rates[4]],  /* TC */
						           [0.0,  0.0,  0.0,  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0]]  /* TG */
						           ;



		R= R+R';
		s = 0.0;
		index = 1;
		Q = R * diag_matrix(freqs);
		for (i in 1:4) {
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
   real random_time_MRCA_below_time_origin(real value, real time_origin_in_oldest_population_time, vector sorted_coal_times_in_oldest_population_time,
                   int[,] topology,int[] number_of_tips_below){

                    real result;
                    int number_candidate_MRCAs=0;
                    int idx_left;
                    int idx_right;
                    int idx_left_below;
                    int idx_right_below;
                    int idx_MRCA;
                    int current_pos=1;
                    int idx_first_coal_time_above_torigin=get_index_coal_time_above_MRCA_in_oldest_population_time_less_than( time_origin_in_oldest_population_time,
                                                          sorted_coal_times_in_oldest_population_time);
                    //vector [idx_first_coal_time_above_torigin-1] coal_times_below_torigin=
                    //segment(sorted_coal_times_in_oldest_population_time, 1, idx_first_coal_time_above_torigin-1);
                    int indexes_nodes_below_torigin[idx_first_coal_time_above_torigin-1]  =topology[1:(idx_first_coal_time_above_torigin-1),1];
                    vector [idx_first_coal_time_above_torigin-1] indexes_candidate_MRCAs;
                    vector [idx_first_coal_time_above_torigin-1] coal_time_candidate_MRCAs;

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

                               if (number_of_tips_below[idx_left]==1){
                                  number_candidate_MRCAs=number_candidate_MRCAs+1;
                                  indexes_candidate_MRCAs[current_pos]=idx_left;
                                  coal_time_candidate_MRCAs[current_pos] =  0.0;
                                  current_pos=current_pos+1;

                               }
                            }
                          if (idx_right_below!=-1 ){
                            number_candidate_MRCAs=number_candidate_MRCAs+1;
                            indexes_candidate_MRCAs[current_pos]=idx_right ;
                            coal_time_candidate_MRCAs[current_pos] = sorted_coal_times_in_oldest_population_time[idx_right_below];
                            current_pos=current_pos+1;
                           }
                           else{

                               if (number_of_tips_below[idx_right]==1){
                                  number_candidate_MRCAs=number_candidate_MRCAs+1;
                                  indexes_candidate_MRCAs[current_pos]=idx_right;
                                  coal_time_candidate_MRCAs[current_pos] =  0.0;
                                  current_pos=current_pos+1;

                               }
                            }
                       }
            // print("number_candidate_MRCAs=", number_candidate_MRCAs);
            // print("indexes_candidate_MRCAs=", indexes_candidate_MRCAs);
            // print("coal_time_candidate_MRCAs=", coal_time_candidate_MRCAs);

            if (number_candidate_MRCAs >0)
               result =- log(number_candidate_MRCAs);
            else
              reject("no candidate MRCAs!");

            //print("result=", result);
            return(result);
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
         real numerator;
         real denominator;
         real logh;
         int order_father_pop;
         real temp;
         real result;
         vector[N] ordered_torigins_in_model_time_oldest_population=sort_asc(torigins_in_model_time_oldest_population_par);
         order_father_pop = find((x_father/pop_sizes_proportion[index_oldest_population]) * father_torigin, ordered_torigins_in_model_time_oldest_population );
         if (order != -1 && order_father_pop!=-1 && order< order_father_pop){

                for(i in (order+1):N){
                  int pos = find((pop_sizes_proportion[index_oldest_population]/pop_sizes_proportion[i] )*torigins_in_model_time_oldest_population_par[i], torigins_in_model_time);
                  denominator = denominator + log( pop_sizes_proportion[pos]);
                  temp= torigin * x / pop_sizes_proportion[pos];
                  logh=log_h(temp, torigins_in_model_time[pos], deltas[pos], K);
                  denominator = denominator  + temp;
                }
         }

         numerator = numerator + log( x_father);
         temp= torigin * x / x_father;
         logh=log_h(temp, father_torigin, father_delta, K);
         numerator = numerator  + temp;

         result= numerator / denominator;
        return(result);
  }
   real log_likelihood_subpopulation(int order, int sample_size, real K, real delta, real torigin, int[] positions,
                                     vector sorted_coal_time_with_previous_torigins_in_model_time, vector previous_torigins_in_model_time){

    real log_lik = 0.0;
    real current_time;
    int  alive_cells;
    real term_only_after_first_coal_event;
    real zero =0.0;
    int currentCoalescentEventInThisEpoch=-1;
    vector[rows(sorted_coal_time_with_previous_torigins_in_model_time)+1] padded_coalescent_times =
                                                  append_row(zero, sorted_coal_time_with_previous_torigins_in_model_time);
    real last_event_time_before_migration=0.0;


   //print("previous_torigins_in_model_time=",previous_torigins_in_model_time);

    alive_cells = sample_size;//alive_cells goes from sample_size downto 2
    for(j in 1:rows(padded_coalescent_times)){
      current_time = padded_coalescent_times[j];

      if (fabs(current_time - torigin)<0.00000001)
        return log_lik;

      if (current_time < 0.0) continue;//those are not real inmigrants

      if (order > 1 && rows(previous_torigins_in_model_time) > 0 && find(current_time, previous_torigins_in_model_time)!= -1){//= -1

        if (alive_cells >1){
            log_lik+= log_lik_no_coalescent_between_times( last_event_time_before_migration,  current_time,  alive_cells,  torigin,  delta,  K);
            alive_cells= alive_cells+1;
            currentCoalescentEventInThisEpoch=1;
            last_event_time_before_migration= current_time;
            }
      }
      else{
         //alive_cells =(sample_size-j+1);//alive_cells goes from sample_size downto 2
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
      return(log_lik);
    }
    vector structured_coalescent_rng(vector deltas,vector torigins_in_model_time, vector torigins_in_model_time_oldest_population,
                                   vector pop_sizes_proportion, int[] sample_sizes,  real K ){
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

    // print("inside structured_coalescent_rng=",sample_sizes);
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
  real structured_coalescent_lpdf(vector coal_event_times_model_time_oldest_population, vector deltas, vector torigins_in_model_time,
                              vector sorted_torigins_in_model_time_oldest_population,  vector pop_sizes_proportion, int[] sample_sizes,  real K){

    real log_lik=0.0;
    int cum_sample_size=0;
    int cum_number_coal=0;
    int first_position_coal_event=1;
    int N= rows(deltas);
    int pos_next_torigin = 0;


    cum_sample_size=0;

    for(i in 1:num_elements(sample_sizes)){
         int order=i;
         int number_events= sample_sizes[i]-1;
         int positions[sample_sizes[i]] = rep_array(1,sample_sizes[i] );
         vector[number_events] current_coal_event_times_with_immigrants;

         if (i==1){
           current_coal_event_times_with_immigrants[1:number_events] = (coal_event_times_model_time_oldest_population[ 1:number_events]) ;
           log_lik= log_lik+ log_likelihood_subpopulation(i, sample_sizes[i], K, deltas[i], torigins_in_model_time[i], positions,
                                      current_coal_event_times_with_immigrants,
                                      rep_vector(1,1));
         // print("i=1 log_lik=",log_lik );
         }

    }
    return(log_lik);
  }

  int get_sample_size(vector sorted_coal_times, real MRCA_time_in_model_time_oldest_population, int[] number_of_tips_below, int[,] topology){
      int result=0;
      int pos = find(MRCA_time_in_model_time_oldest_population, sorted_coal_times);
      //print("pos=", pos);
      //print("MRCA_time_in_model_time_oldest_population=", MRCA_time_in_model_time_oldest_population);

      if (pos >0 && pos<= num_elements(number_of_tips_below))
         result = number_of_tips_below[topology[pos,1]];
      return(result);
  }
  vector  compute_probab_vector(int total_sample_size, int number_candidate_MRCAs, vector candidate_time_MRCAs_in_model_time,
                                vector candidate_time_MRCAs_in_model_time_oldest_population,
                                matrix Q_coeff_hypoexponential, real delta, real torigin_in_model_time, int[] number_of_tips_below,
                                vector sorted_coal_times_in_oldest_population_time,  int[,] topology) {
      real cum_sum=0.0;
      vector[number_candidate_MRCAs] result;
     // print("candidate_time_MRCAs_in_model_time_oldest_population=",candidate_time_MRCAs_in_model_time_oldest_population);

      for(j in 1:number_candidate_MRCAs){
        int sample_size= get_sample_size(sorted_coal_times_in_oldest_population_time,  candidate_time_MRCAs_in_model_time_oldest_population[j],
                         number_of_tips_below, topology);
        int  size_sub_matrix= total_sample_size-sample_size+1;
        matrix [(sample_size-1),(sample_size-1)]  Q;
        real current_time = candidate_time_MRCAs_in_model_time[j];
        // print("current_time=", current_time);
        // print("sample_size=", sample_size);
        // print("size_sub_matrix=", size_sub_matrix);
        Q = Q_coeff_hypoexponential[size_sub_matrix:(total_sample_size-1), size_sub_matrix:(total_sample_size-1)];
        //print("Q=", Q);
        Q = matrix_exp(Q * current_time );
        //print("Q exp=", Q);
        //result[j] = 1 * Q[1,(sample_size-1)];
        result[j] =  Q[1,(sample_size-1)];
        //print("result[j]=", result[j]);
        cum_sum = cum_sum + result[j];
       // print("cum_sum=", cum_sum);
     }
    result = (1.0/cum_sum)*result;
    return(result);
  }
  vector compute_branch_lengths_from_coal_times(int N, int[] map_internal_node_topology_row, int total_sample_size, int[] sample_sizes, vector coal_times_in_model_time_oldest_population, int[,] topology,
              vector torigins_in_model_time_oldest_population){

    vector[2*total_sample_size-2] branch_lengths= rep_vector(0,2*total_sample_size-2);
    vector [total_sample_size-1] only_coal_times=rep_vector(0,total_sample_size-1);
    //vector [rows(coal_times_in_model_time_oldest_population)] sorted_coal_times = sort_asc(coal_times_in_model_time_oldest_population);
    int pos=1;
    int left;
    int right;
    int idx_left;
    int idx_right;
    int cum_sample_size=0;
    int cum_number_coal=0;
    int first_position_coal_event=1;
    int number_inmigrants[N]=rep_array(0,N);
    int pos_next_torigin = 0;
    int is_torigin=1;
    int pos_torigin;
    int current_pop=1;
    only_coal_times = coal_times_in_model_time_oldest_population;
    for(i in 1:N){
         int j=1;
         while(pos<=(total_sample_size-1) && coal_times_in_model_time_oldest_population[j]<torigins_in_model_time_oldest_population[i])
         {
             if (find(coal_times_in_model_time_oldest_population[j], torigins_in_model_time_oldest_population)==-1)
             {
                // print("pos=", pos);
                 only_coal_times[pos]=coal_times_in_model_time_oldest_population[j];
                 //print("only_coal_times[pos]=", only_coal_times[pos]);
                 pos=pos+1;
             }
             j=j+1;
         }
    }

    pos=1;
    for(i in 1:(total_sample_size-1))
    {
        left=topology[i,2];
        right=topology[i,3];

        if (left <= total_sample_size){
           //branch_lengths[pos]= coal_times_in_model_time_oldest_population[current_pop,j]
            branch_lengths[pos] = only_coal_times[i];
             pos+=1;
        }
        else{
           idx_left= find_integer(left,  topology[1:(total_sample_size-1-1),1]);
          // idx_left= map_internal_node_topology_row[left];
           if (idx_left!=-1 && only_coal_times[i] -only_coal_times[idx_left] >0){
               branch_lengths[pos] = only_coal_times[i]-only_coal_times[idx_left];
               pos+=1;
             }
            else{
                reject("the branch lengths must be positive");

            }

        }
        if (right <= total_sample_size){
           branch_lengths[pos] = only_coal_times[i];
           pos+=1;
        }
        else{
           idx_right= find_integer(right,  topology[1:(total_sample_size-1-1),1]);
           //idx_right= map_internal_node_topology_row[right];
           //print("idx_right=", idx_right);
           if (idx_right!=-1 && only_coal_times[i] -only_coal_times[idx_right] >0){
              branch_lengths[pos] = only_coal_times[i]-only_coal_times[idx_right];
              pos+=1;
           }
           else{
                reject("the branch lengths must be positive");

            }
        }
    }
    return(branch_lengths);
  }
  int is_in(int pos,int[] pos_var) {
   
    for (p in 1:(size(pos_var))) {
       if (pos_var[p]==pos) {
       // can return immediately, as soon as find a match
          return 1;
       } 
    }
    return 0;
  }
}
data{
  real<lower=0.0, upper=2.0> K;//parameter for the model M_K
  int<lower=1> N;//number of populations
  int<lower=2> total_sample_size;
  int <lower=0> L;// alignment length
  int<lower=0,upper=16> genotype_tipdata[total_sample_size,L];
  int<lower=1> site_weights[L];
  int<lower=0,upper=2*total_sample_size> topology[total_sample_size-1,3];
  int<lower=0,upper=total_sample_size> tip_association[total_sample_size];

  matrix[2*total_sample_size-2, total_sample_size-1 ] coal_times_to_branch_lengths;
  real<lower=0.0> seq_amp_error;
  real<lower=0.0> allele_dropout_error;

}
transformed data{
    int number_branches = 2*total_sample_size-2;
    int number_of_tips_below[2*total_sample_size-1];

    int number_of_coalescent_times_below[2*total_sample_size-1] ;
    #vector[16] genotype_tip_partials[2*total_sample_size,L];
    matrix[16,L+16] genotype_tip_partials[2*total_sample_size];//the last 16 columns are needed to correct for the ASC 
    int tip_index;

    matrix[total_sample_size-1, 2*total_sample_size-1] indexes_nodes_below= rep_matrix(0, total_sample_size-1,2*total_sample_size-1 );
    int map_internal_node_topology_row[2*total_sample_size-1]= rep_array(0,2*total_sample_size-1);
    //matrix[number_branches, total_sample_size-1 ] coal_times_to_branch_lengths= rep_matrix(0, number_branches, total_sample_size-1);
    vector [3*total_sample_size-4] w;
    int v[3*total_sample_size-4];
    int u[number_branches+1];

    for(i in 1:(total_sample_size)){
       number_of_tips_below[i]=1;
       number_of_coalescent_times_below[i]=0;
    }
      for(i in 1:(total_sample_size-1)){
       number_of_tips_below[topology[i,1]]=number_of_tips_below[topology[i,2]]+ number_of_tips_below[topology[i,3]];
       number_of_coalescent_times_below[topology[i,1]]=number_of_coalescent_times_below[topology[i,2]]+ number_of_coalescent_times_below[topology[i,3]]+1;
       map_internal_node_topology_row[topology[i,1]]=i; 
        coal_times_to_branch_lengths[2*i-1, i]=1.0;
        coal_times_to_branch_lengths[2*i, i]=1.0;
}

 
    for(i in 1:(total_sample_size-1)){
        int parent_node_idx = topology[i,1];
        int left_node_idx = topology[i,2];
        int right_node_idx = topology[i,3];
        int current_row = i;
        number_of_tips_below[parent_node_idx] = number_of_tips_below[left_node_idx]+ number_of_tips_below[right_node_idx];
        number_of_coalescent_times_below[parent_node_idx] = number_of_coalescent_times_below[left_node_idx]+ number_of_coalescent_times_below[right_node_idx]+1;
        
       
    }

 for(i in 1:(total_sample_size-1))
{
        if(topology[i,2] >total_sample_size){
           coal_times_to_branch_lengths[2*i-1,  map_internal_node_topology_row[topology[i,2]]]= -1.0;
        }
        if(topology[i,3] >total_sample_size){
          
           coal_times_to_branch_lengths[2*i,  map_internal_node_topology_row[topology[i,3]]]= -1.0;
        }
    }



    w=csr_extract_w(coal_times_to_branch_lengths);
    v=csr_extract_v(coal_times_to_branch_lengths);
    u=csr_extract_u(coal_times_to_branch_lengths);

    for( n in 1:total_sample_size ) {//rows
          genotype_tip_partials[n]= rep_matrix(0.0, 16, L+16);
          tip_index = tip_association[n];
          for( i in 1:L ) {//columns
                 //print("i=", i);
                 //print("genotype_tipdata[tip_index,i]=",genotype_tipdata[tip_index,i]);
                 if (genotype_tipdata[tip_index,i]!= 0){
                   genotype_tip_partials[n][1:16, i ] = tip_clv(genotype_tipdata[tip_index,i], seq_amp_error, allele_dropout_error);
                   //print("tip_clv=",genotype_tip_partials[n][1:16, i ]);
                   //print("sum=",sum(genotype_tip_partials[n][1:16, i ]));
                 }
                 else{
                   genotype_tip_partials[n][1:16, i ]  = rep_vector(1.0, 16);
                 }
                   
              
          }
          //print("n=",n);
          //print("genotype_tip_partials[n]=",genotype_tip_partials[n]);
          for( i in (L+1):(L+16) ){

             genotype_tip_partials[n][1:16, i ]  = rep_vector(0.0, 16);
             genotype_tip_partials[n][i-L, i ]  = 1.0;

          }
          //print("n=",n);
          //print("genotype_tip_partials[n]=",genotype_tip_partials[n]);
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
  real<lower=0.0> hyper_parameters_exponential;
  //real<lower=0.0> hyper_parameters_exponential_theta;
  real<lower=0.01> gamma;
  real<lower=0.0> torigins_in_model_time_oldest_population;
  real<lower=0.0> theta;

  simplex[total_sample_size] simplex1;
  //simplex[16] frequencies_genotypes;
  //simplex[6] rates;
}
transformed parameters {
  //declarations
  positive_ordered[total_sample_size-1] coal_event_times_model_time_oldest_population;
  coal_event_times_model_time_oldest_population=head(cumulative_sum(simplex1),total_sample_size-1) * torigins_in_model_time_oldest_population;
}
model{
  real b;
  real a;
  real c;
  real exp_2_3;
  real exp_4_3;
  real a_divide_by_3;
  real asc_correction;
  real sum_w;
  #vector[16] left;
  #vector[16] right;
  matrix[16,L+16] left;
  matrix[16,L+16] right;
  matrix[16,L+16] partials[L+16]; //the last  column are needed to correct for the ASC
  #vector[16] partials[total_sample_size,L];  // partial probabilities for the S tips and S-1 internal nodes
  matrix[16,16] p_matrices[number_branches]; // finite-time transition matrices for each branch
  vector[number_branches] branch_lengths;
  vector[total_sample_size-1] J_diag;
  vector[number_branches] J_diag_theta;
  row_vector[16] all_ones = to_row_vector(rep_vector(1, 16));

  hyper_parameters_exponential ~ gamma(0.001,0.001);
  gamma ~ exponential(hyper_parameters_exponential);
  theta ~ exponential(1);

  torigins_in_model_time_oldest_population~conditionalDensityTOrigin(gamma, total_sample_size);


  J_diag = rep_vector(torigins_in_model_time_oldest_population, total_sample_size-1);

  target += sum(log(J_diag));


  coal_event_times_model_time_oldest_population~structured_coalescent([gamma]', [torigins_in_model_time_oldest_population]',
                    [torigins_in_model_time_oldest_population]',
                    [1.0]', {total_sample_size}, K);

 
  branch_lengths =   csr_matrix_times_vector(number_branches, total_sample_size-1, w,  v, u, coal_event_times_model_time_oldest_population);



 for( k in 1:number_branches){
      exp_4_3 = exp(-4.0*branch_lengths[k]*theta/3.0);
      exp_2_3 = exp(-2.0*branch_lengths[k]*theta/3.0);
      b = (1.0/16.0)*exp_4_3*(exp_2_3 - 1)*(exp_2_3 -1);
      p_matrices[k] = rep_matrix(b, 16, 16);
      c = (1.0/16.0)*(1-3*exp_4_3+2*exp_2_3);
 
     for(i in 1:16){
       
          if (i== 5 || i== 6 || i== 7 || i== 11 || i== 12 || i== 13){
            p_matrices[k][i,1]= c;
           }
       
          if (i== 5 || i== 8 || i== 9 || i== 11 || i== 14 || i== 15){
             p_matrices[k][i,2]= c;

           }
       
          if (i== 6 || i== 8 || i== 10 || i== 12 || i== 14 || i== 16){
             p_matrices[k][i,3]= c;
           }
       
          if (i== 7 || i== 9 || i== 10 || i== 13 || i== 15 || i== 16){
             p_matrices[k][i,4]= c;
           }
     
          if (i== 1 || i== 2 || i== 6 || i== 7 || i== 14 || i== 15){
             p_matrices[k][i,5]= c;
           }
        
          if (i== 1 || i== 3 || i== 5 || i== 7 || i== 8 || i== 16){
             p_matrices[k][i,6]= c;
           
           }
         
          if (i== 1 || i== 4 || i== 5 || i== 6 || i== 9 || i== 10){
             p_matrices[k][i,7]= c;
            
           }
    
          if (i== 2 || i== 3 || i== 6 || i== 9 || i== 11 || i== 16){
             p_matrices[k][i,8]= c;
           }
       
          if (i== 2 || i== 4 || i== 7 || i== 8 || i== 10 || i== 11){
             p_matrices[k][i,9]= c;
  
           }
       
          if (i== 3 || i== 4 || i== 7 || i== 9 || i== 12 || i== 14){
             p_matrices[k][i,10]= c;
           }
         
          if (i== 1 || i== 2 || i== 8 || i== 9 || i== 12 || i== 13){
             p_matrices[k][i,11]= c;
            
           }
      
          if (i== 1 || i== 3 || i== 10 || i== 11 || i== 13 || i== 14){
             p_matrices[k][i,12]= c;
        
           }
        
          if (i== 1 || i== 4 || i== 11 || i== 12 || i== 15 || i== 16){
           p_matrices[k][i,13]= c;
           
           }
       
          if (i== 2 || i== 3 || i== 5 || i== 10 || i== 12 || i== 15){
             p_matrices[k][i,14]= c;
           
           }
       
          if (i== 2 || i== 4 || i== 5 || i== 13 || i== 14 || i== 16){
             p_matrices[k][i,15]= c;
            
           }
     
          if (i== 3 || i== 4 || i== 6 || i== 8 || i== 13 || i== 15){
             p_matrices[k][i,16]= c;
           }
        p_matrices[k][i,i]= (1.0/16.0)*exp_4_3*(exp_2_3 +3)*(exp_2_3 +3);
      }
    }


		for( n in 1:(total_sample_size-1) ) {//rows

          left  = (topology[n,2] > total_sample_size)? partials[topology[n,2]-total_sample_size] :  genotype_tip_partials[topology[n,2]];
          right = (topology[n,3] > total_sample_size)? partials[topology[n,3]-total_sample_size] :  genotype_tip_partials[topology[n,3]];
          
          partials[topology[n,1]-total_sample_size] = (p_matrices[2*n-1]*left) .* (p_matrices[2*n]*right);
        }
       
        	partials[2*total_sample_size-total_sample_size] = partials[topology[total_sample_size-1,1]-total_sample_size] .* rep_matrix(1.0/16.0, 16, L+16);
       
        sum_w = sum(to_row_vector(site_weights));
      
         asc_correction = -1.0* sum_w * log(1.0-sum(all_ones*block(partials[total_sample_size], 1, L+1, 16, 1)));
      

        target += sum(to_row_vector(site_weights).*log(all_ones*block(partials[total_sample_size], 1, 1, 16, L)))+ asc_correction;
  
 
}
